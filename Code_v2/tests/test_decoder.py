"""Independent regression checks for the Paulimer integration.
Run: python -m unittest discover -s tests -v
"""
import itertools
from types import SimpleNamespace
import unittest
import numpy as np
from pauli_compiler.decoding.paulimer_backend import derive_decoder
from pauli_compiler.pipeline import compile_operational_measurement
from pauli_compiler.circuits.model import ExtractedCircuit, CircuitOperation


def fixture(pauli, ops):
    sign = -1 if pauli.startswith('-') else 1
    word = pauli.lstrip('+-')
    n = len(word)
    return SimpleNamespace(
        measurement=SimpleNamespace(pauli=word, sign=sign),
        circuit=ExtractedCircuit(tuple(CircuitOperation(name, tuple(wires)) for name,wires in ops),
                                 n, tuple(range(n)), tuple(range(n)), 1))


class PaulimerDecoderTests(unittest.TestCase):
    def test_redundant_measurement_and_record_validation(self):
        decoder = derive_decoder(fixture('Z', [('M', [0]), ('M', [0])]))
        self.assertEqual(decoder.outcome_bit([0, 0]), 0)
        self.assertEqual(decoder.outcome_bit([1, 1]), 1)
        self.assertTrue(decoder.checks)
        with self.assertRaises(ValueError):
            decoder.outcome_bit([0, 1])
        with self.assertRaises(ValueError):
            decoder.outcome_bit([0])
        with self.assertRaises(ValueError):
            decoder.outcome_bit([2, 2])

    def test_signed_y_transpose_and_constant(self):
        decoder = derive_decoder(fixture('-Y', [('MY', [0])]))
        self.assertEqual(decoder.outcome_bit([0]), 1)
        self.assertEqual(decoder.outcome_bit([1]), 0)
        self.assertFalse(any(e.mask for e in decoder.x_corrections + decoder.z_corrections))
        decoder = derive_decoder(fixture('Z', [('X', [0]), ('M', [0]), ('X', [0])]))
        self.assertEqual(decoder.outcome_bit([0]), 1)
        self.assertEqual(decoder.outcome_bit([1]), 0)

    def test_overmeasurement_rejected(self):
        # Has correct ZZ outcome parity but destroys coherence inside its eigenspaces.
        with self.assertRaises(ValueError):
            derive_decoder(fixture('ZZ', [('M', [0]), ('M', [1])]))

    def test_wrong_observable_rejected(self):
        with self.assertRaises(ValueError):
            derive_decoder(fixture('X', [('M', [0])]))

    def test_nonpauli_output_correction_rejected(self):
        with self.assertRaises(ValueError):
            derive_decoder(fixture('Z', [('M', [0]), ('H', [0])]))

    def test_reset_of_input_rejected(self):
        with self.assertRaises(ValueError):
            derive_decoder(fixture('Z', [('R', [0]), ('M', [0])]))

    def test_constant_feedback(self):
        decoder = derive_decoder(fixture('Z', [('M', [0]), ('X', [0])]))
        self.assertTrue(decoder.x_corrections[0].constant)
        for bit in (0, 1):
            self.assertEqual(decoder.outcome_bit([bit]), bit)

    def test_exhaustive_dense_instrument(self):
        # Independent, non-stabilizer input and dense matrix calculation.
        # Check each branch's unnormalized output, including Born probabilities,
        # and then sum over records for each global measurement outcome.
        letters = {'I':np.eye(2), 'X':np.array([[0,1],[1,0]]),
                   'Y':np.array([[0,-1j],[1j,0]]), 'Z':np.diag([1,-1])}
        rng = np.random.default_rng(17)
        for target in ('-Y', '-XY', 'YZ', 'XXXX'):
            with self.subTest(target=target):
                d = compile_operational_measurement(target)
                c = d.compilation.circuit
                n, q = len(c.inputs), c.num_qubits
                psi = rng.normal(size=2**n)+1j*rng.normal(size=2**n)
                psi /= np.linalg.norm(psi)
                initial = np.zeros(2**q, complex)
                for j in range(2**n):
                    physical = sum(((j >> k)&1) << w for k,w in enumerate(c.inputs))
                    initial[physical] = psi[j]

                def tensor(word, wires, count):
                    mapping = dict(zip(wires, word))
                    matrix = np.ones((1,1), complex)
                    for w in reversed(range(count)):
                        matrix = np.kron(matrix, letters[mapping.get(w,'I')])
                    return matrix

                identity = np.eye(2**q)
                expanded = []
                for op in c.operations:
                    if op.name in ('M','MX','MXX','MZZ'):
                        word = {'M':'Z','MX':'X','MXX':'XX','MZZ':'ZZ'}[op.name]
                        expanded.append((True, tensor(word, op.targets, q)))
                    elif op.name == 'R':
                        expanded.append((False, identity))  # fresh zero qubit
                    elif op.name in ('RX','H'):
                        expanded.append((False, (tensor('X',op.targets,q)+tensor('Z',op.targets,q))/np.sqrt(2)))
                    elif op.name in ('S', 'S_DAG'):
                        z = tensor('Z', op.targets, q)
                        phase = 1j if op.name == 'S' else -1j
                        expanded.append((False, (identity+z)/2+phase*(identity-z)/2))
                    elif op.name in ('X','Y','Z'):
                        expanded.append((False, tensor(op.name, op.targets, q)))
                    else:
                        self.fail(f'Unhandled dense test operation: {op.name}')
                rho = [np.zeros((2**n,2**n),complex) for _ in range(2)]
                discarded = [w for w in range(q) if w not in c.outputs]
                p = d.compilation.measurement.sign*tensor(d.compilation.measurement.pauli, range(n), n)
                projected = [(psi+(-1)**r*(p@psi))/2 for r in range(2)]
                alpha_squared = 2.**d.certificate['branch_scalar_squared_power_of_two']
                for record in itertools.product((0,1), repeat=len(d.measurement_operations)):
                    state, cursor = initial.copy(), 0
                    for measured, matrix in expanded:
                        if measured:
                            state = (state+(-1)**record[cursor]*(matrix@state))/2
                            cursor += 1
                        else:
                            state = matrix@state
                    if np.linalg.norm(state) < 1e-12:
                        continue
                    self.assertFalse(any(d.consistency_checks(record)))
                    frame = d.output_pauli_frame(record)
                    frame_word = ''.join('IXYZ'[frame[i]] for i in range(n))
                    state = tensor(frame_word, c.outputs, q)@state
                    # Express the state as output x discarded and trace discarded.
                    reshaped = np.zeros((2**n,2**len(discarded)),complex)
                    for basis, amplitude in enumerate(state):
                        out = sum(((basis>>w)&1)<<j for j,w in enumerate(c.outputs))
                        anc = sum(((basis>>w)&1)<<j for j,w in enumerate(discarded))
                        reshaped[out, anc] = amplitude
                    branch_rho = reshaped@reshaped.conj().T
                    r = d.outcome_bit(record)
                    expected = np.outer(projected[r], projected[r].conj())
                    np.testing.assert_allclose(branch_rho, alpha_squared*expected, atol=1e-12)
                    rho[r] += branch_rho
                for r in range(2):
                    np.testing.assert_allclose(rho[r], np.outer(projected[r],projected[r].conj()),atol=1e-12)


if __name__ == '__main__':
    unittest.main()
