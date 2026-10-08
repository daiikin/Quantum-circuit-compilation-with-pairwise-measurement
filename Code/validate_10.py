"""Ten reproducible end-to-end checks using stabilizer Choi states.

Checks the full conditional linear map up to global phase, not merely a
sample of computational-basis inputs. Also checks its normalization using
postselection probability. No exhaustive enumeration of Pauli strings.
"""
import json
import math
from pathlib import Path
import time
import stim
from pauli_measurements import generate_pauli_measurement
from zx_r4_r5 import compile_measurement


def check_compilation(compilation):
    m, circuit = compilation.measurement, compilation.circuit
    n, q = len(m.pauli), circuit.num_qubits
    actual = stim.TableauSimulator()
    actual.set_num_qubits(q+n)
    for i, wire in enumerate(circuit.inputs):
        actual.h(q+i)
        actual.cnot(q+i, wire)
    random_results = 0
    for op in circuit.operations:
        if op.name in ('M', 'MX', 'MXX', 'MZZ'):
            observable = stim.PauliString(q+n)
            basis = 'X' if op.name in ('MX','MXX') else 'Z'
            for wire in op.targets:
                observable[wire] = basis
            expectation = actual.peek_observable_expectation(observable)
            if expectation == -1:
                raise AssertionError('Extracted postselection is impossible')
            random_results += int(expectation == 0)
            actual.postselect_observable(observable)
        else:
            one = stim.Circuit()
            one.append(op.name, op.targets)
            actual.do_circuit(one)
    # Independently specified ideal Pauli projector on a maximally
    # entangled state. Its normalized Choi state identifies the linear map.
    expected = stim.TableauSimulator()
    expected.set_num_qubits(2*n)
    for i in range(n):
        expected.h(n+i)
        expected.cnot(n+i, i)
    observable = stim.PauliString(m.pauli + 'I'*n)
    observable.sign = m.sign
    expected.postselect_observable(observable, desired_value=(m.outcome == -1))
    # Extend each output/reference stabilizer by identity on physical ancillas.
    # Agreement on 2n independent stabilizers proves equality of the pure
    # retained state and absence of residual entanglement with ancillas.
    mapping = tuple(circuit.outputs) + tuple(q+i for i in range(n))
    for generator in expected.canonical_stabilizers():
        lifted = stim.PauliString(q+n)
        lifted.sign = generator.sign
        for i, dest in enumerate(mapping):
            lifted[dest] = generator[i]
        if actual.peek_observable_expectation(lifted) != 1:
            raise AssertionError('Extracted Choi state differs from ideal measurement')
    # Ideal nontrivial Pauli projector has Choi success probability 1/2.
    # Compare in log domain to avoid underflow for large diagrams.
    log2_scaled_probability = 2*math.log2(abs(circuit.branch_scale))-random_results
    if not math.isclose(log2_scaled_probability, -1, abs_tol=1e-9):
        raise AssertionError(f'Scalar mismatch: log2(probability)={log2_scaled_probability}')
    return {'random_postselections': random_results,
            'choi_generators_checked': 2*n,
            'scaled_choi_probability': .5}


def main():
    folder = Path('generated_10')
    folder.mkdir(exist_ok=True)
    rows = []
    # Covers small/base cases, both parity classes, expansion and splitting,
    # identity wires, signs, and both outcomes. Exactly ten circuits.
    for i, weight in enumerate((1,2,3,4,5,6,7,9,13,20)):
        start = time.monotonic()
        m = generate_pauli_measurement(weight, seed=20261008+i,
                                      total_qubits=weight + (i % 3),
                                      outcome=(-1 if i % 2 else 1),
                                      sign=(-1 if i % 3 == 0 else 1))
        c = compile_measurement(m, check_each=(weight <= 7))
        checks = check_compilation(c)
        # For small diagrams, independently contract the reduced ZX tensor too.
        if len(m.pauli) <= 4:
            import pyzx as zx
            import numpy as np
            from pauli_measurements import expected_projector
            g = c.reduced.graph
            mat = zx.tensor.tensor_to_matrix(zx.tensor.tensorfy(g), len(g.inputs()), len(g.outputs()))
            np.testing.assert_allclose(mat, expected_projector(('-' if m.sign<0 else '')+m.pauli, m.outcome), atol=1e-12)
        stem = f'{i+1:02d}_weight_{weight}'
        (folder/(stem+'.stim')).write_text(str(c.circuit.to_stim()))
        metadata = c.circuit.metadata()
        metadata.update({'pauli': m.pauli, 'sign': m.sign, 'outcome': m.outcome,
                         'rewrite_trace': c.reduced.trace})
        (folder/(stem+'.json')).write_text(json.dumps(metadata, indent=2))
        row = {'weight': weight, 'pauli': ('-' if m.sign<0 else '+')+m.pauli,
               'outcome': m.outcome, 'data_qubits': len(m.pauli),
               'physical_qubits': c.circuit.num_qubits,
               'operations': len(c.circuit.operations),
               'rewrites': len(c.reduced.trace), 'passed': True,
               'seconds': round(time.monotonic()-start, 3), **checks}
        rows.append(row)
        print(json.dumps(row), flush=True)
    (folder/'validation.json').write_text(json.dumps(rows, indent=2))
    print('All 10 circuits passed.', flush=True)


if __name__ == '__main__':
    main()
