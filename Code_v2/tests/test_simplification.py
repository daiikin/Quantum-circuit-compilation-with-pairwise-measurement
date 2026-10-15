"""Run with python -m unittest -v tests.test_simplification."""
import unittest
from pauli_compiler.circuits.model import CircuitOperation, ExtractedCircuit
from pauli_compiler.circuits.simplification import remove_redundant_measurements
from pauli_compiler.pipeline import compile_measurement
from pauli_compiler.decoding.paulimer_backend import derive_decoder
from tests.validate_decoder import check_shot


def circuit(*ops):
    return ExtractedCircuit(tuple(CircuitOperation(name, tuple(targets))
                                  for name, targets in ops), 4, (0,1), (0,1), 2)


class SimplificationTests(unittest.TestCase):
    def test_same_pair_reversed_and_disjoint_gate(self):
        original = circuit(('MZZ',(0,1)), ('H',(3,)), ('MZZ',(1,0)))
        result = remove_redundant_measurements(original)
        self.assertEqual(len(result.operations),2)
        self.assertEqual(len(original.operations),3)
        self.assertEqual(result.branch_scale,original.branch_scale)

    def test_different_partner_kept(self):
        original = circuit(('MZZ',(0,1)), ('MZZ',(0,2)))
        self.assertEqual(remove_redundant_measurements(original),original)

    def test_commuting_measurement_keeps_result(self):
        original = circuit(('MZZ',(0,1)), ('MXX',(0,1)), ('MZZ',(0,1)))
        self.assertEqual(len(remove_redundant_measurements(original).operations),2)

    def test_noncommuting_measurement_invalidates(self):
        original = circuit(('MZZ',(0,1)), ('MX',(0,)), ('MZZ',(0,1)))
        self.assertEqual(remove_redundant_measurements(original),original)

    def test_gate_or_reset_on_either_wire_invalidates(self):
        for name in ('H','X','R','RX'):
            for wire in (0,1):
                original = circuit(('MZZ',(0,1)), (name,(wire,)), ('MZZ',(0,1)))
                self.assertEqual(remove_redundant_measurements(original),original)

    def test_idempotent(self):
        original = circuit(('M',(0,)),('M',(0,)),('M',(0,)))
        result = remove_redundant_measurements(original)
        self.assertEqual(len(result.operations),1)
        self.assertEqual(remove_redundant_measurements(result),result)

    def test_compilation_and_feedback(self):
        for word in ('XXXX','-XYZXYZ','YZZYXZXXY'):
            raw = compile_measurement(word,simplify_measurements=False)
            optimized = compile_measurement(word)
            self.assertLess(len(optimized.circuit.operations),len(raw.circuit.operations))
            decoder = derive_decoder(optimized)
            for seed in range(8):
                check_shot(decoder,seed,bell=seed%2==0,feedback=seed%3==0)


if __name__ == '__main__':
    unittest.main()
