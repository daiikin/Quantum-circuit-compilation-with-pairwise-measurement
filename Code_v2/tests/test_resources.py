import unittest
from types import SimpleNamespace as NS
from pauli_compiler.analysis.scheduling import TimingModel, schedule, parity_ready, feedback_schedule
from pauli_compiler.circuits.model import CircuitOperation as Op
from pauli_compiler.decoding.affine import Affine

class SchedulingTests(unittest.TestCase):
    def test_published_six_step_gadget(self):
        # Grans-Samuelsson et al., Fig. 2 and Eq. (1), A=4 B=5 C=6.
        layers = [ [('MX',(4,))],
                   [('MZZ',(0,4)),('M',(5,)),('MX',(6,))],
                   [('MXX',(4,5)),('MZZ',(6,1))],
                   [('MZZ',(2,4)),('MXX',(5,6))],
                   [('MX',(4,)),('M',(5,)),('MZZ',(6,3))],
                   [('MX',(6,))] ]
        ops = [Op(n,q) for layer in layers for n,q in layer]
        ready, records, entries = schedule(ops,7)
        self.assertEqual(max(ready),6)
        self.assertEqual(len(records),12)
        self.assertEqual([e['record_index'] for e in entries],list(range(12)))

    def test_shared_qubit_is_exclusive_even_when_commuting(self):
        ready, _, entries = schedule([Op('MZZ',(0,1)),Op('MZZ',(1,2)),Op('MZZ',(3,4))],5)
        self.assertEqual([e['start'] for e in entries],[0,1,0])
        self.assertEqual(max(ready),2)

    def test_gate_dependency_and_durations(self):
        ops=[Op('MXX',(0,1)),Op('H',(1,)),Op('MZZ',(1,2))]
        self.assertEqual(max(schedule(ops,3)[0]),3)
        self.assertEqual(max(schedule(ops,3,TimingModel(pair_measurement=5,single_gate=2))[0]),12)
        self.assertEqual(max(schedule(ops,3,TimingModel(single_gate=0))[0]),2)

    def test_parity_feedback_and_classical_latency(self):
        expr=Affine(sum(1 << (i+1) for i in range(4)))
        timing=TimingModel(xor_time=2,feedforward_latency=3,feedback_gate=1)
        self.assertEqual(parity_ready(expr,[1,2,3,4],timing),8)
        decoder=NS(compilation=NS(circuit=NS(outputs=(0,))), x_corrections=[expr],z_corrections=[Affine()])
        ready, entries=feedback_schedule(decoder,[5],[1,2,3,4],timing)
        self.assertEqual(ready,[12]); self.assertEqual(len(entries),1)
        self.assertEqual(feedback_schedule(decoder,[5],[1,2,3,4],TimingModel())[0],[6])

    def test_invalid_model_or_operation(self):
        for value in (-1,float('nan'),float('inf'),True):
            with self.assertRaises(ValueError): TimingModel(pair_measurement=value)
        with self.assertRaises(ValueError): schedule([Op('UNKNOWN',(0,))],1)

if __name__ == '__main__': unittest.main()
