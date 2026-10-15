"""Physical operations, boundary-wire mappings, and raw Stim export."""
from dataclasses import dataclass

@dataclass
class CircuitOperation:
    name: str
    targets: tuple


@dataclass
class ExtractedCircuit:
    """Normalized physical operations plus explicit postselection metadata.

    graph_map = branch_scale * physical_postselected_map, with external
    boundaries ordered by inputs/outputs. Every measurement result must be
    zero. This is not an adaptive all-outcomes implementation.
    """
    operations: tuple
    num_qubits: int
    inputs: tuple
    outputs: tuple
    branch_scale: complex

    def to_stim(self):
        import stim
        circuit = stim.Circuit()
        for op in self.operations:
            if op.name in ('MXX', 'MZZ'):
                f = stim.target_x if op.name == 'MXX' else stim.target_z
                circuit.append('MPP', [f(op.targets[0]), stim.target_combiner(), f(op.targets[1])])
            else:
                circuit.append(op.name, op.targets)
        return circuit

    def metadata(self):
        return {'num_qubits': self.num_qubits, 'input_wires': self.inputs,
                'output_wires': self.outputs,
                'branch_scale': [self.branch_scale.real, self.branch_scale.imag],
                'postselect_measurement_records_to_zero': list(range(self.to_stim().num_measurements)),
                'semantics': 'ZX map = branch_scale * all-zero postselected circuit map'}
