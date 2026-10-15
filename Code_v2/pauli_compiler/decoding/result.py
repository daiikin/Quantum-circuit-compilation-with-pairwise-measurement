"""Decoder result, record evaluation, metadata, and physical feedback export."""
from dataclasses import dataclass
from .affine import Affine

@dataclass
class MeasurementDecoder:
    compilation: object
    outcome: Affine
    x_corrections: tuple
    z_corrections: tuple
    checks: tuple
    random_record_indices: tuple
    measurement_operations: tuple
    certificate: dict

    def _record(self, record):
        record = tuple(record)
        if len(record) != len(self.measurement_operations):
            raise ValueError(f'Expected exactly {len(self.measurement_operations)} measurement bits')
        if any(b not in (0, 1) for b in record):
            raise ValueError('Record entries must be bits')
        return record

    def consistency_checks(self, record):
        return tuple(e.evaluate(self._record(record)) for e in self.checks)

    def outcome_bit(self, record, *, validate=True):
        record = self._record(record)
        if validate and any(self.consistency_checks(record)):
            raise ValueError('Record violates noiseless consistency checks; noise decoding is not implemented')
        return self.outcome.evaluate(record)

    def outcome_eigenvalue(self, record, *, validate=True):
        return (-1)**self.outcome_bit(record, validate=validate)

    def output_pauli_frame(self, record, *, validate=True):
        import stim
        record = self._record(record)
        if validate:
            self.outcome_bit(record)
        letters = ['I', 'X', 'Z', 'Y']
        return stim.PauliString(''.join(
            letters[x.evaluate(record) + 2*z.evaluate(record)]
            for x, z in zip(self.x_corrections, self.z_corrections)))

    def corrected_stim(self):
        """Raw circuit followed by output feedback; raw record indices unchanged."""
        import stim
        circuit = self.compilation.circuit.to_stim()
        count = len(self.measurement_operations)
        for wire, x, z in zip(self.compilation.circuit.outputs,
                              self.x_corrections, self.z_corrections):
            for gate, controlled_gate, expr in [('X', 'CX', x), ('Z', 'CZ', z)]:
                if expr.constant:
                    circuit.append(gate, [wire])
                for index in expr.indices:
                    circuit.append(controlled_gate, [stim.target_rec(index-count), wire])
        return circuit

    def formulas(self):
        lines = [f'Measured Pauli: {self.certificate["signed_pauli"]}',
                 f'Outcome bit: r = {self.outcome}', 'Measured eigenvalue: (-1)^r']
        for j, (wire, x, z) in enumerate(zip(self.compilation.circuit.outputs,
                                             self.x_corrections, self.z_corrections)):
            lines.append(f'D{j} (physical p{wire}): X^({x}), Z^({z})')
        lines += [f'Noiseless check: {check} = 0' for check in self.checks]
        return '\n'.join(lines)

    def metadata(self):
        c = self.compilation.circuit
        return {'backend': 'paulimer', 'outcome': self.outcome.as_dict(),
                'x_corrections': [x.as_dict() for x in self.x_corrections],
                'z_corrections': [z.as_dict() for z in self.z_corrections],
                'checks': [e.as_dict() for e in self.checks],
                'random_record_indices': self.random_record_indices,
                'measurement_operations': self.measurement_operations,
                'input_wires': c.inputs, 'output_wires': c.outputs,
                'certificate': self.certificate}
