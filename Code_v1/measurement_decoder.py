"""Operational Pauli-measurement decoding using Microsoft QDK-EC.

Paulimer OutcomeCompleteSimulation supplies all-branch state/outcome algebra;
Binar supplies GF(2) linear solves. No local stabilizer simulator is used.
Stim is used only for exporting circuits/frames (and independent tests).

Convention: b_i=0 denotes +1. Affine.mask bit 0 is the constant; bit i+1
is b_i. Formulas are valid on noiseless records satisfying checks.
Supports Clifford gates, fresh |0>/|+> preparations, and one/two-body Pauli
measurements from the existing ZX extractor. Resets of used wires are rejected.
"""
from dataclasses import dataclass
from importlib.metadata import version

import paulimer as pm
from binar import BitMatrix, BitVector, EchelonForm


@dataclass(frozen=True)
class Affine:
    mask: int = 0

    def __xor__(self, other):
        return Affine(self.mask ^ other.mask)

    @property
    def constant(self):
        return self.mask & 1

    @property
    def indices(self):
        return tuple(i for i in range(self.mask.bit_length()-1)
                     if (self.mask >> (i+1)) & 1)

    def evaluate(self, record):
        value = self.constant
        for i in self.indices:
            value ^= int(record[i])
        return value

    def __str__(self):
        return ' xor '.join((['1'] if self.constant else []) +
                            [f'b{i}' for i in self.indices]) or '0'

    def as_dict(self):
        return {'constant': self.constant, 'record_indices': self.indices,
                'formula': str(self)}


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


_GATES = {name: getattr(pm.UnitaryOpcode, opcode) for name, opcode in {
    'I': 'I', 'X': 'X', 'Y': 'Y', 'Z': 'Z', 'H': 'Hadamard',
    'S': 'SqrtZ', 'S_DAG': 'SqrtZInv',
    'SQRT_X': 'SqrtX', 'SQRT_X_DAG': 'SqrtXInv',
    'SQRT_Y': 'SqrtY', 'SQRT_Y_DAG': 'SqrtYInv',
    'CX': 'ControlledX', 'CZ': 'ControlledZ', 'SWAP': 'Swap',
}.items()}
_MEASUREMENTS = {'M': 'Z', 'MX': 'X', 'MY': 'Y',
                 'MXX': 'XX', 'MYY': 'YY', 'MZZ': 'ZZ'}


def _pauli(word, wires, sign=1):
    result = pm.SparsePauli(' '.join(f'{letter}{wire}' for letter, wire in zip(word, wires)
                                   if letter != 'I'))
    return -result if sign == -1 else result


def _simulate(compilation):
    """Build the Choi experiment; preparations add no artificial record bits."""
    c = compilation.circuit
    n, q = len(c.inputs), c.num_qubits
    if len(c.outputs) != n or n != len(compilation.measurement.pauli):
        raise ValueError('Expected one input and output for each target Pauli factor')
    for wires in (c.inputs, c.outputs):
        if len(set(wires)) != n or any(w < 0 or w >= q for w in wires):
            raise ValueError('Invalid or duplicate boundary wire')
    sim = pm.OutcomeCompleteSimulation(q+n)
    for j, wire in enumerate(c.inputs):
        sim.apply_unitary(pm.UnitaryOpcode.Hadamard, [q+j])
        sim.apply_unitary(pm.UnitaryOpcode.ControlledX, [q+j, wire])
    used = set(c.inputs)
    records = []
    for op_index, op in enumerate(c.operations):
        targets = tuple(op.targets)
        if not targets or len(set(targets)) != len(targets) or any(w < 0 or w >= q for w in targets):
            raise ValueError(f'Invalid targets at operation {op_index}')
        if op.name in ('R', 'RX'):
            if len(targets) != 1 or targets[0] in used:
                raise ValueError('Only fresh ancilla preparations are supported, not resets of used wires')
            if op.name == 'RX':
                sim.apply_unitary(pm.UnitaryOpcode.Hadamard, list(targets))
        else:
            if not set(targets) <= used:
                raise ValueError('Every ancilla must be explicitly prepared before use')
            if op.name in _MEASUREMENTS:
                word = _MEASUREMENTS[op.name]
                if len(targets) != len(word):
                    raise ValueError(f'Wrong measurement arity: {op.name}')
                index = sim.measure(_pauli(word, targets))
                records.append({'index': index, 'name': op.name,
                                'physical_wires': targets, 'operation_index': op_index})
            elif op.name in _GATES:
                arity = 2 if op.name in ('CX', 'CZ', 'SWAP') else 1
                if len(targets) != arity:
                    raise ValueError(f'Wrong gate arity: {op.name}')
                sim.apply_unitary(_GATES[op.name], list(targets))
            else:
                raise ValueError(f'Unsupported operation {op.name!r}')
        used.update(targets)
    if not set(c.outputs) <= used:
        raise ValueError('Output wire was never prepared or supplied as input')
    return sim, tuple(records)


class _Expectations:
    """Read affine eigenvalue signs from Paulimer's R|A u> representation.

    R^dagger Q R = (-1)^c Z^v implies sign(Q) = c + v^T A u.
    An X/Y in the preimage means Q is not a definite branch stabilizer.
    Masks here use free random bits u, not physical record bits b.
    """
    def __init__(self, sim):
        self.clifford = sim.clifford
        self.rows = [sum(1 << (j+1) for j in row)
                     for row in sim.sign_matrix.sparse_rows()]

    def sign(self, observable):
        preimage = self.clifford.preimage_of(observable)
        if preimage.phase not in (1, -1):
            raise ValueError('Expected a Hermitian observable')
        result = int(preimage.phase == -1)
        for i, letter in enumerate(preimage.characters):
            if letter in 'XY':
                raise ValueError('Circuit does not preserve the required Choi stabilizer; '
                                 'it is not the target measurement up to output Paulis')
            if letter == 'Z':
                result ^= self.rows[i]
        return result


class _Records:
    """Convert free-variable parities to physical record parities using Binar."""
    def __init__(self, sim):
        matrix = sim.outcome_matrix
        self.count = sim.outcome_count
        self.k = sim.random_outcome_count
        shift = sim.outcome_shift
        echelon = EchelonForm(matrix)
        self.free = []
        for j in range(self.k):
            coefficients = echelon.transpose_solve(BitVector([i == j for i in range(self.k)]))
            if coefficients is None:
                raise ValueError('Paulimer outcome matrix has deficient column rank')
            mask = 0
            for i, enabled in enumerate(coefficients):
                if enabled:
                    mask ^= (1 << (i+1)) ^ int(shift[i])
            self.free.append(mask)
        checks = []
        for i, row in enumerate(matrix.sparse_rows()):
            free_mask = int(shift[i]) ^ sum(1 << (j+1) for j in row)
            relation = self.to_record(free_mask).mask ^ (1 << (i+1))
            if relation:
                checks.append(Affine(relation))
        self.checks = tuple(checks)

    def to_record(self, free_mask):
        result = free_mask & 1
        for j in range(self.k):
            if free_mask & (1 << (j+1)):
                result ^= self.free[j]
        return Affine(result)


def derive_decoder(compilation):
    """Derive and verify a complete noiseless Pauli measurement with Paulimer.

    Requires an extracted Clifford circuit plus its intended signed Pauli.
    Fails if the circuit measures extra information or needs corrections that
    cannot be placed solely on the designated output qubits. No enumeration
    of records, and no fallback to the former in-house symbolic simulator.
    """
    c, measurement = compilation.circuit, compilation.measurement
    n, q = len(c.inputs), c.num_qubits
    if measurement.sign not in (-1, 1) or any(p not in 'IXYZ' for p in measurement.pauli):
        raise ValueError('Invalid target Hermitian Pauli')
    if not any(p != 'I' for p in measurement.pauli):
        raise ValueError('Target Pauli must have positive weight')
    sim, operations = _simulate(compilation)
    records = _Records(sim)
    expected_actual = _Expectations(sim)
    # The reference carries P^T, including the minus sign from each Y.
    transpose_sign = measurement.sign * (-1)**measurement.pauli.count('Y')
    reference_p = _pauli(measurement.pauli, range(q, q+n), transpose_sign)
    outcome_free = expected_actual.sign(reference_p)
    if outcome_free >> 1 == 0:
        raise ValueError('Target outcome is not a balanced nonconstant Choi parity')
    outcome = records.to_record(outcome_free)

    # Generate the complete ideal Choi stabilizer basis with Paulimer too.
    ideal = pm.OutcomeCompleteSimulation(2*n)
    for j in range(n):
        ideal.apply_unitary(pm.UnitaryOpcode.Hadamard, [n+j])
        ideal.apply_unitary(pm.UnitaryOpcode.ControlledX, [n+j, j])
    ideal.measure(_pauli(measurement.pauli, range(n), measurement.sign))
    if ideal.random_outcome_count != 1 or ideal.outcome_shift[0] or not ideal.outcome_matrix[0, 0]:
        raise ValueError('Unexpected ideal measurement representation')
    mapping = tuple(c.outputs) + tuple(range(q, q+n))
    coefficients, rhs, targets = [], [], []
    ideal_clifford, ideal_signs = ideal.clifford, ideal.sign_matrix
    for i in range(2*n):
        generator = ideal_clifford.image_z(i)
        target = _pauli(generator.characters, mapping, -1 if generator.phase == -1 else 1)
        target_sign = outcome_free if ideal_signs[i, 0] else 0
        difference = expected_actual.sign(target) ^ target_sign
        # Output C=X^x Z^z flips this generator by z_Q.x + x_Q.z.
        word = generator.characters[:n]
        coefficients.append([p in 'ZY' for p in word] + [p in 'XY' for p in word])
        rhs.append(difference)
        targets.append((target, target_sign))
    system = EchelonForm(BitMatrix(coefficients))
    corrections_free = [0]*(2*n)
    # One decomposition reused for the constant and every free random bit.
    for bit in range(records.k+1):
        solution = system.solve(BitVector([bool(mask & (1 << bit)) for mask in rhs]))
        if solution is None:
            raise ValueError('No output-only Pauli correction realizes the target measurement')
        for j, enabled in enumerate(solution):
            if enabled:
                corrections_free[j] ^= 1 << bit
    corrections = tuple(records.to_record(mask) for mask in corrections_free)
    xs, zs = corrections[:n], corrections[n:]

    # Apply the derived feedback through Paulimer, then recheck ALL 2n ideal
    # generators. Their purity also certifies decoupling from discarded ancillas.
    for wire, x, z in zip(c.outputs, xs, zs):
        for letter, expr in [('X', x), ('Z', z)]:
            pauli = _pauli(letter, [wire])
            if expr.constant:
                sim.apply_pauli(pauli)
            if expr.indices:
                sim.apply_conditional_pauli(pauli, list(expr.indices), parity=True)
    verified = _Expectations(sim)
    if any(verified.sign(target) != sign for target, sign in targets):
        raise ValueError('Paulimer verification of derived feedback failed')
    certificate = {
        'backend': 'paulimer.OutcomeCompleteSimulation',
        'paulimer_version': version('paulimer'), 'binar_version': version('binar'),
        'signed_pauli': ('-' if measurement.sign < 0 else '+') + measurement.pauli,
        'verified_choi_generators': 2*n,
        'all_allowed_records_verified': True,
        'discarded_ancillas_decoupled': True,
        'random_bits': records.k,
        'choi_record_probability_power_of_two': -records.k,
        'branch_scalar_squared_power_of_two': 1-records.k,
        'records_per_decoded_outcome_power_of_two': records.k-1,
        'instrument_normalization_verified': True,
        'scope': 'Noiseless measurement instrument; not a fault-tolerant noise decoder',
    }
    return MeasurementDecoder(compilation, outcome, xs, zs, records.checks,
                              tuple(i for i, bit in enumerate(sim.random_outcome_indicator) if bit),
                              operations, certificate)


def compile_operational_measurement(pauli, outcome=1, *, check_each=False):
    from zx_r4_r5 import compile_measurement
    return derive_decoder(compile_measurement(pauli, outcome, check_each=check_each))
