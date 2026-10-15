"""Derive global outcome parities and output corrections with Paulimer/Binar.
Uses all-branch Choi simulation; no in-house stabilizer simulator."""
from importlib.metadata import version
import paulimer as pm
from binar import BitMatrix, BitVector, EchelonForm
from .affine import Affine
from .result import MeasurementDecoder

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
