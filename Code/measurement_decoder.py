"""Derive affine outcome/Pauli-frame formulas for extracted Clifford circuits.

No outcome enumeration: symbolic stabilizer signs are XOR bitmasks.
Bit 0 is the constant; bit i+1 is measurement result b_i (0 means +1).
The derivation certifies the noiseless corrected measurement instrument on
all allowed records using its Choi stabilizers and branch probabilities.
Consistency checks flag forbidden noiseless records; they are not a noisy
QEC decoder. Requires stim and the existing compiler modules.
"""
from dataclasses import dataclass
import stim


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
        return tuple(i for i in range(self.mask.bit_length()-1) if (self.mask >> (i+1)) & 1)

    def evaluate(self, record):
        value = self.constant
        for i in self.indices:
            value ^= int(record[i])
        return value

    def __str__(self):
        terms = (['1'] if self.constant else []) + [f'b{i}' for i in self.indices]
        return ' xor '.join(terms) or '0'

    def as_dict(self):
        return {'constant': self.constant, 'record_indices': self.indices, 'formula': str(self)}


def _positive(pauli, expr=0):
    p = pauli.copy()
    if p.sign not in (1, -1):
        raise ValueError('Expected Hermitian Pauli')
    expr ^= int(p.sign == -1)
    p.sign = 1
    return p, expr


def _multiply(a, b):
    return _positive(a[0]*b[0], a[1]^b[1])


def _vector(p):
    n = len(p)
    v = 0
    for i in range(n):
        if p[i] in (1, 2):
            v |= 1 << i
        if p[i] in (2, 3):
            v |= 1 << (n+i)
    return v


def _basis(rows):
    basis = {}
    for row in rows:
        v = _vector(row[0])
        while v:
            k = v.bit_length()-1
            if k not in basis:
                basis[k] = (v, row)
                break
            w, b = basis[k]
            v ^= w
            row = _multiply(row, b)
    return basis


def _expectation(basis, observable):
    row = _positive(observable)
    v = _vector(row[0])
    while v:
        k = v.bit_length()-1
        if k not in basis:
            raise ValueError('Requested Pauli is not in the stabilizer group')
        w,b = basis[k]
        v ^= w
        row = _multiply(row, b)
    return row[1]


class _SymbolicState:
    def __init__(self, n):
        self.n = n
        self.rows = []
        for q in range(n):
            p = stim.PauliString(n)
            p[q] = 'Z'
            self.rows.append((p, 0))

    def gate(self, name, targets):
        c = stim.Circuit()
        c.append(name, targets)
        self.rows = [_positive(p.after(c), expr) for p,expr in self.rows]

    def measure(self, observable, expression):
        anti = [i for i,(p,_) in enumerate(self.rows) if not p.commutes(observable)]
        if not anti:
            return False, expression ^ _expectation(_basis(self.rows), observable)
        pivot = self.rows[anti[0]]
        for i in anti[1:]:
            self.rows[i] = _multiply(self.rows[i], pivot)
        self.rows[anti[0]] = _positive(observable, expression)
        return True, 0


def _solve(equations, variables, reduce_expression):
    """GF(2) solve with affine RHS, setting free correction variables to zero."""
    pivots = {}
    for coefficients, rhs in equations:
        rhs = reduce_expression(rhs)
        while coefficients:
            k = coefficients.bit_length()-1
            if k not in pivots:
                pivots[k] = (coefficients, rhs)
                break
            v,e = pivots[k]
            coefficients ^= v
            rhs ^= e
        if not coefficients and reduce_expression(rhs):
            raise ValueError('No output-only Pauli correction solves the instrument')
    solution = [0]*variables
    for k,(v,rhs) in sorted(pivots.items()):
        for j in range(k):
            if (v >> j) & 1:
                rhs ^= solution[j]
        solution[k] = reduce_expression(rhs)
    return solution


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
        if any(b not in (0,1) for b in record):
            raise ValueError('Record entries must be bits')
        return record

    def consistency_checks(self, record):
        record = self._record(record)
        return tuple(e.evaluate(record) for e in self.checks)

    def outcome_bit(self, record, *, validate=True):
        record = self._record(record)
        if validate and any(self.consistency_checks(record)):
            raise ValueError('Record violates noiseless consistency checks; noise decoding is not implemented')
        return self.outcome.evaluate(record)

    def outcome_eigenvalue(self, record, *, validate=True):
        return 1-2*self.outcome_bit(record, validate=validate)

    def output_pauli_frame(self, record, *, validate=True):
        record = self._record(record)
        if validate and any(self.consistency_checks(record)):
            raise ValueError('Record violates noiseless consistency checks')
        p = stim.PauliString(len(self.x_corrections))
        for j,(x,z) in enumerate(zip(self.x_corrections,self.z_corrections)):
            xb,zb = x.evaluate(record),z.evaluate(record)
            p[j] = ('I','X','Z','Y')[xb+2*zb]
        return p  # Data OUTPUT order, not physical wire order; global phase omitted.

    def corrected_stim(self):
        """Physical circuit with real measurement-controlled feedback; no postselection.

        Do not also apply output_pauli_frame after running this circuit.
        Records still need outcome_bit() for the final signed-Pauli result.
        """
        c = self.compilation.circuit.to_stim()
        count = c.num_measurements
        for physical,x,z in zip(self.compilation.circuit.outputs,self.x_corrections,self.z_corrections):
            for gate,expr in (('X',x),('Z',z)):
                if expr.constant:
                    c.append(gate,[physical])
                for i in expr.indices:
                    c.append('C'+gate,[stim.target_rec(i-count),physical])
        return c

    def formulas(self):
        c = self.compilation.circuit
        lines = [f'r = {self.outcome}', 'Measured eigenvalue = (-1)^r']
        for j,(wire,x,z) in enumerate(zip(c.outputs,self.x_corrections,self.z_corrections)):
            lines.append(f'Output data {j} (physical p{wire}): X^({x}) Z^({z})')
        lines += [f'Check: {expr} = 0' for expr in self.checks]
        return '\n'.join(lines)

    def metadata(self):
        return {'outcome':self.outcome.as_dict(),
                'corrections':[{'data_output':j,'physical_wire':w,'x':x.as_dict(),'z':z.as_dict()}
                               for j,(w,x,z) in enumerate(zip(self.compilation.circuit.outputs,
                                                            self.x_corrections,self.z_corrections))],
                'checks':[e.as_dict() for e in self.checks],
                'measurements':self.measurement_operations,
                'input_wires':self.compilation.circuit.inputs,
                'output_wires':self.compilation.circuit.outputs,
                'certificate':self.certificate,
                'semantics':'Unpostselected corrected measurement of the signed Pauli; noise-free'}


def derive_decoder(compilation):
    """Derive and algebraically certify XOR output/correction formulas.

    Raises if output-only Pauli corrections are insufficient. Works with
    either compiled outcome branch: its stored outcome affects the compiled
    gates, while the decoder returns the ACTUAL eigenvalue of the signed P.
    """
    m,c = compilation.measurement,compilation.circuit
    n,q = len(m.pauli),c.num_qubits
    state = _SymbolicState(q+n)
    for j,wire in enumerate(c.inputs):
        state.gate('H',[q+j])
        state.gate('CX',[q+j,wire])
    measurements, constraints, random_indices = [], [], []
    substitutions = {}

    def reduce_expression(expr):
        for i,replacement in sorted(substitutions.items(),reverse=True):
            if expr & (1 << (i+1)):
                expr ^= (1 << (i+1)) ^ replacement
        return expr

    touched = set(c.inputs)
    for op in c.operations:
        if op.name in ('R','RX'):
            if any(w in touched for w in op.targets):
                raise ValueError('Only fresh-wire preparations are supported')
            if op.name == 'RX':
                state.gate('H',op.targets)
        elif op.name in ('M','MX','MXX','MZZ'):
            i = len(measurements)
            obs = stim.PauliString(q+n)
            for w in op.targets:
                obs[w] = 'X' if op.name in ('MX','MXX') else 'Z'
            random, constraint = state.measure(obs,1 << (i+1))
            measurements.append({'index':i,'name':op.name,'physical_wires':op.targets})
            if random:
                random_indices.append(i)
            else:
                constraint = reduce_expression(constraint)
                constraints.append(Affine(constraint))
                substitutions[i] = constraint ^ (1 << (i+1))
        else:
            state.gate(op.name,op.targets)
        touched.update(op.targets)

    basis = _basis(state.rows)
    # Reference-only P^T encodes the measured INPUT observable, unaffected
    # by any output correction. Include the transpose sign of every Y.
    ref = stim.PauliString(q+n)
    ref.sign = m.sign * (-1 if m.pauli.count('Y') % 2 else 1)
    for j,letter in enumerate(m.pauli):
        ref[q+j] = letter
    outcome = reduce_expression(_expectation(basis,ref))
    if not (outcome >> 1):
        raise ValueError('The extracted circuit does not reveal both Pauli outcomes')

    ideal = _SymbolicState(2*n)
    for j in range(n):
        ideal.gate('H',[n+j])
        ideal.gate('CX',[n+j,j])
    target = stim.PauliString(m.pauli+'I'*n)
    target.sign = m.sign
    random,_ = ideal.measure(target,outcome)
    if not random:
        raise AssertionError('Nontrivial Pauli must be random on a Bell input')
    mapping = tuple(c.outputs) + tuple(q+j for j in range(n))
    equations = []
    for p,expr in ideal.rows:
        lifted = stim.PauliString(q+n)
        for j,wire in enumerate(mapping):
            lifted[wire] = p[j]
        actual_expr = _expectation(basis,lifted)
        coefficients = 0
        for j in range(n):
            if p[j] in (2,3):  # Z/Y anticommutes with an X correction
                coefficients |= 1 << j
            if p[j] in (1,2):  # X/Y anticommutes with a Z correction
                coefficients |= 1 << (n+j)
        equations.append((coefficients,actual_expr ^ expr))
    solution = _solve(equations,2*n,reduce_expression)
    for coefficients,rhs in equations:
        for j,value in enumerate(solution):
            if (coefficients >> j) & 1:
                rhs ^= value
        if reduce_expression(rhs):
            raise AssertionError('Correction identity failed symbolic verification')
    # k free outcome bits => each Choi branch probability 2^-k. Nonconstant
    # affine r divides all valid records equally into the two Pauli outcomes.
    # Corrected branch K is proportional to projector with |alpha|^2=2^(1-k).
    # Each outcome has 2^(k-1) records, so the grouped instrument is normalized.
    k = len(random_indices)
    certificate = {'symbolic_choi_generators':2*n,'free_record_bits':k,
                   'valid_records_log2':k,'records_per_outcome_log2':k-1,
                   'branch_projector_scale_squared_log2':1-k,
                   'all_allowed_records_certified':True,
                   'grouped_instrument_normalized':True,
                   'signed_pauli':('-' if m.sign<0 else '+')+m.pauli}
    return MeasurementDecoder(compilation,Affine(outcome),
                              tuple(Affine(x) for x in solution[:n]),
                              tuple(Affine(z) for z in solution[n:]),
                              tuple(constraints),tuple(random_indices),
                              tuple(measurements),certificate)


def compile_operational_measurement(pauli, *, check_each=False):
    """Compile and derive a measurement accepting both outcomes (no postselection)."""
    from zx_r4_r5 import compile_measurement
    return derive_decoder(compile_measurement(pauli,check_each=check_each))
