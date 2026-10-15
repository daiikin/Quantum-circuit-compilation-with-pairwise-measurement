"""Conservative elimination of repeated noiseless Pauli measurements."""
from .model import ExtractedCircuit

def remove_redundant_measurements(circuit):
    """Remove repeated positive Pauli measurements without changing the instrument.

    A measurement is redundant if the same observable was already measured and
    every intervening operation preserves its eigenspaces. Commuting Pauli
    measurements and operations on disjoint wires preserve the cached result.
    For safety, every other operation touching the observable invalidates it
    (including resets and Clifford gates, even ones that might commute).

    XX/YY/ZZ target order is irrelevant; different partners are different
    observables. Noncommuting measurements invalidate cached results.
    Input operations are unchanged. Output measurement records are renumbered:
    call this BEFORE derive_decoder, never reuse an earlier decoder's formulas.
    No noise/fault-detection equivalence is claimed for repeated measurements.
    """
    axes = {'M': 'Z', 'MX': 'X', 'MY': 'Y',
            'MXX': 'X', 'MYY': 'Y', 'MZZ': 'Z'}
    known = set()
    operations = []

    def commute(left, right):
        right = dict(right)
        return sum(q in right and axis != right[q] for q, axis in left) % 2 == 0

    for op in circuit.operations:
        if op.name in axes:
            arity = 2 if op.name in ('MXX', 'MYY', 'MZZ') else 1
            if len(op.targets) != arity or len(set(op.targets)) != arity:
                raise ValueError(f'Invalid targets for {op.name}: {op.targets}')
            observable = tuple(sorted((q, axes[op.name]) for q in op.targets))
            if observable in known:
                continue
            known = {previous for previous in known if commute(previous, observable)}
            known.add(observable)
        else:
            touched = set(op.targets)
            known = {previous for previous in known
                     if not any(q in touched for q, _ in previous)}
        operations.append(op)
    return ExtractedCircuit(tuple(operations), circuit.num_qubits,
                            circuit.inputs, circuit.outputs, circuit.branch_scale)
