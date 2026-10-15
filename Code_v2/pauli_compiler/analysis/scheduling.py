"""Fixed-wire-order ASAP scheduling; see docs/DEPTH.md for the cost model."""
from dataclasses import asdict, dataclass
import math

MEASUREMENTS = frozenset(('M', 'MX', 'MY', 'MXX', 'MYY', 'MZZ'))
PREPARATIONS = frozenset(('R', 'RX', 'RY'))
GATES = frozenset(('H', 'S', 'S_DAG', 'X', 'Y', 'Z', 'CX', 'CZ', 'SWAP'))

@dataclass(frozen=True)
class TimingModel:
    """Durations in caller-chosen units; measurement includes result acquisition.

    xor_time models each level of a balanced binary parity tree, started when
    all its inputs arrive. Constants are folded into the control predicate.
    No resource conflicts beyond shared qubits; no classical processor limit.
    """
    preparation: float = 1
    single_measurement: float = 1
    pair_measurement: float = 1
    single_gate: float = 1
    two_qubit_gate: float = 1
    feedback_gate: float = 1
    xor_time: float = 0
    feedforward_latency: float = 0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f'{name} must be finite and nonnegative')

    def duration(self, operation):
        name = operation.name
        if name in PREPARATIONS:
            return self.preparation
        if name in MEASUREMENTS:
            return self.single_measurement if len(operation.targets) == 1 else self.pair_measurement
        if name in GATES:
            return self.single_gate if len(operation.targets) == 1 else self.two_qubit_gate
        raise ValueError(f'Unsupported timed operation: {name}')


def schedule(operations, num_qubits, timing=None):
    """Earliest schedule preserving order on every wire, including commuting ops.

    Each extracted measurement operation produces exactly one record bit.
    Returned starts/ends are time boundaries (the first layer is [0, 1]).
    """
    timing = timing or TimingModel()
    ready = [0] * num_qubits
    records, entries = [], []
    for index, op in enumerate(operations):
        if not op.targets or len(set(op.targets)) != len(op.targets) or any(q < 0 or q >= num_qubits for q in op.targets):
            raise ValueError('Invalid operation targets')
        start = max(ready[q] for q in op.targets)
        end = start + timing.duration(op)
        record = len(records) if op.name in MEASUREMENTS else None
        if record is not None:
            records.append(end)
        entries.append(dict(operation_index=index, name=op.name, targets=list(op.targets),
                            start=start, end=end, record_index=record))
        for q in op.targets:
            ready[q] = end
    return ready, records, entries


def parity_ready(expression, records, timing):
    indices = expression.indices
    if any(i >= len(records) for i in indices):
        raise ValueError('Correction references a missing measurement record')
    levels = (len(indices)-1).bit_length() if indices else 0
    return max((records[i] for i in indices), default=0) + levels * timing.xor_time


def feedback_schedule(decoder, ready, records, timing):
    """One parity-controlled X then Z per output, with no global final barrier."""
    ready = list(ready)
    entries = []
    for q, x, z in zip(decoder.compilation.circuit.outputs, decoder.x_corrections, decoder.z_corrections):
        for name, expression in (('X', x), ('Z', z)):
            if not expression.mask:
                continue
            control = parity_ready(expression, records, timing)
            if expression.indices:
                control += timing.feedforward_latency
            start = max(ready[q], control)
            ready[q] = start + timing.feedback_gate
            entries.append(dict(name=name, target=q, start=start, end=ready[q],
                                control=expression.as_dict()))
    return ready, entries
