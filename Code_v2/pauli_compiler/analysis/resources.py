"""Resource counts with explicit scheduling assumptions; see docs/DEPTH.md."""
from dataclasses import asdict
import random
from statistics import mean
from .scheduling import TimingModel, MEASUREMENTS, schedule, parity_ready, feedback_schedule


def measure_circuit_resources(decoder, *, timing=None):
    """Structural depths are unitless; latency fields use the supplied timing.

    Compatibility depth_with_feedback now uses aggregated parity controls.
    depth_literal_stim_feedback retains the previous per-term export metric.
    """
    timing = timing or TimingModel()
    c = decoder.compilation.circuit
    unit = TimingModel()
    ready, records, _ = schedule(c.operations, c.num_qubits, unit)
    if len(records) != len(decoder.measurement_operations):
        raise ValueError('Measurement record indexing mismatch')
    unit_depth = max(ready, default=0)
    corrected, feedback = feedback_schedule(decoder, ready, records, unit)
    literal = list(ready)
    literal_count = 0
    for q, x, z in zip(c.outputs, decoder.x_corrections, decoder.z_corrections):
        for expr in (x, z):
            if expr.constant:
                literal[q] += 1
                literal_count += 1
            for i in expr.indices:
                literal[q] = max(literal[q], records[i]) + 1
                literal_count += 1
    # Projected critical-path metrics: preserve every dependency, change costs.
    measurement_model = TimingModel(preparation=0, single_gate=0, two_qubit_gate=0)
    pair_model = TimingModel(preparation=0, single_measurement=0, single_gate=0, two_qubit_gate=0)
    md = max(schedule(c.operations, c.num_qubits, measurement_model)[0], default=0)
    pd = max(schedule(c.operations, c.num_qubits, pair_model)[0], default=0)
    timed, timed_records, entries = schedule(c.operations, c.num_qubits, timing)
    timed_corrected, timed_feedback = feedback_schedule(decoder, timed, timed_records, timing)
    outcome = parity_ready(decoder.outcome, timed_records, timing)
    frame_ready = max((parity_ready(e, timed_records, timing)
                       for e in (*decoder.x_corrections, *decoder.z_corrections)), default=0)
    m = decoder.compilation.measurement
    return {
        'pauli': ('-' if m.sign < 0 else '+') + m.pauli,
        'weight': len(m.support), 'data_qubits': len(c.inputs),
        'physical_qubits': c.num_qubits, 'ancilla_qubits': c.num_qubits-len(c.inputs),
        'quantum_depth': unit_depth, 'measurement_depth': md, 'pair_measurement_depth': pd,
        'depth_pauli_frame': unit_depth,
        'depth_with_feedback': max(corrected, default=0),
        'depth_literal_stim_feedback': max(literal, default=0),
        'outcome_ready_depth': parity_ready(decoder.outcome, records, unit),
        'quantum_latency': max(timed, default=0), 'outcome_ready_time': outcome,
        'pauli_frame_ready_time': frame_ready,
        'latency_pauli_frame': max(max(timed, default=0), outcome, frame_ready),
        'latency_with_feedback': max(max(timed_corrected, default=0), outcome),
        'pair_measurements': sum(op.name in ('MXX', 'MYY', 'MZZ') for op in c.operations),
        'single_measurements': sum(op.name in ('M', 'MX', 'MY') for op in c.operations),
        'quantum_operations': len(c.operations),
        'feedback_operations': len(feedback), 'literal_feedback_operations': literal_count,
        'timing_model': asdict(timing), 'schedule': entries, 'feedback_schedule': timed_feedback,
        'scheduling_assumptions': 'Fixed wire order; disjoint supports parallel; native pair measurements; '
            'unlimited classical parallelism; no routing or shared-control constraints. '
            'Feedback aggregates parity into X then Z. See docs/DEPTH.md.',
    }


def pauli_resources(pauli=None, *, weight=None, seed=None, timing=None):
    """Compile a specific signed string OR generate a random positive-weight one.

    pauli_resources('XYZX')
    pauli_resources(weight=20, seed=42)
    """
    from ..zx.spiders import generate_pauli_measurement
    from ..pipeline import compile_measurement
    from ..decoding.paulimer_backend import derive_decoder
    if (pauli is None)==(weight is None):
        raise ValueError('Specify either pauli or weight, not both')
    m = pauli if pauli is not None else generate_pauli_measurement(weight,seed=seed)
    return measure_circuit_resources(derive_decoder(compile_measurement(m)), timing=timing)


def benchmark_resources(*, weights=range(4,11), samples_per_weight=10, seed=20261009, timing=None):
    """Uniformly sample distinct full-support X/Y/Z strings for each weight.

    Positive sign, reference branch +1, no identity padding. The seed fixes
    the sample; compilation/decoding validates each circuit symbolically.
    """
    if not isinstance(samples_per_weight,int) or samples_per_weight<1:
        raise ValueError('samples_per_weight must be positive')
    rng=random.Random(seed)
    rows=[]
    for weight in weights:
        if not isinstance(weight,int) or weight<1 or samples_per_weight>3**weight:
            raise ValueError('Invalid weight or too many distinct requested strings')
        words=[]
        seen=set()
        while len(words)<samples_per_weight:
            word=''.join(rng.choice('XYZ') for _ in range(weight))
            if word not in seen:
                words.append(word); seen.add(word)
        for j,word in enumerate(words):
            row=pauli_resources(word, timing=timing)
            row.update(sample=j+1,seed=seed)
            rows.append(row)
    summary=[]
    for weight in sorted({r['weight'] for r in rows}):
        group=[r for r in rows if r['weight']==weight]
        s={'weight':weight,'samples':len(group)}
        for key in ('physical_qubits','ancilla_qubits','depth_pauli_frame','depth_with_feedback','measurement_depth','pair_measurement_depth','quantum_depth','depth_literal_stim_feedback','latency_pauli_frame','latency_with_feedback','pair_measurements'):
            values=[r[key] for r in group]
            s[key]={'min':min(values),'mean':mean(values),'max':max(values)}
        summary.append(s)
    return rows,summary


