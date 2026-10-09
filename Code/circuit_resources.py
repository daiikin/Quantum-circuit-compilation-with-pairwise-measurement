"""Resource counts for the Pauli measurement compiler.

Unit-time ASAP depth, unlimited parallel operations on disjoint qubits.
Prep, Clifford gate, single/pair measurement and controlled Pauli each cost 1.
Classical XOR and communication cost 0. No connectivity or qubit reuse assumed.
Feedback depth counts the actual primitive corrections exported by the decoder,
not a collapsed single Pauli-frame update. Run this file for the 70-circuit study.
"""
import json
import random
from pathlib import Path
from statistics import mean

MEASUREMENTS = frozenset(('M','MX','MXX','MZZ'))


def _schedule_raw(operations, num_qubits):
    ready = [0]*num_qubits
    record_ready = []
    for op in operations:
        end = max(ready[q] for q in op.targets)+1
        for q in op.targets:
            ready[q] = end
        if op.name in MEASUREMENTS:
            record_ready.append(end)
    return ready,record_ready


def measure_circuit_resources(decoder):
    """Measure an existing decoded circuit; returns JSON-serializable metrics.

    depth_pauli_frame is raw quantum depth (corrections tracked in software).
    depth_with_feedback schedules every X/Z correction term actually emitted
    by corrected_stim(), subject to qubit AND measurement-record dependencies.
    Neither depth is a globally optimized lower bound or physical runtime.
    """
    c = decoder.compilation.circuit
    ready,record_ready = _schedule_raw(c.operations,c.num_qubits)
    if len(record_ready)!=len(decoder.measurement_operations):
        raise ValueError('Measurement record indexing mismatch')
    raw_depth = max(ready,default=0)
    feedback_operations = 0
    # Keep the same order as corrected_stim(): X terms, then Z terms per output.
    for physical,x,z in zip(c.outputs,decoder.x_corrections,decoder.z_corrections):
        for expr in (x,z):
            if expr.constant:
                ready[physical] += 1
                feedback_operations += 1
            for index in expr.indices:
                ready[physical] = max(ready[physical],record_ready[index])+1
                feedback_operations += 1
    outcome_ready = max((record_ready[i] for i in decoder.outcome.indices),default=0)
    m=decoder.compilation.measurement
    return {
        'pauli':('-' if m.sign<0 else '+')+m.pauli,
        'weight':len(m.support),
        'data_qubits':len(c.inputs),
        'physical_qubits':c.num_qubits,
        'ancilla_qubits':c.num_qubits-len(c.inputs),
        'depth_pauli_frame':max(raw_depth,outcome_ready),
        'depth_with_feedback':max(max(ready,default=0),outcome_ready),
        'outcome_ready_depth':outcome_ready,
        'pair_measurements':sum(op.name in ('MXX','MZZ') for op in c.operations),
        'single_measurements':sum(op.name in ('M','MX') for op in c.operations),
        'quantum_operations':len(c.operations),
        'feedback_operations':feedback_operations,
    }


def pauli_resources(pauli=None, *, weight=None, seed=None):
    """Compile a specific signed string OR generate a random positive-weight one.

    pauli_resources('XYZX')
    pauli_resources(weight=20, seed=42)
    """
    from zx_r4_r5 import generate_pauli_measurement,compile_measurement
    from measurement_decoder import derive_decoder
    if (pauli is None)==(weight is None):
        raise ValueError('Specify either pauli or weight, not both')
    m = pauli if pauli is not None else generate_pauli_measurement(weight,seed=seed)
    return measure_circuit_resources(derive_decoder(compile_measurement(m)))


def benchmark_resources(*, weights=range(4,11), samples_per_weight=10, seed=20261009):
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
            row=pauli_resources(word)
            row.update(sample=j+1,seed=seed)
            rows.append(row)
    summary=[]
    for weight in sorted({r['weight'] for r in rows}):
        group=[r for r in rows if r['weight']==weight]
        s={'weight':weight,'samples':len(group)}
        for key in ('physical_qubits','ancilla_qubits','depth_pauli_frame','depth_with_feedback','pair_measurements'):
            values=[r[key] for r in group]
            s[key]={'min':min(values),'mean':mean(values),'max':max(values)}
        summary.append(s)
    return rows,summary


if __name__=='__main__':
    rows,summary=benchmark_resources()
    Path('resource_results.json').write_text(json.dumps(rows,indent=2))
    Path('resource_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
