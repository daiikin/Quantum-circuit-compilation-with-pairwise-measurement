# Depth and timing of pairwise measurements

## Literature and scope

- Grans-Samuelsson et al., *Improved Pairwise Measurement-Based Surface Code*, Quantum 8, 1429 (2024), Section 2, Eq. (1), Figs. 2–4: https://quantum-journal.org/papers/q-2024-08-02-1429/pdf/ . The single gadget takes six steps; shared boundary measurements allow repeated execution in 4r+2 steps, hence period four. The scheduling test transcribes its twelve-measurement gadget and recovers depth six. Section 6 discusses additional hardware conflicts.
- Gidney, *A Pair Measurement Surface Code on Pentagons*, Quantum 7, 1156 (2023), Section 2, Figs. 2 and 6: https://quantum-journal.org/papers/q-2023-10-25-1156/pdf/ . The six additional layers per syndrome round refer to partially pipelined X/Z stabilizer measurements across the code, not the latency of an arbitrary isolated Pauli gadget.

These papers motivate explicit native-operation schedules. Our compiler schedules one isolated gadget; it does not implement repeated-round pipelining, hardware routing, or the fault-tolerant surface-code constructions of these papers.

## Exact scheduling rule implemented

For operation v with support Q(v), duration d(v), and previous operation p(q) on each occupied wire, compute

    start(v) = max({end(p(q)): q in Q(v)}, default=0)
    end(v) = start(v) + d(v)
    depth = max_v end(v).

This is the weighted longest path of the fixed wire-order dependency graph. With unit durations it is the minimum number of layers for this fixed order under the stated exclusive-qubit model. Disjoint operations may overlap. Sharing even one qubit prevents overlap, including commuting measurements. A native pair measurement occupies BOTH qubits for its full duration; it is not automatically decomposed into CNOTs. Gates/resets retain dependencies even in metrics assigning them zero cost. We do not search over commuting reorderings or alternate ZX rewrites, so this is not the globally smallest equivalent circuit depth.

There are no global barriers between list entries. Renderer columns are for readability and do not define depth. Exported Stim files contain no timing annotations; resources.json contains the actual weighted schedule with original operation and measurement-record indices.

## Metrics

| Field | Meaning |
|---|---|
| quantum_depth / depth_pauli_frame | Unit-cost depth of preparation, gates, and single/pair measurements; no physical corrections |
| measurement_depth | Critical-path measurement count: single/pair measurements cost 1; preparation and gates cost 0 |
| pair_measurement_depth | Critical-path pair count: pair measurements cost 1, everything else costs 0 |
| depth_with_feedback | Unit-cost depth after parity aggregation, at most one X followed by one Z per output |
| depth_literal_stim_feedback | Unit-cost depth of the actual per-record correction gates in corrected_stim() |
| quantum_latency | Completion of quantum operations using TimingModel durations |
| outcome_ready_time | Global measurement parity becomes available |
| pauli_frame_ready_time | All X/Z correction parity values become available |
| latency_pauli_frame | Quantum operation, outcome, and software-frame completion |
| latency_with_feedback | Quantum circuit and aggregated physical correction completion, including outcome availability |

The two measurement-only metrics are projections of the dependency graph, NOT full physical latency when gates or preparation take time. For a measurement-only gadget, quantum_depth and measurement_depth agree. Do not compare pair_measurement_depth with papers that count single-qubit measurements too. Ancillas are not reused.

## Classical computation and corrections

If a correction is X^(b0 xor b3 xor b7), three record-controlled X gates are algebraically equivalent to computing one parity and applying one controlled X. The latter requires one quantum correction slot, not three. We report both implementations; the Stim export remains the literal sequence and its exact depth has its own field. Aggregated feedback is an assumed controller implementation, not a rewritten Stim instruction.

For k nonconstant inputs, the model uses a balanced binary XOR tree with ceil(log2(k)) levels, starting when all inputs are available. Its completion time is max(input-ready times) + ceil(log2(k))*xor_time. This is a specified implementable tree model, not an optimum for unequal arrival times. Constant offsets are folded into the control predicate. Unlimited parallel classical processors and fan-out are assumed. Each conditional correction additionally waits feedforward_latency, then occupies its qubit for feedback_gate. X precedes Z on each output; different outputs overlap. This is a worst-case scheduled slot cost, independent of whether a particular record suppresses a correction. Corrections may start when their output wire and required records are ready; there is no artificial barrier waiting for unrelated ancillas.

Defaults: all quantum durations 1, XOR and communication 0. No claims of hardware runtime follow without hardware durations and conflict constraints. Pauli-frame tracking needs no extra quantum gate, but its classical completion time is included in latency_pauli_frame.

## Usage

```python
from pauli_compiler import run_pipeline, TimingModel

result = run_pipeline(weight=9, seed=54, render=False,
    timing=TimingModel(preparation=1, single_measurement=3,
                       pair_measurement=5, single_gate=0.1,
                       feedback_gate=0.1, xor_time=0.05,
                       feedforward_latency=0.2))
print(result.resources['quantum_depth'])       # dimensionless, unit costs
print(result.resources['measurement_depth']) # measurement-only projection
print(result.resources['latency_with_feedback']) # chosen time units
```

`measure_circuit_resources`, `pauli_resources` and `benchmark_resources` also accept `timing=TimingModel(...)`. Custom timing changes latency fields, not structural depths. Durations must be finite and nonnegative; unsupported operations raise errors. The model has no adaptive measurements inside the extracted circuit; feedback is terminal. Additional bus, connectivity, and crosstalk conflicts are not modeled.

## Migration from the previous version

The earlier raw ASAP depth was already correct for its unit-cost assumptions. Previously depth_with_feedback counted every XOR term as a separate gate. That value is now explicitly depth_literal_stim_feedback; depth_with_feedback uses aggregated X/Z predicates. Likewise feedback_operations now counts aggregated gates; literal_feedback_operations retains the per-term count. Existing benchmarks should be regenerated, not silently mixed with old results. Run `python -m examples.benchmark` for ten Paulis of each weight 4–10.
