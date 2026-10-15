# Resource analysis

`resources.py` counts resources from decoded circuits; it neither changes the
ZX rewrites nor minimizes correction formulas.

```python
from pauli_compiler.analysis.resources import pauli_resources, benchmark_resources

print(pauli_resources('XYZXYZ'))
print(pauli_resources(weight=9, seed=54))
rows, summary = benchmark_resources(weights=range(4,11), samples_per_weight=10)
```

To analyze an existing decoder without recompilation:

```python
from pauli_compiler import compile_operational_measurement
from pauli_compiler.analysis.resources import measure_circuit_resources

decoder = compile_operational_measurement('XYZXYZ')
metrics = measure_circuit_resources(decoder)
```

## Definitions

| Metric | Meaning |
|---|---|
| `data_qubits` | Number of logical inputs, including identity positions |
| `physical_qubits` | Number of extracted physical wires |
| `ancilla_qubits` | Physical wires minus data inputs; no reuse |
| `depth_pauli_frame` | Earliest-finish raw circuit depth, with corrections tracked in software |
| `depth_with_feedback` | Depth with one aggregated X and one Z predicate per output |
| `outcome_ready_depth` | Earliest time when all bits in the decoded outcome are available |
| `pair_measurements` / `single_measurements` | Counts of the respective physical operations |

By default, each preparation, Clifford gate, measurement and conditional Pauli takes one
unit. Disjoint operations may run in parallel. XOR computation and classical
communication are assigned zero duration. Physical feedback waits for the
measurement bits on which it depends. These are logical resource estimates,
not hardware timings or globally minimal depths.

The benchmark samples distinct positive-sign, full-support X/Y/Z strings within
each weight. `rows` contains each circuit's metrics; `summary` contains min, mean
and max grouped by weight. Corrections are not unique, so physical-feedback
counts may differ between valid decoders.

For JSON output and a graph, run from the project root:

```sh
python -m examples.benchmark --min-weight 4 --max-weight 10 --samples 10
```

## Updated depth model

See [Depth and timing](../../docs/DEPTH.md) for definitions, literature, timing parameters, and the change to aggregated feedback. `run_pipeline(..., timing=TimingModel(...))` accepts explicit hardware durations. Structural depths and timed latency are reported separately.
