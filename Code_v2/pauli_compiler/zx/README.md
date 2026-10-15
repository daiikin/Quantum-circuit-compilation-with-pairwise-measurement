# ZX diagrams and rewrites

This folder handles diagrams and directed MC-flow paths. It does not produce
classical outcome formulas or append quantum feedback.

## Files

| File | What it does |
|---|---|
| `spiders.py` | Parse signed Pauli strings and construct a normalized fixed-outcome projector as a ZX star |
| `rewrites.py` | Local R4/R5 templates, plus `rewrite_with_flow` for updating two paths through a matched spider |
| `flow.py` | Validate path coverage and causal ordering, including the explicit Hadamard-box representation |
| `preparation.py` | Add the leaf/gadget structure needed for MC flow on the original measurement star |
| `reduction.py` | Apply R4/R5, higher-degree splitting and leaf expansion while updating paths |

## Construct a diagram

```python
from pauli_compiler.zx.spiders import pauli_measurement, generate_pauli_measurement

specific = pauli_measurement('-XYZIX', outcome=+1)
random = generate_pauli_measurement(6, total_qubits=8, seed=54)
print(random.pauli, random.support)
graph = random.original
```

`outcome` is +1 or -1 for the signed Pauli and specifies the mathematical projector
`(I + outcome*P)/2`. It does not predetermine a physical run's measurement result.
The graph scalar preserves normalization.

`measurement.rewritten` is an optional direct R4/R5 preview. It is **None for
weights other than 4 and 5**. Always check before drawing it. It is not the input
to the operational pipeline.

## Prepare and reduce

```python
from pauli_compiler.zx.preparation import prepare_mc_flow
from pauli_compiler.zx.reduction import reduce_degrees
from pauli_compiler.zx.flow import validate_flow

prepared = prepare_mc_flow(random)
reduced = reduce_degrees(prepared.graph, prepared.paths)
validate_flow(reduced.graph, reduced.paths)
print(reduced.trace)
```

The preparation routine expects the original star and metadata produced by
`spiders.py`, not an arbitrary hand-edited ZX graph. Preparation handles odd
weight by adding a leaf and converts Hadamard edges to explicit H boxes.
Reduction preserves a well-covered MC flow and reduces high-degree zero-phase
spiders to degree at most three. `check_each=True` validates each intermediate
reduction step; it costs additional time.

Functions copy graphs before rewrites rather than modifying the supplied graph.
Use the **returned** graph and paths together; paths refer to vertex IDs in that
graph. A flow layout changes drawing coordinates, not the reduction rules.

For manual R4/R5 calls, ports must be distinct neighbors in the intended cyclic
order, incident edges ordinary, and the spider phase zero. Invalid matches are
rejected. MC-flow validation uses a transitive closure intended for modest
research examples, not a claim of optimized full-code scaling.
