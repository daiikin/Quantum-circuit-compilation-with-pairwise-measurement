# Package and pipeline

For the one-call public interface use `from pauli_compiler import run_pipeline`.
See [the API guide](../docs/PIPELINE_API.md). `workflow.py` coordinates decoding,
resources, rendering and exports. `export.py` contains file-writing only.
The lower-level functions below remain available for experiments.

`pipeline.py` orchestrates the stages; it contains no rewrite or simulator algorithm.
Start reading here to understand the complete workflow, then follow the function
calls into the relevant subpackage.

## Data passed between stages

| Step | Function | Result | Useful fields |
|---|---|---|---|
| 1 | `pauli_measurement` / `generate_pauli_measurement` | `Measurement` | `original`, `pauli`, `sign`, `support` |
| 2 | `prepare_mc_flow` | `PreparedMeasurement` | `graph`, `paths`, `order` |
| 3 | `reduce_degrees` | `ReducedDiagram` | `graph`, `paths`, `order`, `trace` |
| 4 | `extract_circuit` | `ExtractedCircuit` | `operations`, `inputs`, `outputs`, `branch_scale` |
| 5 | `derive_decoder` | `MeasurementDecoder` | `outcome`, `x_corrections`, `z_corrections`, `checks` |

`compile_measurement(pauli, outcome=1, check_each=False)` performs steps 1–4, then removes redundant measurements by default.
Use `simplify_measurements=False` to keep the literal extraction.
Its `Compilation` result stores `measurement`, `prepared`, `reduced`, `circuit` (simplified), and
`raw_circuit` (literal extraction).
Passing an existing `Measurement` uses its stored outcome.

`compile_operational_measurement(...)` also performs step 5 and returns a decoder.
The compilation remains accessible as `decoder.compilation`.

## Follow each stage explicitly

```python
from pauli_compiler.zx.spiders import pauli_measurement
from pauli_compiler.zx.preparation import prepare_mc_flow
from pauli_compiler.zx.reduction import reduce_degrees
from pauli_compiler.circuits.extraction import extract_circuit
from pauli_compiler.pipeline import Compilation
from pauli_compiler.decoding.paulimer_backend import derive_decoder

measurement = pauli_measurement('-XYZXYZ')
prepared = prepare_mc_flow(measurement)
reduced = reduce_degrees(prepared.graph, prepared.paths)
circuit = extract_circuit(reduced.graph, reduced.paths)
compilation = Compilation(measurement, prepared, reduced, circuit)
decoder = derive_decoder(compilation)
print(decoder.formulas())
```

Construction, preparation, reduction, extraction and decoding are independent of
Matplotlib rendering. Visualizations consume their results. Paulimer is loaded
when decoding is requested, so the ZX-only stages do not invoke a simulator.
There are no `sys.path` modifications or imports from the previous flat files.
