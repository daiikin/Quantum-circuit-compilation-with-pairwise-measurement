# Moving from the flat Code directory

Extract the new project into a fresh directory first. Its contents can live under
your existing `Code` directory, but preserve the `pauli_compiler`, `examples`,
`tests`, and `docs` folders. Keep a backup of your existing scripts.

The uploaded redundancy-removal pass is retained in `circuits/simplification.py`.
The new recommended entry point is `from pauli_compiler import run_pipeline`;
see [PIPELINE_API.md](PIPELINE_API.md). Root `main.py` is now an editable function
call, while `python -m examples.measurement` supplies the CLI.

Run from the directory containing the new `main.py`:

```sh
python -m pip install -r requirements.txt
python main.py
```

The new package does not import any old flat module. There are no compatibility
wrappers: there is one implementation per stage. Update imports in personal
scripts using this table.

| Old location | New location |
|---|---|
| `zx_r4_r5.generate_pauli_measurement` | `pauli_compiler.generate_pauli_measurement` |
| `zx_r4_r5.compile_measurement` | `pauli_compiler.compile_measurement` |
| `pauli_measurements.pauli_measurement` | `pauli_compiler.zx.spiders.pauli_measurement` |
| `zx_r4_r5.r4`, `r5`, `rewrite_with_flow` | `pauli_compiler.zx.rewrites` |
| `zx_r4_r5.mc_order` | `pauli_compiler.zx.flow.mc_order` |
| `mc_preparation.validate_flow` | `pauli_compiler.zx.flow.validate_flow` |
| `mc_preparation.prepare_mc_flow` | `pauli_compiler.zx.preparation.prepare_mc_flow` |
| `mc_preparation.apply_prepared_r4` | `pauli_compiler.zx.preparation.apply_prepared_r4` |
| `zx_r4_r5.reduce_degrees` | `pauli_compiler.zx.reduction.reduce_degrees` |
| `zx_r4_r5.remove_redundant_measurements` | `pauli_compiler.circuits.simplification.remove_redundant_measurements` |
| `zx_r4_r5.extract_circuit` | `pauli_compiler.circuits.extraction.extract_circuit` |
| `zx_r4_r5.CircuitOperation`, `ExtractedCircuit` | `pauli_compiler.circuits.model` |
| `zx_r4_r5.Compilation` | `pauli_compiler.pipeline.Compilation` |
| `measurement_decoder.derive_decoder` | `pauli_compiler.decoding.paulimer_backend.derive_decoder` |
| `measurement_decoder.compile_operational_measurement` | `pauli_compiler.compile_operational_measurement` |
| `measurement_decoder.Affine` | `pauli_compiler.decoding.affine.Affine` |
| `measurement_decoder.MeasurementDecoder` | `pauli_compiler.decoding.result.MeasurementDecoder` |
| `flow_layout.draw_flow`, `layout_flow` | `pauli_compiler.visualization.flow` |
| `circuit_renderer.draw_measurement_circuit` | `pauli_compiler.visualization.circuit.draw_measurement_circuit` |
| `circuit_resources` functions | `pauli_compiler.analysis.resources` |

Replacement imports for your previously shared main.py:

```python
import matplotlib.pyplot as plt
import pyzx as zx
from pauli_compiler import generate_pauli_measurement, compile_measurement
from pauli_compiler.decoding.paulimer_backend import derive_decoder
from pauli_compiler.visualization.flow import draw_flow
from pauli_compiler.visualization.circuit import draw_measurement_circuit
from pauli_compiler.analysis.resources import benchmark_resources
```

Compilation/decoder field names and methods are unchanged. Continue checking
`measurement.rewritten is not None` before rendering the optional R4/R5 preview.
Use `result.reduced.graph` for general weights.

The earlier generator/preparation CLI demos and dense validation helpers have
been replaced by centralized examples and tests. The small-system projector
oracle is now `tests.helpers.expected_projector`. Generated results live under
`outputs/`, separate from the source code.

Do not use `from code...` or `from Code...`. Import `pauli_compiler` while running
from this project root. No `sys.path` edits are needed. Run the root launcher or
`python -m examples.measurement`, not individual package files by path.
