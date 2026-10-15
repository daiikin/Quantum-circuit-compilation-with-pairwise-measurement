# Public end-to-end API

```python
from pauli_compiler import run_pipeline
```

`run_pipeline` is the recommended entry point for `main.py`, a notebook or a
larger experiment. It calls the existing stages rather than reimplementing them.

## Minimal main.py

```python
from pauli_compiler import run_pipeline


def main():
    result = run_pipeline(
        weight=9,
        seed=54,
        simplify_measurements=True,
        output_dir='outputs/my_run',
        overwrite=True,
    )
    print(result.decoder.formulas())
    print(result.resources)
    result.close_figures()


if __name__ == '__main__':
    main()
```

For a fixed observable replace `weight=9, seed=54` with `pauli='-XYZXYZ'`.
The string can also be the first positional argument:

```python
result = run_pipeline('-XYZXYZ', render=False)
```

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `pauli` | `None` | Explicit signed string containing I/X/Y/Z; mutually exclusive with weight |
| `weight` | `None` | Positive number of nonidentity factors for random generation |
| `seed` | `None` | Integer random seed; supply it for reproducibility |
| `total_qubits` | `None` | Random input width, at least weight; adds identity positions |
| `sign` | `None` | Random observable sign, +1 or -1; omitted means +1 |
| `outcome` | `1` | Initial ZX branch eigenvalue, +1 or -1; not a forced run outcome |
| `simplify_measurements` | `True` | Remove repeated noiseless Pauli measurements before decoding |
| `check_each` | `False` | Validate MC flow after each reduction step |
| `calculate_resources` | `True` | Calculate resources for the simplified, decoded circuit |
| `render` | `True` | Create all available stage figures |
| `show` | `False` | Also display figures using Matplotlib; requires render=True |
| `labels` | `False` | Draw ZX vertex labels |
| `show_checks` | `True` | Include noiseless record constraints in the circuit figure |
| `show_record_legend` | `False` | Include every record bit's operation/wire mapping |
| `output_dir` | `None` | Write artifacts here; omitted means no filesystem output |
| `image_format` | `'svg'` | `'svg'` or `'png'`; SVG is preferable for large circuits |
| `overwrite` | `False` | Allow reuse of a nonempty output directory |

When a specific Pauli is supplied, do not also supply `seed`, `total_qubits`, or
`sign`. Put its sign and identity factors directly in the string. Contradictory
or invalid inputs raise errors instead of silently choosing an interpretation.

There is no implicit default observable: one of pauli/weight is required.
The illustrative defaults belong to `main.py`, not the API.

## Order of operations

1. Build the normalized Pauli-projector spider representation.
2. Prepare MC flow from the original diagram.
3. Reduce high-degree spiders, updating directed paths.
4. Extract the ordered physical operations and boundary mappings.
5. Optionally eliminate redundant measurements.
6. Derive global outcome parities and output corrections using Paulimer/Binar.
7. Verify the corrected measurement instrument across all allowed records.
8. Optionally calculate resource counts.
9. Optionally render initial, direct-rewrite, prepared, flow-layout and circuit figures.
10. Optionally export and display results.

The direct R4/R5 preview exists only for weights 4 and 5. It is a side
visualization; operational compilation starts from the original diagram.
Redundancy removal acts on the extracted circuit, not on the saved ZX graph.

## Return value: PipelineResult

| Attribute | Contents |
|---|---|
| `measurement` | Pauli string, sign, support, original diagram and optional direct rewrite |
| `compilation` | Measurement, prepared/reduced diagrams, paths and extracted circuits |
| `raw_circuit` | Unsimplified `ExtractedCircuit` |
| `simplified_circuit` | `ExtractedCircuit` used by Paulimer, without feedback |
| `decoder` | Outcome/correction formulas, record checks, metadata and evaluation methods |
| `corrected_circuit` | Stim circuit including physical output feedback |
| `statistics` | Raw/simplified operation and measurement counts, including number removed |
| `resources` | Resource metrics, or None if disabled |
| `parameters` | Requested run settings |
| `figures` | Dictionary of Matplotlib figures, or empty if rendering is disabled |
| `files` | Dictionary of absolute saved-file Paths, or empty if export is disabled |

`result.close_figures()` closes only this run's figures and empties `figures`.
Call it after saving or displaying them in a batch. This function does not delete
files. `show=False` does not change the application's global Matplotlib backend;
for headless use you can set `MPLBACKEND=Agg`. The CLI selects Agg when not showing.

## Output files

| File | Meaning |
|---|---|
| `01_initial.svg` | Initial fixed-branch ZX diagram |
| `02_direct_rewrite.svg` | Optional direct R4/R5 preview |
| `03_prepared.svg` | MC-flow-prepared diagram |
| `04_reduced.svg` | Reduced diagram with flow layout |
| `05_circuit.svg` | Simplified operational circuit with feedback |
| `raw_extracted.stim` | Literal extracted operations before simplification |
| `before_feedback.stim` | Simplified operations used for decoding |
| `measurement.stim` | Simplified circuit plus output feedback |
| `formulas.txt` | Global outcome and correction formulas |
| `measurement.json` | Decoder metadata and physical wire/record mapping |
| `resources.json` | Optional resource counts |
| `run.json` | Parameters, counts, rewrite trace and artifact filenames |

Figure extension changes to `.png` when requested. Formulas refer to
`before_feedback.stim` / `measurement.stim`, not to the old record indices in
`raw_extracted.stim`. The corrected Stim file does not embed an observable label;
evaluate the returned decoder formula, or read the formula in measurement.json.

By default the output directory must be new or empty. With overwrite=True,
only the workflow's fixed artifact names are replaced or cleared; unrelated
files survive. Optional files from earlier runs (such as the direct preview)
are removed to avoid stale results. Use distinct output directories to preserve
several experiments.

## Evaluate a physical record

```python
import stim
from pauli_compiler import run_pipeline

result = run_pipeline('XXXX', render=False)
sim = stim.TableauSimulator(seed=42)
sim.do_circuit(result.corrected_circuit)
record = sim.current_measurement_record()
print(result.decoder.outcome_eigenvalue(record))
print(result.compilation.circuit.outputs)
```

This demonstration starts data in |0...0>; the symbolic certificate covers
arbitrary data inputs. The Bell reference qubits used for analysis are not
present in the physical circuit.

## Updated depth model

See [DEPTH.md](DEPTH.md) for definitions, literature, timing parameters, and the change to aggregated feedback. `run_pipeline(..., timing=TimingModel(...))` accepts explicit hardware durations. Structural depths and timed latency are reported separately.
