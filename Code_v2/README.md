# Pairwise Pauli measurement compilation

This project builds a ZX diagram for a signed Pauli measurement, prepares and
reduces its MC flow, extracts a physical circuit, removes redundant measurements,
and uses Microsoft Paulimer to derive the global outcome and output corrections.

## Run the project

Open a terminal in this directory, containing `main.py` and `pauli_compiler/`:

```sh
python -m pip install -r requirements.txt
python main.py
```

On Windows, `py -3.13` can replace `python`. Edit the arguments in `main.py` to
choose a different Pauli or configuration. The default example uses weight 9,
seed 54, and writes diagrams and results to `outputs/demo/`.

The package is named `pauli_compiler`, not `code` (a Python standard library
module). The enclosing project folder may still be named `Code`.

## One function for all stages

```python
from pauli_compiler import run_pipeline

# Reproducible random Pauli:
result = run_pipeline(weight=9, seed=54, output_dir='outputs/random')

# Or an explicitly specified signed Pauli:
other = run_pipeline('-XYZXYZ', output_dir='outputs/specific')

print(result.decoder.formulas())
print(result.statistics)
print(result.resources)
print(result.corrected_circuit)
result.close_figures()
other.close_figures()
```

Specify either `pauli` or `weight`, never both. `seed`, `total_qubits`, and `sign`
are random-generation options. Use `render=False` to skip figure creation,
`output_dir=None` to avoid file writes, and `overwrite=True` to replace a previous
run's known artifact files. The default refuses a nonempty output directory.

See [docs/PIPELINE_API.md](docs/PIPELINE_API.md) for every parameter, return value,
file name, and a minimal example for your own `main.py`.

## Code organization

| Location | Responsibility |
|---|---|
| [pauli_compiler/zx/](pauli_compiler/zx/README.md) | Spider construction, MC-flow preparation/validation, local rewrites, degree reduction |
| [pauli_compiler/circuits/](pauli_compiler/circuits/README.md) | Circuit model, extraction and redundant-measurement elimination |
| [pauli_compiler/decoding/](pauli_compiler/decoding/README.md) | Paulimer analysis, affine outcome formulas, correction and Stim export |
| [pauli_compiler/visualization/](pauli_compiler/visualization/README.md) | Flow layout, ZX stage figures, circuit rendering |
| [pauli_compiler/analysis/](pauli_compiler/analysis/README.md) | Depth, qubit/ancilla counts and benchmarks |
| [pauli_compiler/](pauli_compiler/README.md) | `pipeline.py` connects compilation; `workflow.py` runs everything; `export.py` writes results |
| [examples/](examples/README.md) | Optional command-line workflow and resource benchmark |
| [tests/](tests/README.md) | Unit tests, dense instrument checks and independent Stim validation |
| [docs/](docs/README.md) | API, migration, derivation and validation notes |
| [archive/](archive/README.md) | Previous outputs preserved from the uploaded project |

## Optional command-line interface

`main.py` is deliberately a short editable Python example. For command-line
arguments use:

```sh
python -m examples.measurement --weight 9 --seed 54 --output outputs/cli
python -m examples.measurement --pauli=-XYZXYZ --output outputs/negative
python -m examples.measurement --weight 6 --seed 54 --no-plots --overwrite
python -m examples.benchmark
```

The benchmark remains separate: it is not automatically run for each single
measurement. By default it compiles ten strings at every weight from 4 through 10.

## Semantics

The selected initial `outcome` specifies a mathematical ZX projector branch.
It does not predetermine the measured result. The decoded operational circuit
accepts both global outcomes, without postselecting the record to zero.

A physical bit `b_i=0` denotes eigenvalue +1. The global outcome is an affine XOR
of the simplified circuit's recorded bits. Use
`result.decoder.outcome_eigenvalue(record)` to evaluate it. Records can include
single-qubit readouts as well as pair measurements.

`result.raw_circuit` is the literal extraction. `result.simplified_circuit` is
the sequence sent to Paulimer. `result.corrected_circuit` is a Stim circuit with
physical feedback appended. Alternatively run the simplified circuit and track
`result.decoder.output_pauli_frame(record)` in software. Do not apply both.
Boundary mappings identify physical data outputs; ancillas may be discarded.

This supports the compiler's positive-weight Pauli diagram family and Clifford
extraction rules, not arbitrary ZX diagrams or noise decoding. Redundancy removal
preserves noiseless behavior, not fault-detection performance of repeated checks.

## Tests and migration

```sh
python -m unittest discover -s tests -v
python -m tests.validate_decoder
python -m tests.validate_extraction
```

See [docs/MIGRATION.md](docs/MIGRATION.md) to update older imports. The archive is
a source distribution: it does not contain `.git` or Python bytecode. Copy the
new source files into your existing Git checkout to retain its history; do not
replace your repository's `.git` directory.

## Updated depth model

See [docs/DEPTH.md](docs/DEPTH.md) for definitions, literature, timing parameters, and the change to aggregated feedback. `run_pipeline(..., timing=TimingModel(...))` accepts explicit hardware durations. Structural depths and timed latency are reported separately.
