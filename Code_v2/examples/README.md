# Examples

Run from the project root. `main.py` is the short editable example calling
`run_pipeline`; it does not parse CLI arguments.

```sh
python main.py
```

## Optional CLI: measurement.py

This is a wrapper around the same function, with no duplicate compilation code.
Run it as a module to resolve package imports correctly:

```sh
python -m examples.measurement --weight 9 --seed 54 --output outputs/random
python -m examples.measurement --pauli=-XYZXYZ --output outputs/specific
python -m examples.measurement --weight 6 --seed 54 --no-plots --no-simplify --overwrite
python -m examples.measurement --help
```

Exactly one of --pauli/--weight is required. Other options include --total-qubits,
--sign, --outcome, --check-each, --no-resources, --show, --labels, --hide-checks,
--record-legend, --format svg/png and --overwrite. CLI output defaults to
`outputs/cli`. Random-only parameters are rejected for an explicit Pauli.

## Resource benchmark: benchmark.py

```sh
python -m examples.benchmark
python -m examples.benchmark --min-weight 4 --max-weight 6 --samples 3 --seed 12
```

The default analyzes ten strings at every weight from 4 to 10 and saves
circuits.json, summary.json and resources.svg to `outputs/benchmark`.
It uses the same compiler and Paulimer decoder, with simplification enabled.
It does not render all intermediate diagrams for every benchmark case.
The batch benchmark is deliberately separate from the single-measurement API.
