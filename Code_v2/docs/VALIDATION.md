# Validation of the reorganized uploaded project

The source was rebuilt from the supplied project ZIP, including its Paulimer
backend and its redundant-measurement elimination pass. Historical main scripts,
notes and outputs are preserved separately; active code uses package imports.

Checks completed on 2026-10-15:

| Check | Result |
|---|---|
| Unit and integration tests | 23 passed |
| Operational decoder | 10 circuits through weight 20, 400 independent Stim checks passed |
| Unsimplified extraction | 10 circuits through weight 20, Choi-state and scalar checks passed |
| End-to-end main.py | Weight 9, seed 54; SVGs, circuits, formulas and JSON exported |
| Resource benchmark | 70 circuits, weights 4–10; JSON tables and plot exported |

The public-workflow tests exercise random and explicit input modes, signed
observables, identity padding, disabled simplification/resources/plots, optional
R4/R5 previews, file exports, stale-artifact removal, and invalid arguments.
Independent dense instrument checks cover probabilities and post-measurement
states. Simplification tests include incompatible intervening measurements,
reset/gate barriers and different partners.

The default example has 53 raw measurements and 28 after simplification: 25
redundant measurements removed before decoding. Logical outcome and correction
expressions use the resulting record indices.

Execution environment: Linux, Python 3.12, Paulimer 0.2.9, Binar 0.1.8,
PyZX 0.10.7, Stim 1.16.0. Windows was not executed here. The upstream Windows
wheels were checked during the preceding backend integration.

See [tests/README.md](../tests/README.md) for reproducible commands.

## Depth update (2026-10-15)

28 unit tests passed, including the published six-step gadget, shared-qubit conflicts, gate dependencies, weighted timing, and aggregated parity feedback. The 70-circuit benchmark (10 strings per weight 4–10, seed 20261015) completed with symbolic Paulimer validation. Its current summary and plot are in `depth_results/`. For weight 9, seed 54: unit quantum depth 9, measurement-only depth 6, pair-only depth 5, aggregated-feedback depth 11, literal-export feedback depth 18. These are distinct metrics, not alternative estimates of the same hardware time.
