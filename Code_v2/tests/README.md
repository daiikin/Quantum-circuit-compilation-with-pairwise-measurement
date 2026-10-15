# Correctness tests

Run from the project root with the dependencies in `requirements.txt` installed.
Tests use independent numerical or Stim calculations where possible; Paulimer
still performs the actual outcome/correction derivation.

## Fast unit and regression suite

```sh
python -m unittest discover -s tests -v
```

- `test_rewrites.py`: R4/R5 tensor equality for port permutations and X/Z colors,
  directed flow updates, and invalid-match rejection.
- `test_decoder.py`: deterministic record constraints, target signs and Y
  transposition, constant feedback, and rejection of overmeasurement, a wrong
  observable, non-Pauli correction requirements, and used-wire resets.
- Its dense test enumerates the small circuits for -Y, -XY, YZ and XXXX on
  arbitrary complex input vectors, checking unnormalized branch states and
  summed outcome probabilities.
- `helpers.py`: a dense target-projector oracle for small tests only.

The old unbounded-by-weight tensor self-test loops are not part of the normal
runtime modules. Dense calculations belong here and should stay small.

## Ten operational circuits

```sh
python -m tests.validate_decoder
```

Exactly ten instances, at weights 1,2,3,4,5,6,7,9,13,20. Includes signs, Y factors,
identity padding and both initial branch conventions. Each is verified
symbolically by Paulimer, then independently checked with 32 Bell-input shots
and 8 other stabilizer-input shots: 400 sampled checks total. Alternates exported
physical feedback and separately evaluated Pauli frames.

Writes formulas, corrected circuits and `validation.json` under
`outputs/validation_decoder/`.

## Ten raw extraction cases

```sh
python -m tests.validate_extraction
```

Tests the earlier stage: the fixed all-zero branch of the extracted circuit,
including its scalar normalization. Uses Stim Choi checks and small-system ZX
tensor comparisons. This complements, rather than replaces, the all-outcome
operational tests. Writes results under `outputs/validation_extraction/`.

## Integration smoke checks

```sh
python main.py
python -m examples.benchmark
```

These exercise package imports, SVG output, resource analysis and the end-to-end
example. The default benchmark analyzes 70 circuits. They are bounded checks,
not an exhaustive enumeration over arbitrary diagram sizes.

## Workflow and simplification regressions

`test_workflow.py` checks both input modes, reproducibility, optional plotting/
resource export, file overwrite behavior, and invalid argument combinations.
`test_simplification.py` checks commuting vs noncommuting intervening operations,
different partners, reset barriers, idempotence and post-decoding correctness.
