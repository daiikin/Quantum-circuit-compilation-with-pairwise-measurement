# Paulimer backend for the ZX Pauli-measurement compiler

This update replaces the symbolic stabilizer simulator in `measurement_decoder.py`
with Microsoft QDK-EC's `paulimer.OutcomeCompleteSimulation`. Binary linear systems
are solved by `binar.EchelonForm`, not by a local elimination implementation.

## Install and run

Place the files in your existing `Code` directory (keep a backup of the previous
`measurement_decoder.py` before replacing it). The ZIP also contains the compiler
and renderer versions used for testing, for a self-contained example.

```sh
python -m pip install -r requirements-paulimer.txt
python main_paulimer.py --weight 9 --seed 54
python main_paulimer.py --pauli=-XYZXYZ --output example_negative
```

`main_paulimer.py` writes SVG diagrams, a corrected Stim circuit, a formula listing,
and JSON metadata. Use `--no-plots` to skip figures. Run in a new output directory
for each example if you want to retain earlier results.

Your existing main.py imports stay the same:

```python
from zx_r4_r5 import generate_pauli_measurement, compile_measurement
from measurement_decoder import derive_decoder
from circuit_renderer import draw_measurement_circuit

measurement = generate_pauli_measurement(9, seed=54)
compilation = compile_measurement(measurement)
decoder = derive_decoder(compilation)
print(decoder.formulas())
print(decoder.corrected_stim())
fig = draw_measurement_circuit(decoder)
fig.savefig('circuit.svg', bbox_inches='tight')
```

Continue to guard the optional direct preview:

```python
if measurement.rewritten is not None:
    save_zx(measurement.rewritten, 'Direct R4/R5 preview', '02_direct.svg')
```

For weights outside 4 and 5, use `compilation.reduced.graph` for the general
rewritten diagram. The actual compilation always starts with the original graph.

## What Paulimer now does

1. Simulates the extracted preparation, Clifford and Pauli-measurement operations
   on half of an input Bell state, keeping all branches symbolically.
2. Supplies the output Clifford, stabilizer-sign matrix, outcome matrix, and
   outcome shift. These capture all allowed measurement records.
3. Supplies Clifford preimages to determine the signs of reference observables
   and the target Choi stabilizers.
4. Applies the derived conditional corrections and verifies the full ideal Choi
   stabilizer basis for every allowed record.

The adapter constructs the ideal measurement, maps wire and outcome IDs, and
forms the linear equations. Binar solves those equations. The Python binding
exposes the primitives used here; this is an adapter built on those APIs, not
an invocation of a built-in `zx_to_circuit` or `derive_corrections` function.
There is no fallback to the old in-house stabilizer simulator. Stim is used only
for output format/Pauli frames and for independent validation tests.

The existing ZX rewrites, MC-flow preparation, and extraction are unchanged.
Paulimer consumes the ordered extracted circuit; it does not itself parse ZX.
The current entry point supports the arbitrary-weight Pauli diagrams accepted
by this compiler, not every possible ZX diagram.

## Outputs and conventions

- `decoder.outcome`: affine XOR expression for the signed target Pauli result.
- `decoder.outcome_bit(record)`: evaluate that expression; bit 0 means +1.
- `decoder.outcome_eigenvalue(record)`: return +1 or -1.
- `decoder.x_corrections`, `decoder.z_corrections`: conditional output Paulis,
  indexed in logical output order, not physical-wire order.
- `decoder.output_pauli_frame(record)`: evaluate the output frame as a Stim Pauli.
- `decoder.corrected_stim()`: append physical X/Z feedback to the raw circuit.
- `decoder.checks`: XOR relations that must vanish for noiseless records.
- `decoder.measurement_operations`: map `b0`, `b1`, ... to original operations.
- `decoder.metadata()`: JSON-ready formulas, mappings, and verification metadata.

The raw record includes all single-qubit and pair measurements in chronological
extraction order. The decoder does not blindly XOR all pair outcomes. A final
formula may use single-qubit readouts when required by the actual circuit.
No Bell reference qubits are added to the emitted physical circuit.
The initially selected ZX branch is a construction convention: the emitted
operational circuit accepts both measurement outcomes, with no postselection.
Apply the physical feedback OR track the returned Pauli frame, not both.

Example for this compiler's `XXXX` circuit:

```
r = b4 xor b5
D0 (physical p0): X^(b2 xor b7), Z^(b0)
D1 (physical p4): X^(b2 xor b6 xor b7), Z^(b4)
D2 (physical p2): X^(b7), Z^(b1)
D3 (physical p5): X^(0), Z^(b5)
Noiseless check: b2 xor b3 = 0
```

Correction formulas are not unique. They can differ from the old decoder by
stabilizers or by relations between valid records. Resource counts for physical
feedback can therefore change; this solver does not minimize feedback depth.

## Mathematical extraction

Paulimer represents the final pure state of the data, retained physical ancillas,
and references as `R |A u>`, with uniformly random independent bits `u`.
The measured record is `b = M u + v0`, over GF(2).

For an observable Q, Paulimer computes `R^dagger Q R`. If this is
`(-1)^c Z^v`, its eigenvalue bit is `c + v^T A u`. Any X/Y factor means the
observable is not definite in a branch; a required target observable then fails
verification.

Apply that rule to `P^T` on the reference register. This yields the global outcome
bit. The transpose includes one minus sign for every Y, as well as the target's
explicit sign. Binar solves `M^T a = h` to express any free-bit parity `h^T u`
as `a^T b + a^T v0` in actual measurement records.

A second Paulimer simulation builds the ideal Choi state for the signed target
Pauli measurement. For each of its 2n independent stabilizers, compare actual
and ideal signs. The commutation of output `X^x Z^z` with that stabilizer gives
one GF(2) equation in x,z. Solve for the constant and every free-bit coefficient
using one Binar echelon decomposition. A missing solution is an error.

Finally, apply the resulting feedback through Paulimer and verify all 2n ideal
Choi generators. This also establishes separation from discarded ancillas.
A nonconstant affine decoded outcome is balanced over the k free Choi bits.
Each record has Choi probability 2^-k, so its corrected Kraus map has squared
scalar 2^(1-k) relative to the ideal projector. There are 2^(k-1) records per
outcome; their weights sum to one. Thus the check covers Born probabilities as
well as normalized output states, up to irrelevant branch-global phases.

## Supported scope and validation

Fresh R/RX preparations, single-qubit Clifford gates, CX/CZ/SWAP, and single or
pair Pauli measurements are supported by the adapter. The current renderer has
its original narrower vocabulary, which covers the compiler's emitted gates.
Resets of already-used wires and unsupported operations are rejected explicitly.
The input circuit is nonadaptive; feedback is synthesized at its outputs.
This is a noiseless instrument analysis, not a noise decoder or a fault-distance
proof. A circuit that overmeasures the input is rejected even if it yields the
correct global parity.

Tests performed in Linux with Paulimer 0.2.9 and Binar 0.1.8:

```sh
python -m unittest -v test_paulimer_decoder.py
python test_decoder_10.py
```

- 10 compiled instances at weights 1,2,3,4,5,6,7,9,13,20; signs, identity padding,
  Y factors, and both construction branches included.
- 400 independent Stim checks of corrected Choi states and other stabilizer
  inputs, alternating physical feedback and evaluated Pauli frames.
- Exhaustive dense-matrix branch checks for -Y, -XY, YZ, and XXXX on random
  non-stabilizer inputs, including probability normalization.
- Rejection tests for overmeasurement, wrong observables, a required non-Pauli
  correction, and resets of data; deterministic record constraints and constant
  signs/feedback are tested too.
- 70 resource-analysis cases: 10 random strings at each weight 4 through 10.
- Weight-9 graphical example successfully rendered.

`validation_10.json` and `paulimer_resources.json` contain the run results.

## Upstream references

- https://github.com/microsoft/qdk-ec
- https://github.com/microsoft/qdk-ec/tree/main/paulimer/bindings/python
- https://github.com/microsoft/qdk-ec/blob/main/paulimer/bindings/python/paulimer.pyi
- https://arxiv.org/abs/2309.08676 (Stabilizer circuit verification)

Versions are pinned for reproducibility; an upstream API change may require
updating the adapter. This project does not vendor the QDK-EC packages.
