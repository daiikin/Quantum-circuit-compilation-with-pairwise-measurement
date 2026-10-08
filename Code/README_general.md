# Arbitrary-weight Pauli measurement compiler

Python 3.10+. Tested with PyZX 0.10.7 and Stim 1.16.0.

```bash
pip install pyzx==0.10.7 stim==1.16.0 numpy matplotlib
```

Keep the three implementation modules together. The general degree-reduction
and extraction functions are implemented in the original `zx_r4_r5.py`.

## Specific Pauli or random weight N

```python
from zx_r4_r5 import compile_measurement, generate_pauli_measurement

# Leftmost character is data qubit zero. I positions do not count as weight.
compiled = compile_measurement('-XYZIXYZ', outcome=-1)

# Exactly N non-identity factors, generated reproducibly.
m = generate_pauli_measurement(
    20, seed=123, total_qubits=24, sign=+1, outcome=+1
)
compiled = compile_measurement(m)

print(compiled.circuit.to_stim())
print(compiled.circuit.metadata())
print(compiled.reduced.trace)
```

Any positive integer weight is accepted; practical time/memory limits still
apply. No all-identity observable is accepted. To choose every Pauli factor,
pass a signed I/X/Y/Z string instead of using the random generator.

The result exposes four stages:

- `measurement.original`: normalized Pauli projector ZX graph.
- `prepared.graph`, `.paths`, `.order`: Proposition 5.1 MC-flow preparation.
- `reduced.graph`, `.paths`, `.trace`: Equation (11) decomposition to degree <=3.
- `circuit`: scheduled preparations, local Clifford gates, and weight-one/two
  measurements, with explicit boundary mappings and branch normalization.

For visualization, use `pyzx.draw(compiled.reduced.graph)`. Coordinates are
layout hints, not a gate schedule; the extracted operation order is the schedule.

## What is implemented

`reduce_degrees(graph, paths)` implements:

- r4 and r5 base cases;
- even splitting r_(2m) into two hubs and m four-legged spiders;
- odd splitting r_(2m+1), with the extra uncovered leg on one hub;
- two-leaf expansion for degrees congruent to 2 or 3 modulo 4;
- updates of every affected directed path, including newly introduced paths.

`extract_circuit(graph, paths)` implements the local cases of Proposition 4.11:

- one-legged endpoints -> normalized preparations/destructive measurements;
- degree-two spiders and H boxes -> single-qubit Clifford gates;
- same-colour leaf attachment -> identity;
- opposite-colour leaf attachment -> weight-one Pauli projection;
- equal-colour degree-three pairs -> XX or ZZ measurement;
- unequal-colour degree-three pairs -> CNOT (use `strict=False`).

Extraction requires a well-covered flow and reduced degree. Interaction
spiders and high-degree spiders must be phaseless; single-qubit Clifford
phases belong on degree-two vertices. H boxes must be degree two, phase one.
Unsupported diagrams fail explicitly. This is not extraction for arbitrary
non-Clifford ZX graphs or automatic MC-flow discovery on arbitrary graphs.

The Pauli compiler produces a strict flow and therefore uses no CNOTs during
extraction. The validator checks O1/O2, P1/P2/P3, PS, and PWC. It reconstructs
a partial order; caller-supplied extra scheduling constraints are not retained.
Use `check_each=True` when debugging intermediate degree-reduction steps.

## Postselection and exact semantics — read before executing

The paper's extraction convention uses +1-postselected measurement diagrams.
This implementation follows it. The output is an executable Clifford circuit,
but its equivalence claim is CONDITIONAL on every measurement record being 0.
It is not a deterministic, all-outcomes implementation of the Pauli measurement.
Automatic correction/Pauli-frame recovery for other internal outcomes is not
implemented. Repeated parity measurements are retained, not optimized away.

Let `K_phys` be the map after preparing the specified ancillary wires and
postselecting every measurement to zero. With data inputs and outputs ordered
as given by the metadata:

```
(I + outcome * signed_Pauli)/2 = branch_scale * K_phys
```

The classical eigenvalue requested by `outcome` is encoded in the diagram's
basis corrections; it is NOT the value of an arbitrary single Stim record.
To represent both ideal measurement branches, compile both outcomes separately.

`input_wires` and `output_wires` can differ: preparation moves quantum data
between physical worldlines. Respect these mappings when composing circuits.
Ancilla preparations are explicit in the operation stream. Input wires are
not reset by the extracted circuit. Stim's default |0> initialization alone
does not test an arbitrary input state.

Metadata includes every required postselection record and the complex scalar.
The scalar can be large and success probability can be small; this prototype
uses floating-point scalar bookkeeping. Its MC-flow closure algorithm is also
intended for research-sized examples, not unlimited-scale compilation.

## Exactly ten end-to-end validation circuits

```bash
python validate_10.py
```

The fixed reproducible suite uses weights 1, 2, 3, 4, 5, 6, 7, 9, 13, and 20,
with random X/Y/Z factors and variation in identity positions, sign and outcome.
There is no exhaustive enumeration of Pauli strings. The ten generated `.stim`
files and matching metadata JSONs are included in `generated_10/`.

For each circuit, the check prepares Bell pairs between every data input and
an untouched reference, executes the extracted operations with postselection,
and compares output/reference stabilizers with the independently specified
ideal Pauli projector's Choi state. It also checks the success probability
against `branch_scale`. Thus it checks the full conditional linear map up to
global phase, not just a few input bitstrings. A few small cases additionally
use direct reduced-ZX tensor contraction. It does not re-prove fault-distance
preservation. The signed examples and both outcomes all passed.

The generated family exercises the even-degree decomposition path needed
after preparation. General odd-degree and non-strict CNOT cases are implemented
but are not separately covered by these ten end-to-end tests.

`generated_10/validation.json` contains the recorded results. The older local
rule tests remain optional; running `validate_10.py` is the recommended check.

Source: Rodatz, Poor, Kissinger, *Floquetifying stabiliser codes with
distance-preserving rewrites*, Theorems 3.8/3.9/3.12/3.13,
Propositions 4.11/4.14/5.1:
https://quantum-journal.org/papers/q-2026-09-03-2202/pdf/
