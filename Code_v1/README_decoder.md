# Pauli outcome decoding and feedback

This adds the missing outcome/feedback layer to the existing ZX compiler.
Install Python 3.10+ dependencies:

```bash
python -m pip install pyzx==0.10.7 stim==1.16.0 numpy matplotlib
```

Extract the files together into your Code folder. Run `python main_decoder.py`.
Do not name your own script `pyzx.py`, `stim.py`, or `measurement_decoder.py`.

## Complete, non-postselected measurement

```python
import stim
from measurement_decoder import compile_operational_measurement

decoder = compile_operational_measurement('XYZXY')
print(decoder.formulas())
circuit = decoder.corrected_stim()

sim = stim.TableauSimulator()
# Initialize any desired input state on decoder.compilation.circuit.inputs.
# Without explicit initialization, Stim starts those inputs in |0>.
sim.do_circuit(circuit)
record = sim.current_measurement_record()
print('Pauli eigenvalue:', decoder.outcome_eigenvalue(record))
# Corrected data now lives on decoder.compilation.circuit.outputs.
```

This circuit accepts both Pauli outcomes. No postselection or branch-scale
rescaling is required. The compiled projector branch is used only to build
the quantum structure; the decoder returns the actual eigenvalue observed.

For an existing compilation:

```python
from zx_r4_r5 import generate_pauli_measurement, compile_measurement
from measurement_decoder import derive_decoder

m = generate_pauli_measurement(20, seed=123)
result = compile_measurement(m)
decoder = derive_decoder(result)
```

The original `result.circuit.to_stim()` still returns the RAW circuit without
corrections. Use `decoder.corrected_stim()` to include classical feedback.

## Two ways to handle corrections

1. Execute `decoder.corrected_stim()`. It appends record-controlled CX/CZ
   instructions (classical controls, not extra quantum two-body gates).
   Do NOT apply the returned Pauli frame again.
2. Execute the raw circuit and evaluate `decoder.output_pauli_frame(record)`.
   This PauliString is in DATA OUTPUT order. Apply it using the physical
   mapping in `result.circuit.outputs`, or track it through subsequent gates
   and measurements in your own Pauli-frame machinery.

Both methods use `decoder.outcome_bit(record)` or `.outcome_eigenvalue(record)`.
Each bit b_i labels the i-th measurement RESULT, not the i-th text line in
Stim (an MPP instruction can produce several results). `measurement_operations`
provides an explicit index legend. Pass exactly this gadget's measurement
record, not earlier measurements from a larger surrounding circuit. If the
gadget is composed into a larger circuit, slice the record appropriately.

All results are relative to the SIGNED Pauli passed in. An identity character
preserves that wire and does not count toward weight. The leftmost character
is data qubit 0. A previously compiled `outcome=-1` changes the starting gates
but does not restrict the operational decoder to that outcome.

## Example formulas for the current XXXX circuit

```
r = b4 xor b5
eigenvalue = (-1)^r

data output 0, physical p0: Z^(b0)
data output 1, physical p4: X^(b6) Z^(b4)
data output 2, physical p2: X^(b2) Z^(b1)
data output 3, physical p5: X^(b2 xor b7) Z^(b5)

consistency check: b2 xor b3 = 0
```

These formulas are for this circuit's exact measurement ordering and wire
mapping. They differ from Gidney's figure because the construction and basis
conventions differ. Derive the formulas again after modifying the circuit.
Correction phases are irrelevant global phases. Different valid corrections
may differ by an output stabilizer.

## Derivation and scope

The implementation keeps one stabilizer generator per qubit, using Stim
PauliStrings for unsigned operators and Python integer masks for affine XOR
signs. Clifford conjugation and Pauli-measurement updates propagate those
signs. Deterministic measurement results produce consistency equations.

A reference-only P^T stabilizer of the circuit's symbolic Choi state gives
the outcome parity; Y transpose signs and signed observables are included.
GF(2) equations then find output-only Pauli corrections making all ideal
postmeasurement Choi stabilizers agree for every allowed record. The decoder
raises an error if this is impossible; it never silently assumes corrections
exist. Ancilla resets are supported only as fresh-wire preparations, matching
the current compiler. Classical adaptation within the raw input circuit is
not supported by this derivation.

With k independent random bits on the Bell input, all 2^k valid Choi records
have probability 2^-k. The nonconstant affine outcome splits them into two
equal sets. For each record the corrected map is alpha times the ideal Pauli
projector with |alpha|^2=2^(1-k). Summing over records therefore establishes
the normalized two-outcome instrument, not just equality of normalized states.
Global phases of individual conditional maps do not affect the instrument.

This is a NOISELESS instrument decoder and correction rule. A nonzero
consistency check signals a record forbidden in the ideal circuit, and by
default raises an error. Setting `validate=False` merely evaluates the XOR
formulas; it does not correct noisy records. This does not implement a
fault-tolerant decoder or prove that physically applied feedback preserves
the paper's distance guarantee under faults. A QEC implementation must use
detectors and an appropriate decoder/Pauli-frame policy.

## Validation

```bash
python test_decoder_10.py
```

Exactly ten circuits at weights 1,2,3,4,5,6,7,9,13,20 are tested. Each gets a
symbolic all-allowed-record check and a normalization check. Independently,
Stim samples 32 Bell-input records and 8 other stabilizer-input records per
circuit, including signed/mixed-Y observables, identity wires, both actual
outcomes and both starting compiled branch choices. Half use actual classical
feedback and half apply the decoded Pauli frame after the raw circuit.
All 400 state checks passed. No exhaustive outcome enumeration is performed.

`decoder_examples/` includes the 10 corrected .stim circuits, formulas,
metadata, and validation results. `certificate` is an algebraic check report
from this implementation, not a proof-assistant certificate or a noise study.
