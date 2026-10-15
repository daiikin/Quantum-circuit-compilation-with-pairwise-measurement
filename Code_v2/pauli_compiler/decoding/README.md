# Global outcomes and output corrections with Paulimer

This stage makes the extracted measurement circuit operational: it computes the
XOR of recorded bits that yields the target Pauli result and the output Paulis
that restore the desired post-measurement state.

## Files

| File | Responsibility |
|---|---|
| `affine.py` | `Affine`: store, display and evaluate XOR expressions |
| `result.py` | `MeasurementDecoder`: user-facing formulas, record checks, JSON and Stim export |
| `paulimer_backend.py` | Translate operations into Paulimer, extract outcome signs, solve corrections using Binar, verify all branches |

Paulimer performs stabilizer simulation; Binar performs binary linear solves.
The adapter wires those APIs to the extracted circuit and intended target Pauli.
There is no fallback to the earlier in-house simulator. Stim is used for export
and independently by tests, not to derive corrections.

## Use the decoder

```python
import stim
from pauli_compiler import compile_operational_measurement

decoder = compile_operational_measurement('XXXX')
print(decoder.formulas())

sim = stim.TableauSimulator(seed=42)
sim.do_circuit(decoder.corrected_stim())
record = sim.current_measurement_record()
print('Global outcome:', decoder.outcome_eigenvalue(record))
print('Keep physical outputs:', decoder.compilation.circuit.outputs)
```

This demonstration initializes the data in |0...0>; the symbolic verification
covers arbitrary inputs, including inputs entangled with a reference.

To track a Pauli frame instead of physically applying feedback, execute the raw
circuit and evaluate `decoder.output_pauli_frame(record)`. Do not also apply the
exported physical feedback. The returned Pauli string is ordered by logical
outputs, so map its positions through `decoder.compilation.circuit.outputs`.

## Formulas and records

- `outcome`: one affine XOR formula; `outcome_bit(record)` evaluates it.
- `x_corrections[j]`, `z_corrections[j]`: X/Z exponents for logical output j.
- `measurement_operations`: mapping from b_i to operation name and physical wires.
- `checks`: affine formulas required to equal zero on noiseless records.
- `consistency_checks(record)`: evaluates the checks; nonzero means violation.
- `metadata()`: JSON-serializable formulas, wire mappings and verification details.
- `certificate`: backend versions and symbolic verification facts, not a formal
  proof-assistant certificate or a noise-correction guarantee.

`Affine.mask` uses bit 0 for the constant and bit i+1 for b_i. For example, a mask
of 11 (binary 1011) represents `1 xor b0 xor b2`. Bit 0 of a physical measurement
record means the measured eigenvalue was +1.

Records include every single-qubit and pair measurement in extraction order.
The decoded parity is not necessarily an XOR of pair measurements alone.
The `outcome_bit`/`outcome_eigenvalue` methods validate noiseless constraints by
default; bypassing validation is not a substitute for noise decoding.

## How the backend works

It runs the circuit on halves of mathematical Bell pairs using
`OutcomeCompleteSimulation`. Paulimer represents all allowed records and final
stabilizer signs as binary linear functions of independent random bits. The sign
of P-transpose on the untouched references gives the decoded global outcome.
The transpose's Y signs and the target's explicit +/- sign are included.

A second Paulimer simulation constructs the ideal measurement Choi stabilizers.
Binar solves for output-only X/Z corrections that align their signs. Paulimer
then applies those corrections and checks all ideal generators, including the
reference/output correlations. The normalization check uses the branch counts
and Choi probabilities. See [the derivation](../../docs/DECODING.md).

Fresh R/RX preparations and the compiler's supported Clifford/Pauli operations
are accepted. Used-wire resets and unsupported operations are rejected. So are
circuits that measure extra information or require non-Pauli output corrections.
The input sequence is nonadaptive; generated feedback is applied at the end.
Different valid correction formulas may differ by stabilizers or record
constraints; the solver does not optimize correction depth.
