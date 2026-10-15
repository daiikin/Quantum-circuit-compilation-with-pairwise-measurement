# Physical circuit model and ZX extraction

This stage converts a reduced graph and directed MC-flow paths into chronological
physical operations. Outcome processing and corrections belong to `decoding/`.

## Files and entry points

- `model.py`: `CircuitOperation(name, targets)` and `ExtractedCircuit`.
- `extraction.py`: `extract_circuit(graph, paths, strict=True)`.
- `simplification.py`: `remove_redundant_measurements(circuit)` removes repeated
  observables through commuting measurements or disjoint operations. Any other
  operation touching either wire invalidates the cached result. Different partners
  are different observables. This is conservative and noiseless only.

`compile_measurement` enables simplification by default, before Paulimer.
Its `.raw_circuit` retains the literal extraction; `.circuit` is simplified.
Measurement indices are regenerated, so never reuse formulas from before this
pass. Boundary mappings and branch_scale are preserved.

```python
from pauli_compiler import compile_measurement

compilation = compile_measurement('XYZXYZ')
raw = compilation.circuit
for operation in raw.operations:
    print(operation.name, operation.targets)
print('Physical input wires:', raw.inputs)
print('Physical output wires:', raw.outputs)
print(raw.to_stim())
```

The extractor supports preparation/measurement endpoints, Clifford degree-two
phase gates, Hadamards, identities, pair projectors, and CNOT patterns (when the
appropriate non-strict flow is allowed). Input diagrams must already be reduced;
it rejects high-degree spiders, overlapping reduced paths and causal cycles.

`inputs[j]` is the physical wire initially carrying logical data qubit j.
`outputs[j]` is its physical output location. Other wires are ancillas; input
and output wire lists need not be identical because the extraction can move data.
`num_qubits` counts physical wires with no reuse optimization.

## Raw circuit semantics

For the selected initial ZX branch:

```
ZX map = branch_scale * all-zero postselected physical circuit map
```

`branch_scale` records scalar normalization. `metadata()` explicitly lists the
all-zero reference record. `to_stim()` exports the operations, but does not itself
force any measurement outcome or add feedback. `MXX` and `MZZ` export as joint
Stim `MPP` measurements, not separate single-qubit measurements.

To obtain an operational nondestructive measurement for all valid records:

```python
from pauli_compiler.decoding.paulimer_backend import derive_decoder

decoder = derive_decoder(compilation)
physical_circuit = decoder.corrected_stim()
```

The full pipeline currently targets positive-weight Pauli measurements with equal
numbers of data inputs and outputs. Extraction is not a universal ZX-to-unitary
synthesis algorithm.
