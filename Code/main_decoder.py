"""Run from the folder containing the compiler modules."""
import json
from pathlib import Path
import stim
from measurement_decoder import compile_operational_measurement

decoder = compile_operational_measurement('XXXX')
print(decoder.formulas())
print('\nMeasurement index legend:')
for entry in decoder.measurement_operations:
    print(f"b{entry['index']}: {entry['name']} on {entry['physical_wires']}")

physical_circuit = decoder.corrected_stim()
print('\nCircuit including feedback:\n', physical_circuit)

# Demonstration: data inputs begin in |0...0>. Ancilla preparations occur
# inside the circuit. No intermediate result is forced or discarded.
simulator = stim.TableauSimulator(seed=42)
simulator.do_circuit(physical_circuit)
record = simulator.current_measurement_record()
print('Raw record:', [int(b) for b in record])
print('Measured eigenvalue:', decoder.outcome_eigenvalue(record))
print('Keep these physical output wires:', decoder.compilation.circuit.outputs)

Path('operational_measurement.stim').write_text(str(physical_circuit))
Path('operational_measurement.json').write_text(json.dumps(decoder.metadata(), indent=2))
