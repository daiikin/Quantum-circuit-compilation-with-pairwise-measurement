"""Compose ZX construction, MC-flow preparation, reduction, and extraction."""
from dataclasses import dataclass
from .zx.spiders import pauli_measurement, generate_pauli_measurement
from .zx.preparation import prepare_mc_flow
from .zx.reduction import ReducedDiagram, reduce_degrees
from .circuits.extraction import extract_circuit
from .circuits.model import ExtractedCircuit
from .circuits.simplification import remove_redundant_measurements

@dataclass
class Compilation:
    measurement: object
    prepared: object
    reduced: ReducedDiagram
    circuit: ExtractedCircuit
    raw_circuit: ExtractedCircuit | None = None


def compile_measurement(pauli, outcome=1, *, check_each=False,
                        simplify_measurements=True):
    """Arbitrary positive-weight Pauli -> prepared ZX -> reduced ZX -> circuit.

    pauli may be a signed string or a Measurement object. For an object,
    its stored outcome is used. Returned circuit is explicitly postselected.
    Repeated measurements are removed before decoding by default; disable with
    simplify_measurements=False to inspect the literal extraction. ZX diagrams
    retain their original vertices; only the extracted operation list changes.
    """
    m = pauli_measurement(pauli, outcome) if isinstance(pauli, str) else pauli
    prepared = prepare_mc_flow(m)
    reduced = reduce_degrees(prepared.graph, prepared.paths, check_each=check_each)
    circuit = extract_circuit(reduced.graph, reduced.paths)
    raw_circuit = circuit
    if simplify_measurements:
        circuit = remove_redundant_measurements(circuit)
    return Compilation(m, prepared, reduced, circuit, raw_circuit)


def compile_operational_measurement(pauli, outcome=1, *, check_each=False,
                                    simplify_measurements=True):
    from .decoding.paulimer_backend import derive_decoder
    return derive_decoder(compile_measurement(pauli, outcome, check_each=check_each,
                                               simplify_measurements=simplify_measurements))
