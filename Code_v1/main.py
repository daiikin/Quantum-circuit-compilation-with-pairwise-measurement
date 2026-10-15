import matplotlib.pyplot as plt
import pyzx as zx

from zx_r4_r5 import generate_pauli_measurement, compile_measurement
from measurement_decoder import derive_decoder
from flow_layout import draw_flow
from circuit_renderer import draw_measurement_circuit
from circuit_resources import benchmark_resources


def save_zx(graph, title, filename):
    fig = zx.draw_matplotlib(graph, labels=True)
    fig.suptitle(title)
    fig.savefig(filename, bbox_inches="tight")
    return fig


# Random weight-4 Pauli. Change seed for another reproducible example.
measurement = generate_pauli_measurement(6, seed=54, outcome=+1)
signed_pauli = ("-" if measurement.sign < 0 else "+") + measurement.pauli
print("Measured Pauli:", signed_pauli)

# 1. Initial spider representation.
save_zx(
    measurement.original,
    f"Initial ZX representation: {signed_pauli}",
    "01_initial.svg",
)

# Optional direct R4/R5 visualization.
if measurement.rewritten is not None:
    save_zx(
        measurement.rewritten,
        "Direct R4/R5 rewrite — before MC-flow preparation",
        "02_direct_rewrite.svg",
    )
else:
    print(
        "No direct R4/R5 preview for this weight. "
        "The general rewritten graph is rendered in step 4."
    )
result = compile_measurement(measurement)

# 3. Prepared graph, before degree reduction.
save_zx(
    result.prepared.graph,
    "After MC-flow preparation",
    "03_mc_prepared.svg",
)

# 4. Reduced graph with its directed MC-flow paths laid out horizontally.
fig = draw_flow(
    result.reduced.graph,
    result.reduced.paths,
    labels=False,
    title="MC-flow preparation + rewrites + flow-based layout",
)
fig.savefig("04_reduced_flow_layout.svg", bbox_inches="tight")

# 5. Derive outcome decoding and corrections, then render the circuit.
decoder = derive_decoder(result)
print("\nOutcome and correction formulas:")
print(decoder.formulas())

fig = draw_measurement_circuit(
    decoder,
    title=f"Measurement of {signed_pauli} with classical feedback",
    show_checks=True,
)
fig.savefig("05_circuit_with_corrections.svg", bbox_inches="tight")




_, summary = benchmark_resources(
    weights=range(4, 11),
    samples_per_weight=10,
    seed=20261009,
)

weights = [row["weight"] for row in summary]

fig, (ax_depth, ax_ancillas) = plt.subplots(
    1, 2, figsize=(12, 4.5)
)

for key, label in [
    ("depth_pauli_frame", "Pauli-frame tracking"),
    ("depth_with_feedback", "Physical feedback"),
]:
    means = [row[key]["mean"] for row in summary]
    minima = [row[key]["min"] for row in summary]
    maxima = [row[key]["max"] for row in summary]

    line, = ax_depth.plot(weights, means, "o-", label=label)
    ax_depth.fill_between(
        weights, minima, maxima,
        color=line.get_color(), alpha=0.15,
    )

ax_depth.set_ylabel("Depth (unit-time steps)")
ax_depth.set_title("Mean depth; shading shows min–max")
ax_depth.legend()

ancillas = [row["ancilla_qubits"]["mean"] for row in summary]
ax_ancillas.plot(weights, ancillas, "o-", color="tab:green")
ax_ancillas.set_ylabel("Ancilla qubits")
ax_ancillas.set_title("Ancilla count — no qubit reuse")

for ax in (ax_depth, ax_ancillas):
    ax.set_xlabel("Pauli weight")
    ax.set_xticks(weights)
    ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig("depth_and_ancillas.svg", bbox_inches="tight")
plt.show()

plt.show()