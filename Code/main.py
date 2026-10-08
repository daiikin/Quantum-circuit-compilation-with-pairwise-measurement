import matplotlib.pyplot as plt
from zx_r4_r5 import generate_pauli_measurement, compile_measurement
from flow_layout import draw_flow

m = generate_pauli_measurement(5, seed=42)
result = compile_measurement(m)

fig = draw_flow(
    result.reduced.graph,
    result.reduced.paths,
    labels=False,
)

fig.savefig("zx_reduced_flow.svg", bbox_inches="tight")
plt.show()