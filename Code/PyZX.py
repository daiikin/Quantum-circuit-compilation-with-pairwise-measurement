import pyzx as zx
from pauli_measurements import pauli_measurement, measurement_branches

# Generate a measurement branch and apply r4 or r5.
m4 = pauli_measurement("XYZX", outcome=+1)
m5 = pauli_measurement("XYZXY", outcome=-1)

zx.draw(m4.original)
zx.draw(m4.rewritten)

# Both outcomes; signs and identity positions are supported.
branches = measurement_branches("-XIYZZ")