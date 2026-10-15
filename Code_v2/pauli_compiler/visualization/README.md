# Diagram and circuit rendering

Rendering consumes the compiled objects. It does not generate corrections,
change the measurement protocol, or optimize the circuit.

- `stages.py`: `render_stages(decoder, ...)` creates the complete figure dictionary
  used by run_pipeline, without writing files.
- `flow.py`: `layout_flow` copies a diagram and assigns coordinates along its
  directed paths; `draw_flow` returns a Matplotlib figure.
- `circuit.py`: `draw_measurement_circuit` returns a figure with joint measurement
  boxes, record labels, output corrections, and global-outcome formulas.

## Draw each stage

```python
import matplotlib.pyplot as plt
import pyzx as zx
from pauli_compiler import compile_operational_measurement
from pauli_compiler.visualization.flow import draw_flow
from pauli_compiler.visualization.circuit import draw_measurement_circuit

decoder = compile_operational_measurement('XYZXYZ')
c = decoder.compilation

fig = zx.draw_matplotlib(c.measurement.original, labels=True)
fig.savefig('initial.svg', bbox_inches='tight')
plt.close(fig)

fig = zx.draw_matplotlib(c.prepared.graph, labels=True)
fig.savefig('prepared.svg', bbox_inches='tight')
plt.close(fig)

fig = draw_flow(c.reduced.graph, c.reduced.paths, labels=False)
fig.savefig('reduced.svg', bbox_inches='tight')
plt.close(fig)

fig = draw_measurement_circuit(decoder, show_checks=True)
fig.savefig('circuit.svg', bbox_inches='tight')
plt.close(fig)
```

SVG is preferred for large circuits because it can be zoomed without losing text.
`main.py` runs this sequence and saves the files under its output directory.

## Reading the diagrams

In the flow layout, arrows follow physical worldlines; dashed edges are uncovered
ZX connections. A crossing is not a graph vertex unless a node is drawn there.

In the circuit rendering:

- `p0`, `p1`, ... are physical wire IDs.
- `D0 in/out`, etc. label logical data boundaries; the positions may differ.
- Connected X/X or Z/Z boxes denote one joint pair measurement with one record bit.
- Dashed X/Z boxes are classical feedback corrections applied to the outputs.
- `discard` marks a physical wire that is not a data output.
- The footer gives the decoded outcome, longer correction formulas, and optional
  noiseless record checks. `show_record_legend=True` adds the full record mapping.

The circuit renderer supports the operations emitted by this compiler. The
Paulimer adapter accepts some additional Clifford/measurement operations whose
special symbols are not implemented here. Rendering positions are chosen for
readability; they should not be read as the resource estimator's depth schedule.
