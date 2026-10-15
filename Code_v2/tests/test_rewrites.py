"""Local R4/R5 tensor equality and path-routing regressions."""
import unittest
from itertools import permutations
import numpy as np
import pyzx as zx
from pyzx.utils import VertexType, EdgeType
from pauli_compiler.zx.rewrites import r4, r5, rewrite_with_flow

def example(n=4, color=VertexType.Z):
    """Open test diagram with two input/output pairs and an optional leaf."""
    if n not in (4, 5):
        raise ValueError("n must be 4 or 5")
    g = zx.Graph()
    v = g.add_vertex(ty=color)
    b = [g.add_vertex(ty=VertexType.BOUNDARY) for _ in range(4)]
    # b order: top input, top output, bottom output, bottom input.
    for p in b:
        g.add_edge(g.edge(v, p))
    g.set_inputs((b[0], b[3]))
    g.set_outputs((b[1], b[2]))
    extra = None
    if n == 5:
        extra = g.add_vertex(ty=color)
        g.add_edge(g.edge(v, extra))
    paths = ((b[0], v, b[1]), (b[3], v, b[2]))
    return g, v, b, extra, paths


class Tests(unittest.TestCase):
    def test_open_tensors_all_port_orders_and_colors(self):
        # Independent tensor contraction, including all 24/120 port orders.
        for n in (4, 5):
            for color in (VertexType.Z, VertexType.X):
                g = zx.Graph()
                v = g.add_vertex(ty=color)
                ports = tuple(g.add_vertex(ty=VertexType.BOUNDARY) for _ in range(n))
                for p in ports:
                    g.add_edge(g.edge(v, p))
                g.set_outputs(ports)
                expected = zx.tensor.tensorfy(g)
                for order in permutations(ports):
                    out = (r4 if n == 4 else r5)(g, v, order)
                    np.testing.assert_allclose(zx.tensor.tensorfy(out.graph), expected, atol=1e-12)
                    self.assertTrue(all(out.graph.vertex_degree(w) == 3 for w in out.ring))
                self.assertIn(v, g.vertices())

    def test_flow(self):
        for n in (4, 5):
            for color in (VertexType.Z, VertexType.X):
                g, v, b, extra, paths = example(n, color)
                out = rewrite_with_flow(g, v, (b[0], b[1]), (b[3], b[2]), paths, extra=extra)
                self.assertNotIn(v, out.graph.vertices())
                self.assertIn(v, g.vertices())
                self.assertEqual(len(out.paths), 2)
                np.testing.assert_allclose(zx.tensor.tensorfy(out.graph), zx.tensor.tensorfy(g), atol=1e-12)

    def test_rejections_leave_input_unchanged(self):
        g, v, b, _, paths = example()
        with self.assertRaises(ValueError):
            r5(g, v, b)
        with self.assertRaises(ValueError):
            r4(g, v, [b[0]]*4)
        with self.assertRaises(ValueError):
            rewrite_with_flow(g, v, (b[0], b[2]), (b[3], b[1]), paths)
        g.set_phase(v, 1)
        with self.assertRaises(ValueError):
            r4(g, v, b)
        g.set_phase(v, 0)
        g.set_edge_type(g.edge(v, b[0]), EdgeType.HADAMARD)
        with self.assertRaises(ValueError):
            r4(g, v, b)
        self.assertEqual(g.num_vertices(), 5)
