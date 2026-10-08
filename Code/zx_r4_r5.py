
    
  
"""Distance-preserving r4/r5 spider expansions in PyZX.

Source: Rodatz, Poor, Kissinger, Quantum 2026-09-03-2202,
Theorems 3.8 (r5), 3.9 (r4), and Equation (11).

Install: python -m pip install pyzx numpy
Test:    python zx_r4_r5.py

Only zero-phase X/Z spiders with distinct neighbours and ordinary edges
are accepted. Functions return a deep-copied graph: input IDs are retained,
except for the removed spider. No generic simplification is performed.
Optional flow support uses directed vertex paths and reconstructs a valid
partial order from O1/O2; it does not retain additional scheduling constraints.
This is a local rewrite module, not a measurement-outcome or circuit compiler.
"""

from copy import deepcopy
from dataclasses import dataclass
from itertools import permutations
import unittest

import numpy as np
import pyzx as zx
from pyzx.utils import VertexType, EdgeType


@dataclass
class RewriteResult:
    graph: object
    ring: tuple
    ports: tuple
    paths: tuple | None = None
    order: frozenset | None = None


def _expand(graph, spider, ports, size):
    ports = tuple(ports)
    if spider not in graph.vertices():
        raise ValueError("Spider does not exist")
    if graph.type(spider) not in (VertexType.Z, VertexType.X):
        raise ValueError("Expected an X or Z spider")
    if graph.phase(spider) != 0:
        raise ValueError("These templates require phase zero")
    if len(ports) != size or len(set(ports)) != size:
        raise ValueError(f"Supply {size} distinct ports in cyclic order")
    if set(ports) != set(graph.neighbors(spider)) or graph.vertex_degree(spider) != size:
        raise ValueError("Ports must be exactly the spider's neighbours")
    if any(graph.edge_type(graph.edge(spider, p)) != EdgeType.SIMPLE for p in ports):
        raise ValueError("Incident edges must be ordinary (not Hadamard) edges")
    if spider in tuple(graph.inputs()) + tuple(graph.outputs()):
        raise ValueError("Cannot remove a designated input/output")
    g = deepcopy(graph)
    color = g.type(spider)
    q, row = g.qubit(spider), g.row(spider)
    ring = tuple(g.add_vertex(ty=color, phase=0,
                             qubit=q + np.sin(2*np.pi*i/size),
                             row=row + np.cos(2*np.pi*i/size))
                 for i in range(size))
    g.remove_vertex(spider)
    for i, p in enumerate(ports):
        g.add_edge(g.edge(p, ring[i]), EdgeType.SIMPLE)
        g.add_edge(g.edge(ring[i], ring[(i+1) % size]), EdgeType.SIMPLE)
    return RewriteResult(g, ring, ports)


def r4(graph, spider, ports):
    """Replace a degree-4 spider by a square; ports are in cyclic order."""
    return _expand(graph, spider, ports, 4)


def r5(graph, spider, ports):
    """Replace a degree-5 spider by a pentagon; ports are in cyclic order."""
    return _expand(graph, spider, ports, 5)


def mc_order(graph, paths, *, strict=True):
    """Validate P1/P2/P3/PWC/(PS), construct and validate O1/O2 order.

    Supports ordinary-edge X/Z/boundary graphs. Hadamard edges/boxes must
    be handled by a future extension. Returns reflexive transitive closure.
    Reconstructing order may change a caller's previous schedule.
    """
    paths = tuple(tuple(p) for p in paths)
    vertices = set(graph.vertices())
    if any(graph.type(v) not in (VertexType.X, VertexType.Z, VertexType.BOUNDARY)
           for v in vertices):
        raise ValueError("Flow checker supports X/Z/boundary vertices only")
    if any(graph.edge_type(e) != EdgeType.SIMPLE for e in graph.edges()):
        raise ValueError("Flow checker requires ordinary edges")
    covered, directed = set(), set()
    for p in paths:
        if len(p) < 2 or len(p) != len(set(p)) or not set(p) <= vertices:
            raise ValueError("Paths must be simple and have at least one edge")
        for v in (p[0], p[-1]):
            if graph.type(v) in (VertexType.X, VertexType.Z) and graph.vertex_degree(v) > 1:
                raise ValueError("PWC: path ends on a multi-leg spider")

