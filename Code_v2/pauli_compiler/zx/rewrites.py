"""Local R4/R5 spider rewrites, with optional directed path updates."""
from copy import deepcopy
from dataclasses import dataclass
import numpy as np
from pyzx.utils import VertexType, EdgeType
from .flow import mc_order

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


def rewrite_with_flow(graph, spider, top, bottom, paths, *, extra=None, strict=True):
    """Equation (11) routing, left-to-right.

    top=(top_in, top_out), bottom=(bottom_in, bottom_out) are neighbours
    of the old spider. For r5, extra is the fifth, uncovered neighbour.
    Paths must traverse top_in -> spider -> top_out and
    bottom_in -> spider -> bottom_out. The fifth edge remains uncovered.
    r4 cycle: top_in, top_out, bottom_out, bottom_in.
    r5 cycle: top_in, top_out, bottom_out, extra, bottom_in.
    Validates old/new flow, reconstructs order, and leaves inputs unchanged.
    """
    paths = tuple(tuple(p) for p in paths)
    mc_order(graph, paths, strict=strict)
    ti, to = top
    bi, bo = bottom
    ports = (ti, to, bo, bi) if extra is None else (ti, to, bo, extra, bi)
    result = (r4 if extra is None else r5)(graph, spider, ports)
    ring = result.ring
    routes = {(ti, to): (ring[0], ring[1]),
              (bi, bo): ((ring[3], ring[2]) if extra is None
                         else (ring[4], ring[3], ring[2]))}
    used, updated = set(), []
    for path in paths:
        if spider not in path:
            updated.append(path)
            continue
        i = path.index(spider)
        if i == 0 or i == len(path)-1:
            raise ValueError("The matched spider must be internal to each path")
        pair = (path[i-1], path[i+1])
        if pair not in routes or pair in used:
            raise ValueError("Paths do not match the supplied two directed port pairs")
        used.add(pair)
        updated.append(path[:i] + routes[pair] + path[i+1:])
    if used != set(routes):
        raise ValueError("Both directed port pairs must occur in the paths")
    result.paths = tuple(updated)
    result.order = mc_order(result.graph, result.paths, strict=strict)
    return result
