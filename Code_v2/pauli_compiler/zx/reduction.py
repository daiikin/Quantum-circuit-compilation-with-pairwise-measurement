"""Reduce high-degree spiders while maintaining their directed MC-flow paths."""
from copy import deepcopy
from dataclasses import dataclass
from pyzx.utils import VertexType
from .flow import validate_flow
from .rewrites import r4, r5

@dataclass
class ReducedDiagram:
    graph: object
    paths: tuple
    order: frozenset
    trace: tuple


def reduce_degrees(graph, paths, *, strict=True, check_each=False):
    """Equation (11): reduce all zero-phase high-degree spiders to degree <=3.

    Implements r4, r5, r_2m, r_(2m+1), and the two-leaf expansion for
    degrees 2/3 modulo four. Requires an existing well-covered MC flow.
    The input is untouched. check_each enables costly intermediate checks.
    """
    g = deepcopy(graph)
    paths = [list(p) for p in paths]
    validate_flow(g, paths, strict=strict)
    trace = []
    while True:
        high = [v for v in g.vertices() if g.type(v) in (VertexType.Z, VertexType.X)
                and g.vertex_degree(v) > 3]
        if not high:
            break
        s = high[0]
        n = g.vertex_degree(s)
        if g.phase(s) != 0:
            raise ValueError('High-degree rules require phase-zero spiders')
        color, row, qubit = g.type(s), g.row(s), g.qubit(s)
        traversals = []
        for pi, p in enumerate(paths):
            if s in p:
                j = p.index(s)
                if not 0 < j < len(p)-1:
                    raise ValueError('Expected well-covered flow')
                traversals.append((pi, j, p[j-1], p[j+1]))
        used = {v for _, _, a, b in traversals for v in (a, b)}
        extra = set(g.neighbors(s)) - used
        if len(traversals) != n//2 or len(extra) != n % 2:
            raise ValueError('Unsupported local path coverage')
        if n in (4, 5):
            t, b = traversals
            ports = [t[2], t[3], b[3]]
            if n == 5:
                ports += [next(iter(extra))]
            ports += [b[2]]
            result = (r4 if n == 4 else r5)(g, s, ports)
            ring = result.ring
            g = result.graph
            routes = [(ring[0], ring[1]),
                      (ring[3], ring[2]) if n == 4 else (ring[4], ring[3], ring[2])]
            for (pi, j, _, _), route in zip(traversals, routes):
                paths[pi][j:j+1] = route
            trace.append({'rule': f'r{n}', 'degree': n})
        elif n % 4 in (2, 3):
            # rfuse twice: introduce a complete new worldline through s.
            a = g.add_vertex(ty=color, qubit=qubit-1, row=row-1)
            b = g.add_vertex(ty=color, qubit=qubit-1, row=row+1)
            g.add_edge(g.edge(a, s))
            g.add_edge(g.edge(s, b))
            paths.append([a, s, b])
            trace.append({'rule': 'two_leaf_expansion', 'degree': n})
        else:
            # m even: m/2 incoming and m/2 outgoing degree-four spokes,
            # each attached to both hubs. Pair paths through alternating hubs.
            m = n//2
            hubs = [g.add_vertex(ty=color, qubit=qubit+h-.5, row=row) for h in range(2)]
            incoming = [g.add_vertex(ty=color, qubit=qubit+i, row=row-1) for i in range(m//2)]
            outgoing = [g.add_vertex(ty=color, qubit=qubit+i, row=row+1) for i in range(m//2)]
            for spoke in incoming + outgoing:
                for hub in hubs:
                    g.add_edge(g.edge(spoke, hub))
            for k, (pi, j, a, b) in enumerate(traversals):
                u, v = incoming[k//2], outgoing[k//2]
                g.add_edge(g.edge(a, u))
                g.add_edge(g.edge(v, b))
                paths[pi][j:j+1] = [u, hubs[k % 2], v]
            if extra:
                g.add_edge(g.edge(hubs[1], next(iter(extra))))
            g.remove_vertex(s)
            trace.append({'rule': 'odd_split' if extra else 'even_split', 'degree': n})
        if check_each:
            validate_flow(g, paths, strict=strict)
    paths = tuple(tuple(p) for p in paths)
    return ReducedDiagram(g, paths, validate_flow(g, paths, strict=strict), tuple(trace))
