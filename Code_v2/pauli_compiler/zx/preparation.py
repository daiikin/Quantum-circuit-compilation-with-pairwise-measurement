"""Prepare an original measurement star for MC-flow extraction (Proposition 5.1)."""
from copy import deepcopy
from dataclasses import dataclass
from pyzx.utils import VertexType as V, EdgeType as E
from .flow import validate_flow
from .rewrites import r4

@dataclass
class PreparedMeasurement:
    measurement: object
    graph: object
    center: int
    paths: tuple
    order: frozenset
    gadgets: tuple
    odd_leaf: int | None
    central_path_indices: tuple


def _wire(g, start, end, center):
    """Unique data chain in the generator's original star diagram."""
    path, previous, current = [start], None, start
    while current != end:
        candidates = [v for v in g.neighbors(current) if v not in (previous, center)]
        if len(candidates) != 1:
            raise ValueError('Expected the original measurement-generator topology')
        previous, current = current, candidates[0]
        if current in path:
            raise ValueError('Cycle on a data wire')
        path.append(current)
    return path


def _expand_hadamards(g, paths):
    """Replace normalized H edges by H boxes; correct unnormalized box scalar."""
    for e in list(g.edges()):
        if g.edge_type(e) != E.HADAMARD:
            continue
        a, b = g.edge_st(e)
        h = g.add_vertex(ty=V.H_BOX, phase=1,
                         qubit=(g.qubit(a)+g.qubit(b))/2,
                         row=(g.row(a)+g.row(b))/2)
        g.remove_edge(e)
        g.add_edge(g.edge(a, h))
        g.add_edge(g.edge(h, b))
        # PyZX H-box tensor [[1,1],[1,-1]] lacks 1/sqrt(2).
        g.scalar.add_power(-1)
        found = 0
        for p in paths:
            for i in range(len(p)-1):
                if {p[i], p[i+1]} == {a, b}:
                    p.insert(i+1, h)
                    found += 1
                    break
        if found != 1:
            raise ValueError('Every Hadamard edge must belong to exactly one path')


def prepare_mc_flow(measurement):
    """Prepare ORIGINAL star diagram using Proposition 5.1.

    Even-indexed support wires: add Z gadget/leaf AFTER the data tap.
    Odd-indexed support wires: add Z gadget/leaf BEFORE the data tap.
    Pair an even wire's input with the following odd wire's output through
    the central X spider. Remaining input/output segments end/start at Z
    leaves. For odd weight, route the last input through the center to a
    new one-legged X spider. Indices refer to support order, not all wires.
    """
    g = deepcopy(measurement.original)
    center = measurement.center
    support = measurement.support
    if not support or g.vertex_degree(center) != len(support):
        raise ValueError('Expected an original positive-weight measurement star')
    wires = [_wire(g, a, b, center) for a, b in zip(g.inputs(), g.outputs())]
    paths = [wires[q][:] for q in range(len(wires)) if q not in support]
    entering, leaving, gadgets = {}, {}, []
    for j, (q, tap) in enumerate(zip(support, measurement.ports)):
        wire = wires[q]
        k = wire.index(tap)
        after = j % 2 == 0
        a, b = (tap, wire[k+1]) if after else (wire[k-1], tap)
        edge_type = g.edge_type(g.edge(a, b))
        gadget = g.add_vertex(ty=V.Z, phase=0, qubit=q, row=4.6 if after else 3.4)
        leaf = g.add_vertex(ty=V.Z, phase=0, qubit=q+.35, row=g.row(gadget))
        g.remove_edge(g.edge(a, b))
        g.add_edge(g.edge(a, gadget), E.SIMPLE if after else edge_type)
        g.add_edge(g.edge(gadget, b), edge_type if after else E.SIMPLE)
        g.add_edge(g.edge(gadget, leaf))
        gadgets.append((q, gadget, leaf, 'after' if after else 'before'))
        if after:
            entering[j] = wire[:k+1]
            paths.append([leaf, gadget] + wire[k+1:])
        else:
            paths.append(wire[:k] + [gadget, leaf])
            leaving[j] = wire[k:]
    central_indices = []
    for j in range(0, len(support)-1, 2):
        central_indices.append(len(paths))
        paths.append(entering[j] + [center] + leaving[j+1])
    odd_leaf = None
    if len(support) % 2:
        odd_leaf = g.add_vertex(ty=V.X, phase=0, qubit=-3, row=5)
        g.add_edge(g.edge(center, odd_leaf))
        # Spider fusion is exact: the leaf's sqrt(2) is cancelled by
        # the central X spider's changed arity. No scalar adjustment.
        central_indices.append(len(paths))
        paths.append(entering[len(support)-1] + [center, odd_leaf])
    _expand_hadamards(g, paths)
    paths = tuple(tuple(p) for p in paths)
    order = validate_flow(g, paths)
    return PreparedMeasurement(measurement, g, center, paths, order,
                               tuple(gadgets), odd_leaf, tuple(central_indices))


def apply_prepared_r4(prepared):
    """For weight four only: r4 with flow-aware square routing."""
    if prepared.graph.vertex_degree(prepared.center) != 4:
        raise ValueError('Prepared weight five has degree SIX; r5 is not applicable')
    paths = prepared.paths
    ids = prepared.central_path_indices
    if len(ids) != 2:
        raise ValueError('Expected two paths through the center')
    a, b = (paths[i] for i in ids)
    ia, ib = a.index(prepared.center), b.index(prepared.center)
    ports = (a[ia-1], a[ia+1], b[ib+1], b[ib-1])
    result = r4(prepared.graph, prepared.center, ports)
    v = result.ring
    updated = list(paths)
    updated[ids[0]] = a[:ia] + (v[0], v[1]) + a[ia+1:]
    updated[ids[1]] = b[:ib] + (v[3], v[2]) + b[ib+1:]
    result.paths = tuple(updated)
    result.order = validate_flow(result.graph, result.paths)
    return result
