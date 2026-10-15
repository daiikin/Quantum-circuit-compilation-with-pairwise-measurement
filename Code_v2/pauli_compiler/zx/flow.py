"""Validation of directed MC-flow paths and their causal order."""
from copy import deepcopy
from pyzx.utils import VertexType, EdgeType
from pyzx.utils import VertexType as V

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
        for a, b in zip(p, p[1:]):
            edge = frozenset((a, b))
            if not graph.connected(a, b) or edge in covered:
                raise ValueError("Path edge is missing or used more than once (P3)")
            covered.add(edge)
            directed.add((a, b))
    for v in vertices:
        incident = {frozenset((v, w)) for w in graph.neighbors(v)}
        missing = incident - covered
        if graph.type(v) == VertexType.BOUNDARY:
            if graph.vertex_degree(v) != 1 or missing:
                raise ValueError("P1: boundary must be fully covered")
        elif len(missing) > 1:
            raise ValueError("P2: more than one uncovered spider edge")
    if strict:
        for e in graph.edges():
            a, b = graph.edge_st(e)
            if (graph.vertex_degree(a) >= 3 and graph.vertex_degree(b) >= 3
                    and graph.type(a) != graph.type(b)
                    and frozenset((a, b)) not in covered):
                raise ValueError("PS: uncovered edge joins unlike high-degree spiders")
    relation = set(directed) | {(v, v) for v in vertices}
    for x, y in directed:
        for z in graph.neighbors(y):
            if (z, y) not in directed:
                relation.add((x, z))
    # Transitive closure; modest graphs only. Replace with a DAG algorithm
    # when scaling to complete codes.
    for k in vertices:
        predecessors = {a for a, b in relation if b == k}
        successors = {b for a, b in relation if a == k}
        relation.update((a, b) for a in predecessors for b in successors)
    if any(a != b and (b, a) in relation for a, b in relation):
        raise ValueError("O1/O2 constraints admit no partial order")
    return frozenset(relation)


def validate_flow(graph, paths, *, strict=True):
    """Validate definitions 4.6/4.8/4.9 for X/Z/H-box/boundary graphs.

    H-boxes here are only degree-two phase-1 Hadamard boxes. All edges
    are ordinary. H-box P1 requires both edges on the SAME path.
    A proxy vertex type is used only to reuse O1/O2/order computation;
    the returned graph and its tensor are never changed by validation.
    """
    paths = tuple(tuple(p) for p in paths)
    proxy = deepcopy(graph)
    for v in graph.vertices():
        if graph.type(v) == V.H_BOX:
            if graph.vertex_degree(v) != 2 or graph.phase(v) != 1:
                raise ValueError('Only degree-two phase-1 Hadamard boxes supported')
            if not any(v in p[1:-1] and set((p[p.index(v)-1], p[p.index(v)+1]))
                       == set(graph.neighbors(v)) for p in paths):
                raise ValueError('P1: Hadamard box must be fully covered by one path')
            proxy.set_type(v, V.Z)
        elif graph.type(v) not in (V.Z, V.X, V.BOUNDARY):
            raise ValueError('Unsupported node type')
    return mc_order(proxy, paths, strict=strict)
