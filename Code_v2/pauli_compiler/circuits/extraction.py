"""Extract an ordered physical circuit from a reduced ZX graph and MC flow."""
from pyzx.utils import VertexType
from ..zx.flow import validate_flow
from .model import CircuitOperation, ExtractedCircuit

def extract_circuit(graph, paths, *, strict=True):
    """Proposition 4.11/4.13 extraction of a reduced, well-covered diagram.

    Handles endpoints, H, Clifford X/Z phase gates, same-colour identities,
    weight-one/two Pauli projectors, and CNOT (strict=False). Interacting
    spiders must be phaseless; phases belong on degree-two gate vertices.
    Returns a Stim-exportable circuit with exact scalar and wire mapping.
    """
    from collections import defaultdict
    import heapq
    validate_flow(graph, paths, strict=strict)
    g = graph
    paths = tuple(tuple(p) for p in paths)
    if any(g.type(v) in (VertexType.Z, VertexType.X) and g.vertex_degree(v) > 3
           for v in g.vertices()):
        raise ValueError('Call reduce_degrees before extraction')
    owner = {}
    covered = set()
    for qi, p in enumerate(paths):
        for v in p:
            if v in owner:
                raise ValueError('Reduced paths must be vertex-disjoint')
            owner[v] = qi
        covered.update(frozenset((a,b)) for a,b in zip(p,p[1:]))
    groups, vertex_event, sqrt_power = [], {}, 0

    def event(vertices, ops):
        i = len(groups)
        groups.append(tuple(ops))
        for v in vertices:
            vertex_event[v] = i

    def phase_ops(color, phase, q):
        k = phase*2
        if int(k) != k:
            raise ValueError('Extraction currently supports Clifford phases only')
        k = int(k) % 4
        if k == 0:
            return []
        middle = {1:'S', 2:'Z', 3:'S_DAG'}[k]
        ops = [CircuitOperation(middle, (q,))]
        if color == VertexType.X:
            ops = [CircuitOperation('H', (q,))] + ops + [CircuitOperation('H', (q,))]
        return ops

    for v in g.vertices():
        if v in vertex_event:
            continue
        ty, deg = g.type(v), g.vertex_degree(v)
        if ty == VertexType.BOUNDARY:
            event([v], [])
        elif ty == VertexType.H_BOX:
            event([v], [CircuitOperation('H', (owner[v],))])
            sqrt_power += 1
        elif deg == 2:
            event([v], phase_ops(ty, g.phase(v), owner[v]))
        elif deg == 1 and v in owner:
            if g.phase(v) != 0:
                raise ValueError('Endpoint spiders must be phaseless')
            q = owner[v]
            start = paths[q][0] == v
            # Z leaf = sqrt(2)|+>; X leaf = sqrt(2)|0>.
            name = ('RX' if ty == VertexType.Z else 'R') if start else ('MX' if ty == VertexType.Z else 'M')
            event([v], [CircuitOperation(name, (q,))])
            sqrt_power += 1
        elif deg in (1, 3):
            if deg == 1:
                w = next(iter(g.neighbors(v)))
            else:
                uncovered = [w for w in g.neighbors(v) if frozenset((v,w)) not in covered]
                if len(uncovered) != 1:
                    raise ValueError('Expected one uncovered interaction edge')
                w = uncovered[0]
            if g.phase(v) != 0 or g.phase(w) != 0:
                raise ValueError('Interaction spiders must be phaseless')
            if g.vertex_degree(w) not in (1, 3):
                raise ValueError('Unsupported interaction partner')
            if deg == 1 and g.vertex_degree(w) == 1:
                raise ValueError('Disconnected scalar components are unsupported')
            same = ty == g.type(w)
            if deg == 1 or g.vertex_degree(w) == 1:
                s = w if deg == 1 else v
                ops = [] if same else [CircuitOperation('M' if g.type(s)==VertexType.Z else 'MX', (owner[s],))]
                sqrt_power += int(not same)
            elif same:
                ops = [CircuitOperation('MZZ' if ty==VertexType.Z else 'MXX', (owner[v], owner[w]))]
            else:
                control, target = (v,w) if ty==VertexType.Z else (w,v)
                ops = [CircuitOperation('CX', (owner[control], owner[target]))]
                sqrt_power -= 1
            event([v,w], ops)
        else:
            raise ValueError(f'Unsupported vertex {v}')
    # Interactions are atomic events shared by their two worldlines.
    succ = defaultdict(set)
    indegree = [0]*len(groups)
    for p in paths:
        for a,b in zip(p,p[1:]):
            x,y = vertex_event[a], vertex_event[b]
            if x != y and y not in succ[x]:
                succ[x].add(y)
                indegree[y] += 1
    ready = [i for i,d in enumerate(indegree) if d==0]
    heapq.heapify(ready)
    operations, done = [], 0
    while ready:
        i = heapq.heappop(ready)
        done += 1
        operations.extend(groups[i])
        for j in succ[i]:
            indegree[j] -= 1
            if indegree[j] == 0:
                heapq.heappush(ready,j)
    if done != len(groups):
        raise ValueError('Extraction event graph contains a causal cycle')
    return ExtractedCircuit(tuple(operations), len(paths),
                            tuple(owner[v] for v in g.inputs()),
                            tuple(owner[v] for v in g.outputs()),
                            complex(g.scalar.to_number()) * 2**(sqrt_power/2))
