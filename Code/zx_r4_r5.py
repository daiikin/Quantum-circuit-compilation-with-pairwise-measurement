"""Distance-preserving ZX degree reduction and measurement-circuit extraction.

Source: Rodatz, Poor, Kissinger, Quantum 2026-09-03-2202,
Theorems 3.8 (r5), 3.9 (r4), and Equation (11).

Install: python -m pip install pyzx numpy
Test:    python zx_r4_r5.py

High-degree rules accept zero-phase X/Z spiders with ordinary edges.
Inputs are not mutated. General APIs: reduce_degrees, extract_circuit,
compile_measurement, generate_pauli_measurement. Circuit export requires Stim.
MC-flow checking reconstructs a valid order from O1/O2; additional scheduling
constraints are not retained. Circuits have explicit all-zero postselection
requirements and normalization factors, not adaptive outcome recovery.
See README_general.md and validate_10.py for the full supported pipeline.
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
    from mc_preparation import validate_flow
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


@dataclass
class CircuitOperation:
    name: str
    targets: tuple


@dataclass
class ExtractedCircuit:
    """Normalized physical operations plus explicit postselection metadata.

    graph_map = branch_scale * physical_postselected_map, with external
    boundaries ordered by inputs/outputs. Every measurement result must be
    zero. This is not an adaptive all-outcomes implementation.
    """
    operations: tuple
    num_qubits: int
    inputs: tuple
    outputs: tuple
    branch_scale: complex

    def to_stim(self):
        import stim
        circuit = stim.Circuit()
        for op in self.operations:
            if op.name in ('MXX', 'MZZ'):
                f = stim.target_x if op.name == 'MXX' else stim.target_z
                circuit.append('MPP', [f(op.targets[0]), stim.target_combiner(), f(op.targets[1])])
            else:
                circuit.append(op.name, op.targets)
        return circuit

    def metadata(self):
        return {'num_qubits': self.num_qubits, 'input_wires': self.inputs,
                'output_wires': self.outputs,
                'branch_scale': [self.branch_scale.real, self.branch_scale.imag],
                'postselect_measurement_records_to_zero': list(range(self.to_stim().num_measurements)),
                'semantics': 'ZX map = branch_scale * all-zero postselected circuit map'}


def extract_circuit(graph, paths, *, strict=True):
    """Proposition 4.11/4.13 extraction of a reduced, well-covered diagram.

    Handles endpoints, H, Clifford X/Z phase gates, same-colour identities,
    weight-one/two Pauli projectors, and CNOT (strict=False). Interacting
    spiders must be phaseless; phases belong on degree-two gate vertices.
    Returns a Stim-exportable circuit with exact scalar and wire mapping.
    """
    from mc_preparation import validate_flow
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


@dataclass
class Compilation:
    measurement: object
    prepared: object
    reduced: ReducedDiagram
    circuit: ExtractedCircuit


def compile_measurement(pauli, outcome=1, *, check_each=False):
    """Arbitrary positive-weight Pauli -> prepared ZX -> reduced ZX -> circuit.

    pauli may be a signed string or a Measurement object. For an object,
    its stored outcome is used. Returned circuit is explicitly postselected.
    """
    from pauli_measurements import pauli_measurement
    from mc_preparation import prepare_mc_flow
    m = pauli_measurement(pauli, outcome) if isinstance(pauli, str) else pauli
    prepared = prepare_mc_flow(m)
    reduced = reduce_degrees(prepared.graph, prepared.paths, check_each=check_each)
    circuit = extract_circuit(reduced.graph, reduced.paths)
    return Compilation(m, prepared, reduced, circuit)


def generate_pauli_measurement(weight, **kwargs):
    """Convenience entry point; see pauli_measurements.generate_pauli_measurement."""
    from pauli_measurements import generate_pauli_measurement as generate
    return generate(weight, **kwargs)


if __name__ == '__main__':
    unittest.main(verbosity=2)
