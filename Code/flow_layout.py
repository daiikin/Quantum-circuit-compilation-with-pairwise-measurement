"""Flow-based layout for the compiler's REDUCED ZX diagrams.

Usage:
    fig = draw_flow(result.reduced.graph, result.reduced.paths)
    fig.savefig('zx_flow.svg', bbox_inches='tight')

Only coordinates change; graph structure, phases, scalar and paths do not.
The original graph is never modified. Requires vertex-disjoint reduced paths.
Solid arrowed edges follow flow; dashed edges are uncovered interactions.
Crossings without a spider are not graph connections.
"""
from copy import deepcopy
from collections import defaultdict
import heapq
import math
from pyzx.utils import VertexType as V


def layout_flow(graph, paths, *, wire_spacing=1.5, time_spacing=1.0):
    """Return a graph copy with row=time and qubit=worldline coordinates.

    Interacting spiders share a column. A longest-path schedule is extended
    to avoid placing unrelated vertices/vertical interactions on the same
    column over overlapping wire intervals. Not a depth optimizer.
    """
    paths = tuple(tuple(p) for p in paths)
    if wire_spacing <= 0 or time_spacing <= 0:
        raise ValueError('Spacing must be positive')
    owner, covered = {}, set()
    for q,p in enumerate(paths):
        for v in p:
            if v in owner:
                raise ValueError('Use reduced paths: vertices must belong to one path only')
            owner[v] = q
        for a,b in zip(p,p[1:]):
            if not graph.connected(a,b):
                raise ValueError('Path contains a missing edge')
            covered.add(frozenset((a,b)))
    # An uncovered edge is an atomic event, including a one-legged leaf.
    groups, event = [], {}
    for v in graph.vertices():
        if v in event:
            continue
        neighbors = [w for w in graph.neighbors(v) if frozenset((v,w)) not in covered]
        if len(neighbors) > 1:
            raise ValueError('Expected a reduced, well-covered graph')
        group = (v, neighbors[0]) if neighbors else (v,)
        if any(w in event for w in group):
            raise ValueError('Uncovered edges do not form separate interactions')
        if not any(w in owner for w in group):
            raise ValueError('Disconnected scalar component has no worldline')
        for w in group:
            event[w] = len(groups)
        groups.append(group)
    succ, pred = defaultdict(set), defaultdict(set)
    for p in paths:
        for a,b in zip(p,p[1:]):
            u,v = event[a],event[b]
            if u == v:
                raise ValueError('An interaction cannot connect a wire to itself')
            succ[u].add(v)
            pred[v].add(u)
    indegree = [len(pred[i]) for i in range(len(groups))]
    ready = [i for i,d in enumerate(indegree) if not d]
    heapq.heapify(ready)
    times, occupied = {}, defaultdict(list)
    while ready:
        e = heapq.heappop(ready)
        wires = [owner[v] for v in groups[e] if v in owner]
        lo, hi = min(wires), max(wires)
        t = max((times[p]+1 for p in pred[e]), default=0)
        while any(not (hi < a or lo > b) for a,b in occupied[t]):
            t += 1
        times[e] = t
        occupied[t].append((lo,hi))
        for nxt in sorted(succ[e]):
            indegree[nxt] -= 1
            if indegree[nxt] == 0:
                heapq.heappush(ready,nxt)
    if len(times) != len(groups):
        raise ValueError('Flow has a causal cycle')
    g = deepcopy(graph)
    last = max(times.values(), default=0)+1
    inputs, outputs = set(g.inputs()), set(g.outputs())
    for v in g.vertices():
        t = 0 if v in inputs else last if v in outputs else times[event[v]]
        if v in owner:
            y = owner[v]
        else:
            neighbor = next(iter(g.neighbors(v)))
            y = owner[neighbor]-.32
        g.set_row(v, t*time_spacing)
        g.set_qubit(v, y*wire_spacing)
    return g


def draw_flow(graph, paths, *, labels=False, title='Reduced ZX diagram with MC flow'):
    """Return a matplotlib Figure; pass to savefig() or call plt.show()."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    paths = tuple(tuple(p) for p in paths)
    g = layout_flow(graph, paths)
    xmax = max((g.row(v) for v in g.vertices()), default=1)
    fig, ax = plt.subplots(figsize=(max(9, min(32, xmax*.48)), max(3, len(paths)*.62)))
    covered = {frozenset((a,b)) for p in paths for a,b in zip(p,p[1:])}
    xy = lambda v: (g.row(v), -g.qubit(v))
    # Uncovered vertical edges: dashed to distinguish crossings from vertices.
    for edge in g.edges():
        a,b = g.edge_st(edge)
        if frozenset((a,b)) not in covered:
            x1,y1 = xy(a)
            x2,y2 = xy(b)
            ax.plot([x1,x2], [y1,y2], color='#6b7280', lw=1.3, ls='--', zorder=1)
    for q,p in enumerate(paths):
        for a,b in zip(p,p[1:]):
            ax.annotate('', xy=xy(b), xytext=xy(a),
                        arrowprops=dict(arrowstyle='->', color='#334155', lw=1.1,
                                        shrinkA=5, shrinkB=5), zorder=2)
        ax.text(-.7, -1.5*q, f'p{q}', ha='right', va='center', fontsize=9)
    for v in g.vertices():
        ty = g.type(v)
        color = {V.Z:'#93d6a2', V.X:'#ee9a9a', V.H_BOX:'#f6d66d', V.BOUNDARY:'#334155'}.get(ty,'white')
        x,y = xy(v)
        ax.scatter([x], [y], s=35 if ty==V.BOUNDARY else 145,
                   marker='s' if ty==V.H_BOX else 'o', c=color,
                   edgecolors='#334155', linewidths=.7, zorder=4)
        if ty==V.H_BOX:
            ax.text(x,y,'H',ha='center',va='center',fontsize=7,zorder=5)
        elif g.phase(v) != 0:
            ax.text(x,y+.3, f'{g.phase(v)}π',ha='center',va='bottom',fontsize=8,zorder=5)
        if labels:
            ax.text(x,y-.3,str(v),ha='center',va='top',fontsize=7,zorder=5)
    ax.set_title(title, pad=18)
    ax.set_xlim(-1.6,xmax+.8)
    ax.set_ylim(-1.5*(len(paths)-1)-.9, 1.1)
    ax.set_axis_off()
    ax.legend(handles=[Line2D([0],[0],color='#334155',label='Directed MC-flow path'),
                       Line2D([0],[0],color='#6b7280',ls='--',label='Uncovered edge')],
              loc='upper center',bbox_to_anchor=(.5,-.02),ncol=2,frameon=False)
    fig.tight_layout()
    return fig
