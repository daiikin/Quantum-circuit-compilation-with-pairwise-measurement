"""Weight-4/5 nondestructive Pauli measurement branches and r4/r5 rewrites.

Requires zx_r4_r5.py alongside this file. Leftmost character is qubit 0.
Each diagram is a Kraus operator K_m=(I+m*P)/2, not a classical-output wire.
Use measurement_branches() to construct the complete two-outcome instrument.
Only the local spider rewrite is performed; no MC flow or circuit extraction.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import argparse
import numpy as np
import pyzx as zx
from pyzx.utils import VertexType as V, EdgeType as E
from zx_r4_r5 import r4, r5


@dataclass
class Measurement:
    pauli: str
    sign: int
    outcome: int
    support: tuple
    original: object
    rewritten: object
    center: int
    ports: tuple
    ring: tuple


def parse_pauli(text):
    if not isinstance(text, str):
        raise ValueError('Pauli must be a string')
    text = ''.join(text.upper().split())
    sign = -1 if text.startswith('-') else 1
    if text.startswith(('+', '-')):
        text = text[1:]
    if not text or any(c not in 'IXYZ' for c in text):
        raise ValueError('Use an optional +/- followed by I, X, Y, Z')
    support = tuple(i for i, p in enumerate(text) if p != 'I')
    if len(support) not in (4, 5):
        raise ValueError('Pauli weight (number of non-I factors) must be 4 or 5')
    return text, sign, support


def pauli_measurement(pauli, outcome=1):
    """Build a normalized branch and apply r4 or r5 to its central spider.

    Examples: 'XYZX', '-XYZXY', 'XIYZZI'. outcome is the eigenvalue of
    the SIGNED operator, +1 or -1. Supports arbitrary identity positions.
    Local Clifford basis changes are undone on the output data wires.
    """
    if outcome not in (-1, 1):
        raise ValueError('outcome must be +1 or -1')
    word, sign, support = parse_pauli(pauli)
    g = zx.Graph()
    center = g.add_vertex(ty=V.X, phase=0, qubit=-2, row=4)
    ins, outs, ports = [], [], []
    # X spider has amplitude 2**(1-w/2) on even-parity strings.
    # This scalar makes the map exactly (I + outcome*P)/2.
    g.scalar.add_power(len(support)-2)
    odd_parity = sign*outcome == -1
    for q, p in enumerate(word):
        inp = g.add_vertex(ty=V.BOUNDARY, qubit=q, row=0)
        out = g.add_vertex(ty=V.BOUNDARY, qubit=q, row=8)
        ins.append(inp)
        outs.append(out)
        if p == 'I':
            g.add_edge(g.edge(inp, out))
            continue
        current = inp
        def node(ty, row, phase=0, edge=E.SIMPLE):
            nonlocal current
            v = g.add_vertex(ty=ty, phase=phase, qubit=q, row=row)
            g.add_edge(g.edge(current, v), edge)
            current = v
            return v
        # Chronological gates: X -> H; Y -> S-dagger, H; Z -> identity.
        # Thus B^dagger Z B equals the requested Pauli factor.
        if p == 'Y':
            node(V.Z, 1, Fraction(-1, 2))
        node(V.Z, 2, edge=E.HADAMARD if p in 'XY' else E.SIMPLE)
        # Conjugate the even-parity projector by X on one support wire
        # for the odd branch, keeping the rewritten spider at phase zero.
        if odd_parity and q == support[0]:
            node(V.X, 3, 1)
        tap = node(V.Z, 4)
        ports.append(tap)
        g.add_edge(g.edge(center, tap))
        if odd_parity and q == support[0]:
            node(V.X, 5, 1)
        node(V.Z, 6, edge=E.HADAMARD if p in 'XY' else E.SIMPLE)
        if p == 'Y':
            node(V.Z, 7, Fraction(1, 2))
        g.add_edge(g.edge(current, out))
    g.set_inputs(tuple(ins))
    g.set_outputs(tuple(outs))
    result = (r4 if len(support) == 4 else r5)(g, center, ports)
    # Layout only: spread the replacement ring above the data wires.
    for i, v in enumerate(result.ring):
        angle = 2*np.pi*i/len(ports)
        result.graph.set_row(v, 4 + 1.2*np.cos(angle))
        result.graph.set_qubit(v, -2 + .65*np.sin(angle))
    return Measurement(word, sign, outcome, support, g, result.graph,
                       center, tuple(ports), result.ring)


def measurement_branches(pauli):
    """Return {+1: Measurement, -1: Measurement}; classical labels are keys."""
    return {m: pauli_measurement(pauli, m) for m in (1, -1)}


def expected_projector(pauli, outcome=1):
    word, sign, _ = parse_pauli(pauli)
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    p = np.array([[sign]], dtype=complex)
    for c in word:
        p = np.kron(p, matrices[c])
    return (np.eye(p.shape[0]) + outcome*p)/2


def verify(measurement):
    """Dense small-system check against the exact signed Pauli projector."""
    label = ('-' if measurement.sign < 0 else '') + measurement.pauli
    target = expected_projector(label, measurement.outcome)
    for g in (measurement.original, measurement.rewritten):
        actual = zx.tensor.tensor_to_matrix(zx.tensor.tensorfy(g),
                                            len(g.inputs()), len(g.outputs()))
        np.testing.assert_allclose(actual, target, atol=1e-12)


def self_test():
    count = 0
    for weight in (4, 5):
        for letters in product('XYZ', repeat=weight):
            for m in (1, -1):
                verify(pauli_measurement(''.join(letters), m))
                count += 1
    for label in ('-XYZX', '-XYZXY', 'IXYIZZI', '-IXYIZXYZI'):
        # The last example deliberately has invalid weight six.
        if label == '-IXYIZXYZI':
            try:
                pauli_measurement(label)
            except ValueError:
                continue
            raise AssertionError('Invalid weight was accepted')
        branches = measurement_branches(label)
        for branch in branches.values():
            verify(branch)
            count += 1
        plus, minus = (expected_projector(label, m) for m in (1, -1))
        np.testing.assert_allclose(plus@minus, 0, atol=1e-12)
        np.testing.assert_allclose(plus.conj().T@plus + minus.conj().T@minus,
                                   np.eye(len(plus)), atol=1e-12)
    print(f'Passed {count} branch checks, before AND after rewriting.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pauli', nargs='?', default='XYZX')
    parser.add_argument('--outcome', type=int, choices=(-1, 1), default=1)
    parser.add_argument('--prefix', default='measurement')
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()
    if args.test:
        self_test()
        return
    m = pauli_measurement(args.pauli, args.outcome)
    verify(m)
    for suffix, g in (('before', m.original), ('after', m.rewritten)):
        with open(f'{args.prefix}_{suffix}.json', 'w') as f:
            f.write(g.to_json())
        figure = zx.draw_matplotlib(g, labels=True)
        figure.savefig(f'{args.prefix}_{suffix}.svg', bbox_inches='tight')
    print(f'Weight {len(m.support)}; applied r{len(m.support)}; projector verified.')
    print(f'Saved {args.prefix}_before/after.json and .svg')


if __name__ == '__main__':
    main()
