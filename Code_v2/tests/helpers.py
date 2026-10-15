"""Dense ideal-projector oracle, only for small numerical tests."""
import numpy as np
from pauli_compiler.zx.spiders import parse_pauli

def expected_projector(pauli, outcome=1):
    word, sign, _ = parse_pauli(pauli)
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    p = np.array([[sign]], dtype=complex)
    for c in word:
        p = np.kron(p, matrices[c])
    return (np.eye(p.shape[0]) + outcome*p)/2
