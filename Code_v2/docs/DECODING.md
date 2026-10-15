# Mathematical basis of the Paulimer decoder

## Mathematical extraction

Paulimer represents the final pure state of the data, retained physical ancillas,
and references as `R |A u>`, with uniformly random independent bits `u`.
The measured record is `b = M u + v0`, over GF(2).

For an observable Q, Paulimer computes `R^dagger Q R`. If this is
`(-1)^c Z^v`, its eigenvalue bit is `c + v^T A u`. Any X/Y factor means the
observable is not definite in a branch; a required target observable then fails
verification.

Apply that rule to `P^T` on the reference register. This yields the global outcome
bit. The transpose includes one minus sign for every Y, as well as the target's
explicit sign. Binar solves `M^T a = h` to express any free-bit parity `h^T u`
as `a^T b + a^T v0` in actual measurement records.

A second Paulimer simulation builds the ideal Choi state for the signed target
Pauli measurement. For each of its 2n independent stabilizers, compare actual
and ideal signs. The commutation of output `X^x Z^z` with that stabilizer gives
one GF(2) equation in x,z. Solve for the constant and every free-bit coefficient
using one Binar echelon decomposition. A missing solution is an error.

Finally, apply the resulting feedback through Paulimer and verify all 2n ideal
Choi generators. This also establishes separation from discarded ancillas.
A nonconstant affine decoded outcome is balanced over the k free Choi bits.
Each record has Choi probability 2^-k, so its corrected Kraus map has squared
scalar 2^(1-k) relative to the ideal projector. There are 2^(k-1) records per
outcome; their weights sum to one. Thus the check covers Born probabilities as
well as normalized output states, up to irrelevant branch-global phases.


## References

- https://github.com/microsoft/qdk-ec
- https://arxiv.org/abs/2309.08676
