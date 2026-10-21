"""Rooted general maps on the sphere and on the projective plane by vertices, faces, root-face degree.

Root-edge deletion (x = vertices, y = faces, z = root-face degree), P planar (incl. the vertex map),
Q projective:
    P = x y + z^2 P^2 / y + y z (P(1) - z P) / (1 - z)
    Q = 2 z^2 P Q / y + y z (Q(1) - z Q) / (1 - z) + z^2 d/dz (z P)
(bridge; non-bridge between two faces; on RP^2 also a twisted edge inside the root face of a planar
map, k + 1 positions on a root face of degree k).  Arrays use t3.T3 (mod p).
"""
import numpy as np

from t3 import T3


def divide_1mz(T, A):
    """A / (1 - z) for A vanishing at z = 1 (prefix sums along z)."""
    return np.cumsum(A % T.p, axis=2) % T.p


def maps_sphere_rp2(T):
    p = T.p
    P = T.zero()
    xy = T.mono(1, 1, 0)
    z = T.mono(0, 0, 1)
    for _ in range(T.X + T.Y + 3):
        P1 = T.at1(P)
        zP = T.mul(z, P)
        nonb = divide_1mz(T, (P1 - zP) % p)
        nonb = T.mul(T.mono(0, 1, 1), nonb)
        br = T.mul(T.mono(0, 0, 2), T.mul(P, P))
        brd = np.zeros_like(br)
        brd[:, :-1] = br[:, 1:]          # / y
        newP = (xy + brd + nonb) % p
        if np.array_equal(newP, P):
            break
        P = newP
    # d/dz (z P): coefficient of z^k times (k + 1)
    k = np.arange(T.Z + 1, dtype=np.int64)
    tw = (P * (k + 1)) % p
    tw = T.mul(T.mono(0, 0, 2), tw)
    Q = T.zero()
    twoP = T.mul(T.mono(0, 0, 2), P) * 2 % p
    for _ in range(T.X + T.Y + 3):
        Q1 = T.at1(Q)
        nonb = T.mul(T.mono(0, 1, 1), divide_1mz(T, (Q1 - T.mul(z, Q)) % p))
        br = T.mul(twoP, Q)
        brd = np.zeros_like(br)
        brd[:, :-1] = br[:, 1:]
        newQ = (brd + nonb + tw) % p
        if np.array_equal(newQ, Q):
            break
        Q = newQ
    return P, Q


if __name__ == '__main__':
    import sys
    from t3 import PRIMES24
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    p = PRIMES24[0]
    T = T3(X, X, 2 * X, p)
    P, Q = maps_sphere_rp2(T)
    # by edges: planar e = v + f - 2, projective e = v + f - 1
    pe, qe = {}, {}
    for v in range(X + 1):
        for f in range(X + 1):
            s = int(P[v, f].sum() % p)
            if s:
                pe[v + f - 2] = pe.get(v + f - 2, 0) + s
            s = int(Q[v, f].sum() % p)
            if s:
                qe[v + f - 1] = qe.get(v + f - 1, 0) + s
    print('planar by edges', [pe.get(e, 0) % p for e in range(X)])
    print('A000168        ', [1, 2, 9, 54, 378, 2916, 24054, 208494][:X])
    print('RP2 by edges   ', [qe.get(e, 0) % p for e in range(1, X)])
    print('A007137        ', [1, 10, 98, 982, 10062, 105024, 1112757][:X - 1])
