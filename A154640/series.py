"""Strategy (a) sequences to any order, without enumerating the group.

Matrix growth series S(x) = sum_w x^l(w) M_w of the Coxeter group W (five
generators, all m_ij = 3) in the Descartes-reflection representation.
Inclusion-exclusion over right descent sets (only the finite parabolics
J = {}, {i}, {i,k} occur) gives

    S * A = I,   A = sum_{J finite} (-1)^|J| x^l(w_J) P_J^{-1} M_{w_J},

with P_J = sum_{v in W_J} x^l(v) M_v.  Spheres <-> pairs (j, coset W_jhat w)
(Tits: row j = e_j^T is a face point of the fundamental chamber whose
stabilizer is W_jhat), and S = P_jhat * (sum over minimal left coset reps), so
the sum over minimal reps is A_jhat * S, where A_jhat is the same alternating
sum restricted to subsets of jhat.

  A154640(n) = [x^n] sum_j  (A_jhat S)_{jj} summed with the trivial rep (counts),
  A154641(n) = [x^n] sum_j  e_j^T A_jhat S b0 / s0,
  A154642(n) = same in the tensor-square representation, (e_j x e_j)^T ... (b0 x b0) / t0.
Power series are truncated at order K with exact integers.
"""
import sys
from itertools import combinations
import numpy as np

K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
N = 5


def R(i):
    m = np.eye(N, dtype=object)
    m[i, :] = 1
    m[i, i] = -1
    return m


class PS:
    """Truncated power series with d x d object-integer matrix coefficients."""

    def __init__(self, c):
        self.c = c  # list of K+1 matrices

    @staticmethod
    def const(m, deg=0):
        d = m.shape[0]
        c = [np.zeros((d, d), dtype=object) for _ in range(K + 1)]
        if deg <= K:
            c[deg] = m.copy()
        return PS(c)

    def __add__(self, o):
        return PS([a + b for a, b in zip(self.c, o.c)])

    def __sub__(self, o):
        return PS([a - b for a, b in zip(self.c, o.c)])

    def __mul__(self, o):
        d = self.c[0].shape[0]
        out = [np.zeros((d, d), dtype=object) for _ in range(K + 1)]
        nz = [k for k in range(K + 1) if any(self.c[k].flat)]
        for i in nz:
            for j in range(K + 1 - i):
                out[i + j] = out[i + j] + self.c[i].dot(o.c[j])
        return PS(out)

    def inv(self):
        # requires c[0] = I
        d = self.c[0].shape[0]
        assert (self.c[0] == np.eye(d, dtype=object)).all()
        out = [np.eye(d, dtype=object)] + [None] * K
        for n in range(1, K + 1):
            acc = np.zeros((d, d), dtype=object)
            for k in range(1, n + 1):
                if any(self.c[k].flat):
                    acc = acc + self.c[k].dot(out[n - k])
            out[n] = -acc
        return PS(out)


def run(rep):
    """rep: function mapping a 5x5 integer matrix to its matrix in the representation."""
    gens = [rep(R(i)) for i in range(N)]
    d = gens[0].shape[0]
    I = np.eye(d, dtype=object)

    memo = {}

    def term(J):
        J = tuple(J)
        if J not in memo:
            memo[J] = term_(J)
        return memo[J]

    def term_(J):
        J = list(J)
        if not J:
            return PS.const(I)
        if len(J) == 1:
            g = gens[J[0]]
            P = PS.const(I) + PS.const(g, 1)
            return PS([-c for c in (P.inv() * PS.const(g, 1)).c])
        a, b = (gens[J[0]], gens[J[1]])
        P = (PS.const(I) + PS.const(a + b, 1) + PS.const(a.dot(b) + b.dot(a), 2)
             + PS.const(a.dot(b).dot(a), 3))
        return P.inv() * PS.const(a.dot(b).dot(a), 3)

    def A(S_):
        tot = PS.const(np.zeros((d, d), dtype=object))
        for r in range(3):
            for J in combinations(S_, r):
                tot = tot + term(J)
        return tot

    S = A(range(N)).inv()
    Ahat = [A([k for k in range(N) if k != j]) for j in range(N)]
    return S, Ahat


def kron(m):
    return np.kron(m, m)


b0 = np.array([1, 1, 1, 0, 0], dtype=object)
s0 = 3
t0 = 3

# counts: trivial representation (1x1)
S1, A1 = run(lambda m: np.array([[1]], dtype=object))
growth = [S1.c[k][0, 0] for k in range(K + 1)]
cnt = [sum((A1[j] * S1).c[k][0, 0] for j in range(N)) for k in range(K + 1)]

S, Ah = run(lambda m: m)
bend = [sum(((Ah[j] * S).c[k].dot(b0))[j] for j in range(N)) for k in range(K + 1)]

S2, Ah2 = run(kron)
bb = np.kron(b0, b0)
sq = [sum(((Ah2[j] * S2).c[k].dot(bb))[j * N + j] for j in range(N)) for k in range(K + 1)]

assert all(v % s0 == 0 for v in bend) and all(v % t0 == 0 for v in sq)
print("growth", growth)
print("A154640", cnt)
print("A154641", [v // s0 for v in bend])
print("A154642", [v // t0 for v in sq])
