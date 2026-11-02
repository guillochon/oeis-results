"""Exact rational g.f.s for A154640-A154642 (strategy (a)) from the Coxeter formula of series.py,
computed over Q(x) with sympy's DomainMatrix (no truncation, no fitting).

The squared bends use the symmetric square Sym^2 (dimension 15) of the 5-dim representation.
"""
from itertools import combinations, combinations_with_replacement
import sympy as sp
from sympy.polys.matrices import DomainMatrix

x = sp.symbols("x")
K = sp.QQ.frac_field(x)
N = 5


def R(i):
    m = sp.eye(N)
    m[i, :] = sp.ones(1, N)
    m[i, i] = -1
    return m


def sym2(m):
    """Matrix of X -> m X m^T on symmetric 5x5 X, coordinates X[a][b] (a <= b)."""
    basis = list(combinations_with_replacement(range(N), 2))
    out = sp.zeros(len(basis))
    for c, (a, b) in enumerate(basis):
        E = sp.zeros(N)
        E[a, b] = 1
        E[b, a] = 1
        Y = m * E * m.T
        for r, (p, q) in enumerate(basis):
            out[r, c] = Y[p, q]
    return out, basis


def dm(m):
    return DomainMatrix.from_Matrix(m).convert_to(K)


def gf(rep):
    gens = [dm(rep(R(i))) for i in range(N)]
    d = gens[0].shape[0]
    I = DomainMatrix.eye(d, K)
    X = K.from_sympy(x)

    def term(J):
        if not J:
            return I
        if len(J) == 1:
            g = gens[J[0]]
            P = I + g * X
            return -(P.inv() * g) * X
        a, b = gens[J[0]], gens[J[1]]
        P = I + (a + b) * X + (a * b + b * a) * X ** 2 + (a * b * a) * X ** 3
        return P.inv() * (a * b * a) * X ** 3

    memo = {}

    def A(S_):
        tot = DomainMatrix.zeros((d, d), K)
        for r in range(3):
            for J in combinations(S_, r):
                if J not in memo:
                    memo[J] = term(J)
                tot = tot + memo[J]
        return tot

    S = A(tuple(range(N))).inv()
    return S, [A(tuple(k for k in range(N) if k != j)) for j in range(N)]


def scalar(e):
    return sp.factor(sp.cancel(K.to_sympy(e)))


# counts (trivial representation)
S, Ah = gf(lambda m: sp.Matrix([[1]]))
cnt = sum(((Ah[j] * S).to_Matrix())[0, 0] for j in range(N))
print("A154640 =", sp.factor(sp.cancel(cnt)))

b0 = sp.Matrix([1, 1, 1, 0, 0])
S, Ah = gf(dm_rep := (lambda m: m))
tot = 0
for j in range(N):
    M = (Ah[j] * S).to_Matrix()
    tot += (M * b0)[j]
print("A154641 =", sp.factor(sp.cancel(tot / 3)))

basis = sym2(sp.eye(N))[1]
S, Ah = gf(lambda m: sym2(m)[0])
# squared bend of row j of M is (M b0)_j^2 = (M B M^T)_jj with B = b0 b0^T
w = sp.Matrix([b0[a] * b0[b] for a, b in basis])
pos = {p: k for k, p in enumerate(basis)}
tot = 0
for j in range(N):
    M = (Ah[j] * S).to_Matrix()
    tot += (M * w)[pos[(j, j)]]
print("A154642 =", sp.factor(sp.cancel(tot / 3)))
