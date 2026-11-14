"""A132183 prototype: ideals of the binary majorization lattice M(n) via the split M(n) = M(k) x M(m).

Checks, for small sizes:
  * brute-force ideal counts a(0..7),
  * the split poset (pattern covers with id / sigma maps) is isomorphic to M(m+k),
  * the k = 2 chain formula and the k = 4 pair-table formula,
  * the self-dual variant against A109456.
Run: python attacks/A132183/proto.py
"""
import sys
from itertools import combinations
import numpy as np

A132183 = [2, 3, 5, 10, 27, 119, 1173, 44315, 16175190, 284432730176]
A109456 = [0, 1, 1, 2, 3, 7, 21, 135, 2470, 319124, 1214554343, 1706241214185942]


def leq(a, b):
    """a <= b in M(n); a, b bitmasks (bit i = element i+1, larger element = stronger)."""
    xa = [i for i in range(a.bit_length()) if a >> i & 1][::-1]
    xb = [i for i in range(b.bit_length()) if b >> i & 1][::-1]
    return len(xa) <= len(xb) and all(x <= y for x, y in zip(xa, xb))


def drop_max(a):
    return a & ~(1 << (a.bit_length() - 1)) if a else 0


class Lat:
    """L = J(M(m)): all ideals, Hasse edges grouped by removed element, sigma and dual maps."""

    def __init__(self, m):
        self.m = m
        n = 1 << m
        elems = sorted(range(n), key=lambda a: (sum(i + 1 for i in range(m) if a >> i & 1), a))
        pos = {a: i for i, a in enumerate(elems)}          # linear extension position
        self.elems, self.pos, self.n = elems, pos, n
        below = [0] * n                                    # strict down-set as bitmask over positions
        for i, a in enumerate(elems):
            for j, b in enumerate(elems[:i]):
                if leq(b, a):
                    below[i] |= 1 << j
        self.below = below
        # enumerate ideals (bitmask over positions)
        ideals = []

        def rec(i, cur):
            if i == n:
                ideals.append(cur)
                return
            rec(i + 1, cur)
            if below[i] & ~cur == 0:
                rec(i + 1, cur | 1 << i)
        sys.setrecursionlimit(10000)
        rec(0, 0)
        ideals.sort()
        self.ideals = ideals
        self.idx = {I: t for t, I in enumerate(ideals)}
        self.N = len(ideals)
        # Hasse edges (hi, lo) with hi = lo + p, grouped by p
        above = [0] * n
        for i in range(n):
            for j in range(i):
                if below[i] >> j & 1:
                    above[j] |= 1 << i
        self.edges = [[] for _ in range(n)]
        for t, I in enumerate(ideals):
            for p in range(n):
                if I >> p & 1 and above[p] & I == 0:
                    self.edges[p].append((t, self.idx[I & ~(1 << p)]))
        self.edges = [np.array(e, dtype=np.int64).reshape(-1, 2) for e in self.edges]
        # sigma(I) = {A : drop_max(A) in I}
        dm = [pos[drop_max(a)] for a in elems]
        self.sig = np.array([self.idx[sum(1 << i for i in range(n) if I >> dm[i] & 1)] for I in ideals])
        # dual(I) = {complement(A) : A not in I}
        comp = [pos[(n - 1) ^ a] for a in elems]
        self.dual = np.array([self.idx[sum(1 << comp[i] for i in range(n) if not I >> i & 1)] for I in ideals])

    # zeta transforms act on axis 0 of an object array
    def up(self, f):            # g(K) = sum_{I >= K} f(I)
        f = f.copy()
        for p in reversed(range(self.n)):
            e = self.edges[p]
            if len(e):
                f[e[:, 1]] += f[e[:, 0]]
        return f

    def down(self, f):          # g(J) = sum_{I <= J} f(I)
        f = f.copy()
        for p in range(self.n):
            e = self.edges[p]
            if len(e):
                f[e[:, 0]] += f[e[:, 1]]
        return f

    def push_sig(self, f):      # (sigma_* f)(K) = sum_{I : sigma(I) = K} f(I)
        g = np.zeros_like(f)
        np.add.at(g, self.sig, f)
        return g

    def contains(self):         # C[a, b] = 1 iff ideal a contains ideal b
        N = self.N
        C = np.zeros((N, N), dtype=object)
        for a in range(N):
            for b in range(N):
                C[a, b] = int(self.ideals[b] & ~self.ideals[a] == 0)
        return C


def ones(N):
    return np.array([1] * N, dtype=object)


def k2(L):
    """a(m+2) and the self-dual count for m+2 variables."""
    u = L.up(ones(L.N))[L.sig]                       # u(I01) = #{I00 >= sigma(I01)}
    w = L.down(L.push_sig(ones(L.N)))                # w(I10) = #{I11 : sigma(I11) <= I10}
    assert all(w[i] == u[L.dual[i]] for i in range(L.N))
    total = int((u * L.down(w)).sum())
    sd = sum(int(u[j]) for j in range(L.N) if L.ideals[L.dual[j]] & ~L.ideals[j] == 0)
    return total, sd


def k4(L):
    """a(m+4) and the self-dual count, via the pair table B(I_14, I_23)."""
    N = L.N
    C = L.contains()
    u = L.up(ones(N))[L.sig]                         # (empty, 1 | I_1)
    v = L.up(u)                                      # v(I_2) = sum_{I_1 >= I_2} u(I_1)
    G = np.zeros((N, N), dtype=object)
    G[np.arange(N), np.arange(N)] = v
    G = L.up(L.up(G).T).T                            # G(a, c) = sum_{I_2 >= a, c} v(I_2)
    T = G[:, L.sig]                                  # T(I_3, I_12)
    T2 = L.up(T.T).T                                 # T2(I_3, I_13) = sum_{I_12 >= I_13}
    T2 = T2 * C[:, L.sig]                            # I_3 >= sigma(I_13)
    T3 = L.up(T2)                                    # T3(I_4, I_13) = sum_{I_3 >= I_4}
    T4 = L.up(T3)[L.sig, :].T                        # T4(I_13, I_14) = sum_{I_4 >= sigma(I_14)}
    T4 = T4 * C                                      # I_13 >= I_14
    B = L.up(T4).T                                   # B(I_14, I_23) = sum_{I_13 >= I_14, I_23}
    Bd = B[L.dual, :][:, L.dual].T                   # Bd(e, f) = B(dual f, dual e)
    total = int((B * Bd).sum())
    sd = sum(int(B[e, L.dual[e]]) for e in range(N))
    return total, sd


def brute(n):
    return Lat(n).N


def check_split(m, k):
    """The poset on M(k) x M(m) generated by pattern covers equals M(m+k)."""
    def gen(a, p):                                   # all (b, q) covering-or-equal images
        out = []
        for i in range(k):                           # shift inside the pattern: id map
            if p >> i & 1 and i + 1 < k and not p >> (i + 1) & 1:
                out.append((a, p ^ (3 << i)))
        if not p & 1 and a >> (m - 1) & 1:           # shift m -> m+1: sigma-type cover
            out.append((a & ~(1 << (m - 1)), p | 1))
        return out
    full = lambda a, p: a | p << m
    # relation claimed: (A, P) <= (B, P') for a pattern cover iff A <= B (id) or drop_max(A) <= B (sigma)
    for p in range(1 << k):
        for q in range(1 << k):
            if p == q:
                continue
            cov_id = any(p >> i & 1 and i + 1 < k and not p >> (i + 1) & 1 and q == p ^ (3 << i) for i in range(k))
            cov_sg = (not p & 1) and q == p | 1
            for a in range(1 << m):
                for b in range(1 << m):
                    real = leq(full(a, p), full(b, q))
                    if cov_id:
                        assert real == leq(a, b)
                    elif cov_sg:
                        assert real == leq(drop_max(a), b)
    return True


if __name__ == "__main__":
    print("brute:", [brute(n) for n in range(7)])
    for m, k in [(2, 2), (3, 2), (3, 3), (2, 4), (3, 4)]:
        assert check_split(m, k)
    print("split relations ok")
    for m in range(0, 7):
        t, sd = k2(Lat(m))
        print(f"k=2 m={m}: a({m+2}) = {t} [{A132183[m+2]}]  selfdual({m+2}) = {sd} [{A109456[m+2]}]")
        assert t == A132183[m + 2] and sd == A109456[m + 2]
    for m in range(0, 6):
        t, sd = k4(Lat(m))
        ref = A132183[m + 4] if m + 4 < len(A132183) else None
        print(f"k=4 m={m}: a({m+4}) = {t} [{ref}]  selfdual({m+4}) = {sd} [{A109456[m+4]}]")
        assert ref is None or t == ref
        assert sd == A109456[m + 4]
