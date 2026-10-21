"""Rooted nonseparable maps and rooted c-nets by vertices, faces and root-face degree (mod p).

Arrays a[i, j, k] (numpy int64, entries < p < 2^24) = coefficient of x^i y^j z^k, rectangular
truncation i <= X, j <= Y, k <= Z.

Nonseparable maps M (>= 2 edges), z = degree of the root face (right of the root edge). Deleting the
root edge leaves a chain of k >= 1 blocks (single edges or nonseparable maps rooted at their edge
leaving the chain's start on the right side); the chain's next attachment vertex sits at distance
j = 1 .. d_B - 1 along the block's right side:
    F(z) = y^2 z sum_{k>=1} x^(1-k) b(z)^k,   b(z) = x^2 z + (z F(1) - F(z)) / (y (1 - z)).
Networks D(z) = F / (x^2 y^2 z) (z marks the length of the right side), Tutte decomposition
    D = z + S + P + H,  S = sum_{k>=2} x^(k-1) A(z)^k (A = D - S),
    P = sum_{k>=2} y^(k-1) B(1)^(k-1) B(z) (B = D - P),
    H = sum_C x^(v-2) y^(f-2) D(1)^(e-d) D(z)^(d-1)   (C rooted c-net, d = its root-face degree)
is inverted for the c-net counts c(v, f, d).
"""
import numpy as np

from fastmul import mulmod

PRIMES24 = [16777213, 16777199, 16777183, 16777153, 16777141, 16777139, 16777127, 16777121, 16777099, 16777049,
            16777027, 16776989]   # the 12 largest primes below 2^24, in decreasing order


class T3:
    def __init__(self, X, Y, Z, p):
        self.X, self.Y, self.Z, self.p = X, Y, Z, p
        self.shape = (X + 1, Y + 1, Z + 1)

    def zero(self):
        return np.zeros(self.shape, dtype=np.int64)

    def mono(self, i, j, k, c=1):
        a = self.zero()
        if i <= self.X and j <= self.Y and k <= self.Z:
            a[i, j, k] = c % self.p
        return a

    def mul(self, a, b):
        return mulmod(a, b, self.p, self.shape)

    def inv(self, a):
        """1/a, a[0,0,0] invertible."""
        p = self.p
        g = self.mono(0, 0, 0, pow(int(a[0, 0, 0]), p - 2, p))
        for _ in range(int(np.log2(self.X + self.Y + self.Z + 2)) + 2):
            g = self.mul(g, (2 * self.mono(0, 0, 0) - self.mul(a, g)) % p)
        return g

    def at1(self, a):
        """a(x, y, 1) as a z-constant array."""
        out = self.zero()
        out[:, :, 0] = a.sum(axis=2) % self.p
        return out


def conv2(a, b, Y, Z, p):
    """2D truncated convolution of (Y+1, Z+1) arrays, entries < 2^24 (exact in int64)."""
    out = np.zeros((Y + 1, Z + 1), dtype=np.int64)
    ra = [j for j in range(Y + 1) if a[j].any()]
    rb = [j for j in range(Y + 1) if b[j].any()]
    for j in ra:
        aj = a[j]
        for l in rb:
            if j + l > Y:
                break
            out[j + l] += np.convolve(aj, b[l])[:Z + 1]
        out %= p  # keep sums small: each convolve sums <= Z+1 products < 2^48
    return out


def nonsep(T):
    """F[i, j, k] = rooted nonseparable maps (>= 2 edges) with i vertices, j faces, root face degree k."""
    p = T.p
    F = T.zero()
    xz2 = T.mono(2, 0, 1)
    one = T.mono(0, 0, 0)
    g = one                    # 1 / (1 - b/x), refined by one Newton step per iteration (warm start)
    for _ in range(2 * (T.X + T.Y) + 8):
        # (z F(1) - F(z)) / (1 - z): prefix sums along z
        P = -F
        P[:, :, 1:] += T.at1(F)[:, :, :-1]
        Q = np.cumsum(P % p, axis=2) % p
        # divide by y: shift down in j (Q has no j = 0 terms since F has y^2 factors)
        Qy = np.zeros_like(Q)
        Qy[:, :-1, :] = Q[:, 1:, :]
        b = (xz2 + Qy) % p
        # b / x: shift down in i (b has x^2 factors)
        bx = np.zeros_like(b)
        bx[:-1] = b[1:]
        # F = y^2 z x * (b/x) / (1 - b/x)
        den = (one - bx) % p
        g = T.mul(g, (2 * one - T.mul(den, g)) % p)
        G = T.mul(bx, g)
        newF = np.zeros_like(G)
        newF[1:, 2:, 1:] = G[:-1, :-2, :-1]
        if np.array_equal(newF, F) and np.array_equal(T.mul(den, g), one):
            return F
        F = newF
    raise RuntimeError('nonsep did not converge')


def cnets_by_rootface(T):
    """c[i, j, k] = rooted c-nets with i vertices, j faces, root face of degree k."""
    p = T.p
    F = nonsep(T)
    # D(z) = F / (x^2 y^2 z), indexed by internal vertices / internal faces / right side length
    D = T.zero()
    D[:-2, :-2, :-1] = F[2:, 2:, 1:]
    x, y, z = T.mono(1, 0, 0), T.mono(0, 1, 0), T.mono(0, 0, 1)
    one = T.mono(0, 0, 0)
    # S = x A^2 / (1 - x A) with A = D - S gives D = A / (1 - x A), so S = x D^2 / (1 + x D); likewise
    # P = y B(1) B / (1 - y B(1)) with B = D - P gives B(1) = D(1) / (1 + y D(1)), P = D y D(1) / (1 + y D(1)).
    D1 = T.at1(D)
    S = T.mul(T.mul(x, T.mul(D, D)), T.inv((one + T.mul(x, D)) % p))
    yD1 = T.mul(y, D1)
    Pn = T.mul(T.mul(D, yD1), T.inv((one + yD1) % p))
    H = (D - z - S - Pn) % p
    # invert H = sum c x^(v-2) y^(f-2) D1^(v+f-2-d) D^(d-1)
    c = T.zero()
    R = H.copy()
    D1pow = [one]
    Dpow = [one]
    for _ in range(T.X + T.Y + 2):
        D1pow.append(T.mul(D1pow[-1], D1))
    for _ in range(T.Z + 1):
        Dpow.append(T.mul(Dpow[-1], D))
    terms = {}                      # D1^(v+f-2-d) D^(d-1) depends only on (v + f, d)
    for tot in range(T.X + T.Y + 1):
        for a in range(min(tot, T.X) + 1):
            bdeg = tot - a
            if bdeg > T.Y:
                continue
            for k in range(T.Z):
                r = int(R[a, bdeg, k])
                if not r:
                    continue
                d = k + 1
                v, f = a + 2, bdeg + 2
                if v <= T.X and f <= T.Y and d <= T.Z:
                    c[v, f, d] = r
                key = (v + f - 2 - d, d - 1)
                if key not in terms:
                    terms[key] = T.mul(D1pow[key[0]], Dpow[key[1]])
                term = terms[key]
                sh = np.zeros_like(term)
                sh[a:, bdeg:, :] = term[:T.X + 1 - a, :T.Y + 1 - bdeg, :]
                R = (R - r * sh) % p
    return c, F


if __name__ == '__main__':
    import sys
    from math import factorial
    X, Y = int(sys.argv[1]), int(sys.argv[2])
    p = PRIMES24[0]
    T = T3(X, Y, X, p)
    c, F = cnets_by_rootface(T)

    def Bf(i, j):
        return factorial(2 * i + j - 2) * factorial(2 * j + i - 2) // (
            factorial(i) * factorial(j) * factorial(2 * i - 1) * factorial(2 * j - 1))

    bad = 0
    for v in range(2, X + 1):
        for f in range(2, Y + 1):
            if int(F[v, f].sum() % p) != Bf(v - 1, f - 1) % p:
                bad += 1
    print('nonsep vs Brown-Tutte: bad', bad)
    from pipeline import networks, rooted_cnets
    from ms import Ctx
    ctx = Ctx(X + Y, p)
    m = rooted_cnets(ctx)
    bad = 0
    for v in range(4, X + 1):
        for f in range(4, Y + 1):
            if v + f - 4 <= X + Y and int(c[v, f].sum() % p) != m[v - 2, f - 2] % p:
                bad += 1
                print('cnet mismatch', v, f, int(c[v, f].sum() % p), m[v - 2, f - 2])
    print('c-nets summed over d vs rooted: bad', bad)
    print('K4', c[4, 4], 'cube rooted by face deg', c[8, 6, :6], 'octa', c[6, 8, :6])
