"""Truncated bivariate power series with exact integer/rational coefficients.

A series is a list of lists a[i][j] = coefficient of x^i y^j, kept for i + j <= N.
"""
from fractions import Fraction


class BS:
    __slots__ = ("N", "c")

    def __init__(self, N, c=None):
        self.N = N
        self.c = c if c is not None else [[0] * (N + 1 - i) for i in range(N + 1)]

    @classmethod
    def const(cls, N, v):
        s = cls(N)
        s.c[0][0] = v
        return s

    @classmethod
    def mono(cls, N, i, j, v=1):
        s = cls(N)
        if i + j <= N:
            s.c[i][j] = v
        return s

    def copy(self):
        return BS(self.N, [r[:] for r in self.c])

    def __getitem__(self, ij):
        i, j = ij
        if i < 0 or j < 0 or i + j > self.N:
            return 0
        return self.c[i][j]

    def __setitem__(self, ij, v):
        i, j = ij
        self.c[i][j] = v

    def _bin(self, o, f):
        if not isinstance(o, BS):
            o = BS.const(self.N, o)
        return BS(self.N, [[f(a, b) for a, b in zip(r, s)] for r, s in zip(self.c, o.c)])

    def __add__(self, o):
        return self._bin(o, lambda a, b: a + b)

    __radd__ = __add__

    def __sub__(self, o):
        return self._bin(o, lambda a, b: a - b)

    def __rsub__(self, o):
        return BS.const(self.N, o) - self

    def __neg__(self):
        return self * -1

    def __mul__(self, o):
        N = self.N
        if not isinstance(o, BS):
            return BS(N, [[a * o for a in r] for r in self.c])
        out = BS(N)
        oc = o.c
        for i, r in enumerate(self.c):
            for j, a in enumerate(r):
                if a == 0:
                    continue
                for k in range(N - i - j + 1):
                    row = oc[k]
                    dst = out.c[i + k]
                    for l in range(N - i - j - k + 1):
                        b = row[l]
                        if b:
                            dst[j + l] += a * b
        return out

    __rmul__ = __mul__

    def inv(self):
        """1/self, requires self[0,0] invertible (+-1 or rational)."""
        N = self.N
        a0 = self.c[0][0]
        out = BS(N)
        # solve out*self = 1 degree by degree
        for d in range(N + 1):
            for i in range(d + 1):
                j = d - i
                s = 1 if d == 0 else 0
                for k in range(i + 1):
                    for l in range(j + 1):
                        if (k, l) != (0, 0):
                            b = self[k, l]
                            if b:
                                s -= b * out.c[i - k][j - l]
                out.c[i][j] = s / a0 if isinstance(a0, Fraction) or a0 not in (1, -1) else s * a0
        return out

    def scale(self, p, q):
        """f(x^p, y^q), truncated."""
        out = BS(self.N)
        for i, r in enumerate(self.c):
            for j, a in enumerate(r):
                if a and p * i + q * j <= self.N:
                    out.c[p * i][q * j] = a
        return out

    def shift(self, p, q):
        """x^p y^q f."""
        out = BS(self.N)
        for i, r in enumerate(self.c):
            for j, a in enumerate(r):
                if a and i + p + j + q <= self.N:
                    out.c[i + p][j + q] = a
        return out

    def items(self):
        for i, r in enumerate(self.c):
            for j, a in enumerate(r):
                if a:
                    yield i, j, a


def compose_xy(m, X, Y):
    """m(X(x,y), Y(x,y)) where X, Y have zero constant term."""
    N = m.N
    out = BS(N)
    # powers of X and Y
    Xp = [BS.const(N, 1)]
    Yp = [BS.const(N, 1)]
    for _ in range(N):
        Xp.append(Xp[-1] * X)
        Yp.append(Yp[-1] * Y)
    for i, j, a in m.items():
        out = out + (Xp[i] * Yp[j]) * a
    return out
