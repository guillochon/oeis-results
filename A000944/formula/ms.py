"""Truncated bivariate power series over Z/p (p < 2^26) on numpy int64 arrays.

a[i, j] = coefficient of x^i y^j, kept for i + j <= N (entries beyond are zero).
Rationals enter only as scalars and are mapped with modular inverses; final integers are
recovered by CRT over several primes (crt below).
"""
from fractions import Fraction

import numpy as np

from fastmul import mulmod

PRIMES = [67108859, 67108837, 67108819, 67108777, 67108763, 67108757, 67108753, 67108747, 67108739, 67108729,
          67108721, 67108709]   # the 12 largest primes below 2^26, in decreasing order


class Ctx:
    def __init__(self, N, p):
        self.N = N
        self.p = p
        i, j = np.indices((N + 1, N + 1))
        self.mask = (i + j <= N)

    def s(self, v):
        """Scalar (int or Fraction) mod p."""
        p = self.p
        if isinstance(v, Fraction):
            return v.numerator % p * pow(v.denominator % p, p - 2, p) % p
        return int(v) % p

    def zero(self):
        return MS(self, np.zeros((self.N + 1, self.N + 1), dtype=np.int64))

    def const(self, v):
        z = self.zero()
        z.a[0, 0] = self.s(v)
        return z

    def mono(self, i, j, v=1):
        z = self.zero()
        if i + j <= self.N:
            z.a[i, j] = self.s(v)
        return z


class MS:
    __slots__ = ('ctx', 'a')

    def __init__(self, ctx, a):
        self.ctx = ctx
        self.a = a

    def _lift(self, o):
        return o if isinstance(o, MS) else self.ctx.const(o)

    def __add__(self, o):
        return MS(self.ctx, (self.a + self._lift(o).a) % self.ctx.p)

    __radd__ = __add__

    def __sub__(self, o):
        return MS(self.ctx, (self.a - self._lift(o).a) % self.ctx.p)

    def __rsub__(self, o):
        return self._lift(o) - self

    def __neg__(self):
        return MS(self.ctx, (-self.a) % self.ctx.p)

    def __mul__(self, o):
        ctx = self.ctx
        N, p = ctx.N, ctx.p
        if not isinstance(o, MS):
            return MS(ctx, self.a * ctx.s(o) % p)
        out = mulmod(self.a, o.a, p)
        out[~ctx.mask] = 0
        return MS(ctx, out)

    __rmul__ = __mul__

    def inv(self):
        """1/self via Newton iteration; requires a nonzero constant term."""
        ctx = self.ctx
        p = ctx.p
        c0 = int(self.a[0, 0])
        g = ctx.const(pow(c0, p - 2, p))
        prec = 1
        while prec <= ctx.N:
            prec *= 2
            g = g * (2 - self * g)
        return g

    def scale(self, s, t):
        """f(x^s, y^t)."""
        ctx = self.ctx
        out = ctx.zero()
        for i, j in zip(*np.nonzero(self.a)):
            if s * i + t * j <= ctx.N:
                out.a[s * i, t * j] = self.a[i, j]
        return out

    def shift(self, s, t):
        """x^s y^t f (s, t may be negative if the low part vanishes)."""
        ctx = self.ctx
        N = ctx.N
        out = np.zeros_like(self.a)
        src = self.a
        i0, j0 = max(0, -s), max(0, -t)
        i1, j1 = min(N + 1, N + 1 - s), min(N + 1, N + 1 - t)
        if i0 < i1 and j0 < j1:
            out[i0 + s:i1 + s, j0 + t:j1 + t] = src[i0:i1, j0:j1]
        out[~ctx.mask] = 0
        return MS(ctx, out)

    def __getitem__(self, ij):
        i, j = ij
        if i < 0 or j < 0 or i + j > self.ctx.N:
            return 0
        return int(self.a[i, j])

    def __setitem__(self, ij, v):
        self.a[ij] = self.ctx.s(v)

    def items(self):
        for i, j in zip(*np.nonzero(self.a)):
            yield int(i), int(j), int(self.a[i, j])

    def iszero(self):
        return not self.a.any()

    def __eq__(self, o):
        return np.array_equal(self.a, o.a)


def crt(residues, primes):
    """Symmetric CRT lift of residues to the integer of least absolute value."""
    x, m = 0, 1
    for r, p in zip(residues, primes):
        t = (r - x) * pow(m, -1, p) % p
        x += m * t
        m *= p
    return x - m if x > m // 2 else x


def frac_mod(fr, p):
    return Fraction(fr).numerator % p * pow(Fraction(fr).denominator % p, p - 2, p) % p


def no_throttle():
    """Opt this process out of Windows EcoQoS (otherwise background jobs land on the E-cores)."""
    import ctypes
    import sys
    if sys.platform != 'win32':
        return

    class PPTS(ctypes.Structure):
        _fields_ = [('Version', ctypes.c_ulong), ('ControlMask', ctypes.c_ulong), ('StateMask', ctypes.c_ulong)]

    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    k32.SetProcessInformation.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_ulong]
    k32.SetProcessInformation.restype = ctypes.c_int
    s = PPTS(1, 1, 0)  # PROCESS_POWER_THROTTLING_EXECUTION_SPEED, state off
    ok = k32.SetProcessInformation(k32.GetCurrentProcess(), 4, ctypes.byref(s), ctypes.sizeof(s))
    if not ok:
        print('warning: SetProcessInformation failed', k32.GetLastError(), file=sys.stderr)
