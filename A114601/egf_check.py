"""Check the closed-form EGF of proofs.md, section 5, against c(k) from formula.py.

C_A(x) = (t - t^2/2 - 2x)/(4x) + x/2,  C_D(x) = U(2x)/2 + t^2/8 - x^2/2 - 16x^3/3! - 208x^4/4!,
t = T(2x), U(z) = -log(1-T(z))/2 - T(z)/2 - T(z)^2/4, T the tree function.
Exact power series with rational coefficients.

    python attacks/A114601/egf_check.py [K]
"""
import sys
from fractions import Fraction
from math import factorial

sys.path.insert(0, __import__("os").path.dirname(__file__))
from formula import c_A, c_D, connected  # noqa: E402

K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
N = K + 2


def mul(a, b):
    out = [Fraction(0)] * N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[: N - i]):
                out[i + j] += x * y
    return out


def scale_arg(a, c):  # a(c z)
    return [x * Fraction(c) ** i for i, x in enumerate(a)]


T = [Fraction(0)] + [Fraction(m ** (m - 1), factorial(m)) for m in range(1, N)]
# -log(1 - T) = sum_j T^j / j
neglog = [Fraction(0)] * N
p = [Fraction(1)] + [Fraction(0)] * (N - 1)
for j in range(1, N):
    p = mul(p, T)
    neglog = [x + y / j for x, y in zip(neglog, p)]
T2 = mul(T, T)
U = [x / 2 - y / 2 - z / 4 for x, y, z in zip(neglog, T, T2)]

t = scale_arg(T, 2)
t2 = mul(t, t)
num = [a - b / 2 for a, b in zip(t, t2)]
num[1] -= 2  # - 2x
# divide by 4x: shift down by one
CA = [num[i + 1] / 4 for i in range(N - 1)] + [Fraction(0)]
CA[1] += Fraction(1, 2)
CD = [a / 2 + b / 8 for a, b in zip(scale_arg(U, 2), t2)]
CD[2] -= Fraction(1, 2)
CD[3] -= Fraction(16, factorial(3))
CD[4] -= Fraction(208, factorial(4))

bad = 0
for k in range(1, K + 1):
    ca, cd = CA[k] * factorial(k), CD[k] * factorial(k)
    ok = ca == c_A(k) and cd == c_D(k)
    bad += not ok
    if k <= 10 or not ok:
        print(f"k={k}: EGF gives c_A = {ca}, c_D = {cd};  formula c_A = {c_A(k)}, c_D = {c_D(k)}  {'OK' if ok else 'MISMATCH'}")
print(f"closed-form EGF vs formula for k = 1..{K}: {'all match' if bad == 0 else f'{bad} mismatches'}")
