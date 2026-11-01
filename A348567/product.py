"""a(n) by the eigenvalue product, plus checks of the sqrt conjecture.

a(n) = prod_{1<=h<=k<=n, h,k odd} (4 - 2cos(h pi/(n+1)) - 2cos(k pi/(n+1))).
The D4-invariant subspace of the grid Laplacian is spanned by
f_hk(i,j) = s_h(i)s_k(j) + s_k(i)s_h(j), s_h(i) = sin(h i pi/(n+1)), h <= k odd;
det of the symmetrized reduced Laplacian (Florescu et al., Prop. 2.10) is the
product of its eigenvalues there.

Usage: python product.py N [bfile]   (cross-checks folded.py exactly for n <= 40)
"""
import sys

from mpmath import cos, log10, mp, nint, pi

from folded import bareiss, folded


def a(n):
    od = range(1, n + 1, 2)
    m = len(od)
    mp.dps = int(m * (m + 1) / 2 * 0.61) + 60
    c = [2 * cos(h * pi / (n + 1)) for h in od]
    P = mp.mpf(1)
    for x in range(m):
        for y in range(x, m):
            P *= 4 - c[x] - c[y]
    v = nint(P)
    assert abs(P - v) < mp.mpf(10) ** -30, n
    return int(v)


def b(n):
    """A348566(n, n): Z/2 x Z/2 symmetric recurrents (all ordered odd pairs)."""
    od = range(1, n + 1, 2)
    mp.dps = int(len(od) ** 2 * 0.61) + 60
    P = mp.mpf(1)
    for h in od:
        for k in od:
            P *= 4 - 2 * cos(h * pi / (n + 1)) - 2 * cos(k * pi / (n + 1))
    return int(nint(P))


def a030451(n):
    return n // 2 if n % 2 == 0 else n // 2 + 2


if __name__ == "__main__":
    N = int(sys.argv[1])
    vals = [1] + [a(n) for n in range(1, N + 1)]
    for n in range(1, min(N, 40) + 1):
        assert vals[n] == bareiss(folded(n)), n
    for n in range(N + 1):
        assert vals[n] ** 2 == b(n) * 2 ** a030451(n), n
    print("ok: folded det n<=%d, sqrt conjecture n<=%d" % (min(N, 40), N))
    for n in range(16):
        print(n, vals[n])
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as f:
            for n, v in enumerate(vals):
                f.write("%d %d\n" % (n, v))
        print("wrote", sys.argv[2], "max digits", len(str(vals[-1])))
