"""A001930 (unlabeled topologies) and A091073 (unlabeled transitive relations) from the Burnside
fixed-point tables of attacks/A000112.

A preorder on n points is a poset Q on its k equivalence classes plus a class size w(q) >= 1 for each q.
Unlabeled preorders are S_k-orbits of pairs (Q labeled on [k], w), so by Burnside

    a(n) = sum_{k<=n} sum_{lam |- k} fix_k(lam)/z_lam * [x^n] prod_{parts l of lam} c(x^l),

with fix_k(lam) = #posets on [k] invariant under a permutation of cycle type lam (attacks/A000112) and
c(x) = x/(1-x): a weighting fixed by the permutation is constant on cycles.
Transitive relations: a class is an irreflexive point or a reflexive class of size >= 1, c(x) = x + x/(1-x).
Posets: c(x) = x.  Labeled versions (c -> exponential) are checked against A000798 / A006905.

Usage: python topologies.py [N]   (default 17)
"""
import sys, os
from fractions import Fraction
from math import gcd, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "A000112")

A001035 = [1, 1, 3, 19, 219, 4231, 130023, 6129859, 431723379, 44511042511, 6611065248783, 1396281677105899,
           414864951055853499, 171850728381587059351, 98484324257128207032183, 77567171020440688353049939,
           83480529785490157813844256579, 122152541250295322862941281269151,
           241939392597201176602897820148085023]
A000112 = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284, 46749427, 1104891746, 33823827452,
           1338193159771, 68275077901156, 4483130665195087, 378180373015320248]
A001930 = [1, 1, 3, 9, 33, 139, 718, 4535, 35979, 363083, 4717687, 79501654, 1744252509, 49872339897,
           1856792610995, 89847422244493, 5637294117525695]
A091073 = [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036, 93810648902, 4461086120602,
           274339212258846, 21775814889230580, 2226876304576948549]
A000798 = [1, 1, 4, 29, 355, 6942, 209527, 9535241, 642779354, 63260289423, 8977053873043]
A006905 = [1, 2, 13, 171, 3994, 154303, 9415189, 878222530, 122207703623, 24890747921947, 7307450299510288]


def load_fix():
    tab = {}
    for name in ("fix_small.txt", "stats.txt"):
        for line in open(os.path.join(SRC, name)):
            f = line.split()
            if not f or f[0].startswith("#"):
                continue
            key = (int(f[0]), tuple(sorted(map(int, f[1].split(",")), reverse=True)))
            v = int(f[2])
            assert tab.get(key, v) == v, key
            tab[key] = v
    return tab


def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for p in partitions(n - k, k):
            yield (k,) + p


def interacting(lam):
    nt = [c for c in lam if c > 1]
    return tuple(sorted((a for i, a in enumerate(nt) if any(j != i and gcd(a, b) > 1 for j, b in enumerate(nt))),
                        reverse=True))


def z(lam):
    r = 1
    for c in set(lam):
        m = lam.count(c)
        r *= c ** m * factorial(m)
    return r


def fix(lam, tab):
    t = interacting(lam)
    return A001035[len(lam)] if not t else tab[(len(lam), t)]


def mul(a, b, N):
    c = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                c[i + j] += x * b[j]
    return c


def cycle_series(l, kind, N):
    """c(x^l) truncated at x^N."""
    s = [0] * (N + 1)
    for m in range(l, N + 1, l):
        s[m] = 1 if kind != "poset" else 0
    if l <= N:
        s[l] += {"poset": 1, "top": 0, "trans": 1}[kind]
    return s


def count(kind, N, tab):
    """a(0..N); also the triangle T[n][k] by number of classes k."""
    tot = [[Fraction(0)] * (N + 1) for _ in range(N + 1)]  # tot[k][n]
    tot[0][0] = Fraction(1)
    cs = {l: cycle_series(l, kind, N) for l in range(1, N + 1)}
    for k in range(1, N + 1):
        for lam in partitions(k):
            w = [1] + [0] * N
            for l in lam:
                w = mul(w, cs[l], N)
            f = Fraction(fix(lam, tab), z(lam))
            for n in range(k, N + 1):
                if w[n]:
                    tot[k][n] += f * w[n]
    tri = [[tot[k][n] for k in range(n + 1)] for n in range(N + 1)]
    for row in tri:
        assert all(x.denominator == 1 for x in row)
    tri = [[int(x) for x in row] for row in tri]
    return [sum(r) for r in tri], tri


def labeled(kind, N):
    """Labeled check: sum_k P(k)/k! * n![x^n] c(x)^k with the exponential class series."""
    out = []
    for n in range(N + 1):
        c = [Fraction(0)] + [Fraction({"top": 1, "trans": 2 if m == 1 else 1}[kind], factorial(m)) for m in range(1, n + 1)]
        p = [Fraction(1)] + [Fraction(0)] * n
        s = Fraction(0)
        for k in range(n + 1):
            s += Fraction(A001035[k], factorial(k)) * p[n]
            p = mul(p, c, n)
        out.append(int(s * factorial(n)))
    return out


def inv_euler(a):
    """b with a = Euler transform of b (a[0] = 1)."""
    n = len(a) - 1
    c = [0] * (n + 1)
    for m in range(1, n + 1):
        c[m] = m * a[m] - sum(c[k] * a[m - k] for k in range(1, m))
    b = [0] * (n + 1)
    for m in range(1, n + 1):
        b[m] = c[m] - sum(d * b[d] for d in range(1, m) if m % d == 0)
        assert b[m] % m == 0
        b[m] //= m
    return b


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 17
    tab = load_fix()
    assert labeled("top", 10) == A000798 and labeled("trans", 10) == A006905
    pos, _ = count("poset", N, tab)
    assert pos == A000112[: N + 1], pos
    top, ttri = count("top", N, tab)
    tra, rtri = count("trans", N, tab)
    for name, got, ref in (("A001930", top, A001930), ("A091073", tra, A091073)):
        m = min(len(ref), N + 1)
        assert got[:m] == ref[:m], (name, got)
        print(f"{name}: reproduces a(0)-a({m - 1})")
        for n in range(m, N + 1):
            print(f"  NEW {name}({n}) = {got[n]}")
    print("A001930:", top)
    print("A091073:", tra)
    print("connected topologies (A001928):", [1] + inv_euler(top)[1:])
    print("connected transitive relations:", [1] + inv_euler(tra)[1:])
    print("topologies with up to n points (A326898):", [sum(top[: i + 1]) for i in range(N + 1)])
    print("topologies by number of classes (rows n, k=0..n):")
    for r in ttri[: min(N, 9) + 1]:
        print("  ", r)
    import json
    json.dump({"A001930": top, "A091073": tra, "top_by_classes": ttri, "trans_by_classes": rtri},
              open(os.path.join(HERE, "results.json"), "w"))
