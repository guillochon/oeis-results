"""Prime posets (A354693) and friends from the cycle index series Z_P of the species of posets, degree <= N.

Z_P = sum_k sum_{lam |- k} fix_k(lam)/z_lam p_lam   (fix tables of attacks/A000112).

Modular decomposition: a poset on >= 2 points is exactly one of
  parallel  (disconnected: a set of >= 2 connected posets),
  series    (ordinal sum of >= 2 posets that are not ordinal sums),
  prime[P]  (a prime poset Q on >= 4 points with a poset substituted for each point).
So with C = connected posets (1 + P = E o C) and S = posets that are not ordinal sums (P = S/(1-S)),
  Pr o P = C + S - P - X,     Z_Pr = (Z_C + Z_S - Z_P - p_1) o Z_P^{<-1>}   (plethystic inverse).

Symmetric functions are dicts {partition (weakly decreasing tuple): Fraction} in the power-sum basis.
"""
import sys, json, os
from fractions import Fraction
from math import factorial
from topologies import load_fix, partitions, fix, z, A000112, inv_euler

N = int(sys.argv[1]) if len(sys.argv) > 1 else 17


def deg(lam):
    return sum(lam)


def add(a, b, c=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + c * v
        if r[k] == 0:
            del r[k]
    return r


def mul(a, b, n=N):
    """Product truncated at degree n."""
    r = {}
    bs = sorted(b.items(), key=lambda kv: deg(kv[0]))
    for la, va in a.items():
        da = deg(la)
        for lb, vb in bs:
            if da + deg(lb) > n:
                break
            k = tuple(sorted(la + lb, reverse=True))
            r[k] = r.get(k, 0) + va * vb
    return {k: v for k, v in r.items() if v != 0}


def pk(g, k, n=N):
    """p_k[g]: every p_j -> p_{jk}."""
    return {tuple(j * k for j in lam): v for lam, v in g.items() if deg(lam) * k <= n}


def plethysm(f, g, n=N):
    """f[g], g without constant term, truncated at degree n."""
    pkg = {}
    powers = {}

    def power(k, m):  # (p_k[g])^m
        if (k, m) not in powers:
            if k not in pkg:
                pkg[k] = pk(g, k, n)
            powers[(k, m)] = pkg[k] if m == 1 else mul(power(k, m - 1), pkg[k], n)
        return powers[(k, m)]

    # group the terms of f in a trie over (part, multiplicity) from the largest part down
    out = {}
    memo = {}

    def prod(key):  # key: tuple of (part, mult)
        if key not in memo:
            if len(key) == 1:
                memo[key] = power(*key[0])
            else:
                memo[key] = mul(prod(key[:-1]), power(*key[-1]), n)
        return memo[key]

    for lam, v in f.items():
        if not lam:
            out[()] = out.get((), 0) + v
            continue
        key = tuple((c, lam.count(c)) for c in sorted(set(lam), reverse=True))
        # cheap lower bound on the degree: g starts in degree >= 1
        if deg(lam) > n:
            continue
        out = add(out, prod(key), v)
    return out


def homog(f, d):
    return {k: v for k, v in f.items() if deg(k) == d}


def count_unlabeled(f, n=N):
    """coefficients of x^d after p_i -> x^i."""
    c = [Fraction(0)] * (n + 1)
    for lam, v in f.items():
        c[deg(lam)] += v
    assert all(x.denominator == 1 for x in c), c
    return [int(x) for x in c]


def count_labeled(f, n=N):
    c = [0] * (n + 1)
    for lam, v in f.items():
        if all(x == 1 for x in lam):
            c[len(lam)] = int(v * factorial(len(lam)))
    return c


def series_log1p(f, n=N):
    """log(1 + f) for f without constant term."""
    out, pw = {}, {(): Fraction(1)}
    for m in range(1, n + 1):
        pw = mul(pw, f, n)
        if not pw:
            break
        out = add(out, pw, Fraction((-1) ** (m + 1), m))
    return out


def mobius(k):
    r, p, m = 1, 2, k
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r


def connected(f, n=N):
    """Z_C with 1 + f = exp(sum_k p_k[Z_C]/k)."""
    L = series_log1p(f, n)
    out = {}
    for k in range(1, n + 1):
        if mobius(k):
            out = add(out, pk(L, k, n), Fraction(mobius(k), k))
    return out


def geometric_inverse(f, n=N):
    """S with f = S/(1-S), i.e. S = f/(1+f)."""
    out, pw = {}, {(): Fraction(1)}
    for m in range(1, n + 1):
        pw = mul(pw, f, n)
        if not pw:
            break
        out = add(out, pw, (-1) ** (m + 1))
    return out


def right_divide(T, g, n=N, lo=1):
    """F with F[g] = T (g = p_1 + higher), solved degree by degree."""
    F, R = {}, dict(T)
    for d in range(lo, n + 1):
        h = homog(R, d)
        if not h:
            continue
        F = add(F, h)
        R = add(R, plethysm(h, g, n), -1)
        assert not homog(R, d)
    assert not R
    return F


if __name__ == "__main__":
    tab = load_fix()
    ZP = {}
    for k in range(1, N + 1):
        for lam in partitions(k):
            ZP[lam] = Fraction(fix(lam, tab), z(lam))
    assert count_unlabeled(ZP)[1:] == A000112[1 : N + 1]
    ZC = connected(ZP)
    ZS = geometric_inverse(ZP)
    T = add(add(add(ZC, ZS), ZP, -1), {(1,): Fraction(1)}, -1)
    assert all(deg(k) >= 4 for k in T), sorted(T, key=deg)[:3]
    ZPr = right_divide(T, ZP, lo=4)
    res = {
        "connected": count_unlabeled(ZC),
        "not_ordinal_sum": count_unlabeled(ZS),
        "prime": count_unlabeled(ZPr),
        "prime_labeled": count_labeled(ZPr),
        "connected_labeled": count_labeled(ZC),
        "not_ordinal_sum_labeled": count_labeled(ZS),
    }
    for k, v in res.items():
        print(k, v)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "prime_results.json"), "w"))
    json.dump({"N": N, "ZPr": [[list(k), str(v)] for k, v in sorted(ZPr.items())]},
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zprime.json"), "w"))
