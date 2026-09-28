"""Brute-force check of the split identity for A006455 (naturally labeled posets).

  a(n+k) = sum over natural Q on [n], natural R on [k] of i(R^op x Q)

where i() counts order ideals (down-sets). Also checks the one-step form
a(n+1) = sum_Q i(Q). Posets are stored as tuples of down-set bitmasks (strict).
"""
import sys
from functools import lru_cache

KNOWN = [1, 1, 2, 7, 40, 357, 4824, 96428, 2800472, 116473461]


@lru_cache(None)
def natural(n):
    """All naturally labeled posets on [n]: element j's strict down-set is an ideal of [j]."""
    if n == 0:
        return [()]
    out = []
    for q in natural(n - 1):
        for I in ideals(q):
            out.append(q + (I,))
    return out


@lru_cache(None)
def ideals(q):
    """Order ideals of poset q (tuple of strict down-set masks) as bitmasks."""
    res = [0]
    for j, d in enumerate(q):  # j's down-set only uses smaller labels
        res += [m | (1 << j) for m in res if m & d == d]
    return tuple(res)


def hom_count(r, q):
    """#order-preserving maps R -> J(Q) = #ideals of R^op x Q."""
    J = ideals(q)
    # R natural: assign phi(r_j) in label order; need phi(r_i) subset phi(r_j) for i in down(r_j)
    def rec(j, phi):
        if j == len(r):
            return 1
        lo = 0
        for i in range(j):
            if r[j] >> i & 1:
                lo |= phi[i]
        return sum(rec(j + 1, phi + (I,)) for I in J if I & lo == lo)
    return rec(0, ())


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for m in range(N + 1):
        assert len(natural(m)) == KNOWN[m], m
        if m:
            assert sum(len(ideals(q)) for q in natural(m - 1)) == KNOWN[m]
        for n in range(m + 1):
            k = m - n
            s = sum(hom_count(r, q) for q in natural(n) for r in natural(k))
            assert s == KNOWN[m], (n, k, s)
        print(f"a({m}) = {KNOWN[m]}: all splits n+k={m} agree")


if __name__ == "__main__":
    main()
