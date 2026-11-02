"""Brute-force reference counter for meanders crossing two perpendicular roads.

Model ("pizza"): the axes are four rays W, S, E, N from the origin; the river is a
simple curve from infinity in quadrant Q3 (-,-) to infinity in a chosen quadrant.
Its crossings on each ray are ordered by distance from O.  Inside each quadrant the
river pieces form a non-crossing perfect matching of the quadrant's boundary points
read in linear order: ray A outward, rim endpoints, ray B inward (the O corner holds
no point, so cyclic = linear).  An isotopy class is exactly a choice of crossing
counts and one matching per quadrant whose union is a single path start -> end.
"""
import sys
from functools import lru_cache
from itertools import product

RAYS = "WSEN"                      # counterclockwise from W
QUADS = [("W", "S"), ("S", "E"), ("E", "N"), ("N", "W")]   # Q3, Q4, Q1, Q2


@lru_cache(None)
def ncm(m):
    """All non-crossing perfect matchings of 0..m-1, as tuples of partner indices."""
    if m == 0:
        return [()]
    if m % 2:
        return []
    out = []
    for j in range(1, m, 2):           # 0 pairs with j
        for a in ncm(j - 1):
            for b in ncm(m - j - 1):
                p = [0] * m
                p[0], p[j] = j, 0
                for i, x in enumerate(a):
                    p[1 + i] = 1 + x
                for i, x in enumerate(b):
                    p[j + 1 + i] = j + 1 + x
                out.append(tuple(p))
    return out


def count(n, end_quad):
    """end_quad: 0..3 index into QUADS (0 = Q3, same as start)."""
    total = 0
    for ks in product(range(n + 1), repeat=4):
        if sum(ks) != n:
            continue
        k = dict(zip(RAYS, ks))
        seqs = []
        for qi, (a, b) in enumerate(QUADS):
            seq = [(a, i) for i in range(k[a])]
            rim = []
            if qi == 0:
                rim.append("start")
            if qi == end_quad:
                rim.append("end")
            seqs.append((seq, rim, [(b, i) for i in reversed(range(k[b]))]))
        # rim orderings: if start and end share Q3, both orders count
        rim_opts = [[r] if len(r) < 2 else [r, r[::-1]] for _, r, _ in seqs]
        for rims in product(*rim_opts):
            lines = [s + r + t for (s, _, t), r in zip(seqs, rims)]
            if any(len(l) % 2 for l in lines):
                continue
            for ms in product(*[ncm(len(l)) for l in lines]):
                adj = {}
                for l, m in zip(lines, ms):
                    for i, j in enumerate(m):
                        adj.setdefault(l[i], []).append(l[j])
                # walk from start
                prev, cur, seen = None, "start", 1
                while cur != "end":
                    nxt = [x for x in adj[cur] if x != prev] if prev is not None else adj[cur]
                    if len(adj[cur]) == 2 and adj[cur][0] == adj[cur][1]:
                        nxt = [adj[cur][0]]
                    prev, cur = cur, nxt[0]
                    seen += 1
                if seen == n + 2:
                    total += 1
    return total


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(N + 1):
        per = [count(n, q) for q in range(4)]
        a906 = per[2] if n % 2 == 0 else per[1]
        print(n, "A076906", a906, "A076907", sum(per), "per-quadrant(Q3,Q4,Q1,Q2)", per)
