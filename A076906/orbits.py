"""Brute-force check of the undirected / up-to-symmetry relatives (A076875, A077550,
A077556) from explicit curves.

A curve is determined by its start quadrant and its crossing sequence
((ray, rank), ...), rank = position along the ray counted from O: consecutive
arcs alternate sides at each crossing, so every arc's quadrant follows.
Rays W,S,E,N = 0..3 counterclockwise; quadrants Q3,Q4,Q1,Q2 = 0..3 (quadrant q
lies between rays q and q+1).
"""
import sys
from itertools import product
from brute import ncm

RAYS = range(4)


def curves(n, sq, eq):
    """Directed curves with n crossings from quadrant sq to quadrant eq."""
    for ks in product(range(n + 1), repeat=4):
        if sum(ks) != n:
            continue
        seqs = []
        for q in range(4):
            a, b = q, (q + 1) % 4
            rim = (["s"] if q == sq else []) + (["e"] if q == eq else [])
            seqs.append(([(a, i) for i in range(ks[a])], rim, [(b, i) for i in reversed(range(ks[b]))]))
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
                path, prev, cur = [], None, "s"
                while cur != "e":
                    nb = adj[cur]
                    nxt = nb[0] if prev is None or nb[0] != prev or len(nb) == 1 else nb[1]
                    if len(nb) == 2 and nb[0] == nb[1]:
                        nxt = nb[0]
                    prev, cur = cur, nxt
                    if cur != "e":
                        path.append(cur)
                if len(path) == n:
                    yield (sq, eq, tuple(path))


def sym(c, g):
    """g = (rot, refl): optional reflection in the x axis, then rotation."""
    sq, eq, path = c
    rot, refl = g
    fr = (lambda r: (-r) % 4) if refl else (lambda r: r)
    fq = (lambda q: (-q - 1) % 4) if refl else (lambda q: q)
    R = lambda r: (fr(r) + rot) % 4
    Q = lambda q: (fq(q) + rot) % 4
    return (Q(sq), Q(eq), tuple((R(r), i) for r, i in path))


def rev(c):
    sq, eq, path = c
    return (eq, sq, path[::-1])


def first_axis(c):
    return "x" if c[2][0][0] in (0, 2) else "y"


if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    G = [(r, f) for r in range(4) for f in (0, 1)]
    for n in range(1, M + 1):
        allc = [c for sq in range(4) for eq in range(4) for c in curves(n, sq, eq)]
        und = {min(c, rev(c)) for c in allc}
        orbits = {min(min(sym(c, g), rev(sym(c, g))) for g in G) for c in allc}
        # A076875: undirected curves with an end in Q3
        a875 = sum(1 for u in und if 0 in (u[0], u[1]))
        # A077556 candidates: undirected curves with an end in Q3 whose path, read
        # from (some) Q3 end, crosses the x axis first
        def x_from_q3(u):
            reads = [c for c in (u, rev(u)) if c[0] == 0]
            return any(first_axis(c) == "x" for c in reads)
        a556 = sum(1 for u in und if 0 in (u[0], u[1]) and x_from_q3(u))
        print(n, "A076875", a875, "A077550", len(orbits), "A077556?", a556, flush=True)
