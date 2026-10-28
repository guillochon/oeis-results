"""Soundness of the DP count for a zipping mode, independent of the interior of G.

Claim: if every locally plausible 2x2 patch of tile types is safe for the mode's frozen boundary, then every motif
counted by the DP (a single a-b path with that boundary) is globally self-avoiding.

Induction: G is locally self-avoiding (DP). If G^k is, then G^(k+1) = G^k o G can only self-contact inside the
2x2 block of copies of G around some interior vertex of G^k; that 2x2 patch of G^k is locally plausible (below),
hence safe, and a single copy never self-contacts. Copies only touch along supertile edges, which only perimeter
tiles reach, so safety depends on the boundary alone.

Locally plausible quad (true of every 2x2 patch of a self-avoiding Hilbert-style path, McKenna Lemma 3.6):
the four segments are distinct edges, no lattice point has two arrivals or two departures, and two edge-adjacent
tiles have opposite parity iff one of them owns their shared edge (a tile's parity says which side of the path
it lies on).

usage: python sound.py N [N ...]
"""
import sys
from itertools import product

from motif import boundary, locally_self_avoiding, modes, segment, transform


def plausible(q):
    cells = [(0, 0, q[0]), (1, 0, q[1]), (0, 1, q[2]), (1, 1, q[3])]
    if not locally_self_avoiding(cells):
        return False
    edge = {(x, y): (lambda s, e: (min(s, e), max(s, e)))(*segment(x, y, t)) for x, y, t in cells}
    par = {(x, y): t // 4 for x, y, t in cells}
    for a, b, shared in [((0, 0), (1, 0), ((1, 0), (1, 1))), ((0, 1), (1, 1), ((1, 1), (1, 2))),
                         ((0, 0), (0, 1), ((0, 1), (1, 1))), ((1, 0), (1, 1), ((1, 1), (2, 1)))]:
        owned = edge[a] == shared or edge[b] == shared
        if (par[a] != par[b]) != owned:
            return False
    return True


def safe(n, bd, q):
    """Do four copies of the boundary, placed as the quad q, avoid each other?"""
    cells = []
    G = [0] * (n * n)
    for (x, y), t in bd.items():
        G[y * n + x] = t
    for k, (ox, oy) in enumerate([(0, 0), (n, 0), (0, n), (n, n)]):
        S = transform(n, tuple(G), q[k])
        # transform moves cells around; map the perimeter cells through it by transforming a marker grid
        for y in range(n):
            for x in range(n):
                if x in (0, n - 1) or y in (0, n - 1):          # the perimeter maps onto the perimeter
                    cells.append((ox + x, oy + y, S[y * n + x]))
    return locally_self_avoiding(cells)


def main():
    quads = [q for q in product(range(8), repeat=4) if plausible(q)]
    print(f"{len(quads)} locally plausible quads of 4096")
    for n in map(int, sys.argv[1:]):
        for m in modes(n):
            bd = boundary(n, m)
            unsafe = [q for q in quads if not safe(n, bd, q)]
            print(f"n={n} mode={m}: {len(unsafe)} unsafe" + (f", e.g. {unsafe[:3]}" if unsafe else ""))


if __name__ == "__main__":
    main()
