"""Check motifs produced by the DP (rust enum/sample output) with the exact G^k test in motif.py.

Each input line gives, per cell (row-major from the bottom), the edge its tile owns. The path is walked from
a = (0, 0) to give every tile its direction (parity); then we check: the walk is one path a -> b using all n^2
tiles, the perimeter matches the frozen boundary including parities, and G is globally self-avoiding.

usage: python verify.py N MODE < motifs.txt
"""
import sys

from motif import as_tuple, boundary, globally_self_avoiding, segment

EDGE = {"B": 0, "R": 1, "T": 2, "L": 3}


def orient(n, word):
    """Edge letters -> tile types (with parity) by walking the path from a."""
    owner = {}
    for idx, ch in enumerate(word):
        x, y = idx % n, idx // n
        s, e = segment(x, y, EDGE[ch])          # forward direction; parity fixed below
        owner[(min(s, e), max(s, e))] = (x, y, EDGE[ch])
    adj = {}
    for k in owner:
        for v in k:
            adj.setdefault(v, []).append(k)
    grid, prev, v = {}, None, (0, 0)
    while True:
        nxt = [k for k in adj[v] if k != prev]
        if not nxt:
            break
        assert len(nxt) == 1, "branching"
        k = nxt[0]
        x, y, o = owner[k]
        w = k[1] if k[0] == v else k[0]
        f = 0 if segment(x, y, o) == (v, w) else 1
        grid[(x, y)] = o + 4 * f
        prev, v = k, w
    assert v == (n, 0) and len(grid) == n * n, "not a single a-b path through all tiles"
    return grid


def main():
    n, mode = int(sys.argv[1]), sys.argv[2]
    bd = boundary(n, mode)
    good = bad = 0
    for line in sys.stdin:
        word = line.strip()
        if not word:
            continue
        grid = orient(n, word)
        ok = all(grid[c] == t for c, t in bd.items()) and globally_self_avoiding(*as_tuple(n, grid))
        good += ok
        bad += not ok
    print(f"n={n} mode={mode}: {good} globally self-avoiding, {bad} not")


if __name__ == "__main__":
    main()
