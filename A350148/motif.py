"""Hilbert-style square-filling motifs (McKenna, ECA 2:2 (2022) #S2R13; OEIS A350148).

Conventions (McKenna, Def. 1.1-1.8). A tile at cell (x, y) of the n x n grid (lattice points 0..n) has one
distinguished edge: orientation o in B, R, T, L (0..3) and parity f (0 = Forward = counterclockwise around the tile,
1 = Backward). A motif is a thread of n^2 tiles whose segments form a path from a = (0, 0) to b = (n, 0).

Composition J o H puts into every tile of J a copy of H rotated by the tile's orientation (k * 90 degrees ccw)
and, for a Backward tile, reversed (same shape, all directions and parities flipped).

`globally_self_avoiding(G)` decides "G^k is locally self-avoiding for every k" exactly: every self-contact in G^k
lies inside some 2 x 2 block of level-1 supertiles, and the 2 x 2 supertile patterns that occur at any level are the
closure of G's own interior-vertex patterns under substitution (finite: at most 8^4).
"""
from functools import lru_cache

B, R, T, L = 0, 1, 2, 3
NAMES = "BRTLbrtl"            # type = o + 4*f


def name_to_type(c):
    return NAMES.index(c)


def segment(x, y, t):
    """Directed segment (start, end) of the tile of type t at cell (x, y)."""
    o, f = t % 4, t // 4
    s, e = [((x, y), (x + 1, y)), ((x + 1, y), (x + 1, y + 1)),
             ((x + 1, y + 1), (x, y + 1)), ((x, y + 1), (x, y))][o]
    return (e, s) if f else (s, e)


def from_thread(word):
    """Build the cell -> type grid from a thread word such as 'llBLtlRb...' starting at a = (0, 0)."""
    n = round(len(word) ** 0.5)
    assert n * n == len(word)
    grid = {}
    p = (0, 0)
    for c in word:
        t = name_to_type(c)
        o, f = t % 4, t // 4
        # the cell is fixed by the segment's start point, orientation and direction
        cands = [(x, y) for x in (p[0] - 1, p[0]) for y in (p[1] - 1, p[1]) if segment(x, y, t)[0] == p]
        assert len(cands) == 1, (c, p, cands)
        cell = cands[0]
        assert cell not in grid and 0 <= cell[0] < n and 0 <= cell[1] < n, (c, cell)
        grid[cell] = t
        p = segment(*cell, t)[1]
    assert p == (n, 0)
    return n, grid


def as_tuple(n, grid):
    return n, tuple(grid[(x, y)] for y in range(n) for x in range(n))


def transform(n, G, t):
    """G (tuple, row-major from the bottom) placed into a tile of type t: n x n tuple of sub-tile types."""
    o, f = t % 4, t // 4
    cur = list(G)
    for _ in range(o):                        # rotate 90 degrees ccw: (x, y) -> (n-1-y, x), orientation + 1
        nxt = [0] * (n * n)
        for y in range(n):
            for x in range(n):
                s = cur[y * n + x]
                nx, ny = n - 1 - y, x
                nxt[ny * n + nx] = (s % 4 + 1) % 4 + 4 * (s // 4)
        cur = nxt
    if f:
        cur = [s % 4 + 4 * (1 - s // 4) for s in cur]
    return tuple(cur)


def locally_self_avoiding(cells):
    """cells: iterable of (x, y, type). No lattice edge used twice; at most one in and one out per vertex."""
    edges, outs, ins = set(), set(), set()
    for x, y, t in cells:
        s, e = segment(x, y, t)
        k = (min(s, e), max(s, e))
        if k in edges or s in outs or e in ins:
            return False
        edges.add(k)
        outs.add(s)
        ins.add(e)
    return True


def globally_self_avoiding(n, G):
    """Exact test that G^k is locally self-avoiding for all k >= 1 (G a tuple from as_tuple)."""
    subs = [transform(n, G, t) for t in range(8)]

    def block(p):
        cells = []
        for q, (ox, oy) in enumerate([(0, 0), (n, 0), (0, n), (n, n)]):
            S = subs[p[q]]
            for y in range(n):
                for x in range(n):
                    cells.append((ox + x, oy + y, S[y * n + x]))
        return cells

    def patches(width, get):
        out = set()
        for j in range(1, width):
            for i in range(1, width):
                out.add((get(i - 1, j - 1), get(i, j - 1), get(i - 1, j), get(i, j)))
        return out

    if not locally_self_avoiding((x, y, G[y * n + x]) for y in range(n) for x in range(n)):
        return False
    todo = list(patches(n, lambda x, y: G[y * n + x]))
    seen = set(todo)
    while todo:
        p = todo.pop()
        cells = block(p)
        if not locally_self_avoiding(cells):
            return False
        grid = {(x, y): t for x, y, t in cells}
        for q in patches(2 * n, lambda x, y: grid[(x, y)]):
            if q not in seen:
                seen.add(q)
                todo.append(q)
    return True


def mode_bits(n, mode):
    """Palindromic zipper bits (length n) from mode bits (the first (n+1)/2, as a string like '00101')."""
    m = mode
    assert len(m) == (n + 1) // 2 and m[-1] == "1"
    return m + m[-2::-1]


def modes(n):
    """All Fibbinary zipping modes for odd n >= 5 (McKenna Thm 3.1): mode bits '00' + fibbinary + '01'."""
    if n == 1:
        return ["1"]
    h = (n + 1) // 2
    out = []
    for v in range(1 << h):
        s = format(v, f"0{h}b")
        if s.startswith("00") and s.endswith("1") and "11" not in s:
            full = s + s[-2::-1]
            if "11" not in full and full[-2:] == "00" and s[-2] == "0":
                out.append(s)
    return sorted(out)


def boundary(n, mode):
    """Frozen perimeter tiles of a left-handed motif with the given zipping mode.

    Each side, walked in the path's flow direction (left side up, top right, right side down: clockwise;
    bottom left to right: counterclockwise), starts on the outer track and swaps track at every 1 bit; the tile
    with a 1 carries the perpendicular segment. The left side owns both left corners, the top owns the top right
    corner and the bottom owns the bottom right corner. Checked against the E motif and by the counts.
    """
    z = mode_bits(n, mode)
    g = {}
    # (cells in flow order, outer, inner, switch-in, switch-out), types as names
    sides = [
        ([(0, y) for y in range(n)], "l", "R", "B", "T"),                     # left, upward
        ([(x, n - 1) for x in range(n)], "t", "B", "L", "R"),                 # top, rightward
    ]
    for cells, outer, inner, sin, sout in sides:
        on_outer = True
        for k, c in enumerate(cells):
            if z[k] == "1":
                t = sin if on_outer else sout
                on_outer = not on_outer
            else:
                t = outer if on_outer else inner
            g.setdefault(c, name_to_type(t))       # first side to claim a corner keeps it
    # bottom, rightward (counterclockwise): ends on the outer track; walk it backwards from b.
    on_outer = True
    for k, x in enumerate(range(n - 1, -1, -1)):
        if z[k] == "1":
            t = "l" if not on_outer else "r"           # walking backwards, swap roles of in/out
            on_outer = not on_outer
        else:
            t = "B" if on_outer else "t"
        if (x, 0) not in g:
            g[(x, 0)] = name_to_type(t)
    # right, downward; both its corners are already owned
    on_outer = True
    for k, y in enumerate(range(n - 1, -1, -1)):
        if z[k] == "1":
            t = "T" if on_outer else "B"
            on_outer = not on_outer
        else:
            t = "r" if on_outer else "L"
        g.setdefault((n - 1, y), name_to_type(t))
    return g


E_WORD = "llBLtlRbRRtLtBBrTbrBLbrBB"   # McKenna's E motif for T_5 (Fig. 12)


if __name__ == "__main__":
    n, grid = from_thread(E_WORD)
    G = as_tuple(n, grid)[1]
    print("E globally self-avoiding:", globally_self_avoiding(n, G))
    bd = boundary(5, "001")
    bad = {c: (NAMES[t], NAMES[grid[c]]) for c, t in bd.items() if grid[c] != t}
    print("E boundary matches rule:", not bad, bad, len(bd))
    n3, g3 = from_thread("lBlRtLLrB")             # McKenna's J (Fig. 9): locally but not globally self-avoiding
    print("J globally self-avoiding (expect False):", globally_self_avoiding(n3, as_tuple(n3, g3)[1]))
    for m in [5, 7, 9, 11, 13, 15]:
        print(m, modes(m))
