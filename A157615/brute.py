"""Exhaustive search for A157615: longest self-avoiding path on an n X n board whose unit steps
alternate horizontal / vertical. Prints a(n) and one optimal path. Feasible for n <= 6 or so.
Usage: python brute.py [max_n]
"""
import sys

sys.setrecursionlimit(10000)


def solve(n):
    best = [0, None]
    seen = [[False] * n for _ in range(n)]
    path = []

    def dfs(x, y, horiz):  # horiz: the next step must be horizontal
        if len(path) > best[0]:
            best[0], best[1] = len(path), list(path)
        moves = [(1, 0), (-1, 0)] if horiz else [(0, 1), (0, -1)]
        for dx, dy in moves:
            u, v = x + dx, y + dy
            if 0 <= u < n and 0 <= v < n and not seen[u][v]:
                seen[u][v] = True; path.append((u, v))
                dfs(u, v, not horiz)
                seen[u][v] = False; path.pop()

    # symmetry: the two mirror flips keep steps horizontal/vertical, so the start can be taken in
    # the top-left quadrant; the transpose swaps H and V, so the first step can be taken horizontal
    for x in range((n + 1) // 2):
        for y in range((n + 1) // 2):
            seen[x][y] = True; path.append((x, y))
            dfs(x, y, True)
            seen[x][y] = False; path.pop()
    return best


def draw(n, p):
    g = [["." for _ in range(n)] for _ in range(n)]
    for i, (x, y) in enumerate(p):
        g[y][x] = "S" if i == 0 else "E" if i == len(p) - 1 else "#"
    return "\n".join(" ".join(r) for r in g)


if __name__ == "__main__":
    known = [1, 4, 7, 14, 19, 32, 39, 58, 67]
    for n in range(1, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 6):
        L, p = solve(n)
        print(f"n={n} a(n)={L} OEIS={known[n-1]} {'ok' if L == known[n-1] else 'MISMATCH'}")
        print(draw(n, p), "\n")
