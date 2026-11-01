"""Brute-force count of fully symmetric recurrent sandpiles on the n x n grid.

Every cell has degree 4 (boundary cells send their missing edges to the sink).
A symmetric configuration is constant on D4 orbits, so enumerate orbit values
in 0..3 and test recurrence with Dhar's burning algorithm: c is recurrent iff
burning from the sink reaches every cell, where a cell burns once c(v) is at
least the number of its unburnt grid neighbours.
"""
import itertools
import sys


def orbits(n):
    rep = {}
    for i in range(n):
        for j in range(n):
            imgs = []
            for a, b in ((i, j), (j, i)):
                for x, y in ((a, b), (n - 1 - a, b), (a, n - 1 - b), (n - 1 - a, n - 1 - b)):
                    imgs.append((x, y))
            rep[(i, j)] = min(imgs)
    reps = sorted(set(rep.values()))
    idx = {r: t for t, r in enumerate(reps)}
    return [[idx[rep[(i, j)]] for j in range(n)] for i in range(n)], len(reps)


def recurrent(c, n):
    nb = [[[(x, y) for x, y in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
            if 0 <= x < n and 0 <= y < n] for j in range(n)] for i in range(n)]
    unburnt = [[len(nb[i][j]) for j in range(n)] for i in range(n)]
    burnt = [[False] * n for _ in range(n)]
    stack = [(i, j) for i in range(n) for j in range(n) if c[i][j] >= unburnt[i][j]]
    count = 0
    while stack:
        i, j = stack.pop()
        if burnt[i][j]:
            continue
        burnt[i][j] = True
        count += 1
        for x, y in nb[i][j]:
            if not burnt[x][y]:
                unburnt[x][y] -= 1
                if c[x][y] >= unburnt[x][y]:
                    stack.append((x, y))
    return count == n * n


def count(n):
    orb, m = orbits(n)
    total = 0
    for vals in itertools.product(range(4), repeat=m):
        c = [[vals[orb[i][j]] for j in range(n)] for i in range(n)]
        total += recurrent(c, n)
    return total


if __name__ == "__main__":
    for n in range(1, int(sys.argv[1]) + 1):
        print(n, count(n), flush=True)
