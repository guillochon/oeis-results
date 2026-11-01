"""Exact a(n) as det of the D4-folded reduced Laplacian (Bareiss, integers).

Rows/columns are D4 orbits of cells.  Lf[O, O'] = sum over cells w in O' of
L[v, w] for a fixed v in O, i.e. the grid Laplacian (degree 4, sink-adjacent
boundary) acting on functions constant on orbits.
"""
import sys

from brute import orbits


def folded(n):
    orb, m = orbits(n)
    rep = {}
    for i in range(n):
        for j in range(n):
            rep.setdefault(orb[i][j], (i, j))
    L = [[0] * m for _ in range(m)]
    for o, (i, j) in rep.items():
        L[o][o] += 4
        for x, y in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
            if 0 <= x < n and 0 <= y < n:
                L[o][orb[x][y]] -= 1
    return L


def bareiss(M):
    M = [row[:] for row in M]
    m, sign, prev = len(M), 1, 1
    for k in range(m - 1):
        if M[k][k] == 0:
            for r in range(k + 1, m):
                if M[r][k]:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, m):
            for j in range(k + 1, m):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1] if m else 1


if __name__ == "__main__":
    for n in range(0, int(sys.argv[1]) + 1):
        print(n, bareiss(folded(n)) if n else 1, flush=True)
