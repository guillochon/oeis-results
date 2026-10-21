"""Rooted nonseparable maps and rooted c-nets (3-connected maps) by vertices and faces.

Conventions (series in x, y, truncated at total degree N):
  B(i,j) = rooted nonseparable maps with i+1 vertices, j+1 faces (i+j edges), Brown-Tutte formula.
  Networks: N = M - root for a rooted nonseparable map M with >= 2 edges; weight x^(v-2) y^(f-2),
  i.e. x marks non-pole vertices, y marks faces other than the two at the root edge.
  D = 1 + S + P + H  (trivial edge, series, parallel, h-networks with a c-net core).
  S = x (D-S)^2 / (1 - x (D-S)),  P = y (D-P)^2 / (1 - y (D-P)),  H = D * m(xD, yD),
  m(X,Y) = sum_c X^(v-2) Y^(f-2) over rooted c-nets.
Running this file prints rooted c-nets by edges and checks them against A000287.
"""
import sys
from math import factorial

from series import BS, compose_xy


def B(i, j):
    return factorial(2 * i + j - 2) * factorial(2 * j + i - 2) // (
        factorial(i) * factorial(j) * factorial(2 * i - 1) * factorial(2 * j - 1))


def networks(N):
    D = BS(N)
    for i in range(1, N + 2):
        for j in range(1, N + 2 - i + 1):
            if (i - 1) + (j - 1) <= N:
                D[i - 1, j - 1] = B(i, j)
    x = BS.mono(N, 1, 0)
    y = BS.mono(N, 0, 1)
    # S = x A^2/(1 - x A) with A = D - S: iterate to the fixed point (each pass gains a degree)
    S = BS(N)
    P = BS(N)
    for _ in range(N + 1):
        A = D - S
        S = x * A * A * (1 - x * A).inv()
        A = D - P
        P = y * A * A * (1 - y * A).inv()
    H = D - 1 - S - P
    return D, S, P, H


def invert_xy(G, X, Y):
    """Find m with m(X, Y) = G where X = x*(1+...), Y = y*(1+...)."""
    N = G.N
    m = BS(N)
    cur = BS(N)  # m(X,Y) so far
    Xp = [BS.const(N, 1)]
    Yp = [BS.const(N, 1)]
    for _ in range(N):
        Xp.append(Xp[-1] * X)
        Yp.append(Yp[-1] * Y)
    for d in range(N + 1):
        for i in range(d + 1):
            j = d - i
            r = G[i, j] - cur[i, j]
            if r:
                m[i, j] = r
                cur = cur + (Xp[i] * Yp[j]) * r
    return m


def rooted_cnets(N):
    """Return (m, D, S, P, H): m[a,b] = rooted c-nets with a+2 vertices and b+2 faces."""
    D, S, P, H = networks(N)
    x = BS.mono(N, 1, 0)
    y = BS.mono(N, 0, 1)
    G = H * D.inv()
    m = invert_xy(G, x * D, y * D)
    return m, D, S, P, H


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    m, *_ = rooted_cnets(N)
    by_e = {}
    for a, b, c in m.items():
        e = a + b + 2
        by_e[e] = by_e.get(e, 0) + c
    ref = {}
    for line in open("../data/b000287.txt"):
        p = line.split()
        if len(p) == 2 and not line.startswith("#"):
            ref[int(p[0])] = int(p[1])
    ok = True
    for e in sorted(by_e):
        if e <= N + 2:
            good = ref.get(e) == by_e[e]
            ok &= good
            print(e, by_e[e], "ok" if good else f"MISMATCH {ref.get(e)}")
    print("all ok" if ok else "FAIL")


if __name__ == "__main__":
    main()
