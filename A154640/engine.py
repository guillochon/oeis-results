"""Exact engine for Mallows's 3D Apollonian packing sequences (A154638-A154645).

A configuration is M*C where M is a product of the Descartes reflections R_i
(5x5 integer matrices) and C is the starting five-sphere configuration in
augmented curvature-center coordinates (invertible).  So a sphere is identified
exactly by an integer row vector v of M, and its bend is v . b0 with
b0 = (1, 1, 1, 0, 0): Mallows's start (three unit spheres between two planes).

Strategy (a): quintuples are the group elements w; generation = length l(w).
              A sphere's generation = least l(w) over configurations containing it.
Strategy (b): a quintuple (set of 5 mutually tangent spheres) exists from the
              generation of its youngest sphere; reflecting a quintuple of
              generation g gives spheres of generation g + 1 (if new).
"""
from collections import defaultdict
import sys

N = 5
B0 = (1, 1, 1, 0, 0)


def reflect(rows, i):
    """Rows of R_i M given the rows of M: row i -> sum of others - row i."""
    s = [sum(c) for c in zip(*rows)]
    r = list(rows)
    r[i] = tuple(s[k] - 2 * rows[i][k] for k in range(N))
    return tuple(r)


def bend(v):
    return sum(a * b for a, b in zip(v, B0))


IDENT = tuple(tuple(int(i == j) for j in range(N)) for i in range(N))


def strategy_a(K):
    """Group BFS by length; spheres by first appearance."""
    seen_el = {IDENT}
    frontier = [IDENT]
    sphere_gen = {r: 0 for r in IDENT}
    elems = [1]
    out = [(5, sum(bend(r) for r in IDENT), sum(bend(r) ** 2 for r in IDENT))]
    for k in range(1, K + 1):
        nxt = []
        for M in frontier:
            for i in range(N):
                M2 = reflect(M, i)
                if M2 not in seen_el:
                    seen_el.add(M2)
                    nxt.append(M2)
        cnt = s = t = 0
        for M in nxt:
            for r in M:
                if r not in sphere_gen:
                    sphere_gen[r] = k
                    b = bend(r)
                    cnt += 1; s += b; t += b * b
        elems.append(len(nxt))
        out.append((cnt, s, t))
        # elements of length k-1 are no longer needed for membership tests:
        # a reflection changes length by exactly 1 in a Coxeter group
        if k >= 2:
            for M in frontier_prev:
                seen_el.discard(M)
        frontier_prev = frontier
        frontier = nxt
    return elems, out, sphere_gen


def strategy_b(K):
    """Level-by-level closure over unordered quintuples (frozensets of rows)."""
    sphere_gen = {r: 0 for r in IDENT}
    seen = [set(), {frozenset(IDENT)}]   # quintuple keys of levels g-1, g
    level = [IDENT]
    quints = []
    out = [(5, sum(bend(r) for r in IDENT), sum(bend(r) ** 2 for r in IDENT))]
    for g in range(0, K):
        # close level g: reflecting a gen-g quintuple may give an old sphere,
        # i.e. another gen-g quintuple (Mallows's "extra" quintuples)
        nxt, nxt_seen = [], set()
        stack = list(level)
        n_level = len(level)
        cnt = s = t = 0
        while stack:
            Q = stack.pop()
            for i in range(N):
                Q2 = reflect(Q, i)
                key = frozenset(Q2)
                if key in seen[0] or key in seen[1] or key in nxt_seen:
                    continue
                r = Q2[i]
                if r not in sphere_gen:
                    sphere_gen[r] = g + 1
                    b = bend(r)
                    cnt += 1; s += b; t += b * b
                if max(sphere_gen[x] for x in Q2) <= g:
                    seen[1].add(key)
                    stack.append(Q2)
                    n_level += 1
                else:
                    nxt_seen.add(key)
                    nxt.append(Q2)
        quints.append(n_level)
        out.append((cnt, s, t))
        seen = [seen[1], nxt_seen]
        level = nxt
    return quints, out, sphere_gen


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    elems, out, _ = strategy_a(K)
    s0, t0 = out[0][1], out[0][2]
    print("elements", elems)
    print("spheres(a)", [c for c, _, _ in out])
    print("s/s0(a)", [s / s0 for _, s, _ in out])
    print("t/t0(a)", [t / t0 for _, _, t in out])
    quints, out, _ = strategy_b(K)
    print("quintuples(b)", quints)
    print("spheres(b)", [c for c, _, _ in out])
    print("s/s0(b)", [s / s0 for _, s, _ in out])
    print("t/t0(b)", [t / t0 for _, _, t in out])
