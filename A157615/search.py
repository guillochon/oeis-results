"""Heuristic search for long alternating self-avoiding paths on the n X n board (A157615).

Randomized depth-first search with Warnsdorff ordering: try the next square with the fewest onward
moves first, break ties randomly, and restart with a new random start after a node budget. Every
path found is validated with path.check and saved to certificates/n<n>.json.

For odd n a path of length n^2 - 2n + 4 proves a(n) = n^2 - 2n + 4, because the four-colouring
argument (README) shows no path is longer. For even n a path is only a lower bound.

    python search.py 13 25          # n = 13..25, target = conjectured value
    python search.py 19 19 --seconds 600
"""
import argparse
import json
import random
import sys
import time
from pathlib import Path

from path import check, conjectured, draw

HERE = Path(__file__).resolve().parent
sys.setrecursionlimit(1_000_000)


def search(n, target, seconds, seed, budget=200_000):
    rng = random.Random(seed)
    t0 = time.time()
    tries = 0
    while time.time() - t0 < seconds:
        tries += 1
        seen = [[False] * n for _ in range(n)]
        x, y = rng.randrange(n), rng.randrange(n)
        seen[x][y] = True
        path = [(x, y)]
        nodes = 0

        def free(u, v):
            return 0 <= u < n and 0 <= v < n and not seen[u][v]

        def moves(u, v, horiz):
            return [(u + 1, v), (u - 1, v)] if horiz else [(u, v + 1), (u, v - 1)]

        def dfs(u, v, horiz):
            nonlocal nodes
            nodes += 1
            if len(path) >= target:
                return True
            if nodes > budget:
                return False
            cand = [c for c in moves(u, v, horiz) if free(*c)]
            cand.sort(key=lambda c: (sum(free(*d) for d in moves(c[0], c[1], not horiz)), rng.random()))
            for a, b in cand:
                seen[a][b] = True
                path.append((a, b))
                if dfs(a, b, not horiz):
                    return True
                seen[a][b] = False
                path.pop()
            return False

        if dfs(x, y, rng.random() < 0.5):
            return path, tries
    return None, tries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lo", type=int)
    ap.add_argument("hi", type=int)
    ap.add_argument("--seconds", type=float, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    (HERE / "certificates").mkdir(exist_ok=True)
    for n in range(a.lo, a.hi + 1):
        target = conjectured(n)
        t0 = time.time()
        p, tries = search(n, target, a.seconds, a.seed + n)
        dt = time.time() - t0
        if p is None:
            print(f"n={n}: no path of length {target} in {dt:.0f} s ({tries} restarts)", flush=True)
            continue
        err = check(n, p)
        assert err is None, err
        kind = "proves a(n) (odd n: matches the upper bound)" if n % 2 else "lower bound (even n)"
        print(f"n={n}: path of length {len(p)} in {dt:.1f} s ({tries} restarts) -> {kind}", flush=True)
        (HERE / "certificates" / f"n{n}.json").write_text(json.dumps(
            {"n": n, "length": len(p), "path": p, "valid": True}), encoding="utf-8")
        if n <= 11:
            print(draw(n, p))


if __name__ == "__main__":
    main()
