"""Burnside engine for A109455: threshold functions of n variables up to permuting variables.

    a(n) = (1/n!) * sum over cycle types lam of |class(lam)| * Fix(lam)

Averaging lemma: a threshold function invariant under a permutation sigma has weights that are
constant on the cycles of sigma. So Fix(lam) = number of threshold functions on the grid
    G = {0..l_1} x ... x {0..l_k}    (l_j = cycle lengths of sigma),
i.e. f(c) = [w.c >= t] with c_j = number of 1s inside cycle j. (sigma = id gives the cube {0,1}^n.)

Counting threshold functions on a grid G (the engine):
  * Flip reduction: reflecting axis j (c_j -> l_j - c_j) maps threshold functions to threshold
    functions, so  #all = sum over POSITIVE functions f (weights >= 0) of 2^(#axes f depends on).
  * Symmetry reduction: axes of equal length can be permuted. Every positive threshold function is
    equivalent to exactly one *regular* one (weights non-increasing within each group of equal
    lengths); its orbit has size prod(m!) / prod(tie-block sizes!).
  * Enumeration: decide f point by point in a linear extension (sum, lex) of the dominance order
        c <= c + e_a            (monotone)
        c <= c + e_i - e_j      (i < j, same group: regularity)
    A point with a predecessor labelled 1 is forced to 1 (no LP needed). Otherwise branch; a
    carried witness (w, t) certifies one branch for free and an LP is solved only for the other.
    Only branch-decided points go into the LP (forced ones are implied).

Usage:
  python A109455/engine.py 6            # a(0..6), all Fix computed by the engine
  python A109455/engine.py 7 --id-known # use A000609(n) for the identity term
  python A109455/engine.py grid 2 1 1 1 1   # count threshold functions on one grid
"""
import itertools
import math
import sys
import threading
import time
from collections import Counter
from multiprocessing import Pool

import numpy as np
from scipy.optimize import linprog

sys.setrecursionlimit(10000)

A000609 = {0: 2, 1: 4, 2: 14, 3: 104, 4: 1882, 5: 94572, 6: 15028134, 7: 8378070864,
           8: 17561539552946, 9: 144130531453121108}
EPS = 1e-9
HEARTBEAT_SECS = 30


class Grid:
    def __init__(self, lengths):
        self.L = tuple(sorted(lengths, reverse=True))  # equal lengths are contiguous
        k = self.k = len(self.L)
        self.groups = [list(g) for _, g in itertools.groupby(range(k), key=lambda a: self.L[a])]
        pts = sorted(itertools.product(*[range(l + 1) for l in self.L]), key=lambda p: (sum(p), p))
        self.pts = pts
        self.P = np.array(pts, dtype=float).reshape(len(pts), k)
        idx = {p: i for i, p in enumerate(pts)}
        self.N = len(pts)
        preds = []
        for p in pts:
            pr = set()
            for a in range(k):
                if p[a] > 0:
                    q = list(p); q[a] -= 1; pr.add(idx[tuple(q)])
            for g in self.groups:
                for x, i in enumerate(g):
                    for j in g[x + 1:]:
                        if p[i] > 0 and p[j] < self.L[j]:
                            q = list(p); q[i] -= 1; q[j] += 1; pr.add(idx[tuple(q)])
            preds.append(sorted(pr))
        self.preds = preds
        # index maps for leaf statistics
        self.step = []  # for each axis a: (from_idx, to_idx) pairs p -> p + e_a
        for a in range(k):
            fr, to = [], []
            for i, p in enumerate(pts):
                if p[a] < self.L[a]:
                    q = list(p); q[a] += 1; fr.append(i); to.append(idx[tuple(q)])
            self.step.append((np.array(fr), np.array(to)))
        self.swaps = {}  # adjacent same-group axes (a, a+1): permutation of point indices
        for g in self.groups:
            for a, b in zip(g, g[1:]):
                perm = []
                for p in pts:
                    q = list(p); q[a], q[b] = q[b], q[a]; perm.append(idx[tuple(q)])
                self.swaps[(a, b)] = np.array(perm)
        # LP static parts: w >= 0, sorted within groups, t free
        rows = []
        for g in self.groups:
            for a, b in zip(g, g[1:]):
                r = np.zeros(k + 1); r[b] = 1; r[a] = -1; rows.append(r)  # w_b - w_a <= 0
        self.A_sort = np.array(rows).reshape(len(rows), k + 1)
        self.bounds = [(0, None)] * k + [(None, None)]


def lp_witness(G, ones, zeros):
    """Find (w, t) with w.x >= t on ones, w.x <= t - 1 on zeros, w >= 0 sorted; or None."""
    k = G.k
    A = [G.A_sort]
    b = [np.zeros(len(G.A_sort))]
    if ones:
        X = G.P[ones]; A.append(np.hstack([-X, np.ones((len(ones), 1))])); b.append(np.zeros(len(ones)))
    if zeros:
        X = G.P[zeros]; A.append(np.hstack([X, -np.ones((len(zeros), 1))])); b.append(-np.ones(len(zeros)))
    c = np.r_[np.ones(k), 0.0]  # minimise total weight: keeps numbers small
    res = linprog(c, A_ub=np.vstack(A), b_ub=np.concatenate(b), bounds=G.bounds, method="highs")
    if res.status != 0:
        return None
    return res.x[:k], res.x[k]


class Counter_:
    def __init__(self, G):
        self.G = G
        self.canonical = 0
        self.total = 0          # number of threshold functions on the grid
        self.positive = 0       # number of positive threshold functions (labelled)
        self.lps = 0
        self.bad = 0
        H = 1
        for g in G.groups:
            H *= math.factorial(len(g))
        self.H = H

    def leaf(self, lab, w, t):
        G = self.G
        f = lab.astype(bool)
        # certify with the witness
        s = G.P @ w - t
        if not (np.all(s[f] >= -1e-7) and np.all(s[~f] < -1e-7)):
            self.bad += 1
        ess = sum(1 for a in range(G.k) if np.any(f[G.step[a][0]] != f[G.step[a][1]]))
        stab = 1
        for g in G.groups:
            run = 1
            for a, b in zip(g, g[1:]):
                if np.array_equal(f, f[G.swaps[(a, b)]]):
                    run += 1
                else:
                    stab *= math.factorial(run); run = 1
            stab *= math.factorial(run)
        orbit = self.H // stab
        self.canonical += 1
        self.positive += orbit
        self.total += orbit << ess

    def run(self, split_depth=None):
        """Full DFS; with split_depth, stop after that many branch decisions and return subtasks."""
        G = self.G
        lab = np.full(G.N, -1, dtype=np.int8)
        wt = lp_witness(G, [], [])
        self.split_depth, self.tasks = split_depth, []
        self._dfs(0, lab, [], [], wt[0], wt[1], 0)
        return self.tasks

    def resume(self, task):
        i, lab, ones, zeros, w, t = task
        self.split_depth = None
        self._dfs(i, lab.copy(), list(ones), list(zeros), w, t, 0)

    def _dfs(self, i, lab, ones, zeros, w, t, depth):
        G = self.G
        if self.split_depth is not None and depth == self.split_depth:
            self.tasks.append((i, lab.copy(), list(ones), list(zeros), w, t))
            return
        while i < G.N and any(lab[q] == 1 for q in G.preds[i]):
            lab[i] = 1
            i += 1
        if i == G.N:
            self.leaf(lab, w, t)
            return
        v = float(G.P[i] @ w - t)
        if v >= EPS:
            free, other = 1, 0
        elif v <= -EPS:
            free, other = 0, 1
        else:
            free = None
        branches = [free, other] if free is not None else [1, 0]
        start = i
        for br in branches:
            if br == free:
                ww, tt = w, t
            else:
                self.lps += 1
                res = lp_witness(G, ones + ([i] if br == 1 else []), zeros + ([i] if br == 0 else []))
                if res is None:
                    continue
                ww, tt = res
            lab[start] = br
            if br == 1:
                self._dfs(start + 1, lab, ones + [start], zeros, ww, tt, depth + 1)
            else:
                self._dfs(start + 1, lab, ones, zeros + [start], ww, tt, depth + 1)
            lab[start:] = -1


def count_grid(lengths):
    G = Grid(lengths)
    C = Counter_(G)
    t0 = time.time()
    C.run()
    return dict(lengths=G.L, total=C.total, positive=C.positive, canonical=C.canonical,
                lps=C.lps, bad=C.bad, secs=round(time.time() - t0, 2))


_WORKER_GRID = {}


def _run_task(args):
    lengths, task = args
    G = _WORKER_GRID.get(lengths)
    if G is None:
        G = _WORKER_GRID[lengths] = Grid(lengths)
    C = Counter_(G)
    C.resume(task)
    return C.total, C.positive, C.canonical, C.lps, C.bad


def count_grid_parallel(lengths, procs, split_depth=12):
    """Split one grid's search tree into subtrees and farm them out to a process pool."""
    G = Grid(lengths)
    C = Counter_(G)
    t0 = time.time()
    tasks = C.run(split_depth=split_depth)  # leaves reached before split_depth are already counted
    done = [0]
    stop = threading.Event()

    def heartbeat():  # periodic status line so long runs can be followed in the console/log
        while not stop.wait(HEARTBEAT_SECS):
            d, el = done[0], time.time() - t0
            eta = f"~{el / d * (len(tasks) - d):.0f}s" if d else "?"
            print(f"    [{time.strftime('%H:%M:%S')}] grid {G.L}: subtrees {d}/{len(tasks)}  "
                  f"canonical so far {C.canonical}  elapsed {el:.0f}s  ETA {eta}", flush=True)

    hb = threading.Thread(target=heartbeat, daemon=True)
    hb.start()
    with Pool(procs) as pool:
        for tot, pos, can, lps, bad in pool.imap_unordered(_run_task, [(G.L, tk) for tk in tasks], chunksize=1):
            C.total += tot; C.positive += pos; C.canonical += can; C.lps += lps; C.bad += bad
            done[0] += 1
    stop.set()
    return dict(lengths=G.L, total=C.total, positive=C.positive, canonical=C.canonical,
                lps=C.lps, bad=C.bad, secs=round(time.time() - t0, 2), tasks=len(tasks))


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def class_size(lam):
    n = sum(lam)
    c = Counter(lam)
    d = 1
    for l, m in c.items():
        d *= l ** m * math.factorial(m)
    return math.factorial(n) // d


def a(n, id_known=False, procs=1, verbose=True):
    lams = [lam for lam in partitions(n)]
    todo = [lam for lam in lams if not (id_known and all(l == 1 for l in lam))]
    if procs > 1:  # classes one after another, each grid's search tree split across the pool
        results = []
        for lam in sorted(todo, key=lambda l: len(l)):
            r = count_grid_parallel(lam, procs)
            print(f"    done {lam}: canonical={r['canonical']} tasks={r['tasks']} {r['secs']}s", flush=True)
            results.append(r)
    else:
        results = [count_grid(lam) for lam in todo]
    fix = {r["lengths"]: r for r in results}
    total = 0
    for lam in lams:
        key = tuple(sorted(lam, reverse=True))
        if id_known and all(l == 1 for l in lam):
            F, info = A000609[n], "(A000609, literature)"
        else:
            r = fix[key]
            F = r["total"]
            info = f"canonical={r['canonical']} lps={r['lps']} bad={r['bad']} {r['secs']}s"
        total += class_size(lam) * F
        if verbose:
            print(f"  {lam}: class {class_size(lam)}, Fix {F}  {info}", flush=True)
    q, rem = divmod(total, math.factorial(n))
    return q, rem


def main():
    args = sys.argv[1:]
    if args and args[0] == "grid":
        nums = [int(x) for x in args[1:] if not x.startswith("--")]
        if "--procs" in args:
            procs = int(args[args.index("--procs") + 1])
            nums = [int(x) for x in args[1:args.index("--procs")]]
            print(count_grid_parallel(nums, procs), flush=True)
        else:
            print(count_grid(nums), flush=True)
        return
    n_max = int(args[0]) if args else 6
    id_known = "--id-known" in args
    procs = int(args[args.index("--procs") + 1]) if "--procs" in args else 1
    n_min = int(args[args.index("--from") + 1]) if "--from" in args else 0
    for n in range(n_min, n_max + 1):
        t0 = time.time()
        print(f"n={n}:", flush=True)
        q, rem = a(n, id_known=id_known, procs=procs)
        print(f"==> a({n}) = {q}   (Burnside remainder {rem}, {time.time() - t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
