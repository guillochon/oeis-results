"""Fast count for matryoshkas with repeated sizes (see repeat_sizes.py for the rules and a slower check).

Burnside: a configuration of identical equal-size pieces is an orbit of configurations of labelled pieces
(D's, B's and T's of each size labelled separately) under prod_j S_{c_j} x S_{o_j} x S_{o_j}. A configuration
fixed by a permutation g is a set of towers permuted by g; a tower has no symmetries, so towers in an L-cycle are
L copies of one tower built from L-cycles of pieces. Hence Fix(g) = prod_L Lab_L, where Lab_L counts labelled
configurations on the L-cycles of g treated as single pieces, each non-base piece weighted L (its alignment with
the L tower copies). Lab_L comes from the free-slot recurrence, run size by size
(largest first); pieces of one size may not use slots created by pieces of the same size.

One DP over sizes then sums over compositions m of N, splits c_j + o_j = m_j, and cycle types of the three
permutations at each size, with weight 1/z for each cycle type.

    python repeat_fast.py 20
"""
import sys
from collections import defaultdict
from functools import lru_cache
from math import factorial


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def z(lam):
    r = 1
    for p in set(lam):
        k = lam.count(p)
        r *= p ** k * factorial(k)
    return r


# Moves for one labelled piece: list of (multiplicity key, consume, pending-add).
# State vectors: stack (s, f, e, t); nest (c, t).
MOVES = {
    'stack': {
        'D': [(None, (0, 0, 0, 0), (1, 0, 1, 0)), (0, (1, 0, 0, 0), (1, 0, 1, 0)),
              (1, (0, 1, 0, 0), (1, 0, 1, 0)), (2, (0, 0, 1, 0), (0, 0, 1, 0))],
        'B': [(None, (0, 0, 0, 0), (0, 1, 0, 0)), (0, (1, 0, 0, 0), (0, 1, 0, 0)),
              (1, (0, 1, 0, 0), (0, 1, 0, 0)), (2, (0, 0, 1, 0), (0, 0, 1, 0))],
        'T': [(None, (0, 0, 0, 0), (1, 0, 0, 1)), (0, (1, 0, 0, 0), (1, 0, 0, 1)),
              (3, (0, 0, 0, 1), (0, 0, 0, 1))],
    },
    'nest': {
        'D': [(None, (0, 0), (1, 0)), (0, (1, 0), (1, 0))],
        'B': [(None, (0, 0), (1, 0)), (0, (1, 0), (1, 0))],
        'T': [(None, (0, 0), (0, 1)), (1, (0, 1), (0, 1))],
    },
}


@lru_cache(None)
def place_group(model, state, nD, nB, nT, L=1):
    """Place nD, nB, nT labelled pieces of one size into slot state; returns {new_state: ways}.
    For L-cycles of pieces, a placement into an existing tower can be aligned with the L copies of that
    tower in L ways, so every non-table move gets a factor L."""
    cur = {(state, tuple([0] * len(state))): 1}
    for kind, n in (('D', nD), ('B', nB), ('T', nT)):
        for _ in range(n):
            nxt = defaultdict(int)
            for (st, pend), w in cur.items():
                for key, cons, add in MOVES[model][kind]:
                    mult = 1 if key is None else st[key] * L
                    if mult <= 0:
                        continue
                    st2 = tuple(a - b for a, b in zip(st, cons))
                    pend2 = tuple(a + b for a, b in zip(pend, add))
                    nxt[st2, pend2] += w * mult
            cur = nxt
    out = defaultdict(int)
    for (st, pend), w in cur.items():
        out[tuple(a + b for a, b in zip(st, pend))] += w
    return dict(out)


def counts(model, N):
    """a(n) for n = 1..N: sum over compositions of n of the number of configurations."""
    zero = (0, 0, 0, 0) if model == 'stack' else (0, 0)
    # per-size choices: m -> list of (weight, {L: (nD, nB, nT)})
    choices = {}
    for m in range(1, N + 1):
        lst = []
        for o in range(m + 1):
            c = m - o
            for lD in partitions(c):
                for lB in partitions(o):
                    for lT in partitions(o):
                        # permutations of this cycle type, over c! o! o! (applied as an exact division)
                        w = (factorial(c) // z(lD)) * (factorial(o) // z(lB)) * (factorial(o) // z(lT))
                        div = factorial(c) * factorial(o) ** 2
                        atoms = defaultdict(lambda: [0, 0, 0])
                        for i, lam in enumerate((lD, lB, lT)):
                            for p in lam:
                                atoms[p][i] += 1
                        lst.append((w, div, {L: tuple(v) for L, v in atoms.items()}))
        choices[m] = lst
    # DP over parts of the composition (sizes from largest to smallest); state = (used, per-L slot states)
    # per-L states stored as a tuple indexed by L = 1..N
    start = tuple([zero] * N)
    # integer weights scaled by K = (N!)^3: each prefix weight K / prod(c! o! o!) is an integer
    K = factorial(N) ** 3
    dp = {start: K}
    by_used = [dict() for _ in range(N + 1)]
    by_used[0] = dp
    for used in range(N):
        for states, w in by_used[used].items():
            for m in range(1, N - used + 1):
                for cw, div, atoms in choices[m]:
                    results = [{states[L - 1]: 1} for L in range(1, N + 1)]
                    for L, (nD, nB, nT) in atoms.items():
                        results[L - 1] = place_group(model, states[L - 1], nD, nB, nT, L)
                    # combine (product over L of independent outcomes)
                    combos = {(): w * cw}
                    for L in range(N):
                        nxt = {}
                        for pre, pw in combos.items():
                            for st, ways in results[L].items():
                                key = pre + (st,)
                                nxt[key] = nxt.get(key, 0) + pw * ways
                        combos = nxt
                    tgt = by_used[used + m]
                    for key, val in combos.items():
                        tgt[key] = tgt.get(key, 0) + val // div
    out = []
    for n in range(1, N + 1):
        # weight sum over all states reached with exactly n dolls... but by_used mixes prefixes; total
        # configurations for n dolls = sum over states in by_used[n]. Cycle lengths L > n never get atoms.
        tot = sum(by_used[n].values())
        assert tot % K == 0, n
        out.append(tot // K)
    return out


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for model in ('nest', 'stack'):
        print(model, counts(model, N), flush=True)
