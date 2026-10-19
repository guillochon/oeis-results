"""A123234 independent check of M_m, the k = 2 term.

M_m = number of ordered 2m-tuples (p_1, ..., p_2m) of permutations of [m] that together cover every
cell (i, j) of the m x m grid exactly twice.  rc_classes.M2 counts them permutation by permutation,
with the grid of remaining capacities as the state.  This program counts the same objects row by row
instead: row i says, for each of the 2m labelled permutations t, which column p_t(i) it uses; every
column must receive exactly two of them, and no t may reuse a column.  The state is the multiset of
used-column sets of the 2m permutations (the number of completions depends only on that multiset).

    python semi_latin.py 5      # prints M_1..M_5 and compares with rc_classes.M2
"""
import sys
from functools import lru_cache

from rc_classes import M2


def count(m):
    full = (1 << m) - 1

    def assignments(used):
        """All ways to give each permutation t a new column, two permutations per column.
        Yields the new used-sets as a tuple (labelled)."""
        t_n = len(used)
        out = []
        new = [0] * t_n
        load = [0] * m

        def rec(t):
            if t == t_n:
                out.append(tuple(new))
                return
            free = full & ~used[t]
            while free:
                j = (free & -free).bit_length() - 1
                free &= free - 1
                if load[j] < 2:
                    load[j] += 1
                    new[t] = used[t] | (1 << j)
                    rec(t + 1)
                    load[j] -= 1

        rec(0)
        return out

    @lru_cache(maxsize=None)
    def f(state, rows_left):
        if rows_left == 0:
            return 1
        tot = 0
        for nxt in assignments(state):
            tot += f(tuple(sorted(nxt)), rows_left - 1)
        return tot

    return f(tuple([0] * (2 * m)), m)


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for m in range(1, top + 1):
        a, b = count(m), M2(m)
        print(m, a, b, "MATCH" if a == b else "MISMATCH", flush=True)
