from collections import Counter
from math import factorial, prod
from functools import cache

def partitions(n, m=None):
    if n == 0: yield (); return
    for p in range(min(n, m or n), 0, -1):
        for r in partitions(n - p, p): yield (p,) + r

def z(lam): return prod(p**k * factorial(k) for p, k in Counter(lam).items())

def a(n, stack):
    # slot moves per piece kind: (multiplicity slot or None for the table, consumed, created)
    if stack:  # state (free heads, open-bottom tower ends, outer closed dolls, outer tops)
        mv = {'D': [(None, (0,0,0,0), (1,0,1,0)), (0, (1,0,0,0), (1,0,1,0)), (1, (0,1,0,0), (1,0,1,0)), (2, (0,0,1,0), (0,0,1,0))],
              'B': [(None, (0,0,0,0), (0,1,0,0)), (0, (1,0,0,0), (0,1,0,0)), (1, (0,1,0,0), (0,1,0,0)), (2, (0,0,1,0), (0,0,1,0))],
              'T': [(None, (0,0,0,0), (1,0,0,1)), (0, (1,0,0,0), (1,0,0,1)), (3, (0,0,0,1), (0,0,0,1))]}
    else:      # state (chains of dolls and bottoms, chains of tops)
        mv = {'D': [(None, (0,0), (1,0)), (0, (1,0), (1,0))], 'B': [(None, (0,0), (1,0)), (0, (1,0), (1,0))],
              'T': [(None, (0,0), (0,1)), (1, (0,1), (0,1))]}
    @cache
    def place(st, cnt, L):  # place labelled L-cycles of one size; non-table moves weigh L
        cur = Counter({(st, (0,) * len(st)): 1})
        for kind, k in zip('DBT', cnt):
            for _ in range(k):
                nxt = Counter()
                for (s, pend), w in cur.items():
                    for key, c, d in mv[kind]:
                        m = 1 if key is None else s[key] * L
                        if m: nxt[tuple(x - y for x, y in zip(s, c)), tuple(x + y for x, y in zip(pend, d))] += w * m
                cur = nxt
        out = Counter()
        for (s, pend), w in cur.items(): out[tuple(x + y for x, y in zip(s, pend))] += w
        return tuple(out.items())
    zero = (0,) * len(mv['D'][0][1])
    K = factorial(n) ** 3
    level = [Counter() for _ in range(n + 1)]
    level[0][(zero,) * n] = K
    for used in range(n):
        for states, w in level[used].items():
            for m in range(1, n - used + 1):
                for o in range(m + 1):
                    for lD in partitions(m - o):
                        for lB in partitions(o):
                            for lT in partitions(o):
                                atoms = Counter()
                                for i, lam in enumerate((lD, lB, lT)):
                                    for p in lam: atoms[p, i] += 1
                                combos = Counter({states: w * factorial(m - o) * factorial(o)**2 // (z(lD) * z(lB) * z(lT))})
                                for L in sorted({p for p, _ in atoms}):
                                    cnt = tuple(atoms[L, i] for i in range(3))
                                    nxt = Counter()
                                    for key, v in combos.items():
                                        for s2, ways in place(key[L - 1], cnt, L):
                                            nxt[key[:L - 1] + (s2,) + key[L:]] += v * ways
                                    combos = nxt
                                for key, v in combos.items():
                                    level[used + m][key] += v // (factorial(m - o) * factorial(o)**2)
    return sum(level[n].values()) // K

if __name__ == '__main__':
    print([a(n, False) for n in range(9)])
    print([a(n, True) for n in range(8)])
