"""Matryoshka halves with flat tops: brute-force enumeration.
unit = (kind, i, inside, on): kind D (closed doll), B (bottom), T (top);
inside/on are tuples holding at most one unit. Slots: inside of D/B/T (nesting), and the flat
top of a D or T (stacking) -- usable only when that piece is not shut in a closed doll or under a top."""
import sys
from functools import lru_cache
from itertools import combinations
ENC = {'inD': True, 'inT': True, 'inB_enc': True, 'inB_free': False, 'on': False, 'table': False}
ACC = {'inD': 'DB', 'inB_enc': 'DB', 'inB_free': 'DB', 'inT': 'T', 'on': 'DBT', 'table': 'DBT'}
def subsets(s):
    s = sorted(s)
    for r in range(len(s) + 1):
        for c in combinations(s, r): yield frozenset(c)
@lru_cache(None)
def units(S, kind, i, enc):
    """all units of given outer kind/index using exactly pieces S (excluding its own)"""
    inslot = {'D': 'inD', 'T': 'inT', 'B': 'inB_enc' if enc else 'inB_free'}[kind]
    out = []
    for X in subsets(S):
        for a in level(X, inslot):
            if enc or kind == 'B':
                if S - X: continue
                out.append((kind, i, a, ()))
            else:
                for b in level(S - X, 'on'): out.append((kind, i, a, b))
    return out
@lru_cache(None)
def level(S, slot):
    if not S: return [()]
    i = max(p[1] for p in S)
    mine = [p for p in S if p[1] == i]
    smaller = frozenset(p for p in S if p[1] < i)
    enc = ENC[slot]
    res = set()
    def finish(us, rem):
        for r in level(rem, slot):
            t = tuple(sorted(us + list(r)))
            if slot != 'table' and len(t) > 1: continue
            if any(u[0] not in ACC[slot] for u in t): continue
            res.add(t)
    if len(mine) == 2:
        for X in subsets(smaller):
            for u in units(X, 'D', i, enc): finish([u], smaller - X)
    def rec(us, rem, todo):
        if not todo: return finish(us, rem)
        k = todo[0][0]
        for X in subsets(rem & smaller):
            for u in units(X, k, i, enc): rec(us + [u], rem - X, todo[1:])
    rec([], smaller, mine)
    return sorted(res)
def a(n):
    return level(frozenset([('B', k) for k in range(1, n + 1)] + [('T', k) for k in range(1, n + 1)]), 'table')
if __name__ == '__main__':
    for n in range(1, int(sys.argv[1]) + 1): print(n, len(a(n)), flush=True)
