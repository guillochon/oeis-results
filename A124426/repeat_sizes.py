"""Matryoshkas with repeated sizes: N dolls, sizes up to shift and gaps (only the order matters), dolls of
equal size identical. A size profile is a composition m = (m_1, ..., m_k) of N (m_j dolls of the j-th size).
Equal sizes never nest or stack on one another.

Two independent counts per profile:
  * polya(m, model): a configuration is a multiset of towers; every tower has strictly decreasing sizes along
    each branch, so it has no symmetries, and the multiset count is the coefficient extraction in
    prod_towers 1/(1 - x^tower).
  * brute(m, model): builds every configuration as a canonical nested tuple and counts distinct ones.
model = 'nest' (A124426 rules) or 'stack' (flat heads).

    python repeat_sizes.py 7          # terms N = 1..7 for both models
"""
import sys
from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb


def compositions(n):
    if n == 0:
        yield ()
        return
    for first in range(1, n + 1):
        for rest in compositions(n - first):
            yield (first,) + rest


# ---------------------------------------------------------------- Polya count
def polya(m, model):
    k = len(m)
    dim = 3 * k  # coordinates (d_j, b_j, t_j) for size j = 1..k (index j-1)

    def ok(v):
        for j in range(k):
            d, b, t = v[3 * j], v[3 * j + 1], v[3 * j + 2]
            if d + b > m[j] or d + t > m[j]:
                return False
        return True

    def mul(p, q):
        r = defaultdict(int)
        for a, x in p.items():
            for b, y in q.items():
                v = tuple(i + j for i, j in zip(a, b))
                if ok(v):
                    r[v] += x * y
        return dict(r)

    def add(p, q):
        r = dict(p)
        for a, x in q.items():
            r[a] = r.get(a, 0) + x
        return r

    zero = tuple([0] * dim)
    one = {zero: 1}

    def var(kind, j):  # j is 1-based size
        v = [0] * dim
        v[3 * (j - 1) + 'DBT'.index(kind)] = 1
        return {tuple(v): 1}

    allowed = {
        'nest': {'D': 'DB', 'B': 'DB', 'T': 'T'},
        'stack': {'D': 'DBT', 'B': 'DB', 'T': 'DBT'},
    }[model]

    @lru_cache(None)
    def enc(kind, i):  # enclosed chain hanging inside a non-enclosed D or T (stacking only)
        if model == 'nest' or kind == 'B':
            return one
        p = one
        for j in range(1, i):
            p = mul(p, add(one, add(var('D', j), var('B', j)) if kind == 'D' else var('T', j)))
        return p

    @lru_cache(None)
    def tower(i, kind):
        cont = one
        for j in range(1, i):
            for y in allowed[kind]:
                cont = add(cont, tower(j, y))
        return mul(mul(var(kind, i), enc(kind, i)), cont)

    towers = {}
    for i in range(1, k + 1):
        for kind in 'DBT':
            towers = add(towers, tower(i, kind))

    F = {zero: 1}
    for tau, c in towers.items():
        if not c:
            continue
        new = defaultdict(int)
        for v, x in F.items():
            n, w = 0, v
            while ok(w):
                new[w] += x * comb(c + n - 1, n)
                n += 1
                w = tuple(a + b for a, b in zip(w, tau))
        F = dict(new)
    total = 0
    for v, x in F.items():
        if all(v[3 * j] + v[3 * j + 1] == m[j] and v[3 * j + 1] == v[3 * j + 2] for j in range(k)):
            total += x
    return total


# ---------------------------------------------------------------- brute force (canonical forms)
# unit = (kind, size, inside, on) with inside/on tuples of at most one unit (on only for 'stack').
ENC = {'inD': True, 'inT': True, 'inB_enc': True, 'inB_free': False, 'on': False, 'table': False}
ACC = {'inD': 'DB', 'inB_enc': 'DB', 'inB_free': 'DB', 'inT': 'T', 'on': 'DBT', 'table': 'DBT'}


def brute(m, model):
    k = len(m)

    def submultisets(ms):
        for c in product(*[range(x + 1) for x in ms]):
            yield c

    def sub(a, b):
        return tuple(x - y for x, y in zip(a, b))

    @lru_cache(None)
    def units(S, kind, i, enc):
        """units with outer piece (kind, i) whose other pieces are exactly the multiset S."""
        if model == 'nest':
            inslot = {'D': 'inD', 'B': 'inD', 'T': 'inT'}[kind]
        else:
            inslot = {'D': 'inD', 'T': 'inT', 'B': 'inB_enc' if enc else 'inB_free'}[kind]
        out = set()
        for X in submultisets(S):
            for a in level(X, inslot):
                rest = sub(S, X)
                if model == 'nest' or enc or kind == 'B':
                    if any(rest):
                        continue
                    out.add((kind, i, a, ()))
                else:
                    for b in level(rest, 'on'):
                        out.add((kind, i, a, b))
        return frozenset(out)

    @lru_cache(None)
    def level(S, slot):
        """all canonical tuples of units using exactly multiset S (3 counts per size: D, B, T) in slot."""
        if not any(S):
            return frozenset([()])
        i = max(j for j in range(k) if S[3 * j] or S[3 * j + 1] or S[3 * j + 2])
        enc = ENC[slot] if model == 'stack' else False
        cap = None if slot == 'table' else 1
        res = set()
        # take one outer piece of size i (largest present); all size-i pieces are outer at this level
        for kind in 'DBT':
            idx = 3 * i + 'DBT'.index(kind)
            if not S[idx] or kind not in ACC[slot]:
                continue
            rem = list(S)
            rem[idx] -= 1
            rem = tuple(rem)
            smaller = tuple(x if j // 3 < i else 0 for j, x in enumerate(rem))
            same = sub(rem, smaller)
            for X in submultisets(smaller):
                for u in units(X, kind, i, enc):
                    for r in level(sub(rem, X), slot):
                        t = tuple(sorted((u,) + r))
                        if cap is not None and len(t) > cap:
                            continue
                        res.add(t)
        return frozenset(res)

    total = 0
    for o in product(*[range(x + 1) for x in m]):
        S = []
        for j in range(k):
            S += [m[j] - o[j], o[j], o[j]]
        total += len(level(tuple(S), 'table'))
    return total


def a(n, model, method=polya):
    return sum(method(c, model) for c in compositions(n))


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for model in ('nest', 'stack'):
        print(model, [a(n, model) for n in range(1, N + 1)], flush=True)
