"""Reproduce Mallows's word rule for strategy (b) and compare it with the geometric generations.

Mallows (Sec. 4): a reduced word of length k containing j sub-words of type 121 is an "extra" word in
generation k - j; the type 12131 counts as one 121, not two (ABA and ACA share the middle A).  We read
j as the largest number of position-disjoint factors x y x, maximized over the reduced expressions of
the group element, and assign the quintuple of element w to generation l(w) - j(w).

The geometric generation of the same quintuple (that of its youngest sphere, engine.strategy_b) is
computed alongside, and spheres are assigned the generation 1 + min(rule generation of a quintuple
that produces them by one reflection), as in Mallows's strategy (b).
"""
import sys
from collections import Counter, defaultdict
from functools import lru_cache

from engine import reflect, IDENT, N, strategy_b

L = int(sys.argv[1]) if len(sys.argv) > 1 else 8

# group elements by BFS, with every reduced expression (words acting as successive left reflections)
length = {IDENT: 0}
words = {IDENT: [()]}
frontier = [IDENT]
for k in range(1, L + 1):
    nxt = {}
    for M in frontier:
        for i in range(N):
            M2 = reflect(M, i)
            if length.get(M2, k) != k:
                continue
            if M2 not in length:
                length[M2] = k
                nxt[M2] = []
            nxt[M2].extend((i,) + w for w in words[M])
    for M2, ws in nxt.items():
        words[M2] = ws
    frontier = list(nxt)


@lru_cache(maxsize=None)
def disjoint_121(w):
    """Largest number of position-disjoint factors x y x in the word w."""
    best = [0] * (len(w) + 1)
    for p in range(len(w) - 1, -1, -1):
        best[p] = best[p + 1]
        if p + 2 < len(w) and w[p] == w[p + 2]:
            best[p] = max(best[p], 1 + best[p + 3])
    return best[0]


rule = {M: length[M] - max(disjoint_121(w) for w in words[M]) for M in length}
_, _, sgen = strategy_b(L)
geo = {M: max(sgen[r] for r in M) for M in length if all(r in sgen for r in M)}

G = L - 3   # rule generations are complete for g <= L - 3 only if no element longer than L can drop to g
q_rule = Counter(g for g in rule.values())
print("quintuples per generation, Mallows's rule:", [q_rule[g] for g in range(G + 1)])
print("quintuples per generation, geometric:     ", [sum(1 for M in geo if geo[M] == g) for g in range(G + 1)])
print("Mallows's table:                           [1, 5, 30, 120, 480, 2070]")
diff = Counter((geo[M], rule[M]) for M in geo if geo[M] != rule[M] and min(geo[M], rule[M]) <= G)
print("elements where the rule and geometry disagree, (geometric gen, rule gen): count", dict(diff))

# spheres by the rule: generation 1 + min over producing quintuples
sph = {}
for M, g in rule.items():
    for r in M:
        if r in IDENT:
            sph[r] = 0
for M, g in sorted(rule.items(), key=lambda t: t[1]):
    for i in range(N):
        r = reflect(M, i)[i]
        if r not in sph or sph[r] > g + 1:
            sph[r] = min(sph.get(r, 99), g + 1)
s_rule = Counter(sph.values())
print("spheres per generation, Mallows's rule:   ", [s_rule[g] for g in range(G + 1)])
print("Mallows's table:                           [5, 5, 20, 90, 330, 1290]")
# quintuples whose rule generation is earlier than their youngest sphere
early = sum(1 for M in geo if rule[M] < geo[M] and rule[M] <= G)
print("quintuples the rule places before their youngest sphere exists:", early)

# Variants that might explain Mallows's 2070: generation taken from a single word instead of the best one.
pairs = Counter()
for M, ws in words.items():
    for g in {length[M] - disjoint_121(w) for w in ws}:
        pairs[g] += 1
print("distinct (element, k - j(word)) pairs per generation:", [pairs[g] for g in range(G + 1)])
mn = Counter(length[M] - min(disjoint_121(w) for w in ws) for M, ws in words.items())
print("j = min over expressions instead of max:          ", [mn[g] for g in range(G + 1)])
