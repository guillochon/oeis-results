"""Sum `selfdual sd` output: self-dual semigroups = sum count/stab, split by idempotents and 3-nilpotency.

    python combine.py out_n8.txt [more files]
"""
import sys
from fractions import Fraction
from collections import defaultdict

tot = defaultdict(Fraction)
extra = defaultdict(Fraction)
seen = set()
for path in sys.argv[1:]:
    for line in open(path):
        f = line.split()
        if not f or f[0] != "class":
            continue
        d = dict(x.split("=") for x in f[3:])
        key = (f[1], d["k"], d["nil3"], d.get("com"), d.get("inv"), d["stab"])
        assert key not in seen, key
        seen.add(key)
        w = Fraction(int(d["count"]), int(d["stab"]))
        tot[(int(d["k"]), int(d["nil3"]))] += w
        for flag in ("com", "inv"):
            if d.get(flag) == "1":
                extra[flag] += w
ks = sorted({k for k, _ in tot})
allsum = sum(tot.values())
assert allsum.denominator == 1
print("self-dual total:", allsum)
print("by idempotents:", [int(tot[(k, 0)] + tot[(k, 1)]) for k in ks])
print("S^3 = 0 (3-nilpotent, zero or trivial):", sum(v for (k, nil), v in tot.items() if nil))
print("not 3-nilpotent:", sum(v for (k, nil), v in tot.items() if not nil))
print("classes seen:", len({s[0] for s in seen}))
for flag, name in (("com", "commutative"), ("inv", "inverse")):
    if flag in extra:
        print(name + ":", extra[flag])
