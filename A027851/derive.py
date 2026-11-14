"""OEIS values that follow from the order-10 self-dual census (`selfdual sdx 10` output).

    python derive.py outx_n10.txt [second run ...]

Self-dual semigroups with S^3 = {0} are not enumerated; they come from the Distler-Mitchell formulas
(A383885, A383886) plus the zero semigroup, and all have exactly one idempotent.
"""
import sys
from fractions import Fraction
from collections import defaultdict

N = 10
A001423_10 = 12418001077381302684   # semigroups up to isomorphism or anti-isomorphism
A383886_10 = 12417282095522918811   # 3-nilpotent, up to equivalence
A383885_10 = 24834563582168716305   # 3-nilpotent, up to isomorphism
A002786_10 = 12417619575092896741   # 1 idempotent, up to equivalence
A002787_10 = 381279977009776        # 2 idempotents, up to equivalence
A005591_10 = 219587421825           # 3 idempotents, up to equivalence
A002788_10 = 7033090                # bands, up to equivalence
A001426_10 = 3518930337             # commutative (Grillet)
A001428_10 = 169163                 # inverse (Malandro)
COMM_S3_ZERO_10 = 3377274622        # commutative with S^3 = 0 (check_commutative.py)


def read(path):
    by_k, extra = defaultdict(Fraction), defaultdict(Fraction)
    seen = set()
    for line in open(path):
        f = line.split()
        if not f or f[0] != "class":
            continue
        d = dict(x.split("=") for x in f[3:])
        key = tuple(f[1:-1])
        assert key not in seen
        seen.add(key)
        assert d["nil3"] == "0"
        w = Fraction(int(d["count"]), int(d["stab"]))
        by_k[int(d["k"])] += w
        for flag in ("com", "inv"):
            if d.get(flag) == "1":
                extra[flag] += w
    assert len({k[0] for k in seen}) <= 42
    row = []
    for k in range(1, N + 1):
        assert by_k[k].denominator == 1, (k, by_k[k])
        row.append(int(by_k[k]))
    return row, {k: v for k, v in extra.items()}


rows = [read(p) for p in sys.argv[1:]]
for r in rows[1:]:
    assert r[0] == rows[0][0], "runs disagree"
row, extra = rows[0]
for r in rows:
    extra.update(r[1])
sd3 = 2 * A383886_10 - A383885_10      # self-dual 3-nilpotent
print("self-dual with S^3 != 0, by idempotents:", row, " total", sum(row))
print("self-dual 3-nilpotent (formula):", sd3, "+ 1 zero semigroup")
row[0] += sd3 + 1
total = sum(row)
if "com" in extra:
    want = A001426_10 - COMM_S3_ZERO_10
    print("check commutative (S^3 != 0):", extra["com"], "expected", want, "OK" if extra["com"] == want else "MISMATCH")
if "inv" in extra:
    print("check inverse:", extra["inv"], "expected", A001428_10, "OK" if extra["inv"] == A001428_10 else "MISMATCH")
print()
print("A029851(10) =", total, "  self-converse semigroups")
print("A027851(10) =", 2 * A001423_10 - total, "  semigroups up to isomorphism")
print("A058118 row 10 =", row)
print("A058119(10) =", row[0], " A058120(10) =", row[1], " A058121(10) =", row[2], " A058122(10) =", row[9])
print("A058109(10) =", 2 * A002786_10 - row[0], "  1 idempotent, up to isomorphism")
print("A058110(10) =", 2 * A002787_10 - row[1], "  2 idempotents")
print("A058111(10) =", 2 * A005591_10 - row[2], "  3 idempotents")
print("A058112(10) =", 2 * A002788_10 - row[9], "  bands")
