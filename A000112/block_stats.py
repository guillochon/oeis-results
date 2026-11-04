"""Turn block-engine statistics (K, V, Ch, Lo at quotient size r) into stats.txt fix values.

    python block_stats.py r K V Ch Lo [--dry]

Two interacting cycles with d = gcd: fix = P(r) + 2(2^d - 2) K(r)       (keys 2,2 4,2 6,2 3,3 4,4)
Three all-even cycles (d = 2 throughout, keys 2,2,2 and 4,2,2):
    fix = 6 [C(r,3) P + 2(r-2) r(r-1) K + 4 (2 r C(r-1,2) V + r(r-1)(r-2) Ch) + 4 r(r-1)(r-2) Lo] / (r)_3
(Σ_v C(deg v, 2) counts Λ-, V- and chain-shaped Hasse 2-paths; Λ = V by duality.)
"""
import sys
from math import comb, gcd

A001035 = [1, 1, 3, 19, 219, 4231, 130023, 6129859, 431723379, 44511042511, 6611065248783,
           1396281677105899, 414864951055853499, 171850728381587059351, 98484324257128207032183,
           77567171020440688353049939, 83480529785490157813844256579]

r = int(sys.argv[1])
K, V, Ch, Lo = (int(x) if x != "-" else None for x in sys.argv[2:6])
P = A001035[r]
out = []
for a, b in [(2, 2), (4, 2), (6, 2), (3, 3), (4, 4)]:
    d = gcd(a, b)
    out.append(f"{r} {a},{b} {P + 2 * (2 ** d - 2) * K}")
if V is not None:
    S = 6 * (comb(r, 3) * P + 2 * (r - 2) * r * (r - 1) * K
             + 4 * (2 * r * comb(r - 1, 2) * V + r * (r - 1) * (r - 2) * Ch)
             + 4 * r * (r - 1) * (r - 2) * Lo)
    f3 = r * (r - 1) * (r - 2)
    assert S % f3 == 0
    out += [f"{r} 2,2,2 {S // f3}", f"{r} 4,2,2 {S // f3}"]
print("\n".join(out))
if "--dry" not in sys.argv:
    with open("stats.txt", "a") as fh:
        fh.write("\n".join(out) + "\n")
