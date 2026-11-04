"""Combine the sliced raw sums of run_stat2.sh into fix values and append them to stats.txt.

    python combine_stats.py R TAG [--dry]

Checks: the slices together hold A000112(R) posets with labeled total A001035(R), and every raw
sum is divisible by the falling factorial (R)_t.
"""
import glob
import sys
from collections import defaultdict

A001035 = [1, 1, 3, 19, 219, 4231, 130023, 6129859, 431723379, 44511042511, 6611065248783,
           1396281677105899, 414864951055853499, 171850728381587059351, 98484324257128207032183,
           77567171020440688353049939]
A000112 = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284, 46749427, 1104891746,
           33823827452, 1338193159771, 68275077901156]

r, tag = int(sys.argv[1]), sys.argv[2]
dry = "--dry" in sys.argv
files = sorted(glob.glob(f"stats_raw/stat{r}_{tag}_p*.txt"))
lab = cnt = 0
raw = defaultdict(int)
for f in files:
    for line in open(f):
        p = line.split()
        if p[0] == "lab":
            lab += int(p[2]); cnt += int(p[3])
        elif p[0] == "raw":
            raw[p[2]] += int(p[3])
assert lab == A001035[r], f"labeled total {lab} != A001035({r}) (slices incomplete?)"
assert cnt == A000112[r], f"{cnt} posets != A000112({r})"
out = []
for key, s in raw.items():
    t = len(key.split(","))
    falling = 1
    for x in range(r - t + 1, r + 1):
        falling *= x
    assert s % falling == 0, key
    out.append(f"{r} {key} {s // falling}")
print(f"{len(files)} slices, {cnt} posets, labeled total OK")
print("\n".join(out))
if not dry:
    with open("stats.txt", "a") as fh:
        fh.write("\n".join(out) + "\n")
