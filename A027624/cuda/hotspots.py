"""Summarise an Nsight Compute SASS+CUDA source CSV (ncu --page source --csv --print-source sass,cuda):
top CUDA lines by warp-stall samples, executed instructions and excessive shared wavefronts.
    python3 hotspots.py SASS.csv"""
import csv, sys
rows = list(csv.reader(open(sys.argv[1], encoding="utf-8", errors="replace")))
hi = next(i for i, r in enumerate(rows) if r and r[0] == "Line No")
hdr = rows[hi]
C = {h: i for i, h in reversed(list(enumerate(hdr)))}
lines = [r for r in rows[hi + 1:] if len(r) > 5 and r[2] == "-"]
def v(r, name):
    try:
        return float(r[C[name]].replace(",", ""))
    except (ValueError, KeyError):
        return 0.0
keys = ["Warp Stall Sampling (All Samples)", "Instructions Executed", "L1 Wavefronts Shared Excessive", "Divergent Branches"]
tot = {k: sum(v(r, k) for r in lines) or 1 for k in keys}
for k in keys[:3]:
    print(f"\n=== top lines by {k} (total {tot[k]:.3e}) ===")
    for r in sorted(lines, key=lambda r: -v(r, k))[:14]:
        if v(r, k) == 0:
            break
        print(f"{100*v(r,k)/tot[k]:5.1f}%  stall {100*v(r,keys[0])/tot[keys[0]]:5.1f}%  inst {100*v(r,keys[1])/tot[keys[1]]:5.1f}%  shx {100*v(r,keys[2])/tot[keys[2]]:5.1f}% | L{r[0]:>4} {r[1].strip()[:95]}")
