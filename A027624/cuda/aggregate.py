"""Combine the per-orbit GPU sums into a(8).

    python attacks/A027624/cuda/aggregate.py [RESULTS.txt] [REPS.bin]

RESULTS.txt lines: rep_index rep_mask orbit_size sum (from `a8 ... --out`). Every one of the
1,228,158 orbit representatives must appear exactly once, with the mask and orbit size from the
reps file; then a(8) = sum over orbits of orbit_size * sum.
"""
import math
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "results" / "a8_sym.txt"
reps_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parent / "rust" / "results" / "ylo_reps_d6.bin"

raw = reps_path.read_bytes()
reps = [struct.unpack_from("<II", raw, 8 * i) for i in range(len(raw) // 8)]
assert len(reps) == 1228158, len(reps)
assert sum(w for _, w in reps) == 2**32

seen = {}
for line in res.read_text().split("\n"):
    if not line.strip():
        continue
    i, mask, w, s = line.split()
    i, mask, w, s = int(i), int(mask), int(w), int(s)
    assert (mask, w) == reps[i], f"rep {i}: mask/weight {mask}/{w} do not match the reps file"
    assert seen.get(i, s) == s, f"rep {i} appears twice with different sums"
    seen[i] = s

missing = len(reps) - len(seen)
total = sum(reps[i][1] * s for i, s in seen.items())
print(f"{len(seen)} of {len(reps)} orbits present ({missing} missing)")
if missing:
    print(f"partial sum so far: {total}")
    sys.exit(1)
print(f"a(8) = {total}")
print(f"     ~ {total:.6e}, {total.bit_length()} bits")
lead = 2 * math.sqrt(math.e) * 2.0**128
print(f"ratio to 2*sqrt(e)*2^128: {total / lead:.4f}  (n = 4..7: 0.88, 1.18, 1.40, 1.29)")
