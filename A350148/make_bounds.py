"""Write bounds/n{n}_{mode}.txt: the frozen perimeter (edge letter per tile, '.' = free), top row first."""
import os
import sys

from motif import boundary, modes

os.makedirs("bounds", exist_ok=True)
for n in map(int, sys.argv[1:] or ["5", "7", "9", "11", "13"]):
    for m in modes(n):
        g = boundary(n, m)
        rows = ["".join("BRTL"[g[(x, y)] % 4] if (x, y) in g else "." for x in range(n)) for y in range(n - 1, -1, -1)]
        with open(f"bounds/n{n}_{m}.txt", "w") as f:
            f.write(f"# n={n} mode={m}\n" + "\n".join(rows) + "\n")
        print(n, m)
