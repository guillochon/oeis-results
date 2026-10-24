"""a(n), n = 3..200, from the exact run: writes results/a_n.txt (b-file style) and fills A_N in the viz page.

n = 3, 4 have narrow plateaus (the f = 1.3 sample falls off them), so they come from the r = sec*(1+1e-4) runs."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "..", "..", "viz", "circle_puzzles", "index.html")

a = {3: 69, 4: 35}
for line in open(os.path.join(HERE, "results", "exact_3_200.txt")):
    m = re.match(r"(\d+) \[([\d, ]+)\] (\d+) (\w+) (\S+)", line)
    if m and int(m[1]) >= 5:
        vals = {int(v) for v in m[2].split(",")}
        assert m[4] == "fit" and m[5] == "check-ok" and len(vals) == 1, line
        a[int(m[1])] = vals.pop()
ns = sorted(a)
assert ns == list(range(3, 201))
with open(os.path.join(HERE, "results", "a_n.txt"), "w") as f:
    f.writelines(f"{n} {a[n]}\n" for n in ns)
if os.path.exists(PAGE):   # the viewer page lives in the oeis-sniper repo only
    page = open(PAGE, encoding="utf8").read()
    page = re.sub(r"const A_N = \[[^\]]*\];", "const A_N = [" + ",".join(str(a[n]) for n in ns) + "];", page)
    open(PAGE, "w", encoding="utf8").write(page)
print(len(ns), "terms; a(200) =", a[200])
