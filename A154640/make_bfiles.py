"""b-files (n = 0..1000) for A154639-A154642 from the proved g.f.s, and DATA strings for all seven entries.

Recurrences (from the denominators; valid once n exceeds the numerator degree):
  A154639: a(n) = 3a(n-1) + 6a(n-2) - 12a(n-3)                                     n >= 4
  A154640: a(n) = 3a(n-1) + 3a(n-2) - 6a(n-3)                                      n >= 4
  A154641: a(n) = 5a(n-1) + 8a(n-2) - 21a(n-3) - 9a(n-4) + 6a(n-5)                 n >= 6
  A154642: a(n) = 10a(n-1) + 7a(n-2) - 85a(n-3) + 24a(n-4) + 27a(n-5) - 6a(n-6)    n >= 7
Initial terms come from the exact series (series.py) and the automaton (reduced_words.py);
every term up to n = 44 (n = 40 for A154639) is checked against those independent lists.
"""
import json
import re
from pathlib import Path

from reduced_words import automaton

HERE = Path(__file__).parent
series = {}
for line in open(HERE / "results" / "series_a_44.txt"):
    k, v = line.split(" ", 1)
    series[k] = [int(t) for t in re.findall(r"-?\d+", v)]
series["A154639"] = automaton(40)

REC = {
    "A154639": ([3, 6, -12], 4),
    "A154640": ([3, 3, -6], 4),
    "A154641": ([5, 8, -21, -9, 6], 6),
    "A154642": ([10, 7, -85, 24, 27, -6], 7),
}
NMAX = 1000
terms = {}
for a, (c, start) in REC.items():
    t = series[a][:start]
    while len(t) <= NMAX:
        n = len(t)
        t.append(sum(ci * t[n - 1 - i] for i, ci in enumerate(c)))
    ref = series[a]
    assert t[:len(ref)] == ref, a
    terms[a] = t
    (HERE / "bfiles").mkdir(exist_ok=True)
    with open(HERE / "bfiles" / f"b{a[1:]}.txt", "w", newline="\n") as f:
        for n, v in enumerate(t):
            f.write(f"{n} {v}\n")

b = json.load(open(HERE / "results" / "stratb_sym_16.json"))
for a in ("A154643", "A154644", "A154645"):
    terms[a] = b[a]


def data(t, limit=260):
    s = ""
    for k, v in enumerate(t):
        nxt = (s + "," if s else "") + str(v)
        if len(nxt) > limit:
            return s, k
        s = nxt
    return s, len(t)


out = {}
for a in sorted(terms):
    s, k = data(terms[a])
    out[a] = {"data": s, "last_n": k - 1, "chars": len(s)}
    print(f"{a}: a(0)..a({k - 1}), {len(s)} chars\n  {s}")
json.dump(out, open(HERE / "results" / "data_strings.json", "w"), indent=1)
