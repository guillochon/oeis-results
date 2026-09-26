"""n = 9 terms of related threshold-function sequences that follow by exact arithmetic.

Inputs (n = 9):
  A000609(9) = 144130531453121108  threshold functions; independently confirmed by our n = 9 identity run
  A000617(9) = 993061484           NP classes; our identity run, agreeing with Kurz (2012) and A000619(9)
  A001532(9) = 175428              self-dual NP classes; OEIS value, reproduced by our identity run
Everything else is from the OEIS entries as of 2026-09. Each relation is first checked on every
known term; then the missing n = 9 terms are derived. Run: python A109455/related_terms.py
"""
from math import comb

A000609 = [2, 4, 14, 104, 1882, 94572, 15028134, 8378070864, 17561539552946, 144130531453121108]
A000615 = [2, 2, 8, 72, 1536, 86080, 14487040, 8274797440, 17494930604032]          # exactly n vars
A002078 = [2, 3, 6, 20, 150, 3287, 244158, 66291591, 68863243522]                   # N-classes, <= n
A002079 = [2, 1, 2, 9, 96, 2690, 226360, 64646855, 68339572672]                     # N-classes, exactly n
A000617 = [2, 3, 5, 10, 27, 119, 1113, 29375, 2730166]          # NP classes, <= n (a(9) omitted: disputed)
A000619 = [2, 1, 2, 5, 17, 92, 994, 28262, 2700791, 990331318]  # NP classes, exactly n
A001529 = [1, 2, 3, 6, 15, 63, 567, 14755, 1366318]             # NPN classes, <= n
A001530 = [1, 1, 1, 3, 9, 48, 504, 14188, 1351563]              # NPN classes, exactly n
A001532 = [0, 1, 1, 2, 3, 7, 21, 135, 2470, 175428]             # self-dual NP classes (offset 1; a(0) := 0)

A000617_9 = 993061484  # corrected value (our identity run; Kurz 2012 Table 2 + 2; A000617(8) + A000619(9))

checks = []


def check(name, cond):
    checks.append((name, cond))
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


def inv_binomial(a, n):
    return sum((-1) ** (n - k) * comb(n, k) * a[k] for k in range(n + 1))


print("Relations on known terms:")
check("A000615 = inverse binomial transform of A000609 (n <= 8)",
      all(inv_binomial(A000609, n) == A000615[n] for n in range(9)))
check("A000615(n) = 2^n * A002079(n) (negate any subset of the n essential variables), n <= 8",
      all(A000615[n] == 2 ** n * A002079[n] for n in range(9)))
check("A002079 = inverse binomial transform of A002078 (n <= 8)",
      all(inv_binomial(A002078, n) == A002079[n] for n in range(9)))
check("A000619(n) = A000617(n) - A000617(n-1), 1 <= n <= 8",
      all(A000619[n] == A000617[n] - A000617[n - 1] for n in range(1, 9)))
check("A001530(n) = A001529(n) - A001529(n-1), 2 <= n <= 8",
      all(A001530[n] == A001529[n] - A001529[n - 1] for n in range(2, 9)))
check("A001529(n) = (A000617(n) + A001532(n)) / 2, n <= 8",
      all(2 * A001529[n] == A000617[n] + A001532[n] for n in range(9)))

print("\nThe disputed A000617(9):")
check("A000617(8) + A000619(9) = 993061484 (A000619(9) is consistent with the corrected value)",
      A000617[8] + A000619[9] == A000617_9)
check("OEIS A000617(9) = 989913346 is NOT consistent with A000619(9)",
      A000617[8] + A000619[9] != 989913346)

print("\nDerived n = 9 terms:")
a615 = inv_binomial(A000609, 9)
q, r = divmod(a615, 2 ** 9)
check("A000615(9) divisible by 2^9 (consistency check on A000609(9))", r == 0)
a2079 = q
a2078 = sum(comb(9, k) * (A002079 + [a2079])[k] for k in range(10))
check("A002078(9) = binomial mean of A000609 as well", a2078 * 2 ** 9 == sum(comb(9, k) * A000609[k] for k in range(10)))
check("A002078(9) equals the identity run's positive count 281814234754247", a2078 == 281814234754247)
two = A000617_9 + A001532[9]
check("A000617(9) + A001532(9) is even", two % 2 == 0)
a1529 = two // 2
a1530 = a1529 - A001529[8]

results = {
    "A000615(9)": a615,
    "A002079(9)": a2079,
    "A002078(9)": a2078,
    "A001529(9)": a1529,
    "A001530(9)": a1530,
}
for k, v in results.items():
    print(f"  {k} = {v}")

assert all(c for _, c in checks), "a check failed"
print(f"\nAll {len(checks)} checks passed.")
