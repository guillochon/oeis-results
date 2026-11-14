"""Combine the residues printed by `regular-ideals k4 m ... prime` (results/k4_m*.txt) by CRT.

a(n) <= a(n-1)^2 (an ideal of M(n) is a pair of ideals of M(n-1)), so for n = 11 the value is
below 4191080125665004647^2 < 2^124 < product of the four primes.
Run: python attacks/A132183/crt.py attacks/A132183/results/k4_m7.txt
"""
import re
import sys

res = {}
for line in open(sys.argv[1]):
    mt = re.search(r"a\((\d+)\) mod (\d+) = (\d+)", line)
    if mt:
        n, p, r = map(int, mt.groups())
        res[p] = r
    mt = re.search(r"selfdual\((\d+)\) = (\d+)", line)
    if mt:
        sd = int(mt.group(2))
x, M = 0, 1
for p, r in res.items():
    x += M * ((r - x) * pow(M, -1, p) % p)
    M *= p
print(f"primes: {len(res)}  modulus bits: {M.bit_length()}")
print(f"a({n}) = {x}  ({x.bit_length()} bits)")
print(f"selfdual({n}) = {sd}")
print(f"a({n}) mod 2^31-1 = {x % 2147483647}")
