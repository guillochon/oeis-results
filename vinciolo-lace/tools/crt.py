"""Recombine residue runs tm_mod_<n>.txt (mod 2^64 and the primes) into a(n) by CRT."""
import sys, re
def read(n):
    L = open(f'tm_mod_{n}.txt').read().split('\n'); res = []
    for k, line in enumerate(L):
        m = re.match(rf'^{n} (\d+) ', line)
        if m: res.append((int(m.group(1)), int(L[k + 2].strip() or 0)))
    return res
def crt(res):
    x, M = 0, 1
    for r, p in res:
        m = p or 2 ** 64
        x += M * ((r - x) * pow(M, -1, m) % m); M *= m
    return x, M
if __name__ == '__main__':
    for n in map(int, sys.argv[1:]):
        res = read(n); x, M = crt(res)
        print(n, len(res), 'moduli', M.bit_length(), 'bits; value', x.bit_length(), 'bits:', x)
