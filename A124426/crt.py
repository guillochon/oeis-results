"""Reassemble residues printed by the Rust counters with the Chinese remainder theorem.

    python crt.py stack stack_residues.txt         # lines "n p r", many primes per n
    python crt.py repeat repeat_nest_residues.txt   # "# primes p1..pP" header, then lines "n r1 .. rP"
Prints "n value" for every n. With --check K the last K primes are held back and must agree.
"""
import sys
from collections import defaultdict


def crt(pairs):
    x, m = 0, 1
    for p, r in pairs:
        t = ((r - x) * pow(m, -1, p)) % p
        x += m * t
        m *= p
    return x, m


def load(kind, path):
    res = defaultdict(list)
    primes = None
    for line in open(path):
        if line.startswith('#'):
            primes = list(map(int, line.split()[2:]))
            continue
        parts = list(map(int, line.split()))
        if kind == 'stack':
            n, p, r = parts
            res[n].append((p, r))
        else:
            n, rs = parts[0], parts[1:]
            res[n] = list(zip(primes, rs))
    return res


def reconstruct(kind, path, check):
    out = {}
    for n, pairs in sorted(load(kind, path).items()):
        use, held = (pairs[:-check], pairs[-check:]) if check else (pairs, [])
        x, m = crt(use)
        # values are nonnegative and assumed < m / 2 ** 8 (checked by the held-back primes)
        for p, r in held:
            if x % p != r:
                raise SystemExit(f"n={n}: reconstruction disagrees with held-back prime {p}")
        if x.bit_length() > m.bit_length() - 8:
            raise SystemExit(f"n={n}: value too close to the modulus ({x.bit_length()} of {m.bit_length()} bits)")
        out[n] = x
    return out


if __name__ == '__main__':
    kind, path = sys.argv[1], sys.argv[2]
    check = int(sys.argv[sys.argv.index('--check') + 1]) if '--check' in sys.argv else 0
    for n, x in reconstruct(kind, path, check).items():
        print(n, x)
