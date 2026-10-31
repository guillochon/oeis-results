"""CRT-combine the per-prime outputs of rust/mtg into exact A000080 terms; the extra prime(s) check the result."""
import sys, glob

PRIMES = [4611686018427387847, 4611686018427387817, 4611686018427387787, 4611686018427387761,
          4611686018427387751]


def crt(res, ps):
    x, m = 0, 1
    for r, p in zip(res, ps):
        t = ((r - x) * pow(m, -1, p)) % p
        x, m = x + m * t, m * p
    return x, m


def main(tag, nuse):
    files = sorted(glob.glob(f'results/{tag}_p*.txt'))
    tab = {}
    for f in files:
        pi = int(f.rsplit('_p', 1)[1].split('.')[0])
        for line in open(f):
            n, v = map(int, line.split())
            tab.setdefault(n, {})[pi] = v
    out = []
    for n in sorted(tab):
        have = sorted(tab[n])
        use = have[:nuse]
        x, m = crt([tab[n][i] for i in use], [PRIMES[i] for i in use])
        assert x < m // 2 ** 20, (n, 'value too close to the CRT modulus')
        for i in have[nuse:]:
            assert x % PRIMES[i] == tab[n][i], (n, 'check prime', i, 'disagrees')
        out.append((n, x))
    return out, len(have) - nuse


if __name__ == '__main__':
    tag = sys.argv[1]
    nuse = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    out, nchk = main(tag, nuse)
    for n, x in out:
        print(n, x)
    print(f'# {len(out)} terms, {nuse} primes for CRT, {nchk} check prime(s) agree', file=sys.stderr)
