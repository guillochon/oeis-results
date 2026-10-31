"""CRT-combine the per-prime outputs of unil.c into exact terms of A035512 (column 1) and A003088 (column 2).
The prime(s) beyond the first NUSE check the result; A035512 is checked against its OEIS b-file and
A003088 against the Python program's b-file (n <= 22).

    python crt.py N NUSE        # reads results/nN_p*.txt, writes results/a003088_nN.txt, a035512_nN.txt
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [int(x) for x in open(os.path.join(HERE, 'primes.txt'))]
# the A035512 b-file (https://oeis.org/A035512/b035512.txt); the check is skipped if it is not found
B035512 = next((f for f in (os.path.join(HERE, 'b035512.txt'),
                            os.path.expanduser('~/oeis-mcp/data/oeisdata/files/A035/b035512.txt')) if os.path.exists(f)), None)


def crt(res, ps):
    x, m = 0, 1
    for r, p in zip(res, ps):
        t = ((r - x) * pow(m, -1, p)) % p
        x, m = x + m * t, m * p
    return x, m


def bfile(path):
    out = {}
    for line in open(path):
        w = line.split()
        if len(w) == 2 and not line.startswith('#'):
            out[int(w[0])] = int(w[1])
    return out


def main(N, nuse):
    tab = {}
    for f in glob.glob(os.path.join(HERE, f'results/n{N}_p*.txt')):
        pi = int(f.rsplit('_p', 1)[1].split('.')[0])
        for line in open(f):
            n, s, a = map(int, line.split())
            tab.setdefault(n, {})[pi] = (s, a)
    have = sorted(tab[N])
    assert len(have) > nuse, f'only {len(have)} primes finished'
    use, chk = have[:nuse], have[nuse:]
    cols = []
    for col in (0, 1):
        vals = []
        for n in range(N + 1):
            x, m = crt([tab[n][i][col] for i in use], [PRIMES[i] for i in use])
            assert x < m >> 20, (n, col, 'value too close to the CRT modulus')
            for i in chk:
                assert x % PRIMES[i] == tab[n][i][col], (n, col, 'check prime', i, 'disagrees')
            vals.append(x)
        cols.append(vals)
    strong, a = cols
    ref = bfile(B035512) if B035512 else {}
    ns = [n for n in range(N + 1) if n in ref]
    assert all(strong[n] == ref[n] for n in ns), 'A035512 mismatch'
    a035512 = f'A035512 matches its b-file for n = {ns[0]}..{ns[-1]}' if ns else 'A035512 b-file not found, check skipped'
    old = bfile(os.path.join(HERE, 'b003088.txt'))
    assert all(a[n] == old[n] for n in old if n <= N), 'A003088 mismatch with b003088.txt'
    for name, vals in (('a003088', a), ('a035512', strong)):
        with open(os.path.join(HERE, f'results/{name}_n{N}.txt'), 'w') as fh:
            fh.writelines(f'{n} {v}\n' for n, v in enumerate(vals))
    print(f'# n <= {N}: {nuse} primes for CRT, {len(chk)} check prime(s) agree; '
          f'{a035512}; A003088 matches b003088.txt for n <= {max(old)}')
    return a


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
