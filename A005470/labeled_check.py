"""Labeled check of the identity sector, mod each prime.

The p1^n coefficient of a cycle index is (labeled count)/n!, and the projection p_k -> 0 (k >= 2) is a
ring map that commutes with d/dp1 and kills plethysms, so the p1^n coefficients of every series in the
tower depend only on the identity sector.  This compares, mod each prime q, against the published
labeled counts (Bodirsky, Groepl, Kang; b-files in the local OEIS data):
  A096330  3-connected  <- kformula rows (identity rows only)
  A096331  2-connected  <- tower blocks B (GPU, same code path as the unlabeled run)
  A096332  connected    <- labeled Robinson: C'(x) = exp(b'(x C'(x))) on the projected B
  A066537  all          <- exp(C)
Usage: python labeled_check.py N NPRIMES
"""
import json
import os
import sys

import numpy as np

import tower_p as T
import tower_gpu as TG

BDIR = os.path.expanduser('~/oeis-mcp/data/oeisdata/files')


def bfile(anum):
    d = {}
    with open(os.path.join(BDIR, anum[:4], 'b' + anum[1:] + '.txt')) as f:
        for line in f:
            s = line.split()
            if len(s) >= 2 and not s[0].startswith('#'):
                d[int(s[0])] = int(s[1])
    return d


def pmul(a, b, N, q):
    r = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                r[i + j] = (r[i + j] + x * b[j]) % q
    return r


def pexp(f, N, q):                      # f[0] == 0
    r = [0] * (N + 1)
    r[0] = 1
    df = [i * f[i] % q for i in range(N + 1)]
    for n in range(1, N + 1):
        r[n] = sum(df[k] * r[n - k] for k in range(1, n + 1)) % q * pow(n, q - 2, q) % q
    return r


def compose(f, g, N, q):                # f(g(x)), g[0] == 0
    r = [0] * (N + 1)
    for c in reversed(f):
        r = pmul(r, g, N, q)
        r[0] = (r[0] + c) % q
    return r


def main():
    N, npr = int(sys.argv[1]), int(sys.argv[2])
    ref = {a: bfile(a) for a in ('A096330', 'A096331', 'A096332', 'A066537')}
    fact = [1] * (N + 1)
    for i in range(1, N + 1):
        fact[i] = fact[i - 1] * i
    bad = 0
    for q in T.PRIMES24[:npr]:
        rows = json.load(open(f'kf_{N}_{q}.json'))
        t3 = [0] * (N + 1)
        for n, vp, eo, w in rows:
            if all(x == 1 for x in vp):
                t3[n] = (t3[n] + w) % q
        ps = TG.PSG(N, q)
        R, Hid, Htau = T.rows_from_zk(rows, q)
        core = TG.CoreG(ps, R, Hid, Htau)
        B, _ = T.blocks(ps, core)
        B = TG.cp.asnumpy(B)
        b = [0] * (N + 1)
        for n in range(1, N + 1):
            b[n] = int(B[ps.idx[(1,) * n]]) % q
        # labeled connected: rooted R(x) = x C'(x) = x exp(b'(R(x)))
        bp = [(k + 1) * b[k + 1] % q for k in range(N)] + [0]
        Rr = [0] * (N + 1)
        Rr[1] = 1
        for _ in range(N + 2):
            e = pexp(compose(bp, Rr, N, q), N, q)
            new = [0] + e[:N]
            if new == Rr:
                break
            Rr = new
        c = [0] + [Rr[n] * pow(n, q - 2, q) % q for n in range(1, N + 1)]
        g = pexp(c, N, q)
        res = {'A096330': t3, 'A096331': b, 'A096332': c, 'A066537': g}
        line = []
        for a, ser in res.items():
            ok = 0
            for n in range(1, N + 1):
                if n in ref[a]:
                    if ser[n] * fact[n] % q == ref[a][n] % q:
                        ok += 1
                    elif ref[a][n] or ser[n]:
                        bad += 1
                        print(f'  MISMATCH {a}({n}) mod {q}', flush=True)
            line.append(f'{a} {ok}/{sum(1 for n in range(1, N + 1) if n in ref[a])}')
        print(q, '  '.join(line), flush=True)
    print('ALL OK' if not bad else f'{bad} MISMATCHES')


if __name__ == '__main__':
    main()
