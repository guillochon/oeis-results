"""Tower from formula rows: python run_tower.py N nprimes  (needs kf_N_<p>.json for the first nprimes
primes, from kformula.py).  Prints A021103, A003094, A005470 for n <= N with a CRT stability check
(the value from one prime fewer must agree) and compares with the OEIS terms in the snapshot."""
import json
import sys

import tower_p as T

KNOWN = {
    'A021103': [0, 0, 0, 1, 3, 9, 44, 294, 2893, 36496, 545808, 9029737, 159563559, 2952794985, 56589742050],
    'A003094': [1, 1, 1, 2, 6, 20, 99, 646, 5974, 71885, 1052805, 17449299, 313372298, 5942258308,
                117243822184, 2385519362107],
    'A005470': [1, 1, 2, 4, 11, 33, 142, 822, 6966, 79853, 1140916, 18681008, 333312451, 6295835195,
                123897775473, 2516143483576],
}


def main():
    N, npr = int(sys.argv[1]), int(sys.argv[2])
    primes = T.PRIMES24[:npr]
    res = []
    for q in primes:
        rows = json.load(open(f'kf_{N}_{q}.json'))
        res.append(T.run_mod(N, q, lambda p, rows=rows: T.rows_from_zk(rows, p)))
    report(N, primes, res)


def report(N, primes, res):
    out = {}
    for name, k in (('A021103', 0), ('A003094', 1), ('A005470', 2)):
        vals = []
        for n in range(N + 1):
            r = [int(x[k][n]) for x in res]
            a = T.crt(r, primes)
            assert a == T.crt(r[:-1], primes[:-1]), (name, n, 'CRT unstable')
            assert a >= 0
            vals.append(a)
        if name == 'A003094':
            vals[0] = 1
        known = KNOWN[name]
        tags = []
        for n, a in enumerate(vals):
            if n < len(known):
                tags.append('ok' if known[n] == a else f'MISMATCH({known[n]})')
            else:
                tags.append('NEW')
        print(name)
        for n, (a, t) in enumerate(zip(vals, tags)):
            print(f'  {n:3d} {a} {t}')
        out[name] = vals
    json.dump(out, open(f'tower_{N}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
