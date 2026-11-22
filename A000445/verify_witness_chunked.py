"""Independent, memory-light check of a cpr SAT witness for large N (numpy, chunked).

Same rule as verify_witness.py: chi is given on the primes in the witness file ('p:v' lines);
every larger prime q gets the free value chi(q) = -w, where w is the smallest value not taken
by chi on 1..floor((N+L-2)/q). chi(n) is then computed per chunk by trial division with the
primes up to sqrt(X), so memory stays at a few hundred MB whatever X is.

Checks that no r in [1, N-1] has chi(r) = chi(r+1) = 0 (L = 2), and with --extend X reports the
first such r in [N, X).

usage: python verify_witness_chunked.py k N witness.txt [--extend X]
"""
import sys

import numpy as np


def main():
    k, N, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    X = int(sys.argv[sys.argv.index("--extend") + 1]) if "--extend" in sys.argv else N
    X = max(X, N)
    given = {}
    for line in open(path):
        line = line.strip()
        if line:
            p, v = line.split(":")
            given[int(p)] = int(v) % k
    P = max(given)
    M = N  # N + L - 2 with L = 2

    # primes up to max(P, sqrt(X)) by a plain sieve
    R = max(P, int(X ** 0.5) + 1)
    s = np.ones(R + 1, bool)
    s[:2] = False
    for i in range(2, int(R ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    primes = np.nonzero(s)[0]
    wp = [int(p) for p in primes if p <= P]
    assert set(wp) == set(given), "witness must list exactly the primes up to its largest one"

    chiP = np.zeros(R + 1, np.int64)  # chi on primes <= R (given, or free for P < q <= R)

    # chi on 1..T for the free rule (T = M // (P+1) is the largest lim a free prime can see)
    T = M // (P + 1) + 1
    small = np.zeros(T + 1, np.int64)
    for n in range(2, T + 1):
        m, c = n, 0
        for p in wp:
            if p * p > m:
                break
            while m % p == 0:
                m //= p
                c += given[p]
        if m > 1:
            c += given[m]
        small[n] = c % k
    # free value for each lim: -(smallest w not in chi(1..lim))
    freev = np.zeros(T + 1, np.int64)
    seen = set()
    for lim in range(0, T + 1):
        if lim >= 1:
            seen.add(int(small[lim]))
        miss = [w for w in range(k) if w not in seen]
        freev[lim] = (-miss[0]) % k if miss else -1  # -1: no free value (must never be used)

    for p in wp:
        chiP[p] = given[p]
    big_primes = primes[primes > P]
    chiP[big_primes] = freev[np.minimum(M // big_primes, T)]
    assert (chiP[big_primes] >= 0).all(), "a prime above the witness has no free value"
    nfree_not1 = int(np.count_nonzero(chiP[big_primes] != 1)) if big_primes.size else 0

    def chi_big(q):  # prime q > R
        lim = np.minimum(M // q, T)
        return freev[lim]

    sq = primes[primes <= int(X ** 0.5) + 1]
    C = 5_000_000
    prev_zero = False
    bad_below_N = None
    first_after = None
    for a in range(1, X + 1, C):
        b = min(a + C, X + 1)  # chunk [a, b)
        rem = np.arange(a, b, dtype=np.int64)
        c = np.zeros(b - a, np.int64)
        for p in sq:
            p = int(p)
            if p * p >= b:
                break
            idx = np.arange((-a) % p, b - a, p)
            while idx.size:
                rem[idx] //= p
                c[idx] += chiP[p]
                idx = idx[rem[idx] % p == 0]
        big = rem > 1
        q = rem[big]
        vals = np.where(q <= R, chiP[np.minimum(q, R)], chi_big(q))
        assert (vals >= 0).all(), "a prime above the witness has no free value"
        c[big] += vals
        z = (c % k) == 0
        if a == 1:
            z[0] = True  # chi(1) = 0
        starts = []
        if prev_zero and z[0]:
            starts.append(a - 1)
        w = np.nonzero(z[:-1] & z[1:])[0]
        starts.extend((a + w).tolist())
        for r in starts:
            if r < N and bad_below_N is None:
                bad_below_N = r
            if r >= N and first_after is None:
                first_after = r
        prev_zero = bool(z[-1])
        print(f"ok to {b - 1}", flush=True)
        if bad_below_N is not None or (first_after is not None):
            break
    print(f"free primes above {P} with chi != 1: {nfree_not1}")
    if bad_below_N is not None:
        print(f"k={k} N={N}: FAIL, zero pair at r = {bad_below_N}")
        return 1
    print(f"k={k} N={N}: no zero pair with r < N: OK")
    if X > N:
        print("first zero pair at r =", first_after)
    return 0


if __name__ == "__main__":
    sys.exit(main())
