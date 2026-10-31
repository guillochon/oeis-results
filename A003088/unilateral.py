"""A003088: unilateral digraphs on n unlabeled nodes, by Burnside over chains of strong components.

Cycle-index series are dicts {partition (sorted tuple, descending): coefficient of p_lambda}, with the
Burnside normalization coefficient = (# structures fixed by a permutation of that cycle type) / z_lambda.
All arithmetic is modulo a Mersenne prime larger than 2^(n^2), so results are exact integers.

1. Digraphs: d_lambda = 2^(g(lambda,lambda) - #cycles) / z_lambda, g(l,m) = sum gcd(a,b) over a in l, b in m
   (orbits of a permutation on ordered pairs; the loops (u,u) form one orbit per cycle).
2. Robinson's source-component identity, with one sign per sigma-orbit of source components:
   E (.) D = 1, where (A (.) B)_nu = sum_{lambda + mu = nu} a_lambda b_mu 2^g(lambda,mu) and
   E = exp(-sum_m Z_S[p_i -> p_im] / m) is the signed set-of-strong-components series.
3. Z_S = sum_m mu(m)/m L[p_i -> p_im], L = -log E: the cycle index of strong digraphs (A035512 = sum of its
   coefficients, A003030(n)/n! = coefficient of p_1^n).
4. Unilateral: the strong components form a chain C_1 -> ... -> C_k (unique order, so every automorphism fixes
   each one); arcs go forward only, any between non-consecutive components, at least one between consecutive.
   a(n) = sum over (lambda_1..lambda_k) of prod c_lambda_i * prod_{i<j} 2^g(lambda_i,lambda_j)
          * prod_i (1 - 2^-g(lambda_i,lambda_i+1)), as a DP over (union so far, last block).

Usage: python unilateral.py N   -> a(0..N) of A003088, plus the A035512 / A003030 / A003029 checks
"""
import sys
import time
from functools import lru_cache
from math import factorial, gcd

MERSENNE = [61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423]


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def merge(a, b):
    return tuple(sorted(a + b, reverse=True))


def submultisets(nu):
    """All (lambda, nu - lambda) splits of the multiset nu."""
    items = sorted(set(nu), reverse=True)
    counts = [nu.count(x) for x in items]

    def rec(i):
        if i == len(items):
            yield (), ()
            return
        for take in range(counts[i] + 1):
            for a, b in rec(i + 1):
                yield (items[i],) * take + a, (items[i],) * (counts[i] - take) + b
    yield from rec(0)


def g(l, m):
    return sum(gcd(a, b) for a in l for b in m)


def z(lam):
    out = 1
    for k in set(lam):
        c = lam.count(k)
        out *= k ** c * factorial(c)
    return out


def mobius(n):
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    return -res if m > 1 else res


def cycle_index(N):
    """(P, parts, Z, strong, strong_lab): Z = strong-digraph cycle index mod the prime P, to degree N."""
    P = (1 << next(e for e in MERSENNE if e > N * N + 8)) - 1
    inv = lambda x: pow(x % P, P - 2, P)
    t0 = time.time()
    parts = {n: list(partitions(n)) for n in range(N + 1)}

    # 1. digraphs
    D = {lam: pow(2, g(lam, lam) - len(lam), P) * inv(z(lam)) % P for n in range(N + 1) for lam in parts[n]}

    # 2. E = D^(.)-1
    E = {(): 1}
    for n in range(1, N + 1):
        for nu in parts[n]:
            s = 0
            for lam, mu in submultisets(nu):
                if mu:
                    s += E[lam] * D[mu] * pow(2, g(lam, mu), P)
            E[nu] = -s % P

    # 3. L = -log E (graded: n L_n = n E_n - sum_{k<n} k L_k E_{n-k}, with E -> -E sign handled below)
    #    log E: n*Lg_n = n*E_n - sum_{k=1}^{n-1} k*Lg_k*E_{n-k}; L = -Lg.
    Lg = {}
    for n in range(1, N + 1):
        acc = {nu: n * E[nu] % P for nu in parts[n]}
        for k in range(1, n):
            for a in parts[k]:
                la = Lg[a]
                if la:
                    for b in parts[n - k]:
                        nu = merge(a, b)
                        acc[nu] = (acc[nu] - k * la * E[b]) % P
        ninv = inv(n)
        for nu in parts[n]:
            Lg[nu] = acc[nu] * ninv % P
    L = {nu: -v % P for nu, v in Lg.items()}

    # Z_S = sum_m mu(m)/m L[m]
    Z = {nu: 0 for n in range(1, N + 1) for nu in parts[n]}
    for m in range(1, N + 1):
        mu_m = mobius(m)
        if not mu_m:
            continue
        w = mu_m * inv(m) % P
        for n in range(1, N // m + 1):
            for lam in parts[n]:
                Z[tuple(m * x for x in lam)] = (Z[tuple(m * x for x in lam)] + w * L[lam]) % P

    def signed(x):
        return x if x <= P // 2 else x - P

    strong = [1] + [signed(sum(Z[nu] for nu in parts[n]) % P) for n in range(1, N + 1)]
    strong_lab = [0] + [signed(Z[(1,) * n] * factorial(n) % P) for n in range(1, N + 1)]
    print(f"# cycle index done in {time.time() - t0:.1f}s", file=sys.stderr)
    return P, parts, Z, strong, strong_lab


def signed_mod(x, P):
    return x if x <= P // 2 else x - P


def direct(N, P, parts, Z):
    """Check of step 4: the sum over sequences of blocks term by term (no DP), for small N."""
    pw = lambda e: pow(2, e, P)

    def rec(left, prev, union):
        if left == 0:
            return 1
        s = 0
        for m in range(1, left + 1):
            for mu in parts[m]:
                if Z[mu]:
                    w = Z[mu] * pw(g(union, mu))
                    if prev is not None:
                        w = w * (1 - inv2 ** g(prev, mu)) % P
                    s += w * rec(left - m, mu, merge(union, mu))
        return s % P
    inv2 = pow(2, P - 2, P)
    return [signed_mod(rec(n, None, ()), P) for n in range(N + 1)]


def main(N):
    t0 = time.time()
    P, parts, Z, strong, strong_lab = cycle_index(N)
    signed = lambda x: signed_mod(x, P)

    # 4. chain DP: F[(Q, lam)] with Q = union of the blocks before the last one, lam = last block.
    pw = lambda e: pow(2, e, P)
    F = {}
    for n in range(1, N + 1):
        for lam in parts[n]:
            if Z[lam]:
                F[((), lam)] = Z[lam]
    total = [0] * (N + 1)
    total[0] = 1
    lab = [0] * (N + 1)  # labeled unilateral (A003029) from the identity-permutation states
    for size in range(1, N + 1):  # process states by |Q| + |lam| = size
        for (Q, lam), v in [(k, v) for k, v in F.items() if sum(k[0]) + sum(k[1]) == size]:
            total[size] = (total[size] + v) % P
            if all(x == 1 for x in Q + lam):
                lab[size] = (lab[size] + v) % P
            Pn = merge(Q, lam)
            for m in range(1, N - size + 1):
                for mu in parts[m]:
                    c = Z[mu]
                    if not c:
                        continue
                    gq, gl = g(Q, mu), g(lam, mu)
                    key = (Pn, mu)
                    F[key] = (F.get(key, 0) + v * c % P * (pw(gq + gl) - pw(gq))) % P
    print(f"# chain DP done in {time.time() - t0:.1f}s, {len(F)} states", file=sys.stderr)
    a = [signed(x) for x in total]
    lab = [1] + [signed(lab[n] * factorial(n) % P) for n in range(1, N + 1)]
    return a, strong, strong_lab, lab


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    a, strong, strong_lab, lab = main(N)
    print("A035512 (strong, unlabeled):", ",".join(map(str, strong)))
    print("A003030 (strong, labeled):  ", ",".join(map(str, strong_lab[1:])))
    print("A003029 (unilateral, lab.): ", ",".join(map(str, lab[1:])))
    print("A003088 (unilateral, unl.): ", ",".join(map(str, a)))
