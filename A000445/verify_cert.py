"""Check an UNSAT proof tree emitted by cpr (CERT=path cpr decide k N).

Claim: no chi: primes -> Z/k with chi(2) in TWOMASK avoids a zero window
chi(r) = chi(r+1) = ... = chi(r+L-1) = 0 for all 1 <= r <= N-1 (L = 2 pairs, L = 3 triples;
header field L, default 2).

Format (preorder):
  B q c      branch on prime q; c children follow, each "V v" + node
  X r        leaf: under the current partial assignment, chi(r) = ... = chi(r+L-1) = 0
Checks:
  * every leaf pair has 1 <= r < N and all prime factors of r, r+1 assigned, both chi = 0;
  * at each branch on q (unassigned), the listed values hit every orbit of the allowed values
    under the units u of Z/k that fix all values assigned so far (chi -> u*chi preserves zero sets
    and the allowed set for chi(2)).
usage: python verify_cert.py cert.txt
"""
import gzip
import sys
from math import gcd


def main():
    path = sys.argv[1]
    f = gzip.open(path, 'rt') if path.endswith('.gz') else open(path)
    head = f.readline().split()
    k, N, twomask = int(head[1]), int(head[3]), int(head[5], 2)
    L = int(head[7]) if len(head) > 7 else 2
    M = N + L - 2
    spf = list(range(M + 1))
    for i in range(2, int(M**0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, M + 1, i):
                if spf[j] == j:
                    spf[j] = i
    units = [u for u in range(1, k) if gcd(u, k) == 1]
    for u in units:
        for v in range(k):
            if (twomask >> v) & 1:
                assert (twomask >> ((u * v) % k)) & 1, 'TWOMASK not unit-invariant'
    toks = iter(f.read().split())
    assign = {}
    stats = {"leaves": 0, "branches": 0}

    def chi(n):
        s = 0
        while n > 1:
            p = spf[n]
            if p not in assign:
                raise AssertionError(f"prime {p} unassigned at leaf")
            s += assign[p]
            n //= p
        return s % k

    # iterative to avoid recursion limits: stack of (q, remaining children)
    def node():
        t = next(toks)
        if t == "X":
            r = int(next(toks))
            assert 1 <= r < N, f"leaf r={r} out of range"
            vals_ = [chi(r + j) for j in range(L)]
            assert all(x == 0 for x in vals_), f"leaf r={r}: chi={vals_} assign={assign}"
            stats["leaves"] += 1
            return
        assert t == "B", t
        q, c = int(next(toks)), int(next(toks))
        assert spf[q] == q and q not in assign, f"bad branch prime {q}"
        stats["branches"] += 1
        G = [u for u in units if all((u * v) % k == v for v in assign.values())]
        allowed = [v for v in range(k) if q != 2 or (twomask >> v) & 1]
        vals = []
        for _ in range(c):
            assert next(toks) == "V"
            v = int(next(toks))
            assert 0 <= v < k
            vals.append(v)
            assign[q] = v
            node()
            del assign[q]
        S = set(vals)
        for v in allowed:
            assert any((u * v) % k in S for u in G), f"branch on {q}: orbit of {v} not covered (G={G}, S={S})"

    sys.setrecursionlimit(100000)
    node()
    rest = list(toks)
    assert not rest, f"trailing tokens: {rest[:5]}"
    print(f"VALID: k={k} L={L} N={N}: every admissible character has a zero window r < {N}. "
          f"branches={stats['branches']} leaves={stats['leaves']}")


if __name__ == "__main__":
    main()
