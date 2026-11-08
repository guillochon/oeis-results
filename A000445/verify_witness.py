"""Independent check of a SAT witness from cpr: chi on small primes (file 'p:v' lines),
completed on larger primes by the free rule (pick chi(q) = -v with v not in chi({1..N//q})),
then verify no r in [1, N-1] has chi(r) = chi(r+1) = 0. Also reports the first zero pair
beyond N if asked (--extend X).

usage: python verify_witness.py k N witness.txt [--extend X]
"""
import sys


def main():
    k, N = int(sys.argv[1]), int(sys.argv[2])
    path = sys.argv[3]
    X = N
    if "--extend" in sys.argv:
        X = int(sys.argv[sys.argv.index("--extend") + 1])
    given = {}
    for line in open(path):
        line = line.strip()
        if line:
            p, v = line.split(":")
            given[int(p)] = int(v) % k
    if k == 8:
        assert given.get(2, 0) % 2 == 0, "chi(2) must be even for p = 1 mod 8"
    M = max(N, X)
    spf = list(range(M + 1))
    for i in range(2, int(M**0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, M + 1, i):
                if spf[j] == j:
                    spf[j] = i
    chi = [0] * (M + 1)
    pval = {}
    for n in range(2, M + 1):
        p = spf[n]
        if p == n:  # prime: value from witness or free rule
            if p in given:
                v = given[p]
            else:
                lim = N // p
                assert p * p > N or lim < k, f"prime {p} not in witness but not free-eligible"
                seen = set(chi[1:lim + 1])
                miss = [w for w in range(k) if w not in seen]
                assert miss, f"prime {p}: chi(1..{lim}) covers Z/{k}"
                v = (-miss[0]) % k
            pval[p] = v
            chi[n] = v
        else:
            chi[n] = (chi[n // p] + pval[p]) % k
    bad = [r for r in range(1, N) if chi[r] == 0 and chi[r + 1] == 0]
    print(f"k={k} N={N}: zero pairs with r<N: {len(bad)}", "OK" if not bad else bad[:5])
    if X > N:
        nxt = next((r for r in range(N, X) if chi[r] == 0 and chi[r + 1] == 0), None)
        print("first zero pair at r =", nxt)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
