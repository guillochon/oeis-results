"""Fast A005162: the transfer matrix of dasasm.py on dense numpy arrays, modulo several primes, combined by CRT.

State array index: bit r (0-based) = running sum of open row r+1; while a column is being processed, the top bit
is the vertical running sum v. See dasasm.py for the model.

Usage: python dasasm_fast.py NMIN NMAX [NPRIMES]
"""
import sys
import time
import numpy as np

def _is_prime(n):
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):  # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _primes(k, top=1 << 62):
    out, q = [], top - 1
    while len(out) < k:
        if _is_prime(q):
            out.append(q)
        q -= 2
    return out


# primes below 2^62, so a sum of two residues fits in uint64
PRIMES = _primes(40)


def addmod(a, b, p):
    s = a + b
    s[s >= p] -= p
    return s


def count_mod(n, p):
    p = np.uint64(p)
    f = np.ones(1, dtype=np.uint64)
    R = 0  # open rows
    for c in range(1, n + 1):
        h = min(c, n + 1 - c)
        g = np.concatenate([f, np.zeros_like(f)])  # v = 0 is the new top bit
        for r in range(h - 1):
            # axes: v, bits above r, bit r, bits below r
            v4 = g.reshape(2, 1 << (R - 1 - r), 2, 1 << r)
            s = addmod(v4[0, :, 0, :], v4[1, :, 1, :], p)
            v4[0, :, 0, :] = s
            v4[1, :, 1, :] = s
        half = g.size // 2
        f = addmod(g[:half], g[half:], p)  # drop v
        if c < n + 1 - c:
            f = np.concatenate([f, f])      # open row c (free bit)
            R += 1
        elif c > n + 1 - c:
            q = f.size // 2
            f = addmod(f[:q], f[q:], p)     # close the top row
            R -= 1
    assert f.size == 1
    return int(f[0])


def crt(residues, primes):
    x, m = 0, 1
    for r, p in zip(residues, primes):
        t = ((r - x) * pow(m, -1, p)) % p
        x += m * t
        m *= p
    return x, m


def digits_needed(n):
    # |DASASM(n)| <= |DASASM(n+1)| and |DASASM(2k+1)| has a product formula: use it as an upper bound
    from dasasm import odd_formula
    k = n // 2
    return len(str(odd_formula(k)))


if __name__ == "__main__":
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    extra = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        need = digits_needed(n) // 18 + 1 + extra   # one or more spare primes as a consistency check
        res = [count_mod(n, p) for p in PRIMES[:need]]
        x, m = crt(res[:-extra], PRIMES[:need - extra])
        ok = all(x % p == r for r, p in zip(res[-extra:], PRIMES[need - extra:need]))
        print(n, x, "spare-prime check ok" if ok else "SPARE PRIME MISMATCH", f"{time.time() - t0:.1f}s", flush=True)
