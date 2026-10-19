"""Drive rust/target/release/dasasm: residues of a(n) modulo enough 31-bit primes, CRT, one spare prime as a check.

Usage: python run_rust.py NMIN NMAX [--w16] [--threads T]   (--w16: u16 storage, primes below 2^16, half the memory)
Appends "n a(n) primes seconds" to results/rust.txt, and compares with results/fast_1_50.txt where available and
with the Behrend-Fischer-Konvalinka product formula for odd n.
"""
import subprocess
import sys
import time
from pathlib import Path

from dasasm import odd_formula

HERE = Path(__file__).resolve().parent


def crt(residues, primes):
    x, m = 0, 1
    for r, p in zip(residues, primes):
        t = ((r - x) * pow(m, -1, p)) % p
        x += m * t
        m *= p
    return x, m


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    extra = sys.argv[3:]
    w16 = "--w16" in extra
    extra = [a for a in extra if a != "--w16"]
    bits = 16 if w16 else 31
    exe = HERE / "rust" / ("target16" if w16 else "target") / "release" / ("dasasm.exe" if sys.platform == "win32" else "dasasm")
    known = {}
    fast = HERE / "results" / "fast_1_50.txt"
    if fast.exists():
        known = {int(l.split()[0]): int(l.split()[1]) for l in fast.read_text().splitlines() if l.strip()}
    for n in range(nmin, nmax + 1):
        bound = odd_formula(n // 2)          # a(n) <= a(2*(n//2)+1)
        k = 1
        while (1 << ((bits - 1) * k)) <= 2 * bound:  # product of k primes > 2^((bits-1)k) > 2*bound
            k += 1
        k += 1                                 # spare prime
        t0 = time.time()
        out = subprocess.run([str(exe), str(n), "0", str(k), *extra], capture_output=True, text=True, check=True).stdout
        rows = [l.split() for l in out.splitlines() if l.strip()]
        ps = [int(r[1]) for r in rows]
        rs = [int(r[2]) for r in rows]
        x, m = crt(rs[:-1], ps[:-1])
        assert x % ps[-1] == rs[-1], f"n={n}: spare prime disagrees"
        notes = []
        if n in known:
            assert known[n] == x, f"n={n}: differs from numpy result"
            notes.append("=numpy")
        if n % 2:
            assert x == odd_formula(n // 2), f"n={n}: differs from odd formula"
            notes.append("=formula")
        line = f"{n} {x} {k} {time.time() - t0:.0f}"
        with open(HERE / "results" / "rust.txt", "a") as f:
            f.write(line + "\n")
        print(line, " ".join(notes), flush=True)


if __name__ == "__main__":
    main()
