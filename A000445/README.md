# A000445: a(8) = 1508324

**Status:** new term found, with a machine-checkable certificate. Not yet submitted.
**Known terms:** a(2)..a(7) = 9, 77, 1224, 7888, 202124, 1649375 (Lehmer, Lehmer, Mills, Dunton, Bierstedt–Mills;
Brillhart–Lehmer–Lehmer 1964). Reble (SeqFan, Dec 2019) showed 1499876 <= a(8) <= 1508324.
**Result:** a(8) = 1508324, so Reble's upper bound is the true value. The extremal pair is

    1508324 = 2^2 * 47 * 71 * 113,    1508325 = 3 * 5^2 * 7 * 13^2 * 17.

**Visual page:** [Consecutive Residues](https://guillochon.github.io/consecutive-residues/) shows the record run as a grid,
the maximal characters for k = 2..8, the first pair for ordinary primes, and statistics of the proof tree.
Source: `viz/A000445/` (`build.py`, then `publish.py`).

## What is computed
For a prime p, let r(p) be the least r >= 1 such that r and r+1 are both n-th power residues mod p.
a(n) is the largest r(p) over all primes that have such a pair. (BLL 1964 write it Λ(n); it is Erdős problem #436's Λ(n,2).)

**Reduction to characters.** Let d = gcd(n, p-1). If d < n, the n-th power residues mod p are the d-th power residues,
so r(p) <= a(d) < a(n). For p ≡ 1 (mod n), fix a primitive root g and put χ(m) = ind_g(m) mod n. Then χ is completely
additive, and m is an n-th power residue exactly when χ(m) = 0. For n = 8, p ≡ 1 (mod 8) makes 2 a quadratic residue,
so χ(2) ∈ {0, 2, 4, 6}. Any choice of values on finitely many primes subject to that one condition is realized
by infinitely many primes p. This follows from Kummer theory over Q(ζ_8) together with Chebotarev; it is the fact
BLL 1964 use, citing Mills, "Characters with preassigned values". Hence

    a(8) = max over χ : primes → Z/8 with χ(2) even, of f(χ) = min{ r : χ(r) = χ(r+1) = 0 }.

Small primes don't change this. If p > R+1, the values mod p on 1..R+1 come from such a χ. If p <= R+1, then r(p) < p <= R+1.

**Decision form.** D(N) asks: is there a χ with no zero pair at any r in [1, N-1]? Then a(8) = max{N : D(N)}.
We show D(1508324) holds (an explicit χ) and D(1508325) fails (a checked proof tree).

## Algorithm (`rust/src/main.rs`, about 400 lines)
- DFS over primes in increasing order. Assigning q fixes χ on every m <= N whose largest prime factor is q.
- Forward checking. Take m whose largest prime Q is unassigned but whose cofactor m/Q^e is fixed. If a neighbor m±1
  is a fixed zero, every value of χ(Q) that would make χ(m) = 0 is struck from dom(Q). An empty domain is a conflict.
- Free primes (makes the SAT side finite). Take q > sqrt(N). Every multiple of q up to N has the form j·q with j <= N/q < q,
  so χ(j) is already fixed. If χ(1..N/q) misses a value v, set χ(q') = -v for q and every larger prime q'.
  Then none of their multiples is a residue, and no new zero pair can appear.
- Symmetry. χ → uχ for a unit u of Z/8 preserves zero sets and the allowed set for χ(2). At each prime the search tries
  one value per orbit of the stabilizer of the current assignment.

This is the BLL 1964 "case vector" proof tree with modern propagation. A run takes about 30 s on one core.

## Checks
1. **Old terms reproduced** by the same program: a(2..7) = 9, 77, 1224, 7888, 202124, 1649375.
2. **Lower bound (SAT witness).** `results/witness_k8_1508324.txt` gives χ on the primes <= 37699; the rest follow the
   free rule. `verify_witness.py` is independent Python. It rebuilds χ on 1..1.6M, confirms there is no zero pair below
   1508324, and reports the first pair at r = 1508324 exactly. The case vector, in BLL's notation R(q) = χ(q):
   - R(q) = 0 for q = 3, 97, 149
   - R(q) = 2 for q = 2, 13, 29, 31, 53, 101, 113, 137, 223, 241, 257, 271, 313, 317, 337, 359, 499, 521, 593, 641, 673, 709,
     883, 1021, 1061, 1361, 1381
   - R(q) = 3 for q = 233, 269, 443, 479, 487, 727, 937, 1171, 1201, 1297, 1301, 1409, 1451, 1543, 1663, 2297, 2381, 2383, 2389,
     2411, 3677, 3821, 4073, 4091, 4177, 4547, 4801, 5407, 5867, 6343, 6733, 6761, 8543, 10271, 15227, 16103, 24049, 24677,
     26029, 32917, 33331
   - R(q) = 1 for every other prime q <= 37699; larger primes are set by the free rule.

   (Like BLL's vectors, this is one maximal vector; others exist.)
3. **Upper bound (UNSAT certificate).** `CERT=... cpr decide 8 1508325` writes the full proof tree,
   `results/cert_k8_1508325.txt.gz` (175,629 branch nodes, 1,229,269 leaves). Every leaf names a pair (r, r+1),
   r < 1508325. Every value struck by forward checking also becomes a leaf with its pair.
   `verify_cert.py` is an independent Python checker. For every leaf it factors r and r+1, checks that all their primes
   are assigned on that path, and checks that χ(r) = χ(r+1) = 0. At every branch it checks that the listed values meet
   every orbit of the allowed values under the current stabilizer. Result: **VALID**. The checker also correctly rejects
   tampered trees (a changed leaf and a dropped branch value).
   The certificates for a(2)..a(6) (`results/cert_k*.txt`) also check VALID.
4. **Consistency with earlier work.** BLL 1964 found a(8) >= 1200744 and Reble 2019 found 1499876 <= a(8) <= 1508324.
   Our D(1499877) is SAT, so the lower bound improves. Our D(1508325) is UNSAT, which agrees with Reble's upper bound.
   The binary search (`results/max_k8.log`) is monotone throughout.

## Reproduce
```
cd rust && cargo build --release
./target/release/cpr max 7 1000000 3000000                         # -> a(7)=1649375
./target/release/cpr decide 8 1508324 ../results/witness.txt          # SAT
CERT=../results/cert.txt ./target/release/cpr decide 8 1508325        # UNSAT + proof tree
python ../verify_witness.py 8 1508324 ../results/witness.txt --extend 1600000
python ../verify_cert.py ../results/cert.txt
```

## About "recent 3SAT advances"
The subexponential worst-case 3SAT bounds (bigospeedrun.com/3sat, Oct 2026) were not useful here. With about 10^4–10^5
variables, 2^{n/sqrt(log n)} is still astronomically large. This instance is easy because of its structure: smooth
consecutive pairs drive the propagation, and the free-prime argument removes the large primes. The whole a(8) refutation
explores 78k nodes.

## Possible next steps
- a(9): BLL 1964 showed a(9) > 10^7. The same code should handle it, since the work scales with the number of
  smooth pairs rather than with N. Expect N around 10^7–10^8, with memory for the arrays of about 30 bytes × N.
- a(10) >= 22458303 (BLL), with constraint R(5) even.
- Erdős #436 lists Λ(k,2) only up to k = 7; a(8) could be reported there.

## Log
- 2026-11-07: read BLL 1964 and the Erdős #436 page (the Reble SeqFan post was unreachable because the archive was down).
  Solver written; a(2..7) reproduced. k=8: SAT at 1499877, UNSAT at 1508325; binary search → a(8) = 1508324.
  SAT witness verified, UNSAT certificate verified (also from the gzip file).
