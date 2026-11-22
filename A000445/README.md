# A000445: a(8) = 1508324, and a(9) >= 550000096

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

## a(9) >= 550000096 (lower bound only)
BLL 1964 showed a(9) > 10^7. For odd k the Kummer extensions have full degree, so every assignment χ : primes → Z/9 is
realized by infinitely many primes p ≡ 1 (mod 9), including χ(3). (We also checked χ(3) empirically: over 11,727 primes
p ≡ 1 (mod 9), its 9 classes are as evenly spread as those of 2 and 11.)

**Early DFS runs (Nov 8).** The DFS above, even with conflict-directed backjumping, seemed to stall for k = 9. A run that
hits its time limit is not evidence of UNSAT: D(25770880) took 57 minutes to come back SAT, after 1-hour runs at nearby N
had timed out. Given enough time it does succeed: see the N = 3·10^8 run below.

**What worked.** Every good witness has χ(q) = 1 for almost all primes, so χ(n) is close to Ω(n) mod 9 and zero pairs
are rare. `ls/beam.py` searches over small exception sets E (prime → value). It computes χ on [1, X] with numpy, finds
the first zero pair (r, r+1), and tries changing χ at each prime dividing r or r+1. It keeps the 6 sets with the latest
first pair. A run with X = 3·10^8, 40 rounds, takes 2.5 hours. The search finds the true a(6) exactly, but only 77% of a(5)
and 17% of a(7), so its results are lower bounds, possibly far from the truth.

**Beam result** (`ls/best_k9.txt`). χ(q) = 1 for every prime except:
- χ = 0: 3, 11, 277, 509, 3631, 5903, 9257, 20773, 25261, 47939, 409729
- χ = 2: 127, 173, 179, 269, 347, 619, 907, 1747
- χ = 3: 769

The first zero pair is

    178708831 = 7^8 * 31,    178708832 = 2^5 * 19 * 41 * 67 * 107,

so a(9) >= 178708831. `ls/verify_chunked.py` checks this independently. It computes χ(n) by trial division over every
n <= 1.8·10^8, in chunks of 10^7, and reports the same first pair (`results/verify_k9_178708831.log`). An earlier
6-exception witness was also checked with `verify_witness.py`, giving no pair below 6·10^7.

**DFS at N = 3·10^8 (Nov 21).** After the memory cut (10 bytes per n, commit 6e756b6), `DEADLINE=43200 cpr decide 9
300000000` ran at below-normal priority and returned **SAT** after 3,541,060 nodes in 6.0 hours (6.1 ms/node;
`results/decide_k9_300000000.log`). The witness (`results/witness_k9_3e8.txt`) sets χ on the 90,900 primes up to
1,171,867, and larger primes follow the free rule. χ(q) = 1 for all of them except:
- χ = 0: 3, 13, 59, 71, 101, 137
- χ = 2: 71 primes from 107 to 713771 (listed in the witness file)
- χ = 3: 379, 2557

`verify_witness_chunked.py` is a new independent numpy checker. It uses the same free rule as `verify_witness.py`, but
works in chunks, so memory stays small. It confirms there is no zero pair below 3·10^8, and with this completion the first
pair is at r = 300000031, a prime above N (`results/verify_k9_3e8.log`). It reproduces the k=8 witness check exactly
(first pair 1508324). So **a(9) >= 300000031**. That is 1.7 times the beam bound, and the beam had again fallen short.
The bound is set only by N: a run at larger N would probably push it further.

**Faster solver and portfolio ladder (Nov 21-22).** Commit 79d15f3 makes the DFS 2.7 times faster at large N (k=9,
N = 10^8, 20,000 nodes: 74.7 s before, 27.2 s after). Each prime's list is now visited in one pass (set χ and check the
windows together, valid because q >= L), the exponent table is replaced by division, and both hot loops use software
prefetch. Node counts for k = 2..8 are unchanged, and the k=8 certificate is byte-identical. A new `portfolio` mode
shares the read-only tables between threads (about 16 bytes per n with 4 threads). Thread 0 runs the complete DFS, so an
UNSAT from it is still a proof. The other threads run randomized Luby restarts, and thread 1 starts from a hint witness.

`overnight.sh` climbed N with 4 threads at below-normal priority, using each witness as the next hint
(`results/overnight_k9.log`):

| N | result | time | nodes (winning thread) | first pair after completion |
|---|---|---|---|---|
| 3.5·10^8 | SAT, hint from 3·10^8 | 104 min | 1,569,515 | 350000086 |
| 4·10^8 | SAT | 84 s | 598,540 | 400000066 |
| 4.5·10^8 | SAT | 98 s | 612,149 | 450000016 |
| 5·10^8 | SAT | 109 s | 625,720 | 500000008 |

Each rung's winner was thread 1, following the hint. Once the hint covered the hard small primes, a rung cost only
about one node per new prime. (A single-thread warm start at 3.3·10^8 had looked stuck below p = 293 after 5 minutes,
so the 104-minute first rung was the hard part.) `verify_witness_chunked.py` checks every witness
(`results/verify_k9_*.log`): no zero pair below N, with the first pair just above N. The ladder ended at 00:08 because it
ran out of rungs, not time. So **a(9) >= 500000008** (`results/witness_k9_500000000.txt`, 145,713 primes up to 1,953,109).
Below 1,171,867 it is close to the 3·10^8 witness: χ = 0 at 3, 7, 71, 137, 293, 3217, and 141 other non-1 values. The
primes above that have only a few multiples below N, and the search set about 41,700 of them to values other than 1. The
free rule gives χ = 1 to every prime above the witness.

**Single-thread ladder (Nov 22).** With `THREADS=1` and a hint, `portfolio` now runs only the hinted restart search
(about 7.75 bytes per n, so 10^9 would fit in 8 GB). From the 5·10^8 witness, N = 5.5·10^8 was SAT in 11 minutes
(1,440,775 nodes; first pair after completion at 550000096, checked by `verify_witness_chunked.py`), so
**a(9) >= 550000096**. N = 6·10^8 then timed out after 8 hours (4.7 million nodes): the hint broke at a small
prime, and the restarts stalled around p = 600. Steps of 5·10^7 are too coarse once the hint stops carrying over.

The exact a(9) is open. It needs a refutation of D(a(9)+1), with a(9) now known to be above 5.5·10^8.

## How hard is an exact a(9)? (scaling runs, 2026-11-21)
Proof tree size of `cpr decide k a(k)+1` (UNSAT), with the SAT search at a(k) for comparison:

| k | a(k)+1 | UNSAT nodes | time | SAT nodes at a(k) |
|---|---|---|---|---|
| 2 | 10 | 2 | 0 s | |
| 3 | 78 | 5 | 0 s | |
| 4 | 1225 | 23 | 0 s | |
| 5 | 7889 | 319 | 0 s | 392 |
| 6 | 202125 | 1004 | 0.03 s | 1175 |
| 7 | 1649376 | 44712 | 8 s | 42156 |
| 8 | 1508325 | 77949 | 14 s | 49499 |

(Committed solver, at the memory-cut commit 6e756b6, which searches identically. The first scaling runs used an older
`cpr.exe` binary that was about 3 times slower and took a few hundred more nodes; the other node counts below come from
that binary.)

At fixed k the tree shrinks slowly as N grows past a(k) (k=8: 78k nodes at 1.5M, 58k at 1.2·10^7; k=7: 45k to 40k
from 1.6M to 6M), while time per node grows roughly linearly with N.

Fits of log(nodes) against k and log a(k) (k=3..8, residual factor about 2) predict only 1-3 million nodes for k=9, which
works out to roughly one day at N = 1.8·10^8. **That extrapolation is contradicted by the k=9 runs.** For k <= 8 the SAT
search at a(k) uses 60-100% of the nodes of the UNSAT tree at a(k)+1, and far below a(k) it is much cheaper (k=8 at
10^6: 20k nodes, 26% of the tree). For k=9, SAT searches at N = 1.5-2.6·10^7, at most 14% of the known lower bound,
already took 0.65-3.9 million nodes. So the real k=9 tree is at least tens of millions of nodes, and likely far more.
The N = 3·10^8 run confirms this. It needed 3.5 million nodes (6.1 ms each) just to find a solution, and a(9) is beyond
that N. A refutation tree above 3·10^8 is likely tens of millions of nodes or more: days to weeks of single-core time
at the very least, with no useful upper estimate.

## Possible next steps
- a(9): continue the hinted ladder past 5·10^8. Memory is about 5 bytes per n shared plus 2.75 per thread, so a 1- or
  2-thread run at 10^9 needs 8-10.5 GB. An exact value needs a far stronger refutation method.
- a(10) >= 22458303 (BLL), with constraint R(5) even.
- Erdős #436 lists Λ(k,2) only up to k = 7; a(8) could be reported there.

## Log
- 2026-11-07: read BLL 1964 and the Erdős #436 page (the Reble SeqFan post was unreachable because the archive was down).
  Solver written; a(2..7) reproduced. k=8: SAT at 1499877, UNSAT at 1508325; binary search → a(8) = 1508324.
  SAT witness verified, UNSAT certificate verified (also from the gzip file).
- 2026-11-08: a(9). DFS with backjumping verified witnesses up to N = 25770880, then stalled. Beam search over exception
  sets: a(9) >= 6·10^7 (verify_witness), then a(9) >= 178708831 (X = 3·10^8, checked by verify_chunked.py).
- 2026-11-21: scaling runs for a(9) cost; solver memory cut from 30 to 10 bytes per n (same search, same k=8 certificate);
  DFS at N = 3·10^8 SAT in 6 h, so a(9) >= 300000031 (checked by verify_witness_chunked.py).
- 2026-11-22: solver 2.7x faster plus a portfolio mode (79d15f3); overnight hinted ladder SAT at 3.5, 4, 4.5 and
  5·10^8, so a(9) >= 500000008 (all four witnesses checked by verify_witness_chunked.py). Single-thread ladder:
  5.5·10^8 SAT in 11 min, so a(9) >= 550000096; 6·10^8 timed out after 8 h.
