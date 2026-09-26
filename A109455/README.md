# A109455 — Equivalence classes of threshold functions under permutations of the variables

**Status:** idea. **The existing a(6) appears to be wrong** (proof below).
**Known terms:** a(0)..a(6) = 2, 4, 10, 34, 178, 1720, 590440 (Knuth, 2005).
**Targets:**
1. Correct a(6). The true value must lie in [20,872, 71,232]. One guess: a typo for **59044**, with an extra 0.
2. New terms a(7) and a(8).

**Literature check (2026-09-25):** nothing found on permutation-only classes. Related sequences with data:
- A000609, labeled threshold functions, known to n = 9;
- A000617, NP-classes, known to n = 9;
- A002078, N-classes;
- A001529, NPN-classes.

## What is being counted
A *threshold function* is f(x) = [w_1·x_1 + … + w_n·x_n ≥ t] for some real weights w and threshold t. Negative weights are allowed. Examples: AND, OR, majority, x_1 ∧ ¬x_2.
We count them up to **renaming variables only** (the group S_n); negating variables is not allowed.

## Why the listed a(6) = 590440 is impossible
- A000617(6) = 1113 is the number of classes up to renaming **and negating** variables. That group is bigger than S_n by a factor of 2^6 = 64.
- Each class of the bigger group splits into at most 64 classes of S_n (index argument).
- Therefore **a(6) ≤ 64 × 1113 = 71,232**.
- Lower bound: a(6) ≥ A000609(6)/6! = 15,028,134/720 ≈ 20,872.
- Sanity check on the other terms: a(5) = 1720 fits [94572/120 ≈ 788, 32 × 119 = 3808]. a(4) and a(3) fit their windows too.

The only escape would be if A109455 and A000617 counted different universes of functions, but both say "threshold functions of n variables". This is worth reporting as an OEIS correction together with a freshly computed correct value.

## How the known terms were probably computed
Knuth probably generated all threshold functions and reduced them by symmetry, which is fine up to n = 6 (15 million functions). For n = 7 there are 8.4×10^9 labeled functions, and for n = 8 there are 1.76×10^13.

## The clever idea: Burnside plus "threshold functions on a grid"
Burnside's lemma: a(n) = (1/n!) Σ_{σ ∈ S_n} Fix(σ), where Fix(σ) = # threshold functions with f∘σ = f. Only the cycle type of σ matters: 15 types for n = 7, 22 for n = 8.

**Key lemma (averaging):** if f∘σ = f and w separates f, then w∘σ also separates f. The set of separating (w,t) is convex, so the average over the cycle group separates f too. Hence **a σ-invariant threshold function always has weights that are constant on each cycle of σ**.

So if σ has cycles of lengths ℓ_1..ℓ_k, a σ-invariant threshold function depends only on the counts c_j = (number of 1s in cycle j) ∈ {0..ℓ_j}. It is a threshold function on the **grid** G = {0..ℓ_1} × … × {0..ℓ_k} ⊂ R^k:

    Fix(σ) = number of linearly separable 0/1 labelings of the grid G_λ   (λ = cycle type)

- σ = id: the grid is the cube {0,1}^n, so Fix = A000609(n), already known.
- Every other σ has k < n dimensions, which is a strictly smaller and easier problem.
- The hardest non-identity class is a single transposition: k = n−1 dimensions with 3·2^(n−2) points.

Counting separable labelings of a point set means counting the regions of a hyperplane arrangement in (w,t)-space, one hyperplane per point. Methods:
- [ ] Incremental region enumeration, adding one point at a time. This is the standard threshold-function enumeration algorithm, adapted to grids.
- [ ] Zaslavsky's theorem: regions = |χ(−1)|, where χ is the characteristic polynomial, computed via the intersection lattice or the Athanasiadis finite-field method.
- [ ] Reuse known results. Grid threshold functions ("threshold functions on k-valued logic / digital halfplanes") have their own literature, and some counts may already exist.

## Plan
1. Compute a(0..5) by Burnside with grid counting. They must match the listed values.
2. Compute a(6) and confirm it lies in [20872, 71232]; is it 59044?
3. As a second method for a(6), enumerate all 15,028,134 threshold functions of 6 variables directly and reduce by S_6.
4. Compute a(7): 15 cycle types, and every non-identity grid has ≤ 6 dimensions.
5. Compute a(8): the transposition class, with a 7-dimensional grid of 192 points, is the main computational job.

## Verification plan
- Two independent methods for a(6) (Burnside and direct enumeration).
- For n = 7, cross-check with an NP-class splitting argument. If the lists of NP-class representatives exist (weighted-game lists by Kurz or Tautenhahn), the number of S_n-classes inside each NP-class is the number of double cosets S_n \ G / Stab_G(f).
- Every Fix(σ) must be ≤ A000609(n), and the result must lie in [A000609(n)/n!, 2^n·A000617(n)].

## Files and how to run

| File | What it does |
|---|---|
| `bruteforce.py` | Generates all threshold functions for n ≤ 6 and counts orbits (method 1). |
| `engine.py` | Burnside plus grid engine in Python (method 2). `python A109455/engine.py 7 --from 7 --procs 12` |
| `rust/` | Compiled port of `engine.py`, meant for a(8), a(9) and beyond. |
| `layman.md` | Plain-English explanation of the sequence, why it matters, and what we found. |
| `logs/` | Run logs. |

Rust engine (run from `A109455/rust`):
```
cargo build --release
./target/release/threshold-burnside 8 --from 0            # everything, identity included (validation)
./target/release/threshold-burnside 9 --id-known          # identity term from A000609(9)
./target/release/threshold-burnside grid 2 1 1 1 1 1 1 1  # one grid only
```
- Results go to `results/n<N>.tsv`, one row per finished cycle type.
- Every finished subtree is also appended to `results/grid_*.tasks`.
- **Rerunning the same command resumes**: finished classes and subtrees are skipped.
- **Progress:** a heartbeat line goes to stderr every 30 s (`--progress SECS`, 0 turns it off). It shows the class number, subtrees done/total, canonical functions found and their rate, LP count, and an ETA. Example:
  `[ 120s] class 21/22 n=8 2-1-1-1-1-1-1 [searching]  subtrees 812/4096 (19.8%)  canonical 1.234e6 (1.1e4/s)  LPs 2.5e6  class 95s  ETA ~380s`
- To run in the background and follow along (PowerShell, from `A109455/rust`):
  ```
  Start-Process -NoNewWindow ./target/release/threshold-burnside.exe -ArgumentList "9 --id-known --threads 14" -RedirectStandardError ../logs/n9.err -RedirectStandardOutput ../logs/n9.log
  Get-Content ../logs/n9.err -Wait -Tail 5
  ```
- The Python engine prints a similar heartbeat every 30 s (`HEARTBEAT_SECS` in `engine.py`).
- **Do not move or rename `engine.py` while a Python run is going.** On Windows, multiprocessing workers re-load the script from its path, and the run crashes when a new pool starts.
- The split into subtrees does not depend on the thread count, so a run can be resumed with a different `--threads`.

Beyond a(9): the identity term A000609(10) (threshold functions of 10 variables) is **not known**, and it is a hard, well-studied problem. So a(10) needs either that number or a new idea.

## Log
- 2026-09-25: Found that the listed a(6) violates the bound a(6) ≤ 2^6·A000617(6) = 71,232. Averaging lemma plus grid reduction sketched.
- 2026-09-25: **Computed a(6) = 38456** with `bruteforce.py` (Burnside, remainder 0). The typo guess 59044 was wrong. Pipeline checks:
  - positive threshold functions match A002078(0..6);
  - all threshold functions match A000609(0..6) (15,028,134 at n = 6), so the generated set is complete;
  - a(0..5) match the listed values;
  - Burnside and direct canonical-form counts agree for n ≤ 5.

  Fix(σ) for n = 6, by cycle type:

  | cycle type | 1^6 | 2·1^4 | 2²·1² | 2³ | 3·1³ | 3·2·1 | 3² | 4·1² | 4·2 | 5·1 | 6 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | Fix | 15028134 | 643854 | 36566 | 2702 | 23418 | 1998 | 174 | 1158 | 150 | 94 | 14 |

  Fix(6-cycle) = 14 = the number of threshold cuts of 7 collinear points, as the grid lemma predicts. 38456 appears nowhere else in the OEIS.
- 2026-09-25: Built `engine.py`: Burnside over cycle types, using the grid lemma, flip and symmetry reduction, and point-by-point search with a carried witness plus LP. On n = 6 it reproduces all 11 brute-force Fix values, so **a(6) = 38456 is confirmed by two independent methods**.
- 2026-09-25: **a(7) = 2490634** (log: `logs/n7.log`, 163 s on 12 processes). Checks:
  - the identity term, computed by the engine, reproduces A000609(7) = 8378070864;
  - its canonical count, 29375, equals A000617(7);
  - the Burnside remainder is 0;
  - the value lies inside [A000609(7)/7!, 2^7·A000617(7)] = [1662315, 3760000];
  - bad = 0, meaning every function was certified by an explicit weight vector.
- 2026-09-25: Rust port (`rust/`) reproduces every Fix value for n ≤ 7. It is about 100× faster than Python after the iterative-deepening split fix.
- 2026-09-25: **a(8) = 550112272** (`logs/n8_rust.log`, `rust/results/n8.tsv`, 433 s on 12 threads). Checks:
  - the identity term reproduces A000609(8) = 17561539552946;
  - the canonical count 2730166 = A000617(8) and the positive count 68863243522 = A002078(8);
  - the Burnside remainder is 0;
  - the value lies inside [435554056, 698922496];
  - bad = 0.

  The Python run crashed after 17/22 classes because `engine.py` was moved mid-run; all 17 of its canonical counts match Rust. Load imbalance seen: one subtree of the identity grid ran alone for about 3 minutes.
- 2026-09-25: Engine performance fixes, found while watching the a(9) run:
  1. **Load balancing:** replaced rayon with a shared work queue plus an explicit DFS stack. Idle workers receive the *shallowest* unexplored branch, and each checkpoint "family" is recorded once all work donated from it is done. Rayon's `join` rarely got stolen here because the heavy branch is usually the second one.
  2. **LP size:** only *maximal* branch-decided zeros go into the LP. q ≤ p for all feasible weights iff the within-group prefix sums of q are ≤ those of p. Deep nodes had hundreds of redundant rows (LPs up to ~700 µs each).

  n = 8 identity grid: 403 s → 95 s on 8 threads, same exact counts (A000609(8)); a(0..7) re-validated. The a(9) run was resumed from its checkpoints with the new binary.
- 2026-09-25: **Root cause of the "only ~3.3 cores busy" slowdowns:** this machine has an i5-12600K (6 P-cores + 4 E-cores). Windows 11 treated the engine, launched in the background, as EcoQoS work and kept it on the 4 E-cores; the P-cores sat idle. The engine now calls `SetProcessInformation(ProcessPowerThrottling)` at startup to opt out. Measured on the n = 8 transposition grid: 3.3 → 7.2 of 8 threads busy, 25 s, same exact counts. Also removed per-LP heap allocations (stack buffers). The a(9) run was resumed with this binary.
