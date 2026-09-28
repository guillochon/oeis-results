# A027624 — Number of independent vertex sets in the n-hypercube graph Q_n

**Status:** **a(7) = 78685477899897082403 computed 2026-09-26** (Rust, 0.75 s; see Results). **a(8) = 1268098993536094508894717661843009268823 computed 2026-09-27** (CUDA, 23.6 h on an RTX 3080; see "a(8): result"). OEIS edit drafted in `oeis_submission.md`.
**Known terms:** a(0)..a(6) = 2, 3, 7, 35, 743, 254475, 19768832143 (7 terms). a(6) is old; there is no b-file.
**Target:** a(7). Sapozhenko's asymptotic gives a(n) ~ 2·√e·2^(2^(n−1)), so a(7) ≈ 6×10^19, which probably **exceeds 2^64 ≈ 1.8×10^19**. Use 128-bit integers.
**Literature check (2026-09-25):** only asymptotic results (Sapozhenko; Galvin [1901.01991]; Jenssen–Perkins [1907.00862]). No exact a(7) found.

## What is being counted
Q_n has 2^n vertices, the binary strings of length n. Two vertices are adjacent when they differ in exactly one bit. An *independent set* is a set of vertices with no two adjacent, and the empty set counts.
Example: Q_2 is a square 00–01–11–10–00. Its independent sets are {}, the 4 single vertices, and the 2 diagonals {00,11} and {01,10}, so a(2) = 7.

## How the known terms were computed
The Maple and Mathematica programs in the entry are plain recursive branching: "take vertex v or not". The cost grows roughly like the answer itself, and a(6) ≈ 2×10^10 was already a stretch that way. At about 6×10^19, a(7) is **roughly 3×10^9 times more work** by branching. That rules out brute force.

## The clever idea: turn Q_7 into a product and count with matrices
**Key fact:** Q_7 = Q_5 × C_4, because Q_2 is the 4-cycle C_4.
Picture Q_7 as 4 copies of Q_5 arranged in a ring, S_0, S_1, S_2, S_3, where copy i is joined to copy i±1 vertex-by-vertex.
An independent set of Q_7 is a choice of independent sets S_0..S_3 of Q_5 such that neighbouring copies are **disjoint**: S_i ∩ S_{i+1} = ∅, with indices taken mod 4.

So with M the 0/1 "disjointness" matrix on the 254475 independent sets of Q_5:

    a(7) = trace(M^4) = Σ_{S,U} (M²[S,U])²,

where M²[S,U] = #{T independent in Q_5 : T ∩ (S∪U) = ∅} = f(complement of S∪U), and
**f(W) = the number of independent sets of Q_5 contained in the vertex set W**.

    a(7) = Σ over independent sets S, U of Q_5 of  f(~(S ∪ U))²

**Validated (2026-09-25):** the identical formula one dimension down (Q_d × C_4 = Q_{d+2}) reproduces a(4) = 743, a(5) = 254475 and a(6) = 19768832143 exactly. It is also checked via Q_d × K_2 = Q_{d+1}, using a(d+1) = Σ_S f(~S).

## Making it fast
1. **Computing f:**
   - *Option A (simple):* store f for all 2^32 subsets W with a subset-sum ("zeta") transform of the indicator of independent sets. That is 2^32 × uint32 = 16 GiB of RAM.
   - *Option B (low memory):* split Q_5 = Q_4 × K_2, so W = (W_lo, W_hi) with 16 bits each. Then f(W) = Σ_{A indep in Q_4, A ⊆ W_lo} g(W_hi \ A), where g is the 2^16-entry zeta table for Q_4. That is at most 743 terms per evaluation and needs no big table.
2. **Symmetry:** Aut(Q_5) has 2^5·5! = 3840 elements (flip any bits, permute coordinates). The summand is invariant when S and U are moved together by the same symmetry. So:
   - sum S only over **orbit representatives**, weighted by orbit size (a few hundred representatives);
   - keep U over all 254475 sets.

   That comes to about 10^8 evaluations of f, which should take minutes to hours on one PC.
3. The C_4 rotation and reflection symmetries give a further small constant factor if needed.

## Other levers (if the above hits a wall)
- [ ] Q_7 = Q_4 × Q_3. Count homomorphisms from Q_3 into the 743-vertex "disjointness graph" of Q_4 independent sets. The same idea applies, organised differently.
- [x] Count by size: track Σ x^|I| to get A354802 (independent sets by size) for n = 7 as a bonus sequence. Done: row 7 in `b354802.txt` (n = 0..135).
- [x] A354082(7) = 7715 reproduced (check).

## Stretch: a(8)
a(8) ≈ 10^39. The same trick needs Q_6 independent sets, of which there are 2×10^10, so the pairs sum is infeasible without a new idea. Possible directions:
- transfer matrices over a Q_4 × Q_4 split;
- inclusion–exclusion over the two sides of the bipartition, which the Sapozhenko/Galvin container method organises around "small defects".

This is genuinely hard and would need real novelty.

## Verification plan
- Reproduce a(0..6) with the same code path, which is already done in Python for the formula.
- Compute a(7) twice: (a) with the full zeta table, (b) with the low-memory split, and ideally (c) with the Q_4 × Q_3 organisation.
- Sanity check against the Jenssen–Perkins refined asymptotic expansion.
- Check parity or modular identities: for example, an involution argument (complementing all bits maps independent sets to independent sets) gives constraints on a(n) mod 2.

## Log
- 2026-09-25: Triage picked this target. Literature check found nothing exact. The trace(M^4) identity was validated on a(4), a(5), a(6) in Python.

## Results (2026-09-26)
**a(7) = 78685477899897082403** (≈ 7.87×10^19, above 2^64, so accumulated in u128).

Program: `rust/` (`hypercube-indep`). It implements exactly the formula above:
- f(W) goes through Q_5 = Q_4 × K_2 (option B), with the lists precomputed:
  - g = zeta table of Q_4 (2^16 entries);
  - for each of the 2^16 possible W_lo, a CSR list of the Q_4 independent sets inside it (5.0×10^6 entries in total).

  f(W) = Σ_{A in list[W_lo]} g[W_hi & ~A].
- S runs over the orbit representatives of Aut(Q_5) (288 orbits), found by union-find over the group generators and weighted by orbit size. U runs over all 254475 sets.
- Work items (S-rep, U-chunk) are handed out through an atomic counter, one worker thread per core. The program opts out of Windows EcoQoS.

| Check | Result |
|---|---|
| `validate`: d = 1..4, each of the full / flips-only / trivial groups; both Q_d×K_2 → a(d+1) and Q_d×C_4 → a(d+2) | all reproduce a(2..6) |
| `spotcheck 5 4000`, `spotcheck 4 4000`: CSR f vs direct recursive count on random W | 0 mismatches |
| d = 5 with the Q_5×K_2 identity | reproduces a(6) = 19768832143 |
| d = 5, full group Aut(Q_5) (288 reps) | 78685477899897082403, 0.75 s |
| d = 5, bit flips only (8394 reps) | 78685477899897082403, 18.5 s |
| d = 5, no symmetry (all 6.5×10^10 pairs) | 78685477899897082403, 627 s |
| `poly 5`: the same sum by set size (independence polynomial of Q_7) | coefficients sum to a(7); value at −1 is **7715 = A354082(7)** (Flippen–Taylor 2022, deletion–contraction: an independent method); rows 2..6 match A354802 (row 6 against its b-file). 4.9 s |

Sanity checks against the asymptotics: a(n) / (2√e·2^(2^(n−1))) = 0.88, 1.18, 1.40, **1.29** for n = 4..7, heading back toward 1 as expected. The Jenssen–Perkins first-order correction 1 + (3d²−3d−2)/(8·2^d) gives 1.17 (d = 6) and 1.12 (d = 7). The leftover ratios, 1.19 and 1.15, shrink.

Novelty: the value is not in our 2026-09-25 OEIS clone (`git grep` over all of `seq/`). Valency search finds only asymptotic papers.

### GPU?
Not needed for a(7): the whole computation is under a second on the CPU. For a(8), see below: no brute-force formulation is in reach even on a GPU, so the obstacle is the algorithm, not the hardware.

### a(8)? A feasible route: union grouping (2026-09-26)
a(8) ≈ 1.2×10^39 ≈ 2^130. That is **above 2^128**, so the sums use 256-bit accumulators.

The pair formula one level up is a(8) = Σ_{S,U ∈ I(Q_6)} f_6(~(S∪U))². Done naively that is about
8.6×10^15 expensive f_6 evaluations even after using symmetry. Two ideas make it tractable
(`rust/src/ygroup.rs`):

1. **Group the pairs by their union Y = S ∪ U.** The number of pairs with a given union has a closed form,
   c(Y) = 3^{#isolated vertices of Q_d[Y]} · 2^{#other components}.
   An isolated vertex can be in S, U or both. A larger component is connected and bipartite, so it
   has exactly two S/U colourings. Therefore
   **a(d+2) = Σ_{Y ⊆ V(Q_d)} c(Y) · f_d(~Y)²**.
2. **Compute whole slices of f at once.** Write Q_6 = Q_4 × C_4 as four 16-bit fibres, Y = (Y0, Y1, Y2, Y3).
   Fibres 0,1 form the lo copy of Q_5 and fibres 3,2 the hi copy.
   - Y_lo = (Y0, Y1) runs over the orbits of Aut(Q_5): **1,228,158 orbits = A000616(5)**.
   - Y3 and Y2 run over all 2^16 values each.
   - For fixed (Y0, Y1, Y3), the map W2 ↦ f_6(W) is one subset-sum transform over 2^16 values:
     f_6(W) = Σ_{t2 ⊆ W2} G(t2), with G(t2) = Σ_{t0 ⊆ W0} f_4(W1∖(t0∪t2))·f_4(W3∖(t0∪t2)).
     This works out to about 6 ns per term.

Validation:
- `ygroup`: the same generic code with d = 3, 4, 5 reproduces a(5), a(6) and a(7) exactly (0.8 s).
  Its orbit counts (6, 22, 402) match A000616.
- `ybench`: at d = 6, the 16-bit f_6 slice matches f_6 evaluated independently through Q_6 = Q_5 × K_2
  (via the a(7) engine's f_5) on 256 random points.
- `ybench`: the 64-bit fibre-layout neighbour function is correct for every vertex, and c(Y)
  agrees with a brute-force count of the S/U labellings on 300 random Y.

**Engineering (2026-09-26), measured:**

| Version | Time per orbit | Full run |
|---|---|---|
| CPU, first version (35 ns/term) | ~30 s/core | ~200 days on 12 threads |
| CPU, fast c(Y) (component merge via byte tables) + vectorised transform (7 ns/term) | ~100 s at 12 threads in parallel | ~36 days on paper; ~4 months in practice (hyperthreading and E-cores) |
| **GPU, RTX 3080** (`cuda/a8.cu`, CUDA 13.3 in WSL) | 0.15 s | 51 h |
| **GPU + lo↔hi swap symmetry** (`--sym`) | 0.111 s (uniform sample of 300 orbits) | ≈ 38 h |
| + lane-indexed tables moved from constant to shared memory (constant memory serialises divergent reads) | 0.100 s | ≈ 34 h |
| **+ per-high-byte 8-bit transform in registers via warp shuffles**  | 0.094 s | ≈ 32 h |
| **+ Nsight Compute tuning** (see below) | **0.067 s** (uniform sample of 300) | **≈ 23 h** |

Nsight Compute findings (2026-09-26, after enabling the GPU performance counters and restarting WSL), and what was done:

| Finding | Fix | Effect |
|---|---|---|
| Occupancy 33%: 34.5 KB of shared memory, 2 blocks per SM | reuse the per-warp buffer for the 35 filtered sums; < 33 KB | 3 blocks per SM, 50% occupancy |
| Filter loop: 25% of instructions (64-bit mask tricks, strided reads) | G stored transposed, 32 + 3-bit masks, rows 32..34 in the same loop | ~9% of instructions |
| 11% of stall samples: warps idle at block end (swap skipping gives warps unequal work) | high bytes handed out dynamically, most work first | – |
| 192-bit product | the block constant 3^base_iso·2^base_c2 is applied once per block, so the per-term multiplier is usually one limb | ~3% |

Together: 0.094 → 0.067 s per orbit (kernel 71 → 51 ms per orbit launch, IPC 2.44 → 2.75).

Tried and rejected, each still exact:
- **No popcount ordering** (conflict-free reads, but skips are per lane): 0.078 s.
- **Popcount order regrouped for distinct low nibbles**: 0.074 s.
- **XOR-swizzled f storage**: no change.
- **Table-free flood fill**: 0.098 s.
- **512 or 384 threads per block**: no change.

Remaining profile:
- the flood fill's random table lookups: about 17% of instructions and half the bank conflicts;
- the register transform's shuffles: about 19%, close to minimal for a dense 256-point transform of 35-bit values;
- the product: about 8%.

Robustness: one run stalled once, with the host spinning and the GPU idle, right after a profiling session. 15 consecutive reference checks afterwards were clean. The run script therefore processes 4000-orbit chunks under a timeout and retries any chunk that stalls.

Tried and rejected: a precomputed 35×256 table per block (`-DZH`). It removes the per-byte work but needs 70 KB of shared memory, so only one block fits per SM: 0.117 s with 1024 threads, 0.142 s with 512.

Profile by ablation (Nsight Compute cannot read the performance counters until they are enabled in the Windows NVIDIA Control Panel), per orbit without `--sym`, before REGZ:
- per-high-byte transform ≈ 0.07 s;
- c(Y) ≈ 0.035 s;
- 192-bit product ≈ 0.015 s;
- bare loop ≈ 0.01 s;
- G ≈ 0.

Switches: `-DSKIP_G/-DSKIP_C/-DSKIP_MUL/-DTERM_MIN/-DSKIP_PERA` (`cuda/timing.sh`, results wrong); variants in `cuda/variants.sh`.

GPU design: one block per (orbit, Y3) slice of 2^16 terms.
- G(t2) is stored as a 35×35 matrix, since every I_4 set is a pair of I_3 bytes.
- Each warp takes one high byte of Y2 at a time. It builds 35 filtered sums (the loop is the same for every lane) and runs an 8-bit transform in shared memory.
- c(Y) is computed with byte-table flood fills.
- Each term is accumulated exactly as 3^iso·(f² ≪ comps) in 192 bits.

Swap symmetry: the summand is invariant under Y_lo ↔ Y_hi and popcount is Aut-invariant. So a term counts 0/1/2 times as |Y_hi| is <, = or > |Y_lo|, and low bytes are visited in popcount order so that whole warps skip together. The weighted total equals the plain total at d = 3, 4, 5 (`ygroup`).

Orbits: `hypercube-indep reps6` enumerates the **1,228,158 = A000616(5)** orbit representatives of subsets of V(Q_5) (a 2^32-bit visited bitmap, 205 s). The orbit sizes sum to 2^32.

Validation of the GPU code: on 12 orbits, including Y_lo = ∅ and Y_lo = everything, the GPU sums equal the CPU sums from `hypercube-indep cpu6`, both plain (`logs/cpu6_reference.txt`) and swap-weighted (`logs/cpu6_reference_sym.txt`). The run script re-checks this before starting.

**Running a(8) natively on Windows** (PowerShell, no WSL at run time; resumable; about 23 h of GPU time; the desktop will be sluggish while it runs):
```
cd C:\Users\guill\oeis-sniper\attacks\A027624\rust
cargo build --release
.\target\release\hypercube-indep.exe a8run      # rerun the same command to resume
python ..\cuda\aggregate.py                     # checks every orbit is present once, prints a(8)
```
`a8run` first re-checks the embedded kernel against the CPU reference sums. It then processes 4000-orbit chunks, each in a child `a8` process that is restarted if it stalls (watchdog, exit code 3) or times out. The log is `logs/a8_gpu.log`.

How it works: `rust/src/gpu.rs` loads `nvcuda.dll` (the CUDA driver API, part of the NVIDIA driver) at run time and launches the kernel from `cuda/a8.cubin` (sm_86 SASS, embedded with `include_bytes!`). No CUDA toolkit is needed on Windows. After changing `a8.cu`, rebuild the cubin in WSL with `cuda/ptx.sh`.

Checks:
- The native build matches both CPU reference sets.
- It gives 0.067 s/orbit on the uniform sample, the same as WSL.
- Its results for orbits 0–39 are byte-identical to the WSL run's.
- Resume after a truncated line and the stall watchdog behave as in WSL.
- The supervisor retries and gives up as intended.

The result files are interchangeable. The WSL route still works: `wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/guill/oeis-sniper/attacks/A027624/cuda/run_a8.sh --fg`.

### a(8): result (2026-09-27)

> **a(8) = 1268098993536094508894717661843009268823** (≈ 1.268×10^39, 130 bits)

The native Windows run (`a8run`, `--sym`) started 2026-09-26 23:33 and finished 2026-09-27 23:08: **23.6 h wall**, 308 chunks of 4000 orbits, no stalls and no retries (`logs/a8_gpu.log`). `cuda/aggregate.py` checks that each of the 1,228,158 orbits appears exactly once with the mask and orbit size from the representatives file, then sums orbit_size × sum. The per-orbit table `results/a8_sym.txt` (56 MB) is kept locally, outside git.

Checks on the value:

| Check | Result |
|---|---|
| Every orbit present exactly once; orbit sizes sum to 2^32 | yes (`aggregate.py`) |
| GPU vs CPU (`cpu6`), 12 chosen orbits incl. Y = ∅ and Y = V, plain and swap-weighted | equal (`logs/cpu6_reference*.txt`, re-checked by `a8run` before the run) |
| GPU vs CPU, **24 orbits chosen at random after the run** (seed 20260927), swap-weighted | all 24 equal (`logs/cpu6_spotcheck_sym.txt`, indices in `logs/cpu6_spotcheck_idx.txt`) |
| Same generic code at d = 3, 4, 5 | reproduces a(5), a(6), a(7) |
| Novelty | the value is not in the 2026-09-25 OEIS clone |

Ratio to Sapozhenko's 2√e·2^(2^(n−1)): 0.88, 1.18, 1.40, 1.29, **1.13** for n = 4..8. After dividing out the Jenssen–Perkins first-order correction 1 + (3n²−3n−2)/(8·2^n) = 1.081 at n = 8, the leftover ratio is 1.05 (it was 1.15 at n = 7, 1.19 at n = 6), so the value sits where the asymptotics say it should.

**Not done:** the fully independent recomputation, a second full run without the swap symmetry (`hypercube-indep.exe a8run --nosym`, output `results/a8_plain.txt`, roughly 1.5 days). Its per-orbit sums are entirely different numbers, but the total must agree. The OEIS draft flags this.

### Files
- `rust/`: the engine. Build with `cargo build --release`. `.cargo/config.toml` sets target-cpu=native.
  - `./target/release/hypercube-indep validate`
  - `./target/release/hypercube-indep run 5 [--group full|flips|none] [--threads T]`
  - `./target/release/hypercube-indep spotcheck 5 4000`
  - `./target/release/hypercube-indep ygroup` (union grouping, d = 3..5) and `ybench [samples] [y3]` (d = 6 checks + cost estimate)
- `logs/`: run logs for a(7) (`a7_*.log`), for the Q_7 independence polynomial (`poly7_full.log`), the a(8) GPU run (`a8_gpu.log`) and the CPU reference sums for the GPU (`cpu6_*.txt`).
- `cuda/`: the a(8) kernel (`a8.cu`), its compiled `a8.cubin` (embedded by the Rust binary), `aggregate.py`, and the WSL build/profile/run scripts.
- `b027624.txt`: A027624 n = 0..8 (optional b-file). `b354802.txt`: A354802 rows 0..7 (n = 0..135), for the A354802 edit.
- `oeis_submission.md`: the drafted edits for A027624 (a(7), a(8)) and A354802 (row 7).
  - `./target/release/hypercube-indep poly 5`: Q_6 and Q_7 polynomials, checked against A354802 / A354082.
