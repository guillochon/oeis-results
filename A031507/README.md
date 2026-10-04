# A031507 / A031508 — Smallest k > 0 such that y² = x³ + k (resp. y² = x³ − k) has rank n

**Status:** **SOLVED, 2026-09-28.** **A031507(7) = 47550317** and **A031508(7) = 56877643**, proved without GRH or BSD
(so also A373795(7) = 47550317). This page is the full account of that result. The sections from "What is being
counted" through "A literature note" are a self-contained explanation: what is counted, the theorem, the method stage by
stage, every check, and how to reproduce it. The run summary is `logs/a7_run.log`. Below that, under "Working notes",
are the original plan, the literature check and the log, kept as they were written.
**Known terms:** a(0)..a(6) in both sequences (7 terms each). DATA was checked on oeis.org on 2026-09-27.
- A031507: 1, 2, 15, 113, 2089, 66265, 1358556
- A031508: 1, 2, 11, 174, 2351, 28279, 975379

All 14 terms are **unconditionally proved**. a(0..5) are from Gebel–Pethő–Zimmer 1998. a(6) comes from Womack's
exhaustive mwrank search over −10⁶ < k < 1.4×10⁶ (Sept–Nov 2001), described on his web page as "without using
conjectures". The terms were last extended in 2001–2003 (Womack); since then there have been only upper-bound comments (Womack, Aranda,
Hasler, Elkies 2024).
**Target:** **a(7) for both sequences.** Womack's 2002 candidates are A031507(7) = 47550317 and A031508(7) = 56877643. Both curves were
known to have rank exactly 7 (mwrank, unconditional for rank ≤ 8, per Womack). What was missing is minimality: every
smaller k must have rank ≠ 7. Womack (2002): "checking slightly over 100 million curves … would take mwrank roughly
one year running on 500 1GHz computers." The values are Womack's; what is new here is the proof that they are the
smallest. This also settles A373795(7) = min of the two = 47550317.

## What is being counted
For each n, a(n) is the least k ≥ 1 such that the Mordell curve E_k: y² = x³ + k (A031507), or y² = x³ − k (A031508), has
Mordell–Weil rank exactly n. Only sixth-power-free k matter, since k and u⁶k give isomorphic curves.
E_k and E_{−27k} are 3-isogenous, so they have the same rank.
- Tiny example: A031507(2) = 15. For k = 1..14 the ranks of y² = x³ + k are 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0
  (A060950; checked with `ellrank`), and y² = x³ + 15 has rank 2, with independent points (1, 4) and (109, 1138)
  (the height-pairing determinant is 3.76).
- A031508(1) = 2: y² = x³ − 1 has rank 0, and y² = x³ − 2 contains (3, 5), a point of infinite order.

## Result
**Theorem.** Every sixth-power-free k with 1 ≤ k < 47550317 has rank(y² = x³ + k) ≤ 6, and every sixth-power-free
k with 1 ≤ k < 56877643 has rank(y² = x³ − k) ≤ 6. The curves y² = x³ + 47550317 and y² = x³ − 56877643 have rank
exactly 7. Hence A031507(7) = 47550317 and A031508(7) = 56877643, and A373795(7) = min of the two = 47550317.
(A non-sixth-power-free k gives the same curve as k/u⁶, which is smaller and already covered.) No conjecture is used; the only unproved ingredient is PARI's BPSW
primality test inside the factoring. BPSW is proven for n < 2⁶⁴ and has no known counterexample above that.

**How the ~1.14×10⁸ curves were eliminated** (both signs; counts are for k below the candidate):

| Stage | Tool | + curves left | − curves left |
|---|---|---:|---:|
| 1 | 3-isogeny ambient bound B ≤ 6 (Rust, 2.5 s for everything) | 32,546 | 39,962 |
| 2 | Cassels' formula for dim Sel^φ − dim Sel^φ̂ (Tamagawa numbers, periods) | 12,355 | 18,943 |
| 4 | explicit Selmer non-members (locally insoluble coverings at primes ≤ 1000) | 121 | 143 |
| 5 | the same, over the whole span of the generators and all primes of S | 120 | 143 |
| 6 | 2-descent with class groups proved by `bnfcertify`: R ≤ 5 for all | **0** | **0** |

1. **Stage 1, the ambient bound.** For E_v: y² = x³ + v, the 3-isogeny φ (kernel (0, ±√v)) and its dual have Selmer
   groups in the minus parts of K(S,3), K = Q(√−3v) and Q(√v), S = {p | 6k}. By Kummer theory and the S-unit and
   class-group sequences, each has F₃-dimension at most s(K) = h₃(K) + [K real or Q(√−3)] + #{p ∈ S split in K}.
   The isogeny formula 3^r = |Im α||Im α̂| / (torsion terms) gives rank ≤ s + s' − δ (δ = 1 iff v or −3v is a square).
   - **Unconditional h₃.** By Hasse, a fundamental D has exactly (3^h₃ − 1)/2 cubic fields of discriminant D. PARI's
     `nflist("S3", …)` (Belabas' enumeration, unconditional) plus `nfdisc` gave every h₃ for |D| ≤ 6.83×10⁸
     (`tools/h3chunk.gp`, `tools/run_h3.sh`; about 20 min on 10 processes). Checks: every count had the form
     (3^h − 1)/2, and the table agrees with `quadclassunit` for all 1,215,866 fundamental |D| ≤ 2×10⁶.
   - This replaces 23 CPU-hours of GRH-conditional `quadclassunit` calls with a table lookup: the whole stage takes
     2.5 s in Rust (`rust/`).
2. **Stage 2, Cassels.** |Sel^φ(E)|/|Sel^φ̂(E')| = |E(Q)[φ]| Ω(E') Πc_p(E') / (|E'(Q)[φ̂]| Ω(E) Πc_p(E)). This fixes
   g = dim Sel^φ − dim Sel^φ̂ exactly, so rank ≤ min(2s_a − g, 2s_b + g) − δ (`tools/stage2.gp`). In every case the
   ratio came out an exact power of 3.
3. **Stages 4–5, witnesses.** Since g is known, one element of the ambient group proved *not* to be in the Selmer group
   lowers the total bound by 2. Candidates are units, generators of split primes, and a³ for 3-torsion classes.
   `bnfinit` (GRH) is only used to *find* them; membership is re-checked with exact ideal arithmetic. Each c = c₁ + c₂√u
   gives the covering Z³ = c₂X³ + 3c₁X²Y + 3uc₂XY² + uc₁Y³, and insolubility over Q_p is decided exactly by a p-adic tree
   search (`tools/witness.gp`). The fast version only descends into roots of the reduction. It returns "soluble" by
   the Weil bound for p > 50, but that answer can never produce a false witness.
4. **Stage 6, certified 2-descent** for the 264 curves the 3-isogeny method could not close (the 121 + 143 stage-4
   survivors; stage 5 closed one of the + curves, but stage 6 was run on all 121 anyway). Their φ-Selmer groups
   apparently fill the ambient bound, which points to nontrivial Sha[3]. `ellrankinit(E)` holds the `bnf` of
   Q[x]/(x³ + v) used by the 2-descent, and `bnfcertify` proves that `bnf` without GRH. `ellrank` on the same
   structure then gives R ∈ {1, 3, 5} for all 264 curves. This was the slow part: 19 + 224 CPU-minutes, one curve
   alone took an hour.
5. **The rank-7 curves.** `ellrank` finds 7 points on each curve. The height-pairing matrix is positive definite
   (regulators 181126.71 and 155212.41, smallest eigenvalue > 1.2, the same at 38 and 115 digits). The unconditional
   bound of stage 2 is exactly 7 for both (s_a = 4, s_b = 3, g = 1): so rank = 7. Points in `logs/a7_run.log`.

**Checks** (logs in `logs/`):
- Stage 2: over all 9,832 curves with k ≤ 5000, the bound is never below `ellrank`'s proven lower bound r. With the
  Cassels ratio inverted the same test fails 3,099 times, so the test has teeth (`validate_stage2_5000.log`).
- Witnesses: the images y + √u of 427 actual rational points are never declared insoluble. On all 1,968 curves with
  k ≤ 1000 the improved bound is never below r (`test_witness_1000.log`).
- Stage 1 vs an independent GP implementation of the Cohen–Pazuki bound (k ≤ 10⁵): ours is ≥ theirs everywhere and
  equal on 2/3 of the curves. The difference is that we keep 3 in S.
- **End-to-end:** the same pipeline reproduces the known, unconditionally proved a(6) of both sequences. It shows rank ≤ 5
  for every k below 1358556 (+) and 975379 (−), in 7 s and 5 s, and confirms the a(6) curves have rank 6 with a
  certified 2-descent (`a6_reproduction.log`). This is independent of Womack's 2001 mwrank search.

**Reproduce** (PARI/GP 2.15 and Rust; about 5 CPU-hours; the intermediate files in `data/` take 1.4 GB):
```bash
bash tools/run_h3.sh 450000000 10000000 10               # unconditional h3 table, |D| <= 4.5e8
SIZEMAX=10000000000 bash tools/run_h3.sh 683000000 2500000 10 450000001   # the rest (bigger PARI stack)
cd rust && cargo build --release && ./target/release/a031507 table 683000000
./target/release/a031507 filter 56877643 7 10             # stage 1, both signs
cd .. && bash tools/pipeline.sh a7plus 1 data/cand_plus_7_56877643.txt 47550317 6
bash tools/pipeline.sh a7minus -1 data/cand_minus_7_56877643.txt 56877643 6
gp -q tools/rank7_certificate.gp                          # the two rank-7 curves
```
`tools/pipeline.sh` runs stages 2–6 for one sign. With target 5 and the a(6) values it reproduces a(6).

**A literature note.** arXiv:math/0403116 is **Elkies–Rogers**, "Elliptic curves x³ + y³ = k of high
rank" (ANTS 2004), not Elkies–Watkins. Its statements that ranks 4 and 5 are proved minimal and ranks 6 and 7 are minimal under
weak BSD + GRH are about the *cube-sum* family E_k: x³ + y³ = k, which is y² = x³ − 432k². They are not about
A031507/A031508. For our sequences rank 6 is already unconditional. The new content is **the first a(7), proved
unconditionally, in both sequences.**

---

# Working notes
Everything below was written while planning and running the computation (from 2026-09-27). The account above
supersedes it where they differ: for example, the plan's `ellrank` stage became stages 2–6, and the h₃ table was built
with PARI's `nflist` rather than Belabas's CUBIC.

## How the known terms were computed
- GPZ 1998: ranks and integral points for |k| ≤ 10⁴ (later up to 10⁵), by 2-descent. This gives a(0..5).
- Womack 2001: mwrank (2-descent) on every k in (−10⁶, 1.4×10⁶); eight weeks on three Athlons. This proves a(6). The
  larger values a(7..12) were *found* by sieving for curves with many small integral points, or by Quer's 3-rank
  construction. None of them was proved minimal.
- The OEIS program (`ellrank` loop) is fine for small n but runs a full 2-descent on every k. At about 20 ms per curve that is
  about 1e8 × 20 ms ≈ 23 CPU-days for a(7). That is possible, but it wastes time on curves that obviously have small rank. Also,
  about 1% of curves have nontrivial Sha[2], so there `ellrank` returns r < R and 2-descent alone cannot decide.

## Clever ideas to pursue

1. **Cheap unconditional filter by descent via 3-isogeny, then 2-descent on the few survivors** (main plan).
   - *Formula* (Cohen–Pazuki 2009, arXiv:0903.4963, with a = 0). Write k = D·b², where D is a fundamental discriminant
     (D = 1 allowed), K = Q(√D), and 2b ∈ Z. Prop. 2.2 gives |Im α|·|Im α̂| = 3^(r+δ), where δ = 1 iff D ∈ {1, −3}.
     Cor. 4.3 and Lemma 4.4 put Im α inside a group of F₃-dimension at most
     **s(k) = h₃(K) + [D > 0 or D = −3] + #{p | 2b : p splits in K}**. (For D = 1, Thm 3.1(3) gives ω(2b).) The dual curve is
     E_{−27k}, with K̂ = Q(√−3D). Hence
     **rank E_k ≤ B3(k) := s(k) + s(−27k) − δ.**
     Here h₃ is the 3-rank of the class group, and the unit term is the dimension of U(K)/U(K)³.
   - This ignores all local conditions, so it bounds the 3-isogeny Selmer bound from above. Consistent with this
     are the Selmer-in-class-group bounds of Jha–Majumdar–Shingavekar (arXiv:2207.12487, Thm 3.18/3.22) and the rank bound
     "#rows + #cols − 2·rank(A) − 1" of Elkies–Rogers §3 for the cube-sum subfamily. For squarefree k, B3 is roughly
     h₃(Q(√k)) + h₃(Q(√−3k)) + 1. So rank 7 needs 3-ranks of about 3 in *both* fields (as in Quer's construction), or many split primes in b.
   - *Pipeline:*
     1. For every sixth-power-free k < 47550317 (+) and < 56877643 (−), compute B3(k), and discard k if B3(k) ≤ 6.
     2. For the survivors, run `ellrank` (2-descent; its upper bound R is unconditional), and discard k if R ≤ 6.
     3. Any k still left gets an exact 3-isogeny Selmer computation (the Cohen–Pazuki local tests, §5–6), `ellrank` with more
        effort, and a point search. The only k that cannot be settled are those with both Sha[2] and Sha[3] large;
        expect none or a handful.
   - *Calibration (2026-09-27, see the log):* 0.8 ms per k. About 2×10⁻⁴ of the k near 4.75×10⁷ have B3 ≥ 7, and each of
     those took `ellrank` 10–35 ms, all with R ≤ 3 in the sample.
   - **Estimated total: about 1.03×10⁸ curves × 0.8 ms ≈ 23 CPU-hours ≈ 2–3 h wall on 12–16 threads.** Stage 2 adds
     at most about 2×10⁴ × 30 ms ≈ 10 min. Compare Womack's estimate of 500 machine-years for mwrank.
   - *The GRH caveat:* `quadclassunit` assumes GRH, and a class group computed too small would make B3 too small. The
     unconditional replacement is Hasse's theorem: a fundamental D has exactly (3^h₃ − 1)/2 cubic fields of
     discriminant D. Belabas's CUBIC program (Davenport–Heilbronn reduction of binary cubic forms, unconditional,
     https://www.math.u-bordeaux.fr/~kbelabas/research/cubic.html, v1.4, builds against libpari) counts all cubic fields
     with |disc| < 10⁷ in about 0.6 s. The range needed here is |D| ≤ 12·5.7×10⁷ ≈ 7×10⁸, which is about 1 minute of counting
     or about 10 minutes with output. Build a bitmap/table of h₃(D) for all fundamental |D| ≤ 7×10⁸, both signs, and B3 becomes
     a lookup plus trial factoring of k. That is *faster* than the GRH version, and fully unconditional. (PARI's own
     `nflist("S3", [a,b])` agrees with `quadclassunit` for |D| ≤ 2×10⁵, but it crashed with a thread-stack error at 10⁶, and it is much slower.)
2. **Rank-7 curves themselves.** Re-verify with `ellrank`: 7 independent points (nonzero regulator from
   `ellheightmatrix`, checked at two precisions) and upper bound R = 7. B3 already equals 7 for both curves,
   so the 3-isogeny bound alone gives rank ≤ 7 once h₃ is certified. Cross-check with mwrank (eclib) if it can be
   installed (fetch on Windows; WSL has no DNS).
3. **Stretch goals.**
   - a(8): the bounds are 1632201497 (+) and 2520963512 (−). That is about 30× more curves, so about 3–4 CPU-weeks with a cached
     h₃ table to 3×10¹⁰; CUBIC to 3×10¹⁰ is about 1 hour. Feasible, but only after a(7) is done.
   - Side target: make the Elkies–Rogers conditional minimality for x³ + y³ = k (ranks 6 and 7: k = 9902523 and
     1144421889) unconditional with the same filter, since E_k is y² = x³ − 432k². Look up the OEIS entry for that family first.

Checklist:
- [x] Pruning via necessary conditions: B3 is a necessary condition for rank ≥ 7
- [x] Related sequences: A060950/A060951 (ranks), A373795 (min over signs), the cube-sum family
- [x] Symmetry: only sixth-power-free k; k ↔ −27k have the same rank (not otherwise used)
- [ ] Parallelize: `parfor` over blocks of k in gp, or a Rust driver using the h₃ table
- [ ] n/a: transfer matrix, generating function, SAT

## Literature check (2026-09-27: Valency corpus, arXiv full text, oeis.org, Womack's page)
**Verdict: nobody has proved a(7) for either sequence; the OEIS lists only a(0..6).**

| Source | What it has | Relevance |
|---|---|---|
| oeis.org/A031507, A031508, A373795 (live, 2026-09-27) | DATA stops at rank 6; a(7) ≤ 47550317 / 56877643 only as comments; no conditional flags | target confirmed |
| Womack, web page (2002, web.archive 2017 snapshot) | green = proved minimal via exhaustive mwrank for −10⁶ < k < 1.4×10⁶ (ranks ≤ 6); rank ≤ 8 curves' ranks proved by mwrank; rank-7 minimality would need about 10⁸ curves, "one year on 500 computers" | source of a(6) and of the a(7) candidates |
| Womack, PhD thesis (Nottingham 2003), Table 3.3 | same tables | could not download (Warwick host refused); not read |
| Gebel–Pethő–Zimmer, Compositio 110 (1998) | ranks for small |k| | a(0..5) |
| Elkies–Rogers, arXiv:math/0403116 (ANTS 2004) | x³ + y³ = k family: ranks ≤ 5 minimal unconditionally, 6 and 7 minimal under weak BSD + GRH; 3-isogeny rank bound via a cubic-residue matrix | **different family** (y² = x³ − 432k²); method is the same idea as ours |
| Cohen–Pazuki, arXiv:0903.4963 (Acta Arith. 140, 2009) | complete 3-isogeny descent for y² = x³ + D(ax + b)², including all local solubility tests | **the formula we use** (Prop. 2.2, Cor. 4.3, Lemma 4.4, §6) |
| Jha–Majumdar–Shingavekar, arXiv:2207.12487 | upper and lower bounds on φ-Selmer of y² = x³ + a in terms of S-class groups (Thm 3.14, 3.18, 3.22) | independent second formula, for cross-checking B3 |
| Chan, arXiv:2211.06062 (IMRN 2024); Shingavekar arXiv:2406.03066; Bhargava–Klagsbrun–Lemke Oliver–Shnidman arXiv:1709.09790 | statistics of 3-isogeny Selmer groups | explains why the filter kills almost everything |
| Bennett–Ghadermarzi, arXiv:1311.7077 | all integral points for |k| ≤ 10⁷ via Thue equations | no ranks, not relevant to minimality |
| Elkies, ANTS-XVI (2024) | a(16) upper bound, rank ≥ 17 exists | upper bounds only |

## Verification plan
- Reproduce a(0..6) of both sequences with the same pipeline (threshold n instead of 7). For k < a(n) show rank ≠ n.
  B3 ≤ n − 1 settles most k; otherwise use the exact rank from `ellrank` (r = R).
- Test the B3 formula: B3 ≥ the `ellrank` lower bound for every k (already true for 0 < |k| ≤ 3000). Compare with the
  JMS Thm 3.18 bound on a random sample. On a sample of about 10⁵ k, check B3 ≥ R_exact whenever R_exact = r.
- Unconditional h₃: CUBIC counts compared with `quadclassunit` on all |D| ≤ 10⁷ (they must agree unless GRH fails).
  Use the CUBIC table in the production run.
- Sixth-power-free normalisation of −27k (when 27 | k, reduce by 3⁶ before computing s(−27k)).
- For every survivor, log k, B3, and the `ellrank` output. The final claim rests on a file listing every k with B3 ≥ 7 and how it was
  closed.
- Independent second run of stage 2 with mwrank (if installable) or with a different `ellrank` effort/seed.

## First milestone
Build CUBIC in WSL (fetch the tarball on Windows) and produce the h₃ table for |D| ≤ 7×10⁸. Then reproduce
**A031507(6) = 1358556 and A031508(6) = 975379** unconditionally (about 2.4×10⁶ curves, a few minutes). After that, run the
full a(7) range in blocks of 10⁶ k with a checkpoint file.

## Log
- 2026-09-27: attack file started. Literature check (see table). An earlier note's citation of "Elkies–Watkins" is actually
  Elkies–Rogers on x³ + y³ = k. For A031507/8 rank 6 is already unconditional (Womack 2001), so the target is a(7).
- 2026-09-27: `selmer3_bound.gp` (B3 from Cohen–Pazuki, h₃ via `quadclassunit`), `calibrate.gp`, `candidates_test.gp`:
  - Record curves, ranks 0..8, both signs: **B3 = rank exactly for all 18** (e.g. B3(47550317) = 7, B3(−56877643) = 7,
    B3(1632201497) = 8). So the bound is tight on the extremal curves.
  - 0 < |k| ≤ 3000, sixth-power free: B3 ≥ the `ellrank` lower bound for every k (0 failures). The slack B3 − r had
    histogram 0:1878, 1:2521, 2:990, 3:439, 4:66, 5:6.
  - 9831 k in [47540000, 47550000], each sign: 7.7–8.5 s (0.8 ms per k). B3 histogram (+): 1:2930, 2:3663, 3:2191,
    4:829, 5:197, 6:20, 7:1. The (−) histogram is almost identical.
  - k ∈ [4.750×10⁷, 4.755×10⁷] with B3 ≥ 7: 9 (+) and 7 (−). `ellrank` settles each in 10–35 ms, with R ∈ {0, 1, 2, 3}.
  - `h3_cubic_test.gp`: cubic-field counts from `nflist` give the same 3-rank as `quadclassunit` for all fundamental
    |D| ≤ 2×10⁵. `nflist("S3", [1, 10^6])` crashed in a worker (thread stack), so use Belabas's CUBIC for the real table.
- 2026-09-27/28: **solved.** The unconditional h₃ table was built from cubic-field counts, to |D| ≤ 6.83×10⁸; the
  4 GB PARI stack cap had to be raised above 4.5×10⁸. The Rust stage-1 filter replaced the planned 23 CPU-hours with
  2.5 s. With the Cassels refinement, the witnesses and the certified 2-descent (for 264 curves), every smaller k is
  shown to have rank ≤ 6. The rank-7 curves have 7 independent points and the unconditional bound 7. The pipeline
  reproduces a(6) of both sequences. Totals: about 1 CPU-hour for the table and stages 1–5, and 4 CPU-hours for
  stage 6.
