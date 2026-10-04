# A109459 — Inequivalent Krom (2SAT) functions under permuting and complementing variables

This folder computes new terms of three OEIS sequences that count Krom (2-SAT) functions of n
variables, all from one run of one program:

| sequence | counts Krom functions | known terms | new terms |
|---|---|---|---|
| **A109457** | labelled (no identification) | a(0)..a(8) | a(9) = 1517417860795700, a(10) = 2026467721823749888 |
| **A109458** | up to permuting the variables | a(0)..a(8) | a(9) = 4787612373, a(10) = 612671534302 |
| **A109459** | up to permuting and complementing the variables | a(0)..a(7) | a(8) = 315080, a(9) = 11073674, a(10) = 667787954 |

The known terms are Knuth's (2005). The sections below explain everything: what is counted, the
structure theorem that reduces the count to "skew posets", how the new terms were computed, every
check, the asymptotics, and how to reproduce. Working notes and the dated log come after.

## What is being counted
A *Krom function* is a Boolean function that can be written as an AND of clauses with at most
2 literals each, a 2-CNF like (x∨¬y)∧(y∨z). Equivalently its set of true points S ⊆ {0,1}^n is
**closed under bitwise majority**. The empty set counts too.
- **A109457** counts all Krom functions of n variables.
- **A109458** counts them up to renaming the variables (the symmetric group S_n, order n!).
- **A109459** counts them up to renaming variables and/or flipping some variables. That is the
  hyperoctahedral group B_n, of order 2^n·n!.

## The structure theorem (proved in the log below and in Lean, and validated by 26 known terms)
Take a satisfiable 2-CNF and close it under all implied 2-clauses (including unit clauses, and
weakening of units). The closed form splits into three layers:
1. **Forced variables** (k of them): value fixed to 0 or 1.
2. **Equivalent variables:** among the rest, literals forced equal form **blocks** (a block is a
   set of literals, e.g. {x, ¬y, z}, up to global negation). Keep one literal per block.
3. **Skew poset:** the remaining implications among the r block-literals and their negations form
   a strict partial order on 2r literals with  l < m  ⟺  ¬m < ¬l  and no literal comparable to its
   own negation.

Conversely every such triple gives a distinct Krom function (proof sketch in the log; Lean proof
described below). Hence
```
A109457(n) = 1 + Σ_{k=0..n} C(n,k) 2^k L(n−k),   L(m) = Σ_{r=0..m} S2(m,r) 2^(m−r) T(r)
A109459(n) = 1 + Σ_{m=0..n} U(m)
A109458(n) = 1 + Σ_{m=0..n} (n−m+1) U'(m)
```
where
- **T(r)** = number of skew posets on r labelled pairs of literals;
- **U(m)** = number of skew posets whose points (pairs) are weighted by block sizes ≥ 1 summing
  to m, counted up to isomorphisms commuting with negation;
- **U'(m)** = the same, but each block-literal is coloured by (number of positive, number of
  negative literals in its block), its negation by the reversed pair, total weight m, up to
  colour-preserving isomorphism.

Where each factor comes from. The leading 1 is the empty function. In the labelled count,
C(n,k) 2^k chooses the k forced variables and their values; S2(m,r) partitions the other m
variables into r blocks, and 2^(m−r) chooses the signs of the literals inside each block relative
to one of them (a block of size s has 2^(s−1) sign patterns). Up to permuting and complementing
(A109459), the forced variables are all alike, so only their number n − m matters, and a block
is described by its size alone. Up to permuting only (A109458), a forced variable still carries
its value, so the forced layer contributes the number of ways to choose how many of the n − m
forced variables are forced to 1, which is n − m + 1; and a block now carries its numbers of
positive and negative literals, which is why U' uses signed colours.

In the language of Bollobás, Brightwell and Leader (2003) the forced variables are the "spine",
the blocks are the "associated pairs", and a function with neither is "elementary": the
elementary functions are exactly the skew posets (their posets on the 2r literals), and their
H(r) is our T(r).

## How the counts were computed
`c/skew.c` (C + nauty 2.8.8, OpenMP) enumerates skew posets up to isomorphism, level by level in
the number r of pairs. Every skew poset on r pairs has a maximal literal y; removing its pair
leaves a skew poset on r − 1 pairs, and y is recovered by choosing B = {l : ¬y < l}, which can
be any up-set of the parent with no complementary pair (proof in the log). Children are
canonicalised with nauty and deduplicated. For every class P it records |Aut(P)|, a subgroup of
B_r; T(r) = Σ_P |B_r|/|Aut(P)|.

The unlabelled counts are Pólya sums over Aut(P) ⊆ B_r: for a signed cycle of length k the
weight is x^k/(1−x^k) (block sizes), or (2t−t²)/(1−t)² at t = x^k for an even-sign cycle and
t²/(1−t²) for an odd-sign cycle (signed colours). Only **connected** skew posets are Pólya-counted;
the rest follow by an Euler transform (multisets of connected components), so huge groups such as
the antichain's B_9 (order 1.86·10^8) are never enumerated. `sequences.py` then applies the
three formulas above.

Level r = 9 takes 20 s on 16 threads (64.5 M nauty calls). Level 10 (573012592 classes,
4.4·10^9 candidates) runs without storage via canonical augmentation (`./skew 10 a`), 20 min 40 s
wall (318 CPU-min); the largest connected automorphism group is S_10 (the "at most one true" structure).

## Results (c/skew.c + sequences.py, 2026-09-26)
| r | skew posets S(r) | connected | labelled T(r) |
|---|---|---|---|
| 0 | 1 | – | 1 |
| 1 | 1 | 1 | 1 |
| 2 | 2 | 1 | 5 |
| 3 | 5 | 3 | 69 |
| 4 | 18 | 12 | 2153 |
| 5 | 90 | 69 | 138057 |
| 6 | 736 | 627 | 17221933 |
| 7 | 9879 | 9035 | 4036310797 |
| 8 | 227348 | 216538 | 1736041818705 |
| 9 | 8839096 | 8599886 | 1348201535646097 |
| 10 | 573012592 | 563918911 | 1869742968781218773 |

U(m), m = 0..10: 1, 1, 3, 8, 31, 151, 1164, 14271, 299449, 10758594, 656714280.
U'(m), m = 0..10: 1, 1, 6, 28, 224, 2460, 43742, 1254942, 59502105, 4664654461, 603158463959.

None of S(r), T(r), U(m), U'(m) is in the OEIS (searched 2026-09-26); they are candidates for
new entries ("skew posets" = closed implication structures of 2-CNFs with no forced or
equivalent variable).

## Verification
- **26 known terms reproduced exactly:** A109457(0..8), A109458(0..8), A109459(0..7). The
  labelled check T(8) ≈ 1.7·10^12 involves every one of the 227348 classes with the right |Aut|.
- **Independent second method for the unlabelled counts** (`c/skewcol.c`): generate the
  *coloured* skew posets directly with nauty vertex colours and count canonical forms — no Pólya
  sums, no Euler transform, no automorphism groups. Agrees with U(m) for all m ≤ 9 (so with
  A109459(8) and A109459(9)) and with U'(m) for all m ≤ 8.
- Internal checks in skew.c: the group closure from nauty's generators has exactly nauty's group
  order; every Pólya sum is divisible by |Aut|; |Aut| divides |B_r|; the exponential formula
  applied to the connected labelled counts reproduces T(r).
- Canonical-augmentation mode (no hash table) reproduces the hash-table results for r ≤ 9
  exactly (S, connected, T, all Pólya sums).
- **Caveat for n = 10:** the n = 10 terms have no second independent computation (the coloured
  check would need ~37 GB); they rest on the internal checks above and on the augmentation mode
  matching the stored mode at r ≤ 9. So the second method covers A109459(8), A109459(9) and
  A109458 up to n = 8; A109458(9) needs U'(9) and A109457 is labelled, so those terms (like all
  n = 10 terms) rest on skew.c with its internal checks and the 26 reproduced known terms.

## Asymptotics
Details, proofs and the exact remaining gap are in `asymptotics.md`. In short:
- **Proved:** A109457(n) = (1 + o(1)) · 2^(n(n+1)/2). Conjectured by Bollobás, Brightwell and
  Leader (2003), proved by Allen (2007), second proof by Ilinca and Kahn (2012). Equivalently
  T(n) = (1 + o(1)) · 2^(n(n+1)/2). The convergence is slow: A109457(n)/2^(n(n+1)/2) is still
  56.2 at n = 10.
- **Almost all Krom functions are unate** (monotone after complementing some variables; Allen),
  but at n = 10 only 1.8 % of them are: the ratio of all to unate Krom functions is 55.70 at
  n = 10 and still growing.
- **Conjecture (supported by n ≤ 10):** almost all Krom functions have trivial stabiliser, so
  ```
  A109459(n) ~ 2^(n(n−1)/2)/n!,      A109458(n) ~ 2^(n(n+1)/2)/n!.
  ```
  By Burnside's lemma for B_n, 2^n n! · A109459(n) = A109457(n) + Σ_{g≠1} Fix(g), where Fix(g)
  is the number of Krom functions fixed by g. The second term is 0.951, 0.577, 0.356, 0.225
  times A109457(n) for n = 7..10 (shrinking by about 0.63 per step). For A109458,
  n! · A109458(n)/A109457(n) = 1.32, 1.22, 1.14, 1.10 for n = 7..10.

## How to reproduce
In WSL (or any Linux) with nauty 2.8.8 built in `$NAUTY`:
```
cd c
gcc -O2 -fopenmp -I$NAUTY skew.c    $NAUTY/nauty.a -lm -o skew
gcc -O2 -fopenmp -I$NAUTY skewcol.c $NAUTY/nauty.a -lm -o skewcol
./skew 9    > ../out9.txt        # r <= 9, hash-table mode, 20 s
./skew 10 a > ../out10a.txt      # r <= 10, last level by canonical augmentation, ~20 min
./skewcol u 9 > ../col_u9.txt    # independent check of U(m), m <= 9
./skewcol s 8 > ../col_s8.txt    # independent check of U'(m), m <= 8
cd ..
python sequences.py out10a.txt   # prints S, T, U, U' and the three sequences, checks known terms
```
`sequences.py` prints "ALL KNOWN TERMS REPRODUCED" and the new terms listed above.

## References
- P. Allen, Almost every 2-SAT function is unate, Israel J. Math. 161 (2007), 311-346.
- B. Bollobás, G. Brightwell and I. Leader, The number of 2-SAT functions, Israel J. Math. 133
  (2003), 45-60.
- L. Ilinca and J. Kahn, On the number of 2-SAT functions, Combin. Probab. Comput. 21 (2012),
  621-636; arXiv:1005.2863.
- D. E. Knuth, the terms of A109457-A109459 for small n (2005); his program `krom-count.ch` is a
  change file for `horn-count.w` on his programs page.

## Lean formalization (2026-09-27)
`lean/OeisLean/A109459.lean` (Lean 4 + Mathlib v4.33.1, standard axioms only) proves the
structure theorem behind the formulas, for an arbitrary finite type of variables:
- `sol_valid_eq` — Krom's theorem: a nonempty majority-closed set equals the solution set of the
  implications (2-clauses) valid on it. Proof: the "agree on any coordinate set" induction.
- `valid_sol_iff`, `exists_sat` — for a skew preorder R (reflexive, transitive, reversed by
  negation, no l ≤ ¬l) the implications valid on Sol R are exactly R, and no literal is forced.
  Proof: the extension lemma `extend` (closed consistent literal sets extend to solutions).
- `kromEquivSkew` — forced-free Krom functions ≃ skew preorders on the literals.
- `fiberEquiv`, `card_krom` — the forced-variable layer: |Krom functions on V| =
  1 + Σ over patterns f : V → Option Bool of |skew preorders on {v | f v = none}|. Grouping the
  patterns by the number k of forced variables gives the C(n,k)·2^k of the A109457 formula.
Not formalized: the Pólya/Euler orbit counting and the enumeration itself.

## Files
- `c/skew.c` — enumerator (C + nauty 2.8.8, OpenMP). Build in WSL, see header.
- `c/skewcol.c` — independent coloured-generation check.
- `sequences.py` — turns `out*.txt` into the three sequences and checks them against the OEIS.
- `out8.txt`, `out9.txt`, `out9a.txt`, `out10a.txt` — raw outputs (r ≤ 8, r ≤ 9, r ≤ 9 augmented, r ≤ 10 augmented).
- `col_u9.txt`, `col_s8.txt` — outputs of the independent check.
- `b109457.txt`, `b109458.txt`, `b109459.txt` — tables n = 0..10 in b-file format (not submitted:
  with 11 terms each, DATA already holds the whole table).
- `oeis_submission.md` — draft edits (one document for all three sequences).
- `layman.md` — plain-English write-up.
- `asymptotics.md` — literature (BBL, Allen, Ilinca–Kahn), the asymmetry conjecture, partial proofs.
- `../../lean/OeisLean/A109459.lean` — Lean proof of the structure theorem.

## Working notes
**Status:** SOLVED for n = 8, 9, 10 (2026-09-26). n = 11 would need ~4·10^10 classes (out of reach without a new idea).
All three sequences are tagged `hard,more`.

Nothing beyond the OEIS terms exists in the literature (checked 2026-09-25/26: only asymptotic
2-SAT counts, Allen; Ilinca–Kahn; Dong–Mani–Zhao; the values are not on the web). Knuth's own
program (`krom-count.ch`, a change file for `horn-count.w` on his programs page) is a direct
backtracking enumeration of median-closed sets; that is how n <= 7 was done in 2005.

Literature, found 2026-09-27: Bollobás–Brightwell–Leader (Israel J. Math. 133, 2003) define
the same decomposition (spine = forced variables, associated pairs = blocks, "elementary"
functions = skew posets, their H(n) = our T(n)), and Allen (2007) / Ilinca–Kahn (2012) prove
A109457(n) ~ 2^(n(n+1)/2). Nobody computed terms. See `asymptotics.md` for the asymptotic
picture, the conjecture A109459(n) ~ 2^(n(n-1)/2)/n!, what we could prove towards it and the
exact gap that remains.

## Log
- 2026-09-25: Definition confirmed by brute force (n ≤ 4). Structure decomposition sketched, and
  a(1), a(2) checked by hand.
- 2026-09-26: Proof details for the structure theorem. Given a skew poset P, an assignment is a
  choice of one literal per pair that is an up-set. For l ≤ m not in P the up-set generated by
  {l, ¬m} contains no complementary pair (each of the three ways it could would force l ≤ ¬l,
  l ≤ m or ¬m ≤ m), and any consistent up-set extends to a full assignment (add any unassigned
  literal's up-set; it cannot meet the negation of the current set, because up-sets are
  closed). So the implication relation of S(P) is exactly P, no variable is forced, and distinct
  P give distinct functions. Every Krom function with no forced variable arises: Schaefer/Krom
  say a majority-closed S is the solution set of the 2-clauses it satisfies.
- 2026-09-26: Generation rule. Every skew poset on r pairs has a maximal literal y; removing its
  pair leaves a skew poset, and y is recovered by choosing B = {l : ¬y < l}, any up-set of the
  parent with no complementary pair. Conversely every such B gives a valid child (y is not forced:
  ↑y = {y}; ¬y is not forced: B ∪ {¬y} is consistent; transitivity ¬b < y < nothing is fine, and
  c < ¬b implies b < ¬c so ¬c ∈ B).
- 2026-09-26: skew.c written and run: r ≤ 9 in 20 s. All 26 known terms match. New:
  A109459(8) = 315080, A109459(9) = 11073674, A109458(9) = 4787612373,
  A109457(9) = 1517417860795700. Coloured-generation check agrees (U to m = 9, U' to m = 8).
  Augmentation mode validated on r = 9; r = 10 done in 20 min: A109459(10) = 667787954,
  A109458(10) = 612671534302, A109457(10) = 2026467721823749888.
