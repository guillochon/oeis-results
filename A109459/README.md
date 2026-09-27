# A109459 — Inequivalent Krom (2SAT) functions under permuting and complementing variables

**Status:** SOLVED for n = 8, 9, 10 (2026-09-26). n = 11 would need ~4·10^10 classes (out of reach without a new idea).
**Known terms:** a(0)..a(7) = 2, 3, 6, 14, 45, 196, 1360, 15631 (Knuth, 2005).
**New terms:** a(8) = 315080, a(9) = 11073674, a(10) = 667787954.
**Bonus (same run):**
- **A109457** (labelled count): a(9) = 1517417860795700, a(10) = 2026467721823749888.
- **A109458** (up to permutations only): a(9) = 4787612373, a(10) = 612671534302.

All three are tagged `hard,more`. Nothing beyond the OEIS terms exists in the literature
(checked 2026-09-25/26: only asymptotic 2-SAT counts, Allen; Ilinca–Kahn; Dong–Mani–Zhao; the
values are not on the web). Knuth's own program (`krom-count.ch`, a change file for
`horn-count.w` on his programs page) is a direct backtracking enumeration of median-closed sets;
that is how n <= 7 was done in 2005.

**Literature, found 2026-09-27:** Bollobás–Brightwell–Leader (Israel J. Math. 133, 2003) define
the same decomposition (spine = forced variables, associated pairs = blocks, "elementary"
functions = skew posets, their H(n) = our T(n)), and Allen (2007) / Ilinca–Kahn (2012) prove
A109457(n) ~ 2^(n(n+1)/2). Nobody computed terms. See `asymptotics.md` for the asymptotic
picture, the conjecture A109459(n) ~ 2^(n(n-1)/2)/n!, what we could prove towards it and the
exact gap that remains.

## What is being counted
A *Krom function* is a Boolean function that can be written as an AND of clauses with at most
2 literals each, a 2-CNF like (x∨¬y)∧(y∨z). Equivalently its set of true points S ⊆ {0,1}^n is
**closed under bitwise majority**. The empty set counts too. Two functions are equivalent if one
becomes the other by renaming variables and/or flipping some variables. That is the
hyperoctahedral group B_n, of order 2^n·n!.

## The structure theorem (proved in the log below, and validated by 26 known terms)
Take a satisfiable 2-CNF and close it under all implied 2-clauses (including unit clauses, and
weakening of units). The closed form splits into three layers:
1. **Forced variables** (k of them): value fixed to 0 or 1.
2. **Equivalent variables:** among the rest, literals forced equal form **blocks** (a block is a
   set of literals, e.g. {x, ¬y, z}, up to global negation). Keep one literal per block.
3. **Skew poset:** the remaining implications among the r block-literals and their negations form
   a strict partial order on 2r literals with  l < m  ⟺  ¬m < ¬l  and no literal comparable to its
   own negation.

Conversely every such triple gives a distinct Krom function (proof sketch in the log). Hence
```
A109457(n) = 1 + Σ_k C(n,k) 2^k L(n−k),   L(m) = Σ_r S2(m,r) 2^(m−r) T(r)
A109459(n) = 1 + Σ_{m≤n} U(m)
A109458(n) = 1 + Σ_{m≤n} (n−m+1) U'(m)
```
where T(r) = number of labelled skew posets on r labelled pairs, U(m) = number of skew posets
with points coloured by block sizes summing to m, up to isomorphisms commuting with negation,
and U'(m) = the same with each block-literal coloured by (number of positive, number of negative
literals in the block), up to isomorphism.

The unlabelled counts are Pólya sums over Aut(P) ⊆ B_r: for a signed cycle of length k the
weight is x^k/(1−x^k) (block sizes), or (2t−t²)/(1−t)² at t = x^k for an even-sign cycle and
t²/(1−t²) for an odd-sign cycle (signed colours). Only **connected** skew posets are Pólya-counted;
the rest follow by an Euler transform (multisets of connected components), so huge groups such as
the antichain's B_9 (order 1.86·10^8) are never enumerated.

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

Level r = 9 takes 20 s on 16 threads (64.5 M nauty calls). Level 10 (573012592 classes,
4.4·10^9 candidates) runs without storage via canonical augmentation (`./skew 10 a`), 20 min 40 s
wall (318 CPU-min); the largest connected automorphism group is S_10 (the "at most one true" structure).
The n = 10 terms have no second independent computation (the coloured check would need ~37 GB);
they rest on the internal checks below and on the augmentation mode matching the stored mode at r ≤ 9.

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
- `b109457.txt`, `b109458.txt`, `b109459.txt` — b-files, n = 0..10.
- `oeis_submission.md` — draft edits (one document for all three sequences).
- `layman.md` — plain-English write-up.
- `asymptotics.md` — literature (BBL, Allen, Ilinca–Kahn), the asymmetry conjecture, partial proofs.
- `../../lean/OeisLean/A109459.lean` — Lean proof of the structure theorem.

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
