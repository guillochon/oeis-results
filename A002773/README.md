# Matroids on 9 and 10 elements: simple and loopless counts (A002773 and family)

**A002773(n)** counts simple matroids (combinatorial geometries) on n points up to isomorphism. The known
terms stopped at a(9) = 376467 (Mayhew–Royle 2008). In 2026 Sunghyeon Jo counted *all* matroids on 10
elements ([github.com/ainta/matroid-count v2.0.0](https://github.com/ainta/matroid-count/tree/v2.0.0); A055545(10) =
3232000741644), by rank, both unlabeled and labeled. That count does not separate simple or loopless matroids.
This page derives them exactly from Jo's totals and the Mayhew–Royle catalogue of matroids on at most 9 elements:

    A002773(10) = 3231997106202     simple matroids on 10 points
    A058718(10) = 3232000358472     loopless matroids on 10 points

Pictures of every simple matroid on up to 8 points, a sample on 9 and the counts by rank for 10:
https://guillochon.github.io/matroid-gallery/

The same data complete about fifteen related entries (table in section 4). No new enumeration of 10-element
matroids is involved: everything follows from Jo's by-rank totals by two exact identities plus a computation over the 9-element
catalogue that takes about 2 minutes.

## 1. Identities

Write U(n, r) for the number of unlabeled matroids on n elements of rank r, and similarly for loopless
(UL) and simple (US). Write L, LL and LS for the labeled counts.

1. **Loops.** A matroid with at least one loop is a matroid on n − 1 elements plus a loop. Deleting a loop is
   a bijection on isomorphism classes, because all loops are alike. So

       UL(n, r) = U(n, r) − U(n − 1, r).

2. **Parallel classes.** A loopless matroid N determines its simplification M, a simple matroid on k ≤ n points.
   It also determines the sizes w_1, …, w_k ≥ 1 of its parallel classes, with sum n and the same rank.
   Conversely (M, w) rebuilds N, and N ≅ N′ iff M ≅ M′ and w′ = w∘g for some g ∈ Aut(M). Hence

       UL(n, r) = Σ_{k ≤ n} Σ_{simple M on k points, rank r}  #{Aut(M)-orbits on compositions of n into k parts}.

   The k = n term is US(n, r), since the only composition is all ones. For n = 10, every other term has
   k ≤ 9 and is in the catalogue:

       US(10, r) = U(10, r) − U(9, r) − Σ_{k ≤ 9} Σ_M #orbits(M, 10).

3. **Labeled.** Choosing the loops, and then the partition into parallel classes, gives

       L(n, r)  = Σ_j C(n, j) · LL(n − j, r),        LL(n, r) = Σ_k S2(n, k) · LS(k, r),

   which invert recursively (S2 is the Stirling number of the second kind). L(10, r) is in Jo's result.
   L(n, r) for n ≤ 9 is Σ n!/|Aut(M)| over the catalogue.

## 2. Computation

* **Catalogue:** Mayhew and Royle, *Matroids on 9 elements*, Zenodo record 6825419, file `matroids09_rankLine`
  (385370 matroids, one rank function per line; SHA-256 `38afc036f42c0ae981354f4cdfd79f2f42b3e57ba1eb5c1ffa9838d33e9399dd`).
* **Jo's data:** `jo_result.json` (the n = 10 totals by rank from the v2.0.0 release) and `jo_cmp.json` (his
  comparison with the catalogue, which gives U(n, r) and L(n, r) for n ≤ 9).
* **Two independent orbit programs** (both give identical output, bit for bit):
  - `orb.c` (Nov 2 2026). Aut(M) by backtracking on the rank function over all subsets, pruned by
    element invariants. Orbits by **Burnside**: the average over Aut(M) of the number of compositions fixed by g
    (the coefficient of x^n in Π_cycles x^ℓ/(1 − x^ℓ)).
  - `orbits_hyp.py` (Nov 3 2026). Aut(M) as permutations preserving the **hyperplanes** (which determine
    the matroid). A strong generating set is found level by level along the point-stabilizer chain, pruned by pairwise
    hyperplane-size profiles. Orbits are counted **explicitly**, by union–find of compositions under the generators,
    with no Burnside sum. It also returns |Aut(M)| as the product of the basic orbit lengths, and so the labeled counts.
* `assemble.py` combines everything, checks it (section 3) and writes `results.json`. `drafts.py` lays out the
  OEIS rows and b-files.

The non-simple loopless part of n = 10 (output row s = 10, by rank 0..10):

    0, 1, 40, 3513, 1616845, 1626070, 5577, 208, 15, 1, 0      (total 3252270)

so A002773(10) = (3232000741644 − 383172) − 3252270 = 3231997106202.

## 3. Checks

* **Reproduction:** every known term of 21 related entries comes out of the same pipeline: A055545, A058718,
  A002773, A058731, A058719, A058693, A336704, A058673, A058712, A058721, A058722, A058715, A058681,
  A100728, A058687, A058669, A058710, A058720, A058716, A058730 and A053534 (`assemble_out.txt`). In particular,
  identity 2 at n ≤ 9 equals the catalogue's own loopless counts, and the labeled counts Σ n!/|Aut| equal Jo's
  independently generated labeled counts for every (n ≤ 9, r).
* **Values at n = 10 already known from other sources:**
  - US(10, 3) = 5249 (A058731) and US(10, 4) = 4884573865 (A058730, Matsumoto–Moriyama–Imai–Bremner 2012).
    These are the two values our subtraction produces that were computed independently before.
  - UL(10, 3) = 8762 (A058719).
  - LL(10, 3) = 16122092568 (A058715).
  - L(10, 2) = 677546 (A058681, closed form).
  - LS(10, 8) = 650761 (A100728).
* **Closed forms at n = 10:**
  - US(10, 9) = 8 and LS(10, 9) = 968 = 2^10 − 1 − 10 − 45: a rank-(n−1) matroid has one circuit, which must
    have size at least 3.
  - US(10, 8) = 104, by duality: the dual is a rank-2 matroid with no coloops and no series pairs. These are
    counted by partitions, and the count gives 1, 4, 11, 22, 40, 66, 104 for n = 4..10.
  - UL(10, 2) = 41 = p(10) − 1 and LL(10, 2) = 115974 = Bell(10) − 1.
  - Every count is non-negative, and the labeled and unlabeled rows are rank-symmetric where they should be
    (all matroids under duality).

## 4. New values

| entry | what | new |
|---|---|---|
| A002773 | simple, unlabeled | a(10) = 3231997106202 |
| A058730 | simple, unlabeled, by rank | row 10: 1, 5249, 4884573865, 3222226143160, 4886374072, 9742, 104, 8, 1 |
| A058718 | loopless, unlabeled | a(10) = 3232000358472 |
| A058716 / A058717 | loopless, unlabeled, by rank | row 9: 0, 1, 29, 950, 189274, 189889, 1217, 79, 8, 1; row 10: 0, 1, 41, 8762, 4886190710, 3222227769230, 4886379649, 9950, 119, 9, 1 |
| A053534 | all, unlabeled, by rank | row 10 (Jo): 1, 10, 128, 10037, 4886380924, 3222227959444, 4886380924, 10037, 128, 10, 1 |
| A058673 | all, labeled | a(9) = 126240373826, a(10) = 11727995799397397461 |
| A058669 | all, labeled, by rank | rows 9 and 10 (b-file) |
| A058687 | all, labeled, rank 3 | a(9) = 78968940, a(10) = 16869404175 |
| A058712 | loopless, labeled | a(9) = 126147335851, a(10) = 11727994537462688545 |
| A058710 / A058711 | loopless, labeled, by rank | rows 9 and 10 (b-files) |
| A058721 | simple, labeled | a(9) = 125807231446, a(10) = 11727988869192909352 |
| A058720 | simple, labeled, by rank | rows 9 and 10 |
| A058722 | simple, labeled, rank 4 | a(9) = 62680769286, a(10) = 17713104354273649 |

Full rows are in `results.json` and `assemble_out.txt`.

## 5. Reproduce

    # catalogue: https://zenodo.org/records/6825419 -> matroids09_rankLine (check the SHA-256 above)
    python orbits_hyp.py matroids09_rankLine 6 > orbits_hyp_out.txt     # ~2 min on 6 cores
    gcc -O2 -o orb orb.c && ./orb matroids09_rankLine                    # same rows, Burnside route (~40 s)
    python assemble.py > assemble_out.txt                                # identities + all checks
    python drafts.py                                                     # OEIS rows, b-files

## Log

- 2026-11-02: lead found while revisiting high-significance entries. `orb.c` and `lab.py` gave the derived values.
  Ranks 3 and 4 matched the published values.
- 2026-11-03: second implementation, `orbits_hyp.py` (hyperplanes, strong generators, explicit orbits, labeled
  counts). It reproduces `orb.c` exactly. Its labeled n ≤ 9 counts match Jo's. `assemble.py` now checks 21 entries.
  Confirmed on oeis.org that A002773 and its siblings were unchanged since the Oct 23 snapshot.
