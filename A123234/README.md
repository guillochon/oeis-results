# A123234: n X n Latin squares up to row and column permutation

**New terms:**

    a(9)  = 41955285698826624
    a(10) = 758072148410166701948352
    a(11) = 487630706661579208919970322218240

**Status (Oct 18 2026):** computed and checked; OEIS edit drafted, not yet submitted.

This page is the explanation linked from the OEIS entry: what is counted, the formula and its derivation, how each
new term was computed, and how it was checked. Everything needed to reproduce the values is in this folder.

## What is counted

A Latin square of order n is an n X n array with entries from {0, ..., n-1} in which every symbol occurs exactly once
in each row and once in each column. Two squares are in the same class when one is obtained from the other by
permuting its rows and permuting its columns; the symbols are not relabelled ("principal isotopy" or
"RC-equivalence"). a(n) is the number of classes. The OEIS entry gives a(1)..a(8), the last of them from McKay,
Meynert and Myrvold (2007).

## The formula

Write L(n) for the number of Latin squares of order n (A002860; known up to n = 11). For k | n, k >= 2, put m = n/k
and let D_k(n) be the number of ways to choose n permutations from the wreath product C_k wr S_m (the centralizer
in S_n of a permutation with m cycles of length k) that are pairwise disjoint, i.e. no two of them send any point
to the same place. Then

    a(n) = L(n)/(n!)^2  +  sum_{k | n, k >= 2}  n! * D_k(n) / (k^m * m!)^2,      m = n/k.

### Derivation

The group S_n x S_n acts on the L(n) Latin squares, a pair (alpha, beta) moving the entry in cell (r, c) to cell
(alpha(r), beta(c)). By Burnside's lemma a(n) = (1/(n!)^2) * sum Fix(alpha, beta), where Fix counts the squares the
pair leaves unchanged. The identity pair fixes all L(n) squares, which gives the first term.

**Squares as permutations.** For each symbol s let sigma_s be the permutation sending row r to the column where s
lies in row r. A square is the same thing as n permutations sigma_0, ..., sigma_(n-1) that are pairwise disjoint
(sigma_s(r) != sigma_t(r) for s != t): disjointness says each cell gets one symbol, and each sigma_s being a bijection
says each symbol appears once per row and once per column. The pair (alpha, beta) turns sigma_s into
beta sigma_s alpha^(-1), so the square is fixed exactly when

    sigma_s alpha = beta sigma_s      for every symbol s.

**1. alpha and beta are conjugate, and either both are the identity or both are fixed-point-free.** The equation
says alpha = sigma_s^(-1) beta sigma_s, so the two have the same cycle type. If alpha fixes a row r, then
beta(sigma_s(r)) = sigma_s(r) for every s; the n columns sigma_s(r) are all different, so beta fixes every column,
beta is the identity, and then so is alpha.

**2. The fixed squares of one pair.** Fix any tau with tau alpha tau^(-1) = beta. The permutations sigma with
sigma alpha = beta sigma are exactly sigma = tau g with g in the centralizer C(alpha), and tau g, tau h are disjoint
exactly when g, h are. So a fixed square is an ordered n-tuple of pairwise disjoint elements of C(alpha), and

    Fix(alpha, beta) = n! * D(C(alpha)),

where D(G) is the number of n-element sets of pairwise disjoint permutations in G. Such a set sends each point to
each point exactly once (n disjoint permutations, n possible images), so it is a sharply transitive subset of G.

**3. Only cycle types k^m contribute.** If alpha is fixed-point-free with cycles of different lengths, let X be the
set of points lying in cycles of one length k. Every element of C(alpha) maps X onto itself, so none sends a point
of X outside X, and no set of elements of C(alpha) can send that point to every point. So D(C(alpha)) = 0 unless
all cycles have the same length k, and then C(alpha) = C_k wr S_m.

**Putting it together.** The class of type k^m has n!/z elements, z = k^m * m!, so the pairs of that type
contribute (n!/z)^2 * n! * D_k(n) to the Burnside sum. Dividing by (n!)^2 gives the formula. (k = n is the n-cycle,
with D = 1: the only choice is the cyclic group itself.)

### Computing D_k(n)

**General k: clique search.** Let D1 be the number of such sets that contain the identity. For any such set S and
any g in S, g^(-1) S is again such a set and contains the identity; this gives n * D = |G| * D1. A set containing the
identity is the identity plus n - 1 pairwise disjoint derangements in G, exactly one of them sending point 0 to t for
each t = 1, ..., n-1, so D1 is a small clique count. This is how D_3(9) = 2612736 and D_5(10) = 14400 are computed.

**k = 2: a matrix count.** C_2 wr S_m acts on m blocks of two points; an element is a permutation pi of the blocks
together with a flip (swap the two points or not) at each block. Take n = 2m elements forming a sharply transitive
set. For blocks i and j, the two points of block i must reach the two points of block j through two different
elements, so exactly two of the pi's send i to j, and those two have opposite flips at i. Conversely, given an
ordered 2m-tuple of block permutations in which every pair (i, j) is used exactly twice, there are 2^(m^2) ways to
choose the flips (one binary choice per pair (i, j)), and every choice gives 2m pairwise disjoint, distinct elements.
Hence

    D_2(2m) = 2^(m^2) * M_m / (2m)!,

where M_m is the number of ordered 2m-tuples of m X m permutation matrices summing to the all-2 matrix
(M_1..M_5 = 1, 6, 900, 12640320, 41376005798400). M_m is computed by dynamic programming over the matrix of remaining
capacities, fixing the first permutation to the identity and multiplying by m! (relabelling columns acts freely).

### Values

| n  | L(n)/(n!)^2      | nonzero terms of the sum (k: D_k(n))   | a(n) |
|----|------------------|-----------------------------------------|------|
| 8  | 6.691018 * 10^10 | 8: 1, 4: 576, 2: 20545536               | 66915816462 |
| 9  | 4.195529 * 10^16 | 9: 1, 3: 2612736                        | 41955285698826624 |
| 10 | 7.580721 * 10^23 | 10: 1, 5: 14400, 2: 382591593086976     | 758072148410166701948352 |
| 11 | 4.876307 * 10^32 | 11: 1                                   | 487630706661579208919970322218240 |

## Checks

1. **Old terms.** The formula reproduces a(1)..a(8) = 1, 1, 1, 4, 16, 1868, 2420400, 66915816462, and the sum is an
   integer for every n <= 11 (asserted in the program; the individual terms are not integers, so this is a real check).
2. **Burnside without the reduction (`fixcount.c`, `burnside_check.py`).** For every pair of cycle types
   (lambda, mu) of S_n other than the identity pair, a direct backtracking search counts the Latin squares fixed by
   one pair (alpha, beta) of those types, using nothing but the condition L(alpha(r), beta(c)) = L(r, c). (Relabelling
   symbols commutes with the action and acts freely, so it counts the fixed squares whose first row is 0..n-1 and
   multiplies by n!.) Summing over all pairs, weighted by class sizes, and adding L(n)/(n!)^2 gives a(n) again:
   - n = 2..9: all pairs searched (899 at n = 9); **a(9) confirmed**. The only pairs that fix anything are the
     types k^m against themselves, with exactly the counts D_k(n) above, confirming steps 1-3 and the D values.
   - n = 10, 11: all pairs searched except (2^5, 2^5) at n = 10, whose roughly 3.8 * 10^14 normalized fixed squares
     are too many to list one by one; that term is taken from the k = 2 formula (check 3). **a(10), a(11) confirmed.**
3. **M_m two ways (`semi_latin.py`).** M_m is recounted row by row (state: the multiset of columns each of the 2m
   permutations has used so far) rather than permutation by permutation; both give M_1..M_5 above, including
   M_5 = 41376005798400 (23 s). The formula D_2(2m) = 2^(m^2) M_m/(2m)! itself is confirmed by the direct search for
   m = 2, 3, 4 (640 at n = 6, 20545536 at n = 8).
4. **Size.** a(n) is dominated by L(n)/(n!)^2, the count if no square had a symmetry; a(n)(n!)^2/L(n) - 1 =
   8.6e-10, 1.2e-10 and 6.8e-28 for n = 9, 10, 11.
5. **Not already known.** None of the three values, nor the D or M sequences, is in the OEIS (search Oct 18 2026).
   The entry's most recent reference (Cambie and Wanless, arXiv:2610.20008, Oct 2026) quotes values only up to n = 7.

The new terms depend on the published values of L(9), L(10) and L(11) (McKay, Meynert and Myrvold 2007; McKay and
Rogoyski 1995; McKay and Wanless 2005), which enter the first term.

## Reproduce

    python rc_classes.py              # a(4)..a(11), about 30 s (the k = 2 term at n = 10 dominates)
    bash run_fixcount.sh all 2 3 4 5 6 7 8 9          # direct Burnside search, seconds (needs gcc)
    SKIP2=skip2 bash run_fixcount.sh all 10 11
    python burnside_check.py logs/fixcount_*.txt      # sums the searches and compares with rc_classes
    python semi_latin.py 5            # M_1..M_5 two ways

Logs of the runs are in `logs/` (`run_fixcount.log` has the timings, `burnside_check.txt` the comparison).

## Files

- `rc_classes.py`: the formula; L(n) table, clique search for D_k, matrix DP for M_m.
- `fixcount.c`, `run_fixcount.sh`, `burnside_check.py`: the direct Burnside check.
- `semi_latin.py`: the second count of M_m.
- `logs/`: outputs of the checks.

---

## Working log

- 2026-09-28: Burnside formula found during round-4 triage; a(9)-a(11) computed with `rc_classes.py`, reproducing
  a(4)-a(8).
- 2026-10-18: derivation written out; direct Burnside search over all cycle-type pairs (`fixcount.c`) and a second
  count of M_m (`semi_latin.py`) added; both agree.
