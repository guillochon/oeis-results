# A006455 — Naturally labeled posets on {1, …, n}

**Status:** new terms found, each confirmed by at least 3 splits: **a(13) = 10252978223661714**, **a(14) = 2223739823604074838**, **a(15) = 657800677264030983838**, **a(16) = 264134770940774037129920**. The formula is proved in Lean. OEIS edit drafted.
**Known terms:** a(0)..a(12) (13 terms). a(10)..a(12) are from Knuth's POSETS program, Dec 2001; nothing new since.
Cheon–Giraudo ([arXiv:2512.17749](https://arxiv.org/abs/2512.17749), 2025, p. 2) say the sequence is known only to n = 12.
**Reached:** a(13)..a(16). a(17) (split 9+8, roughly 5–8 h) and a(18) (9+9, days) are feasible; beyond that needs a new idea (see the log).
Recent papers cite it: Bevan–Cheon–Kitaev 2023,
Cheon et al. 2026, Mazorchuk 2026, Sack 2026.

## What is being counted
Partial orders R on [n] with x R y ⇒ x < y. Equivalently: posets together with a linear extension, up to relabeling.
The same count also gives n×n upper-unitriangular Boolean matrices B with B = B², and naturally labeled posets are the
states of "grow a poset by adding a new maximal element".
Example n = 3 (a(3) = 7): ∅, {1<2}, {1<3}, {2<3}, {1<2, 1<3}, {1<3, 2<3}, the chain.

Notation: J(P) is the lattice of order ideals (down-sets) of P, i(P) = |J(P)|, e(P) is the number of linear
extensions, and w(P) = e(P)/|Aut P| is the number of natural labelings of P up to equality.

## How the known terms were computed
Knuth's POSETS grows naturally labeled posets one element at a time: element n+1 is maximal, and its down-set is any
order ideal of the poset on [n]. That gives **a(n+1) = Σ over natural Q on [n] of i(Q)**, but it visits every one of the
a(n) objects. a(12) needed a(11) ≈ 5.7×10¹¹ visits; a(13) would need 6.5×10¹³, and each further term is about 100× more.

## Clever ideas to pursue

1. **Split identity (proved below; checked numerically in `split_check.py` for every split n + k ≤ 7).**
   **a(n+k) = Σ over natural Q on [n] and natural R on [k] of i(R^op × Q) = Σ over unlabeled P on n and S on k of w(P) w(S) i(S^op × P).**
   *Proof.* Take a natural poset on [n+k] and restrict it to [n] (giving Q) and to {n+1, …, n+k} (giving R, relabeled to
   [k]). No element of [n] lies above a new element. For each new element y, let φ(y) be its down-set inside [n]; φ(y) is an
   ideal of Q. Transitivity is equivalent to φ being order-preserving R → J(Q), and every triple (Q, R, φ) glues back
   to a natural poset. Order-preserving maps R → J(Q) correspond to down-sets of R^op × Q, via S = {(r, q) : q ∈ φ(r)}. ∎
   - Since i(X) = i(X^op), the summand is symmetric under swapping (P, S), so only unordered pairs are needed.
   - k = 1 recovers Knuth's recursion. k = 2 gives the extension count f₂(P) = i(P)² + #{I ⊆ J in J(P)}.
   - **Cost:** the input is only the unlabeled posets on ⌈N/2⌉ points (A000112: 2045 on 7, 16999 on 8, 183231 on 9), plus
     their e(P) and |Aut P|.

   | Term | Split | Unordered pairs | Size of R^op × Q |
   |---|---|---:|---:|
   | a(13) | 7 + 6 | 6.5×10⁵ | 42 |
   | a(14) | 7 + 7 | 2.1×10⁶ | 49 |
   | a(15) | 8 + 7 | 3.5×10⁷ | 56 |
   | a(16) | 8 + 8 | 1.4×10⁸ | 64 |
   | a(17) | 9 + 8 | 3.1×10⁹ | 72 |
   | a(18) | 9 + 9 | 1.7×10¹⁰ | 81 |

   The inner cost is counting down-sets of a product poset of about 50–80 elements. Do it as a DP over the S side: assign φ(s) ∈ J(P)
   in linear-extension order, with a precomputed ⊆-table on J(P) (at most 512 ideals for |P| = 9). Memoize on the values of
   the elements of S that still have unassigned elements above them (their pathwidth). Expect a(13) and a(14) in minutes;
   a(16) is plausible on the GPU; a(17) and a(18) depend on how fast the inner count gets.
2. **Aggregate the inner sum.** For a fixed lattice L = J(P), precompute F_k(L) = Σ_S w(S)·Hom(S, L) directly, for
   example by a DP over "grow S one maximal element at a time, carrying the values φ on the current maximal elements".
   This could replace the pairwise loop entirely.
3. **Fallback: one-step formula over unlabeled posets on N−1 points,** a(N) = Σ_P w(P) i(P). This needs
   generating all posets on 12 points (1.1×10⁹; nauty `genposetg`, already built in WSL at ~/a109459). It reaches a(13) only,
   so use it as an independent check of a(13).

Checklist:
- [x] Symmetry: sum over unlabeled posets weighted by e/|Aut|
- [x] Decomposition: split identity (meet in the middle)
- [ ] Transfer matrix / DP: down-set counting in R^op × Q
- [ ] Related sequences: A000112 (e(P), |Aut P| per poset); A006455 feeds the average-linear-extension ratio in the
      %C line (a(n)·n!/A001035(n)), and A001035 is known to n = 18
- [ ] Parallelize / GPU once a(13) and a(14) are done on the CPU

## Verification plan
- Reproduce a(0..12) with every split n + k of each N (splits of the same N must agree exactly; a strong internal check).
- a(13) via two splits (7 + 6 and 8 + 5) and via the one-step formula over the posets on 12 points (fallback 3).
- Check the per-poset data for each n: Σ_P n!/|Aut P| = A001035(n) (labeled posets), and Σ_P w(P) = a(n).
- Check u64 overflow: a(14) is about 2×10¹⁸ and a(15) is about 5×10²⁰, so accumulate in u128 from the start.

## Literature check (2026-09-27, Valency corpus + arXiv full text + oeis.org + Knuth's site)
**Verdict: a(13) and a(14) appear new.** Nobody has published beyond a(12), and the split identity for *natural*
labelings does not appear anywhere. It is a close cousin of the Erné–Stege decomposition for labeled posets, though,
so the write-up must credit that. **Update:** the bijection itself is Theorem 2.1 of a Campo–Erné 2018 (below).

| Source | What it has | Relevance |
|---|---|---|
| oeis.org/A006455 (live, checked 2026-09-27) | DATA and b-file stop at a(12) | a(13), a(14) not listed |
| Knuth, POSETS0/POSETS (Dec 2001, `programs/posets.w`) | row-by-row DP over upper-triangular Boolean matrices, with memoized weighted lists; "n = 12 at least" | the source of a(10..12); a different method |
| Cheon–Giraudo–Kwon–Lee, [arXiv:2512.17749](https://arxiv.org/abs/2512.17749) v2 (May 2026) | "known up to n = 12"; one-step recursion a(n+1) = Σ_A (number of ideal vectors of A), plus a Burnside version for the unlabeled count | no new terms; no split |
| Cheon et al., [arXiv:2602.04533](https://arxiv.org/abs/2602.04533) (Feb 2026) | matrix approach, indexing via the binary Pascal matrix | no new terms (abstract) |
| Bevan–Cheon–Kitaev, [arXiv:2311.08023](https://arxiv.org/abs/2311.08023) | naturally labelled posets vs 12-34-avoiding permutations | bijective and structural, no terms |
| Ayala, [arXiv:2606.31526](https://arxiv.org/abs/2606.31526) (Jun 2026) | labeled posets A001035(19), via the **Erné–Stege moment reduction** P(N) = C(N+1,2)P(N−1) − Σ(−1)^(N−m) C(N−1−m,2) C(N,m) G(m, N−m), where G(m,k) = Σ_Q (m!/\|Aut Q\|)·d(Q)^k over unlabeled Q, and d is the ideal count; code at github.com/Rafael-Ayala/posets-and-topologies-19 | **closest prior art.** d(Q)^k = Hom(antichain_k, J(Q)) is the antichain case of our Hom(S, J(P)) term, and Erné–Stege moments are the labeled analogue of our split. Does not mention A006455 or natural labelings |
| Erné–Stege, "Counting finite posets and topologies", Order 8 (1991) | the original moment reduction | not in the corpus. **Still to check:** does it contain a natural-labeling version (Σ over S of w(S)·Hom(S, J(Q)))? |
| a Campo–Erné, [arXiv:1802.01419](https://arxiv.org/abs/1802.01419) (2018), Theorem 2.1 | orders on X ∪ Y that restrict to O on X and P on Y, with nothing in Y below anything in X ("generalized vertical sums"), correspond one-to-one to order-preserving maps (Y,P) → D(O); isomorphic to U(O) ⊗ D(P) | **this is our bijection**, so the identity is known. New here: applying it to *natural* labelings (sum over all natural O, P) and the computation. Credit it in the OEIS comment. (Erné–Stege 1991 itself is paywalled and not read.) |
| Mei–Wang, [arXiv:1709.01234](https://arxiv.org/abs/1709.01234) | counting order-preserving maps of posets | general theory; no A006455 |

Credit in the write-up: "the natural-labeling analogue of the Erné–Stege decomposition; the S = antichain terms are
Erné–Stege moments".

## Log
- 2026-09-27: attack file started. Split identity found, proved, and checked by brute force for all splits n + k ≤ 7
  (`split_check.py`, 1 s).
- 2026-09-27: Rust engine (`rust/`, `gen_posets.sh`). Posets from nauty `genposetg` (in WSL); e(P) by DP over ideals,
  |Aut P| by backtracking. `check 9`: Σ w = A006455(n) and Σ n!/|Aut| = A001035(n) for all n ≤ 9.
  Hom(S, J(P)) by a DP over S in label order. The state is one lattice value per *class* of active elements, where a class
  groups the elements with the same set of future coverers (only the join of a class matters). That made a(12) 150× faster
  than keeping one value per active element (74 s → 0.5 s).
  Reproduced a(2..12) with splits 1+1, 2+3, 3+3, 4+4, 4+5, 5+4, 5+5, 5+6, 6+5, 6+6.
  New: a(13) = 10252978223661714 (splits 7+6 in 8 s, 6+7 in 23 s, 8+5 in 5 s);
  a(14) = 2223739823604074838 (splits 7+7 in 412 s, 8+6 in 171 s). Log: `logs/a13_a14.log`.
- 2026-09-27: speedups (`split n k [threads] [flags]`; `logs/speedups.log`):
  - `o`: processing order of the structure side chosen by a shortest path over its ideal lattice, using a cost model;
  - `d`: each pair runs in the cheaper direction (S over J(P), or P^op over J(S^op); i(S^op × P) = i(P × S^op));
  - `s`: unordered pairs when n = k;
  - u64 counters whenever nk < 64.

  | Split | Before | o | d | ods |
  |---|---:|---:|---:|---:|
  | 7+6 → a(13) | 6.5 s | 1.1 s | 1.2 s | **0.27 s** |
  | 7+7 → a(14) | 412 s | 24.7 s | 25.9 s | **1.8 s** |

  With all flags on: 8+6 in 4.4 s, 6+8 in 4.9 s, 8+5 in 0.45 s, 5+8 in 0.70 s. a(12) takes 0.02 s.
  a(13) is now confirmed by 5 splits (7+6, 6+7, 8+5, 5+8, and the unoptimized 7+6); a(14) by 4 (7+7 with and without
  speedups, 8+6, 6+8). As agreed, nothing beyond a(14) was run.
- 2026-09-27: a(15) = 657800677264030983838 from splits 8+7 (71 s) and 7+8 (72 s), `logs/a15_a16.log`. Runs for
  9+6, then 8+8, 9+7 and 7+9 (a(16)) are in progress. (The log also has a header line from an 8+8 run that was stopped
  by mistake and then restarted.) The a Campo–Erné 2018 paper was read: its Theorem 2.1 is the bijection. Drafted
  `layman.md`, `oeis_submission.md` (a(16) placeholders) and `publish.py` (not yet run).
- 2026-09-27: a(15) confirmed a third time (split 9+6, 123 s). a(16) = 264134770940774037129920 from split 8+8
  (784 s, u128 counters, unordered pairs); confirming split 9+7 is running.
- 2026-09-27: **Lean** (`lean/OeisLean/A006455.lean`, standard axioms only, 11 s to check):
  - `equiv`: natural orders on Fin (n+k) ≃ triples (P, Q, f), with f order-preserving from Q to the down-sets of P;
  - `a_add`: A006455(n+k) = Σ_P Σ_Q #{f}, the identity the program evaluates;
  - `card_linExt`: e(s) = |Aut s| · #{natural orders isomorphic to s}, which justifies the weights w = e/|Aut| used for
    the sum over unlabeled posets;
  - `a_three`: a(3) = 7 by `decide +kernel`.

  Not formalized: that genposetg lists each unlabeled poset exactly once, and the program's evaluation of the sum
  (the terms rest on the computation and the split-agreement checks).
- 2026-09-27: **explored idea 2** (removing the sum over pairs). Result: no cheap route found.
  - *Low rank?* The matrix M[P][Q] = Hom(Q, J(P)) over unlabeled posets (`a006455 matrix n k`) has **full rank**:
    16/16 (4×4), 63/63 (5×5), 63 (5×6), 318/318 (6×6), exact mod 2^61 − 1. So the counts cannot be factored through
    fewer statistics than there are posets; this is consistent with Lovász-type results that homomorphism counts
    separate structures. As a side check, wᵀMw reproduced a(8), a(10), a(11), a(12).
  - *Per-lattice DP?* F(L) = Σ_Q w(Q)·Hom(Q, L) needs, as DP state, the ideal lattice of the partial Q with lattice
    labels. Growing P instead (L = J(P + x_I)) gives Hom(Q, J(P + x)) = Σ over up-sets U of Q of
    #{φ: Q → J(P) monotone, φ(U) ⊆ ↑I}. That turns the state into (Q, U) pairs, which is bigger, not smaller.
  - Realistic remaining gains are constant factors: sharing DP prefixes between Q's in a trie (the Q's share their
    first 4–5 elements; perhaps about 2×, and it conflicts with per-Q order optimization), and the GPU (about 10×?).
    These make a(17) and a(18) cheaper but do not reach a(19) or a(20). A moment-style reduction (as in Erné–Stege
    for labeled posets) remains the open research direction.
- 2026-09-27: **explored a moment-style (Erné–Stege) reduction.** Conclusion: it does not beat the split.
  - The natural analogue of the Erné–Stege moments is μ(m, j) = Σ_{Q natural on [m]} d(Q)^j. By the split bijection
    (with Q an antichain) this is the number of natural posets on [m+j] whose top j labels form an antichain. It is cheap:
    linear in the number of posets on m points.
  - Peeling the longest suffix of labels that are maximal elements gives
    a(N) = 1 + Σ_{t≥1} [μ(N−t, t) − μ(N−t−1, t+1)], using Σ_{Q on [n]} d(Q minus its top label)^t = μ(n−1, t+1).
    This is true (checked for N ≤ 8), but it **telescopes** to the trivial a(N) = μ(N−1, 1): no information.
  - Peeling the whole set of maximal elements, as Erné does, breaks down because of the labels: each new maximal
    element may only see the elements with smaller labels, so the weights become products of the ideal counts
    d(Q_{≤c}) of the *prefixes* of Q (a complete homogeneous polynomial in them). That is a sum along the whole growth
    tree of Q, on about N − #max ≈ 3N/4 points: more expensive than the split.
  - Even for labeled posets, the Erné–Stege reduction gains only about 3 points (Ayala: P(19) from posets on ≤ 16 points),
    because it still needs posets on up to N − 3 points. The split already needs only N/2. At N = 19: posets on 16 points
    (4.5×10¹⁵) vs pairs on 9 + 10 points (4.7×10¹¹). So a moment reduction would be a step backwards unless combined
    with the split.
  - Possible hybrid, untested: label N is always maximal, so a(n+k) = Σ_{P, Q′ on k−1} Σ_φ i(P ∪_φ Q′). This absorbs one
    element into the DP (a(17) from 9+7 pairs instead of 9+8, about 8× fewer pairs). But the DP state must also carry
    the membership of active elements in the extra ideal and its running join (a factor of roughly 2^a·|L|), which
    probably costs more than it saves.
- 2026-09-27: a(16) confirmed by splits 9+7 (1808 s) and 7+9 (1752 s), matching 8+8. Every new term now has at least
  3 agreeing splits. b-file `b006455.txt` (n = 0..16); OEIS draft and `layman.md` completed.
