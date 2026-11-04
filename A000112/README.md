# A000112(17): unlabeled posets on 17 points by Burnside over the labeled counts

**a(17) = 378180373015320248.** It follows that A000608(17) = 373626005846478701 (connected posets,
the inverse Euler transform of A000112).

Interactive companion page: https://guillochon.github.io/seventeen-point-posets/ (the method worked out
for N ≤ 6, where a(n) comes from, and every poset up to 7 points).

A000112 stopped at a(16) (Brinkmann and McKay, Order 19 (2002)). The labeled counts P(n) = A001035(n)
reach further, now to n = 19 (Ayala, arXiv:2606.31526, 2026).

This page turns the labeled counts into the unlabeled ones with Burnside's lemma:
- 99.34% of the sum for a(17) is a combination of known A001035 values;
- the remaining 0.66% comes from a few labeled statistics that we compute exactly;
- the same program reproduces a(1)–a(16).

## 1. Burnside over S_N

    N! a(N) = Σ_λ |C_λ| fix(λ).

Here λ runs over the cycle types of S_N, and |C_λ| is the class size. fix(λ) counts the partial orders
on [N] invariant under one fixed permutation g of type λ.

## 2. The quotient theorem

Let g have cycles C_1..C_r of lengths l_1..l_r, and let < be a g-invariant partial order.

1. **Each cycle is an antichain.** If x < g^k x, then x < g^k x < g^{2k} x < … < x.
2. **The cycles carry a quotient order Q.** Write C_i ≺ C_j if some element of C_i is below some element
   of C_j. Antisymmetry follows from 1. Transitivity holds because < is g-invariant.
3. **Decorations.** Index C_u by Z_{l_u}, with g acting as +1. For C_u ≺ C_v, the relation between the two
   cycles is a nonempty union of ⟨g⟩-orbits on C_u × C_v. These orbits are the residues of (t − s) mod
   d_uv = gcd(l_u, l_v). So the relation is a nonempty set S_uv ⊆ Z_{d_uv}, and
   x_u(s) < x_v(t) iff (t − s) mod d_uv ∈ S_uv.
4. **Transitivity through u ≺ w ≺ v.** Let g = gcd(d_uw, d_wv), which divides d_uv. Then every
   c ∈ Z_{d_uv} with c mod g ∈ (S_uw + S_wv) mod g must lie in S_uv. In particular, a middle cycle with
   g = 1 (a fixed point, or any cycle coprime to u or to v) forces S_uv = Z_{d_uv}.

Conversely, every labeled poset Q on the cycles, together with sets S_uv satisfying 4, gives exactly one
invariant order. Hence fix(g) = Σ_Q W(Q), where W(Q) is the number of valid decorations.

- **Coprime classes.** If the lengths are pairwise coprime (fixed points count as length 1), then every
  d = 1 and W ≡ 1, so **fix(g) = P(r) = A001035(r)**. This covers the identity, transpositions, single
  cycles, 2+3 and so on: 99.34% of the sum for a(17).
- **Interacting cycles.** A cycle is interacting if its length is > 1 and shares a factor with another
  such cycle. Only interacting cycles matter; every other cycle behaves like a fixed point. With t
  interacting cycles,

      fix = (1/(r)_t) Σ_{Q labeled on r} Σ_{injective φ: interacting cycles → Q} W(Q, φ),

  and W depends only on the *pattern* of the marked points. The pattern records which marked pairs are
  comparable and which comparable pairs are forced, i.e. have between them an unmarked point or a
  marked cycle of gcd 1.

## 3. How each fix value was computed

Write r for the number of cycles. "Hard" classes are those with interacting cycles.

| r | method | tool |
|---|---|---|
| ≤ 10 | sum over unlabeled posets Q (weight r!/\|Aut Q\|) of the generic decoration count | `burnside`, posets from `genposetg r t` |
| 11, 12, 13 | the same sum, streamed from nauty's `genposetg` (sweeps) | `stat`, `run_stat2.sh`, `combine_stats.py` |
| 14, 15 | four labeled statistics by a three-block split, with no enumeration of 14- or 15-point posets | `blocks1`, `block_stats.py` |

For N ≤ 17, only these values with r > 10 are needed (one line each in `stats.txt`):
- **r = 11:** the sweep gives (2,2), (3,3), (4,4), (2,2,2), (2,2,2,2), (2,2,2,2,2), (2,2,2,2,2,2), (3,3,3)
  and (3,3,2,2). (4,2,2,2) equals (2,2,2,2); see §3.3.
- **r = 12:** the sweep gives (2,2,2,2) and (2,2,2,2,2). The blocks give (2,2), (2,2,2) and the equivalent
  pairs; (2,2) and (2,2,2) were also swept.
- **r = 13:** the sweep gives (2,2,2,2). The blocks give (2,2), (2,2,2), (3,3), (4,4), …; (2,2) and
  (2,2,2) were also swept.
- **r = 14:** the blocks give (2,2), (2,2,2) and the equivalent pairs and triples.
- **r = 15:** the blocks give (2,2).

### 3.1 Two and three interacting cycles: Hasse-diagram statistics

Let E be the number of Hasse edges (covers). For two interacting cycles with d = gcd, W = 2^d − 1 on a
cover and 1 otherwise, so

    Σ_φ W = r(r−1) + 2(2^d − 2) E,   fix(r; pair) = P(r) + 2(2^d − 2) K(r),

where K(r) = #{labeled posets on [r] with 1 ⋖ 2}.

For three 2-cycles, W = 3^{covers in the triple}, with one exception: on a chain a ⋖ b ⋖ c whose interval
(a, c) is exactly {b}, S_ac must contain S_ab + S_bc, so W = 13 instead of 9. Hasse diagrams have no
triangles, which gives

    fix(r; 2,2,2) = 6[C(r,3)P + 2(r−2)·r(r−1)K + 4(2r·C(r−1,2)V + r(r−1)(r−2)Ch)
                    + 4 r(r−1)(r−2) Lo] / (r)_3,

where:
- V = #{1 ⋖ 3, 2 ⋖ 3} (and Λ = V by duality);
- Ch = #{1 ⋖ 2 ⋖ 3};
- Lo = #{1 ⋖ 2 ⋖ 3 with (1,3) = {2}}.

Both formulas are asserted equal to the generic count on every poset with ≤ 10 points.

### 3.2 The three-block split (K, V, Ch, Lo at 14 and 15 points)

**The split.** Take the hub h to be point 2 (point 3 for V). Let P1 = down(h), P3 = up(h), and P2 = the
rest. Then:
- P1 and P1 ∪ P2 are ideals;
- all of P1 lies below h, and below all of P3;
- the poset is exactly the data (P1, P2, P3, φ: P2 → J(P1), ψ: P3 → J(P2)), with φ and ψ
  order-preserving (J = lattice of ideals).

The marks say that 1 (and 2, for V) is maximal in P1, and that 3 is minimal in P3 (Ch, Lo).

**The formula.** Write Φ_m,κ(Y) = Σ_{X labeled on m, marks κ} Hom(X, J(Y)), and use
Hom(X, J(Y)) = Hom(Y, J(X)). Then

    stat(N) = Σ_{b1+b2+s=N−1} multinomial · Σ_{Y unlabeled on b2} (b2!/|Aut Y|) Φ_b1(Y) Φ_s(Y).

Lo also needs that no b ∈ P2 lies strictly between 1 and 3. Let U = {b : 1 ∈ φ(b)}, which is a filter.
The condition is ψ(3) ∩ U = ∅, so Lo is

    Σ_U Σ_{V ideal, V ∩ U = ∅} Φ_{b1−1}(Y^{+U}) Φ_{s−1}(Y^{V+}),

which is one zeta transform over J(Y).

**Data needed, by block size:**
- **≤ 10 points:** held in memory. A Hom count is a dynamic program over one side's plan and the other
  side's ideal lattice; each pair uses the orientation whose dense state is smaller.
- **11 and 12 points:** these occur only next to partners of ≤ 3 points. They come from streamed
  aggregates (`run_agg.sh`, `aggs/`), which store Hom(T, J(Q)) for every Q and every T on ≤ 3 points.
- **13 points:** these occur only with trivial partners, and come from identities checked for small m.
  Let G(m,k) = Σ_{X labeled on m} |J(X)|^k. Then:
  - #{X on m with a given point maximal} = G(m−1,1), and #{X with two given points maximal} = G(m−2,2);
  - Σ_{X on m, 1 maximal} |J(X)| = G(m−1,2) + Σ_X Hom(X, C_3);
  - peeling all maximal elements gives P(N) = Σ_{k≥1} (−1)^{k+1} C(N,k) G(N−k,k), hence G(13,1)
    from P(14).

So K(15) needs nothing beyond 12-point posets.

### 3.3 Four or more 2-cycles (sweeps at 11–13 points)

For t two-cycles, W(M) = 3^{e(M)}, where e counts the covers inside the marked set M. The exception is
when M contains a long unforced pair. With exactly one such pair (a lonely chain), W = 13·3^{e−2}. So

    Σ_M W = Σ_M 3^{e(M)} + Σ_{lonely L} Σ_{M ⊇ L} 4·3^{e(M)−2} + (exact corrections for the few M
            with several long pairs, or one whose interval has ≥ 2 points).

For t = 4, Σ_M 3^{e(M)} is a closed form in E, 2-paths, 3-stars, 3-paths and 4-cycles of the Hasse
diagram (`all2` in main.rs). The corrections are computed exactly and cached by pattern. Every fast path
is asserted equal to the generic count on all posets with ≤ 10 points.

Cycle types in which every pair has gcd 2, such as (4,2,2,2) and (2,2,2,2), behave identically: all
d = 2, so the residue conditions are the same.

## 4. Result and checks

From `logs/burnside17.out` (`a000112 burnside 17`, 4.9 h, mostly spent on the generic r ≤ 10 classes):

    a(1)–a(14)                     [OK]
    a(15) = 68275077901156         [OK]  hard classes 1.33%
    a(16) = 4483130665195087       [OK]  hard classes 0.94%
    a(17) = 378180373015320248     [new] hard classes 0.66%

- **Known terms.** a(1)–a(16) are reproduced exactly. a(15) and a(16) use the same machinery as a(17):
  sweeps at 11–12 points and block statistics at 13–14 points.
- **Divisibility.** The Burnside sum is asserted divisible by N!.
- **Labeled totals.** For every r ≤ 10, Σ r!/|Aut Q| = A001035(r). Every sweep reproduces the exact count
  A000112(R) of posets and the labeled total A001035(R).
- **Values computed two ways agree exactly:**
  - K(11) and K(12): sweep vs. blocks;
  - fix(12; 2,2) and fix(12; 2,2,2): 12-point sweep vs. blocks;
  - fix(13; 2,2) and fix(13; 2,2,2): 13-point sweep vs. blocks;
  - fix(10; 2,2) and fix(10; 2,2,2): 10-point sweep vs. blocks.

  `burnside` also asserts that duplicate keys in stats.txt are equal.
- **Small cases.**
  - K, V, Ch and Lo were brute-forced for N ≤ 6 (K = 1, 5, 49, 839, 23461; Lo = 0, 1, 7, 93, 2119).
  - The quotient theorem was checked against brute force for n ≤ 6 (`burnside_check.py`).
- **Growth.** a(n)/a(n−1) = 39.6, 51.0, 65.7, 84.4 for n = 14..17, with ratios of ratios 1.29, 1.29, 1.28.

Statistics used (all in `stats.txt`):

    K(13) = 21977493194713940999         V(13)  = 3122366604409097065
    Ch(13) = 1765581128430222778         Lo(13) = 829651602074902773
    K(14) = 12254356018534195811701      V(14)  = 1723018854782285348913
    Ch(14) = 947969580026294188976       Lo(14) = 424414419295836649159
    K(15) = 9417709865854114258491335
    fix(13; 2,2,2,2) = 1814990835384517121055   (13-point sweep, 3.4e10 posets)
    fix(12; 2,2,2,2,2) = 21702243686686532739   (12-point sweep, 1.1e9 posets)

A000608(17) is the inverse Euler transform of a(0..17). The same transform reproduces A000608(0..16).

## 5. Reproduce

Requirements: Rust, and nauty 2.8.9 (`genposetg`), which was run in WSL here.

    genposetg n t > ../A006455/data/posets{n}.d6     # n ≤ 10: the posets for the small r
    cd rust && cargo build --release
    bash run_agg.sh 11 6; bash run_agg.sh 12 6        # aggregates for the block engine (~80 min)
    a000112 blocks1 14 K   (and V, Ch, Lo; then 15 K)  # block statistics
    python block_stats.py 14 K V Ch Lo                 # → stats.txt lines
    bash run_stat2.sh 11 "<keys>" 6 s11; python combine_stats.py 11 s11      # sweeps
    bash run_stat2.sh 12 "2,2,2,2/2,2,2,2,2/2,2/2,2,2" 6 s12; python combine_stats.py 12 s12
    bash run_stat2.sh 13 "2,2,2,2/2,2/2,2,2" 6 s13; python combine_stats.py 13 s13
    a000112 burnside 17

Costs on an i5-12600K with 6 workers:
- 13-point sweep (3.4e10 posets, 6 slices): ~10.5 h per slice;
- 12-point sweep: 5.4 h on 2 slices;
- K(15): 11 h and 7.5 GB;
- `burnside 17`: 4.9 h.

`stat` checkpoints every 50M posets (A000112_CKPT). genposetg's output order is deterministic, so a
restarted slice resumes where it stopped.

## Files

- `burnside_check.py`: brute-force checks, n ≤ 6.
- `rust/src/main.rs`: `burnside`, `stat` (generic and fast W, `all2`), automorphism counts.
- `rust/src/blocks.rs`: the three-block engine (`blocks`, `blocks1`, `agg`).
- Sweeps: `run_stat2.sh` + `combine_stats.py` (sliced), `run_stat.sh` (single key).
- Aggregates: `run_agg.sh` + `aggs/`.
- `block_stats.py`: turns K, V, Ch and Lo into stats.txt lines.
- `stats.txt`: one "r lengths value" line per fix value with r > 10. `stats_raw/`: raw sweep sums per
  slice.
- `logs/`: run output, including `burnside16.out` and `burnside17.out`.
- `b000112.txt`, `b000608.txt`: the tables to n = 17.

## Notes

- **Speed.**
  - |Aut| factors out incomparable twins (k twins give k!) and backtracks only the rest.
  - The exact multi-set weights in `all2` are cached by pattern.
  - Together these gave ~50× on 12–13 points.
  - Without choosing the Hom orientation per pair, K(15) needed > 30 GB.
- **a(18)** would need (16; 2,2), that is K(16). That means blocks of 13–14 points next to small partners,
  at the scale of Ayala's computation (github.com/Rafael-Ayala/posets-and-topologies-19).
- **A001930** (unlabeled topologies) is the same construction on preorders with A000798. There, cycles can
  merge into equivalence classes, so the theorem needs redoing.
