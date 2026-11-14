# A001930(17): unlabeled topologies on 17 points, and 29 related entries, from the cycle index of posets

**A001930(17) = 457507251208437773** (topologies, or quasi-orders, on 17 unlabeled points; the entry
stopped at n = 16, Brinkmann and McKay 2005).

**A091073(16) = 292885432498738892585** and **A091073(17) = 49464136530826608408726** (transitive
relations on n unlabeled points; the entry stopped at n = 15).

Nothing new was enumerated. The computation of A000112(17) (unlabeled posets on 17 points, see
[../A000112](../A000112)) produced, as a by-product, the number of posets fixed by a permutation of
every cycle type, for every n ≤ 17. That table is the **cycle index series of the species of posets**
to degree 17. Topologies, transitive relations, prime posets and the other entries below are all
substitutions into that series, so each takes seconds.

`python family.py` (3 s) recomputes everything, checks it against all terms in the OEIS (30 entries,
every known term reproduced), and writes `family.md`, `family.json` and the b-files.

## 1. The cycle index of posets

For a permutation g of [k] with cycle type λ, let fix_k(λ) be the number of partial orders on [k] that
g preserves, and z_λ = Π_i i^{m_i} m_i! (m_i = number of i-cycles). Put

    Z_P = Σ_{k ≥ 1} Σ_{λ ⊢ k} fix_k(λ)/z_λ · p_λ,       p_λ = p_{λ1} p_{λ2} ⋯  (power sums).

- Setting p_i → x^i gives Σ a(k) x^k with a = A000112 (Burnside's lemma).
- The coefficient of p_1^k is A001035(k)/k! (labeled posets).

The values fix_k(λ) are in `../A000112/fix_small.txt` and `../A000112/stats.txt` (classes whose cycle
lengths are pairwise coprime have fix = A001035(number of cycles)); the A000112 page explains how they
were computed and checked. `poset_cycle_index.json` holds Z_P itself, all 1211 coefficients for
k = 1..17.

## 2. Topologies and transitive relations: weights on the points

**Topologies.** A topology on a finite set is the same thing as a quasi-order (reflexive and transitive
relation). Its equivalence classes x ≤ y ≤ x form a poset Q, and the quasi-order is Q together with a
size w(q) ≥ 1 for every point q. So unlabeled quasi-orders on n points are the orbits of S_k on pairs
(Q, w) with Q a poset on [k] and Σ w = n. A pair is fixed by g when Q is fixed and w is constant on the
cycles of g, so by Burnside

    A001930(n) = Σ_{k ≤ n} Σ_{λ ⊢ k} fix_k(λ)/z_λ · [x^n] Π_{parts l of λ} x^l/(1 − x^l).          (1)

In species language: quasi-orders = posets ∘ nonempty sets, and (1) is Z_P with p_i → x^i/(1 − x^i).

For n = 17 the terms of (1) by number of classes k are the row n = 17 of the triangle T(n, k) in
`results.json`. The k = 17 term is A000112(17) (83%), the k = 16 term is the number of 16-point posets
with one marked point (15%), and so on.

**Transitive relations.** In a transitive relation, x ~ y iff x = y or (xRy and yRx) is an
equivalence; a class with two or more points is a reflexive clique, a one-point class is reflexive or
irreflexive, and the classes form a (strict) poset. So the weight of a point becomes x + x/(1 − x):

    A091073(n) = Σ_{k ≤ n} Σ_{λ ⊢ k} fix_k(λ)/z_λ · [x^n] Π_l ( x^l + x^l/(1 − x^l) ).              (2)

**Antisymmetric transitive relations** (A079265): every point is reflexive or irreflexive, weight 2x,
so a(n) = Σ_{λ ⊢ n} fix_n(λ)/z_λ · 2^{number of parts}.

The connected versions (A001928, A296105) are inverse Euler transforms, A326898 is the partial sums of
A001930, and A399304 is the first differences of A091073.

## 3. Prime, irreducible and reduced posets: plethysm

These need the whole series Z_P, not only its specializations. Symmetric functions are kept in the
power-sum basis; f[g] (plethysm) replaces each p_k in f by g with every p_j → p_{jk}; species
composition F ∘ G has cycle index Z_F[Z_G] (`prime.py`).

**Irreducible** (Wright): not an ordinal sum of two nonempty parts. Every poset is uniquely an ordinal
sum of irreducible ones, so P = S/(1 − S) and the g.f. is 2 − 1/B(x), B = g.f. of A000112 (A046907) or
of A001930 (A046911). A disconnected quasi-order is irreducible, which gives A046909.

**Prime** (A046904, A354693): no nontrivial module, i.e. no proper subset T with ≥ 2 points such that
every outside point is below all of T, above all of T, or incomparable to all of T. By modular
decomposition, a poset on ≥ 2 points is exactly one of

- disconnected (a set of ≥ 2 connected posets),
- an ordinal sum (a sequence of ≥ 2 irreducible posets),
- a prime poset on ≥ 4 points with a poset substituted for each point.

With C = connected and S = irreducible posets this reads Pr ∘ P = C + S − P − X, so

    Z_Pr = (Z_C + Z_S − Z_P − p_1) [ Z_P^{⟨−1⟩} ],

where Z_C = Σ_k μ(k)/k · p_k[log(1 + Z_P)], Z_S = Z_P/(1 + Z_P), and the plethystic inverse is solved
degree by degree. The labeled version (A046905) is the same identity on exponential generating functions
and needs only A001035, so it reaches n = 19.

**Reduced** (Bower's A066302–A066305; Bergeron–Labelle–Leroux p. 57): P = R ∘ L₊, a reduced poset with
a nonempty chain substituted for each point; equivalently no two-point chain is a module. The labeled
entry A066302 states the e.g.f. B(x/(1 + x)). Unlabeled: Z_R[p_1/(1 − p_1)] = Z_P. This definition is
inferred: the unlabeled entries carry none, but it reproduces all their terms (7 and 6).

**Prime and composite topologies** (Kreminski, A028849–A028852 and A028855–A028858): a space is
composite if it is homeomorphic to a product of two spaces with ≥ 2 points. The counts follow from
A001930, A000112, A001928 and A000608 by a "multiplicative Euler transform": total(n) = number of
multisets of prime spaces whose sizes multiply to n. This assumes unique factorization. For connected
spaces that is a theorem (Hashimoto 1951 for posets, McKenzie 1971 for reflexive structures). For
disconnected spaces it can fail in general (the classical example (1 + A + A²)(1 + A³) = (1 + A)(1 + A² + A⁴)
has at least 63 points), and I have not proved that it cannot fail below 18 points. All 10 known terms
of each composite entry agree. Treat A028849, A028850, A028855 and A028856 beyond n = 11 as
conditional on that.

## 4. Results

New terms only; `family.md` has the full list, `b*.txt` the tables.

| Entry | What | Was | New |
|---|---|---|---|
| A001930 | topologies, unlabeled | n ≤ 16 | a(17) = 457507251208437773 |
| A001928 | connected topologies | n ≤ 16 | a(17) = 451681979218132943 |
| A326898 | topologies with up to n points | n ≤ 16 | a(17) = 463236301242035204 |
| A091073 | transitive relations | n ≤ 15 | a(16) = 292885432498738892585, a(17) = 49464136530826608408726 |
| A296105 | connected transitive relations | n ≤ 15 | a(16) = 288339872230971802596, a(17) = 48869097696701613134501 |
| A399304 | transitive relations, no isolated point | n ≤ 15 | a(16) = 290658556194161944036, a(17) = 49171251098327869516141 |
| A079265 | antisymmetric transitive relations | n ≤ 15 | a(16) = 276814330775821716921, a(17) = 47194522493307302356312 |
| A354693, A046904 | prime posets | n ≤ 12 (n ≤ 7) | a(13) = 14768822778, a(14) = 668464176022, a(15) = 38242248559266, a(16) = 2766797802270655, a(17) = 253291359664557411 |
| A066304 | reduced posets | n ≤ 6 | a(7)–a(17): 984, 8548, 97870, 1467323, 28625346, 722710229, 23511594539, 982371552249, 52587632218131, 3600161822649381, 314783357690940315 |
| A066305 | connected reduced posets | n ≤ 6 | a(7)–a(17): 824, 7505, 88974, 1366952, 27136223, 693843768, 22785439520, 958795636035, 51603690210047, 3547524286134711, 311181144011713324 |
| A046907 | irreducible posets | n ≤ 16 | a(17) = 369142784853079801 |
| A046911 | irreducible quasi-orders | n ≤ 7 | a(8)–a(17): 24591, 275056, 3860200, 68757519, 1568806975, 46116937192, 1751373984753, 85975310423697, 5451790807552624, 445954681497401168 |
| A046909 | connected irreducible quasi-orders | n ≤ 7 | a(8)–a(17): 18127, 226447, 3398240, 63087678, 1477637487, 44186559430, 1697591820169, 84010188430934, 5357998016607666, 440129409507096338 |
| A381121 | posets covering n points | n ≤ 16 | a(17) = 373697242350125161 |
| A028849–52, A028855–58 | composite and prime topologies (Kreminski) | n ≤ 11 or n ≤ 7 | to n = 17 (see family.md; conditional, §3) |
| A046905 | prime posets, labeled | n ≤ 8 | a(9)–a(19) (family.md) |
| A066302, A066303 | reduced posets, labeled | n ≤ 15 | a(16)–a(19) (family.md) |
| A085628, A001927, A046908, A046906 | labeled, from A001035(19) | n ≤ 18 | a(19) (family.md) |

## 5. Checks

- **Known terms.** `family.py` holds the OEIS terms of all 30 entries (snapshot Oct 23 2026, rechecked
  live Nov 14) and reproduces every one: 17 terms of A001930, 16 of A091073, 12 of A354693, and so on.
  The reproductions use the fixed-point table with many different weights, so they test the individual
  values fix_k(λ), not just the Burnside total.
- **Integrality.** Every coefficient of the triangle T(n, k), of the prime-poset series and of the
  reduced-poset series is asserted to be an integer; the denominators involved go up to 17!.
- **Labeled versions.** Formulas (1) and (2) with exponential weights reproduce A000798 and A006905
  (n ≤ 10). The p_1^n coefficients of Z_Pr and Z_R equal the labeled prime and reduced counts computed
  from A001035 alone (n ≤ 17).
- **An independent a(8).** A046905(8) = 97463520 was added to the OEIS on Nov 7 2026 by another
  contributor; it agrees.
- **What a(17) rests on.** A001930(17) and the other n = 17 values use A000112(17) = 378180373015320248
  (in the OEIS since Nov 6 2026) and the n = 16 fixed-point values. The n = 16 values are tested by
  A000112(16) and A001930(16); the n = 17 values have the checks listed on the A000112 page
  (divisibility by 17!, statistics computed two ways) but no second published source. The n = 16 terms
  of A091073, A296105, A399304 and A079265 do not depend on n = 17 data at all.

## 6. Files

- `topologies.py`: formulas (1) and (2), the triangle by number of classes, labeled checks.
- `prime.py`: symmetric functions in the power-sum basis, plethysm, plethystic inverse; prime posets.
- `family.py`: all 30 entries, the comparison with the OEIS terms, b-files.
- `family.md`, `family.json`: every entry's terms, new ones marked.
- `poset_cycle_index.json`: Z_P to degree 17. `zprime.json`: the cycle index of prime posets.
- `results.json`: A001930, A091073 and their triangles by number of classes.
- `b*.txt`: tables for the entries with new terms and more than 15 terms.

## Notes

- A001930(18) needs A000112(18), which is at the scale of Ayala's computation of A001035(19).
- Not derivable from Z_P: self-dual posets, posets with trivial automorphism group (A399107), counts by
  height or by maximal elements, and A079263 (posets with a marked up-set).
- The triangle T(n, k) of topologies by number of T_0-classes, the non-T_0 topologies
  (0, 0, 1, 4, 17, 76, 400, 2490, …) and the pointed posets (1, 3, 11, 47, 243, 1533, 12038, …) are not in
  the OEIS.
