# A154639–A154645: Mallows's 3D Apollonian sphere packing

Colin Mallows (J. Integer Sequences 12 (2009), article 09.2.1, "Growing Apollonian packings") grew the
three-dimensional Apollonian packing from five mutually tangent spheres and tabulated, for generations 0–6, the
number of spheres added in each generation and the sums of their bends and squared bends, counted in two ways
("strategy (a)" and "strategy (b)"). These became A154639–A154645. Mallows wrote that he could make little
progress and that all four counting problems remained open.

This folder contains:

- **Proved generating functions** for A154639, A154640, A154641 and A154642.
- **Corrections** to A154640, A154642, A154643, A154644 and A154645, all in the terms Mallows computed.
- **Terms to generation 16** for the strategy (b) entries A154643–A154645.
- **A conjectured rational generating function** for A154643, plus a new companion sequence: quintuples per
  generation.

Pictures: [Apollonian Sphere Generations](https://guillochon.github.io/apollonian-spheres/), a three.js viewer of the packing through generation 12.

## Results

| Entry | What a(n) is | Change | Formula |
|---|---|---|---|
| A154639 | Mallows's "reduced words" of length n | a(7)– | g.f. (1+2x−x²+2x³)/(1−3x−6x²+12x³), **proved** |
| A154640 | spheres added in generation n, strategy (a) | **a(6) = 2340, not 3330**; a(7)– | g.f. 5(1−x)(1−x−3x²)/(1−3x−3x²+6x³) = 5·A154638(x)/A162740(x), **proved** |
| A154641 | bend sum of generation n / bend sum of generation 0, (a) | a(7)– | g.f. (1−x)²(1+x)(1−x−3x²)/(1−5x−8x²+21x³+9x⁴−6x⁵), **proved** |
| A154642 | squared-bend sum / initial, (a) | **a(1) = 7, a(2) = 76, not 6, 77**; a(7)– | g.f. (1−x)³(1+x)(1−x−3x²)/(1−10x−7x²+85x³−24x⁴−27x⁵+6x⁶), **proved** |
| A154643 | spheres added in generation n, strategy (b) | **a(5) = 1350, not 1290**; a(6)–a(16) | g.f. 5(1−x)(1−3x)(1+2x+3x²)/(1−3x−3x²−9x³+18x⁴), conjectured (fits a(0)–a(16)) |
| A154644 | bend sum / initial, (b) | **a(5) = 8706, not 8454**; a(6)–a(16) | none found up to order 6 |
| A154645 | squared-bend sum / initial, (b) | **a(1) = 7, a(2) = 76, a(5) = 220338, not 6, 77, 216366**; a(6)–a(16) | none found up to order 6 |
| (new) | Descartes quintuples of generation n, strategy (b) | 1, 5, 30, 120, 480, 1980, …, n ≤ 15 | g.f. (1+2x+12x²+6x³+3x⁴)/(1−3x−3x²−9x³+18x⁴), conjectured |

Terms: `bfiles/b1546xx.txt` (n = 0..1000) for A154639–A154642; `results/stratb_sym_16.json` for strategy (b).

## 1. The setting

Mallows starts from three unit spheres sandwiched between two parallel planes at distance 2 (bends 1, 1, 1, 0, 0).
In 3 dimensions the bends b of a Descartes configuration (five mutually tangent spheres) satisfy
3·Σb² = (Σb)², so any one sphere can be swapped for the other sphere tangent to the remaining four. This is
**Descartes reflection** R_i, a 5×5 integer matrix: the identity except row i, which is (1, …, 1) with −1 on the
diagonal. A configuration after reflections R_{i1}, …, R_{ik} is M·C with M = R_{i1}⋯R_{ik} and C the
starting configuration.

**Spheres are integer vectors.** Let C be the starting configuration in augmented curvature-center coordinates
(b̄, b, b·x, b·y, b·z) (Graham–Lagarias–Mallows–Wilks–Yan), an invertible 5×5 matrix. Then every sphere of the
packing is v·C for a unique integer row vector v, namely a row of some M. Its bend is v₀+v₁+v₂. Two spheres v, w
are tangent exactly when v·G₀·wᵀ = −1, with G₀ = 2I − J the Lorentz Gram matrix of the start (each R_i preserves
it). There is no floating point anywhere in the computation; geometry (centers, radii) is used only for pictures.

The reflections satisfy R_i² = I and (R_iR_j)³ = I. The group is the Coxeter group W with five generators and all
m_ij = 3. Its growth series (distinct elements by length) is A154638, (1+x)(1+x+x²)/(1−3x−3x²+6x³).

**Generations.** Generation 0 is the starting five. The two strategies differ only in which quintuples may be
reflected.

- **Strategy (a)**: the quintuples of generation k are the configurations M_w·C for the group elements w of length
  k. A sphere belongs to the generation in which it first appears.
- **Strategy (b)**: a quintuple (any five mutually tangent spheres) belongs to the generation of its youngest
  sphere. The spheres of generation g+1 are the new spheres obtained by reflecting the quintuples of generation g.
  Already in generation 2, ten pairs of new spheres touch, closing up ten "extra" quintuples (three old spheres and
  two new). Strategy (b) reflects these at once.

This is how Mallows describes the strategies in the paper and in the entries. For strategy (b) he did not compute
the quintuples geometrically. He used a word rule instead: a word of length k with j factors of type 121 is an
"extra" word of generation k − j. That rule agrees with the geometric definition through generation 4 and
differs from generation 5 on (§4).

## 2. Strategy (a): proof of the formulas

**Theorem.** A sphere's strategy (a) generation is the length of the shortest element in a coset W_ĵ·w. Here
W_ĵ is the parabolic subgroup generated by the four reflections other than R_j. Each sphere arises from exactly
one pair (j, W_ĵw).

*Proof.* Write the action on row vectors: a sphere of configuration M_w is e_jᵀM_w. Since R_k = I − e_k·α_k with
α_k = 3e_kᵀ − 1ᵀ, the map φ ↦ φR_k is φ ↦ φ − φ_k·α_k. The Cartan matrix α_k(e_l) = 3δ_kl − 1 has 2 on the
diagonal and −1 off it, and it is nonsingular (eigenvalues 3, 3, 3, 3, −2). So this is exactly the contragredient
geometric (Tits) representation of W, with fundamental chamber {φ : φ_k > 0 for all k}. The starting sphere e_jᵀ
lies in the closed chamber and is fixed by R_k for every k ≠ j. By Tits's theorem (Humphreys, *Reflection Groups
and Coxeter Groups*, §5.13), if φ, ψ lie in the closed chamber and φ·w = ψ, then φ = ψ, and the stabilizer of φ
is the parabolic subgroup generated by the simple reflections fixing it. So e_jᵀM_w = e_kᵀM_{w'} implies j = k
and w' ∈ W_ĵ·w. Conversely, the whole coset gives the same row. The configurations containing the sphere are
the M_u with u in the coset, so its first generation is the minimal length in the coset. ∎

The same argument gives a local rule. Reflecting a negative coordinate of v (v ↦ vR_k) shortens the minimal coset
representative by one, so a sphere's generation is the number of such steps needed to reach a unit vector
(`export_viz.py`, `gen_a`).

**Counting.** The minimal left coset representatives of W_ĵ have growth series W(x)/W_ĵ(x). By symmetry all five
W_ĵ have growth series A162740, (1+x)(1+x+x²)/((1−x)(1−x−3x²)). So

  A154640(x) = 5·A154638(x)/A162740(x) = 5(1−x)(1−x−3x²)/(1−3x−3x²+6x³).

(n = 0 gives 5: the identity in each of the five cosets. These are the five starting spheres, planes included,
as in Mallows's table.)

**Bend sums.** Let S(x) = Σ_w x^ℓ(w)·M_w (a matrix of power series). Inclusion–exclusion over right descent sets
gives the matrix analogue of Steinberg's formula. Only the finite parabolic subgroups J = ∅, {i}, {i,k} occur;
{i,k} gives the symmetric group S₃ of order 6, and any three generators give an infinite affine group.

  S(x)·A(x) = I,   A(x) = Σ_{J finite} (−1)^|J| x^ℓ(w_J) P_J(x)⁻¹ M_{w_J},   P_J(x) = Σ_{v ∈ W_J} x^ℓ(v) M_v.

Proof: w ↦ (u, v) with w = uv, u ∈ W^J, v ∈ W_J is a bijection with lengths adding, so S = S^J·P_J. The
elements with J ⊆ D_R(w) are u·w_J with u ∈ W^J, and Σ_J (−1)^|J| [J ⊆ D_R(w)] = [w = e]. Applying the same
formula inside W_ĵ shows that the sum over minimal left coset representatives is A_ĵ(x)·S(x). Hence

  Σ_k x^k·(bend sum of generation k) = Σ_j e_jᵀ·A_ĵ(x)·S(x)·b₀,   b₀ = (1,1,1,0,0)ᵀ,

and the squared bends use the representation X ↦ M·X·Mᵀ on symmetric matrices, with X = b₀b₀ᵀ and the (j,j)
entries. `exact_gf.py` evaluates these over ℚ(x) with sympy (no truncation) and gets exactly the three g.f.s in
the table. Dividing by s₀ = 3 and t₀ = 3 gives integer ratios for this start.

**Checks.**
- Brute-force group enumeration (`engine.py`, strategy_a) to generation 11 matches all three series.
- The truncated-series implementation (`series.py`, 45 terms) matches all three series.
- Mallows's a(0)–a(6) of A154641 all match, and so do a(3)–a(6) of A154642, including generation 6 where
  A154640 disagrees.

**The corrections in (a).**
- *A154640 a(6) = 2340, not 3330.* Mallows's own generation-6 bend ratio 20460 (A154641) and squared-bend ratio
  707076 (A154642) are reproduced exactly by these 2340 spheres. And 3330 is more than the 2730 group elements of
  length 6 (A154638, which already corrects Mallows's 780 at n = 5 to 810), while each element adds at most one
  sphere that is new to strategy (a).
- *A154642 a(1) = 7, a(2) = 76.* Mallows lists the bends himself: generation 1 is 1, 1, 1, 3, 3 (squares sum to
  21 = 7·3) and generation 2 is 1⁶3⁶4⁶6² (squares sum to 228 = 76·3). His 6 and 77 add up to the same 83 as
  7 + 76.

## 3. A154639: Mallows's "reduced words"

These are words in five letters with no factor V V and no factor V₁V₂⋯V₂m (m ≥ 2) with V₁ = V₃,
V_{2j} = V_{2j+3} for 1 ≤ j ≤ m−2, and V_{2m−2} = V_{2m} (the GLMWY III definition Mallows quotes). Such a
factor is recognized by remembering the last three letters and the parities of the partial matches in progress.
`reduced_words.py` builds this automaton (226 reachable states) and checks it against brute force for n ≤ 8.
Mallows's seven terms are reproduced. Since a(n) = uᵀTⁿv with T of size 226, the residual
a(n) − 3a(n−1) − 6a(n−2) + 12a(n−3) satisfies a linear recurrence of order ≤ 226. It vanishes for 4 ≤ n ≤ 238,
so it vanishes for all n ≥ 4. That proves the g.f. (1+2x−x²+2x³)/(1−3x−6x²+12x³).

## 4. Strategy (b): computation

`c/stratb.c` closes the generations level by level:
1. Reflect every quintuple of generation g.
2. A result whose youngest sphere is not newer than g is another generation-g quintuple (an "extra" one), and it
   is reflected too.
3. A result containing a new sphere is a generation-(g+1) quintuple.

Every quintuple ever found is kept. A counter records any newly found quintuple of generation < g, which would
mean the closure missed it at its own level. It stayed 0 in every run.

`c/stratb_sym.c` runs the same search on orbits of the S₅ symmetry that permutes the five generators. A
permutation matrix P conjugates W to itself, maps spheres to coordinate-permuted spheres and preserves strategy
(b) generations, which do not depend on the bends.
- A sphere class is its sorted coordinate vector. A quintuple class is the least P·M·P⁻¹ (stored as a 128-bit
  fingerprint).
- Orbit sizes are 120/|stabilizer|. Over the arrangements of a class's coordinates the bend v₀+v₁+v₂ sums to
  #arr·(3/5)·Σv, and its square to #arr·(3/10)·(Σv² + (Σv)²).
- This cuts memory about 100-fold. Generation 16 (2.2·10⁹ quintuples in generation 15, 6.1·10⁹ new spheres in
  generation 16) runs in minutes on one core within a few GB (18.3 million quintuple orbits at the last level).

**Checks.**
- `engine.py` (Python, unordered quintuples) and `stratb.c` agree through generation 10.
- `stratb.c` and `stratb_sym.c` agree through generation 12.
- `check_cliques.py` independently finds every set of five mutually tangent spheres among the spheres of
  generation ≤ 7. The test is u·G₀·vᵀ = −1. Sorted by youngest sphere, they give 1, 5, 30, 120, 480, 1980, 7920,
  31860, exactly the quintuple counts of the search. So the search reflects precisely the quintuples the
  definition allows.
- Generations 0–4 reproduce Mallows's quintuples (1, 5, 30, 120, 480), spheres, and both bend ratios (174, 1170;
  1278, 15978).
- At generation 5 Mallows's word rule gives 2070 quintuples and 1290 spheres. But up to generation 5 there are
  only 2616 sets of five mutually tangent spheres in the whole packing, against his total of 2706, so his
  generation-5 values cannot all be quintuples. The geometric counts are 1980 and 1350.

**Pattern.** The quintuple and sphere counts satisfy one order-4 recurrence,
a(n) = 3a(n−1) + 3a(n−2) + 9a(n−3) − 18a(n−4), from n = 5 on, checked on 7 and 8 terms beyond the fit. No proof
yet. The bend sums fit no recurrence of order ≤ 6 with 17 terms.

## 5. Independent verification of every corrected term

Every check above uses the integer sphere vectors. Because these results correct published terms, everything
was recomputed a second way that shares no code or representation with them.

**Plain floating-point geometry** (`verify_geom.py`). Spheres are real (center, radius) pairs and the plates
are planes, from Mallows's start. A Descartes reflection of member E of {A, B, C, D, E} is computed as the
inversion of E in the sphere through the tangency points of A, B, C, D. That sphere is orthogonal to all four,
so the inversion fixes them and swaps E with the other sphere tangent to all four. When both plates are among
the four, it becomes the mirror reflection in a plane. Spheres are identified by rounded bend and center, and
every bend is checked to be within 10⁻⁶ of an integer. Every quintuple the search produces is checked to be
five mutually tangent objects, by distances.

Results (`results/verify_geom_8.txt`): all seven sequences and the quintuple counts agree with the integer
engines for n ≤ 8, in particular:
- A154640(6) = 2340
- A154642(1), A154642(2) = 7, 76
- A154643(5), A154644(5), A154645(5) = 1350, 8706, 220338

**Geometric completeness of strategy (b)** (`verify_cliques_geom.py`). This finds every set of five mutually
tangent objects among the spheres of generation ≤ 7, using only distances (|c₁ − c₂| = r₁ + r₂; sphere–plate
|z| + r = 1). Sorted by youngest member, they number 1, 5, 30, 120, 480, 1980, 7920, 31860, equal to the
quintuples reflected by the search. So the search neither misses a quintuple nor invents one.

**Mallows's own rule** (`mallows_rule.py`). Mallows assigned the quintuple of a reduced word of length k with j
factors of type 121 to generation k − j, counting overlapping factors such as 12131 once. Here j is read as the
largest number of position-disjoint factors x y x over all reduced expressions of the group element. For every
element of length ≤ 9 this rule gives exactly the geometric generation, with no exceptions. So his rule, applied
to every element, gives 1980 quintuples and 1350 spheres in generation 5, the values above. His 2070 and 1290
are a slip in carrying it out, not a different definition. His table has a similar slip in the same row: 780
group elements of length 5 instead of 810, already corrected in A154638. Two other ways of applying the rule
were tried (taking j from a single expression, or the minimum over expressions); neither gives 2070.

This also gives a combinatorial description of strategy (b): the generation of a quintuple is
ℓ(w) − (largest number of disjoint braid factors xyx in a reduced expression of w). It is a candidate route to
proving the order-4 recurrence.

## 6. Reproduce

```
python engine.py 9                 # both strategies by brute force (Python)
python series.py 44                # strategy (a), 45 terms by the Coxeter formula
python exact_gf.py                 # strategy (a), exact g.f.s over Q(x) (about 1 min)
python reduced_words.py 8          # A154639 brute force; automaton in the same file
python check_cliques.py 7          # independent check of strategy (b)
python make_bfiles.py              # b-files and DATA strings
gcc -O2 -o c/stratb_sym c/stratb_sym.c && c/stratb_sym 16 80000000 40000000
gp -q check_pari.gp                # PARI programs vs b-files, n = 0..1000
python export_viz.py 8             # sphere data for the picture page
python verify_geom.py 8 8          # independent float-geometry recount, both strategies (about 1 min)
python verify_cliques_geom.py 7    # geometric 5-clique completeness check (about 1 min)
python mallows_rule.py 9           # Mallows's word rule vs the geometric generations
```

## 7. Working log

- Nov 1 2026: read Mallows (lit/mallows8.tex). Exact integer engine. Strategy (a) reproduces A154641 through
  generation 6 but gives 2340 spheres at generation 6 (3330 in A154640). Coset formula found and proved (Tits
  representation). Strategy (b) reproduces generations 0–4 and differs at 5. Tangency-clique check confirms
  through generation 7. Matrix Steinberg formula gives the bend g.f.s; exact over ℚ(x). A154639 automaton.
  S₅-orbit engine to generation 16. three.js page.
- Nov 1 2026 (later): independent verification. A float-geometry engine (inversion in dual spheres) agrees
  with all seven sequences for n ≤ 8. A geometric 5-clique search confirms the strategy (b) quintuples through
  generation 7. Mallows's own word rule, applied to every element of length ≤ 9, equals the geometric
  generation.
