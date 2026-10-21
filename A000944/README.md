# A000944: polyhedra with n vertices — a(19) to a(40) by formula

**a(19) = 1837856160118250, a(20) = 31696106245572086, a(21) = 552533560690636208**, and further terms
to a(40) (§7, `b000944.txt`). A000944(n) is the
number of combinatorial types of convex polyhedra with n vertices, i.e. 3-connected simple planar graphs
on n vertices. The terms up to a(18) were found by exhaustive generation (Brinkmann and McKay's plantri;
a(18) in 2006). Generation cannot reach n = 19 (about 2·10^15 graphs). Here the new terms are computed
instead by counting: Burnside's lemma over map automorphisms, with each kind of symmetry counted by a
decomposition of 3-connected planar maps. No graph is generated except for the small brute-force checks.

## 1. Setting

By Whitney's theorem a 3-connected planar graph has a unique embedding in the sphere up to
homeomorphism, so polyhedra are the same as **c-nets** (3-connected planar maps) up to all
homeomorphisms of the sphere, orientation-reversing ones included. Let U(v,f) be the number of
polyhedra with v vertices and f faces (A212438), so A000944(n) = Σ_f U(n,f).

Burnside's lemma, written over pairs (map, automorphism): for each geometric **type** τ of map
automorphism let

    W_τ(v,f) = Σ_M #(elements of Aut M of type τ) / |Aut M|,

summed over c-nets M with v vertices and f faces (Aut includes reflections). Then
U(v,f) = Σ_τ W_τ(v,f). A finite-order homeomorphism of the sphere is conjugate to an isometry, so the
types are: the identity; a rotation of order k with two poles, each pole a vertex, a face or (k = 2) an
edge midpoint; a reflection; or a rotatory reflection (no fixed point). Write

    U = S/2 + W_ref + W_rr,

where S = 2·(identity + rotation terms) is the number of polyhedra up to orientation-preserving
homeomorphism (A239893).

All counts below are generating functions in x (vertices) and y (faces), computed modulo primes
p < 2^26 with exact rational coefficients reduced mod p, and combined by the Chinese remainder theorem.

## 2. Identity term: rooted c-nets

W_id = R(v,f)/(4e), e = v + f − 2, where R is the number of rooted c-nets. Rooted nonseparable maps
with i+1 vertices and j+1 faces number (2i+j−2)!(2j+i−2)!/(i! j! (2i−1)! (2j−1)!) (Brown and Tutte).
Removing the root edge gives a *network*, and Tutte's decomposition reads

    D = 1 + S + P + H,   S = x(D−S)²/(1 − x(D−S)),   P = y(D−P)²/(1 − y(D−P)),
    H = D · C(xD, yD),

where D counts networks by inner vertices and inner faces, S and P are series and parallel networks,
H networks with a 3-connected core, and C(X,Y) = Σ R(v,f) X^(v−2) Y^(f−2). Inverting the last equation
gives R(v,f); summed over v + f this is A000287.

## 3. Rotations

For a rotation g of order k the weighted count over nonseparable maps is known from Liskovets and
Walsh's quotient maps (Howroyd's PARI program in A342061 gives every pole type by vertices and faces).
The Tutte decomposition tree of a nonseparable map is g-invariant. On its barycentric subdivision g
fixes a subtree, and #fixed nodes − #fixed edges = 1 turns Burnside's sum into sums over a fixed node
(a c-net, a polygon or a bond, with its virtual edges carrying networks in g-orbits) and over fixed tree
edges.
* k ≥ 3: no tree edge is fixed, the fixed node is unique, and every edge orbit has size k. With
  u = x^k, w = y^k this gives W^nonsep_τ = x^(pv) y^(pf) C_τ(uD(u,w), wD(u,w)) + (polygon or bond term),
  a polynomial-time inversion for the c-net weights C_τ.
* k = 2: fixed tree edges occur (both sides kept or both swapped). Edges on the axis carry networks
  symmetric under the half-turn that swaps the poles and the two sides, refined by the cell (vertex, edge
  or face) at their centre. The resulting system is solved degree by degree for the six pole types.

Together these give S(v,f).

## 4. Reflections

A c-net M with a reflection σ is determined by its **half**: the vertices on one side of the mirror, the
mirror vertices, and each edge crossing the mirror shortened to its endpoint. The half H is a
2-connected plane graph whose outer cycle C' carries labels V (vertex on the mirror), X (endpoint of a
crossing edge) or P (plain). Doubling H across the mirror recovers M. Using the criterion that a
2-connected plane map is 3-connected exactly when every two faces meet in nothing, a vertex or an edge,
M is 3-connected iff

  (a) the cone over C' is 3-connected;
  (b) every inner face of H contains at most one V, or exactly two V joined by a mirror edge on it;
  (c) every inner face meets every *piece* of C' (arc between consecutive V/X labels) in nothing, a
      vertex or an edge;
  (d) there are at least 3 labels V/X.

Then W_ref = ½ Σ_(H, labels) 1/|Aut(H, labels)|. To count halves, H is decomposed along its 2-cuts,
which all lie on C'. The nodes are R (c-nets whose outer face carries the virtual edges, counted by
outer-face degree), S (faces of H touching C' at least three times) and P (a chord of C'). In this
decomposition, (b) and (c) become local conditions on the labels: every pocket contains a V/X label;
a touch point between two C'-edges of the same face is labelled; and the V count is restricted per
face. These are handled by 3×3 transfer matrices, with one correction at the root. Rooting at a labelled
vertex makes Aut act freely. The number A of labels is marked by a variable α, evaluated at α = 0..N
and interpolated, so W_ref = ½ Σ_(A ≥ 3) [α^A](rooted count)/(2A).

c-nets by outer-face degree d come from the same network inversion with z marking the length of one
side of the network. Rooted nonseparable maps by root-face degree satisfy
F(z) = y² z Σ_(k≥1) x^(1−k) b(z)^k, with b(z) = x² z + (z F(1) − F(z))/(y(1−z)) (delete the root edge;
the rest is a chain of blocks).

## 5. Rotatory reflections

Let g be orientation-reversing without fixed points, of (even) order n. Either some power of g is a
reflection, or none is.

**No reflection in ⟨g⟩.** Then g² is a rotation of order n/2 (the identity if n = 2) and g swaps its
two poles. The quotient of the sphere by ⟨g⟩ is the projective plane, with one cone point of order n/2
(at a vertex, a face, or for n = 4 at the end of a semi-edge) when n ≥ 4.
1. *General maps.* Rooted maps on the projective plane, by vertices, faces and a marked vertex or face,
   satisfy a Tutte equation coupled to planar maps:
   P = xy + z²P²/y + yz(P(1) − zP)/(1 − z) (sphere), and
   Q = 2z²PQ/y + yz(Q(1) − zQ)/(1 − z) + z²·∂_z(zP) (projective plane), with z the root-face degree.
   Each rooted quotient with its cone cell lifts to φ(n) pairs (M, g), giving W^gen for every pole type
   (n = 2: Q/(4e) in half variables).
2. *Blocks.* g fixes exactly one block (no vertex is fixed). Its 2e corners form free orbits, each
   decorated by a rooted planar map Γ = P(1)/(xy) hanging in that corner; this is inverted for the
   nonseparable weights.
3. *Tutte tree.* g fixes a unique tree node (a c-net, a polygon or a bond) or, for n = 2, inverts one
   tree edge. Edge orbits have size n, except that for n = 4 with edge poles the two pole edges carry
   half-turn-symmetric networks (the k = 2 rotation series refined by centre cell). Inverting gives the
   c-net weights.

**A reflection in ⟨g⟩** (group C_mh, m = n/2 odd ≥ 3: a mirror plus a perpendicular rotation of order
m). The map is again determined by its mirror half, which now carries a rotation of order m about an
interior point. The half grammar of §4 is reused with Burnside over this rotation: the fixed centre of
the half's decomposition is either an S-node (a face meeting C' in an m-periodic pattern) or an R-node,
a c-net with a rotation whose pole is its outer face, counted by outer-face degree. Those come from the
same grammar in a "cone" mode (a disk plus an apex vertex is a c-net rooted at a vertex) by duality,
recursively in m. The rest of the boundary is one period of labels, counted by transfer matrices.

For odd v only C_mh contributes; for v = 19 that is just the 17-gonal bipyramid (W = 16/68 at f = 34).

## 6. Checks

* **Brute force.** Every polyhedron with n ≤ 11 vertices was generated with nauty (geng | planarg,
  3-connectivity in Python). For each one, every automorphism was classified by type (`brute/`),
  rotatory reflections by order and by whether ⟨g⟩ contains a reflection or else the pole type of g².
  The rooted counts, every rotation type of every order, the half-turn pole types, the reflection
  weights, the antipodal map (30 cells), rotatory reflections of even order ≥ 4 (756 cells) and C_mh
  (48 cells) all agree with the formulas cell by cell. The half-map characterization (a)–(d) was checked
  against a direct 3-connectivity test on all 138105 labelled halves whose cone has at most 8 vertices.
* **Known terms.** The same programs give A000944(4..18) exactly, all 196 entries of A212438 (rows
  4–17; these include U(19, f) for f = 12..17 and U(20, f), U(21, f) for f ≤ 17), all 196 entries of the
  A239893 b-file, A002840 (polyhedra by edges) for all 30 terms up to 35 edges, and A006445 (achiral
  polyhedra by edges) for all 33 terms up to 33 edges. The last two exercise every symmetry type at
  even vertex counts.
* **Duality.** U(v,f) = U(f,v) holds on all 324 pairs inside the 40×76 box, although the half-map counts treat
  vertices and faces differently.
* **Integrality.** S/2 + W_ref + W_rr is a sum of rationals. It comes out integral in every cell, and the
  CRT value is the same from one prime fewer (8 primes below 2^24 for the a(40) run: values below 2^192,
  checked with 2^168).
* **Independent runs.** The boxes 21×38, 25×46, 30×56 and 40×76 were computed separately (different
  truncations); every cell they share agrees. The fast arithmetic of §8 reproduces the arrays of the
  original slow code bit for bit at 21×38.
* **Growth.** 3-connected planar graphs on n vertices grow like c·γ^n·n^(−7/2), γ ≈ 21.049 (Bender,
  Gao and Wormald), and almost all are asymmetric, so a(n)/a(n−1)·(n/(n−1))^(7/2) should approach γ
  with corrections of order 1/n. It rises smoothly: 20.64 (n = 20), 20.79 (25), 20.87 (30), 20.92 (35),
  20.95 (40).

## 7. The new terms, and U(n, f) for n = 19, 20, 21

| n | a(n) |
|---|---|
| 19 | 1837856160118250 |
| 20 | 31696106245572086 |
| 21 | 552533560690636208 |
| 22 | 9725100290645366750 |
| 23 | 172665610387355491030 |
| 24 | 3089897351402967680382 |
| 25 | 55693596773314728386434 |
| 26 | 1010472821318274587081170 |
| 27 | 18444616002584669852951766 |
| 28 | 338559437591583889391837842 |
| 29 | 6246524274375883224870258320 |
| 30 | 115802469478907466726560502772 |
| 31 | 2156392465904152551633439664656 |
| 32 | 40321607164314613909499184810474 |
| 33 | 756885454198595130457093610733207 |
| 34 | 14259300397464793694906262785553986 |
| 35 | 269553677716154334083059055531235023 |
| 36 | 5111909800989602220520842699871144004 |
| 37 | 97236992305707437010775781801607307042 |
| 38 | 1854883285393075863218182649498233922286 |
| 39 | 35478950056586460924968154727569837727537 |
| 40 | 680349069787813513725946584317522497406424 |

By faces:


| f | U(19,f) | U(20,f) | U(21,f) |
|---|---|---|---|
| 12 | 70454 | 7595 | |
| 13 | 7706577 | 2599554 | 527235 |
| 14 | 274542869 | 168992630 | 74424566 |
| 15 | 4939809506 | 4686995652 | 3380569040 |
| 16 | 54271705726 | 73247405678 | 77220397213 |
| 17 | 404008232288 | 741171341224 | 1075323264149 |
| 18 | 2173270387051 | 5270785332349 | 10150757285258 |
| 19 | 8816236975026 | 27777985734886 | 69400639266068 |
| 20 | 27777985734886 | 112504197045246 | 359289661461953 |
| 21 | 69400639266068 | 359289661461953 | 1453200013970915 |
| 22 | 139494234729610 | 921617699752397 | 4697109097593154 |
| 23 | 227766462394472 | 1924083064415544 | 12335533935051363 |
| 24 | 303855809020774 | 3299394559021312 | 26641659035841264 |
| 25 | 331946749588737 | 4674086493547056 | 47728392180472226 |
| 26 | 296607096533750 | 5485782863160364 | 71332393707130745 |
| 27 | 215721437561310 | 5334049720569423 | 89221991800700926 |
| 28 | 126528648782160 | 4284218978307195 | 93471186140228468 |
| 29 | 58961388185063 | 2824592039010974 | 81884917062469287 |
| 30 | 21327809721265 | 1512860451465421 | 59732351730524213 |
| 31 | 5774235892834 | 647880752888956 | 36017376484761941 |
| 32 | 1101217837858 | 216583366005639 | 17751099378124911 |
| 33 | 131959532817 | 54460097581123 | 7032837045047755 |
| 34 | 7475907149 | 9688979580528 | 2185475988514592 |
| 35 | | 1087376999834 | 513025172527074 |
| 36 | | 57896349553 | 85534438582808 |
| 37 | | | 9027039627035 |
| 38 | | | 453382272049 |
| **sum** | **1837856160118250** | **31696106245572086** | **552533560690636208** |

All cells for v ≤ 40, f ≤ 76: `results/U_40_76.json`; sensed counts S(v,f) on that box:
`results/sensed_116.json`; the run's output with all checks: `results/run_a40.log` (the a(21) run:
`results/run_a21.log`, `U_21_38.json`, `sensed_59.json`).

## 8. Reproduce

Python 3 with numpy and networkx (nauty only for the brute-force checks). The checks read OEIS b-files from
`data/` (`b000944.txt`, `b212438.txt`, `b239893.txt`, `b000287.txt`, `b002840.txt`, `b006445.txt`).
From `formula/`:

    python sensed.py 116 7 40 76     # S(v,f) on the box v <= 40, f <= 76, 7 primes
    python run_wref.py 40 76 i       # W_ref mod the i-th prime, i = 0..7 (run in parallel)
    python run_wrr.py 40 76 i        # W_rr mod the i-th prime, i = 0..7
    python combine_full.py 40 76 8   # U, all checks, A000944(n) for n <= 40

On an i5-12600K with all 17 jobs in parallel this takes about 80 minutes (a(30): 18 minutes, a(25): 8,
a(21): 2 per prime). Series products use `fastmul.py`: exact multiplication mod p by floating-point FFT
(operands split into limbs so that every convolution sum stays below 2^44, and each result checked to be
within 0.25 of an integer), operands cropped to the bounding box of their nonzero terms, and the anchor
markers α of the grammars evaluated as one batch. The series/parallel network equations are solved in
closed form (S = xD²/(1 + xD), and P = D·yD(1)/(1 + yD(1)) with the catalytic variable), and the
nonseparable-map equation by fixed-point iteration with a warm-started Newton inverse.
    python check_brute.py 20         # rooted / rotation / half-turn vs brute force (needs brute/stats_*.json)
    python refl.py 11 20             # reflection weights vs brute force
    python rotorefl.py 20; python check_rr.py 20; python test_wrrh.py   # rotatory reflections vs brute force

Brute force (`brute/`): `gen_poly.sh 11 <dir>` (in WSL with nauty), then `symstats.py planar2c_*.g6`;
`refl_check.py 8` checks the half-map characterization.

## 9. Files

* `formula/fastmul.py`: exact FFT multiplication of truncated series mod p.
* `formula/ms.py`, `series.py`: truncated bivariate series mod p (numpy) and exact versions.
* `formula/rooted.py`, `pipeline.py`: networks, rooted c-nets, rotations of order ≥ 3, half-turns.
* `formula/t3.py`: nonseparable maps and c-nets by vertices, faces and root-face degree.
* `formula/refl.py`: the half-map grammar for reflections.
* `formula/rp2.py`, `rotorefl.py`: projective-plane maps; antipodal and even-order rotatory reflections.
* `formula/disk.py`: the generic disk grammar and its rotationally symmetric version (C_mh).
* `formula/sensed.py`, `run_wref.py`, `run_wrr.py`, `combine_full.py`: the a(21) run
  (`combine.py`: the earlier odd-n-only a(19) run).
* `formula/check_*.py`, `test_*.py`: comparisons with brute force.
* `brute/`: nauty generation, automorphism classification, half-map checks.
* `results/`: `run_a19.log`, `U_19_34.json`, `sensed_53.json` (a(19)); `run_a21.log`, `U_21_38.json`,
  `sensed_59.json` (a(20), a(21)); `run_a40.log`, `U_40_76.json`, `sensed_116.json` (to a(40)).
* `b000944.txt`: a(1)–a(40).

## 10. Open

Nothing new is needed for more terms: the cost is polynomial in n (roughly n^4.5 per prime here), plus
more primes as the terms grow.

---

## Working log

* 2026-10-20: brute force n ≤ 11; rooted c-nets; rotations k ≥ 3 (needed Howroyd's F(i,j) read
  correctly: its second term is added also when i = 0); half-turns (digon skeleton weight x²y², d | e
  for the d = 2 terms); sensed counts match A239893. Reflections: characterization brute-checked (needs at
  least 3 labels), first grammar missed P-nodes (chords; e.g. the square pyramid's diagonal mirror),
  then exact on 270 rooted coefficients and on the weights for v ≤ 11. A000944(13) reproduced, then
  a(19) (commit ea1633a).
* 2026-10-20/21: rotatory reflections. Projective-plane Tutte equation (checked against A000168 and
  A007137); antipodal and even orders (the corner decoration is one rooted map Γ = P(1)/(xy), not a
  sequence of pendants); C_mh via the rotational disk grammar (centre face counted once, not n times;
  truncation sizes rounded up). Sensed counts to v + f = 59 needed 6 primes (4 were not CRT-stable at
  26 vertices). a(20) and a(21) with every check passing.
* 2026-10-21: profiling showed the time in many small series products. FFT multiplication, batched
  anchor markers, closed-form series/parallel networks and a warm-started Newton inverse: about 25× per
  prime, bit-identical output. a(22)–a(25) (8 min), a(26)–a(30) (18 min), a(31)–a(40) (80 min). The
  first a(40) combine failed its CRT check: the 8th entry of the prime list (16777111 = 13·1290547) was
  not prime; fixed (the first six entries, used for every earlier result, were correct).
