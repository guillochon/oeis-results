# Asymptotics of A109457 / A109458 / A109459: what is known, what we can prove, what is open

*2026-09-27. Notation: G(n) = A109457(n) (all Krom functions), H(n) = T(n) = number of
labelled skew posets on n points, B_n = hyperoctahedral group, |B_n| = 2^n n!.*

## 1. Literature (found 2026-09-27)

Bollobás–Brightwell–Leader, *The number of 2-SAT functions*, Israel J. Math. 133 (2003) 45–60,
already use exactly our decomposition. In their words a 2-SAT function has a **spine** (the
variables taking only one value on the solution set: our forced variables), **associated**
pairs (x ⇔ y or x ⇔ ¬y on the solution set: our equivalence blocks), and is **elementary** when
both are empty. Elementary functions are in bijection with posets on the 2n literals in which
each pair x, x̄ is incomparable and x < y ⇔ ȳ < x̄ — our skew posets. So their H(n) is our
T(n), and their inequality
```
H(n) ≤ G(n) ≤ 1 + Σ_k H(n−k) C(n,k) (2n−2k+2)^k
```
is the crude form of our exact formula G(n) = 1 + Σ_k C(n,k) 2^k Σ_r S2(n−k,r) 2^(n−k−r) T(r).

Theorem (Allen 2007, conjectured by BBL; second proof Ilinca–Kahn 2012, arXiv:1005.2863):
```
G(n) = (1 + o(1)) 2^(C(n+1,2)),   equivalently   H(n) = (1 + o(1)) 2^(C(n+1,2)),
```
and almost every 2-SAT function is unate (monotone after complementing some variables). The
Ilinca–Kahn proof gives the rate: the number of skew posets whose colored graph is not
"blue-bipartite" (essentially the non-unate ones) is < 2^(−Ω(n)) · 2^(C(n+1,2)).

**These are new formula lines for the OEIS entries** (A109457 has no FORMULA at all).

### How far the data are from the limit
Unate(n) = 2^n Σ_k C(n,k) Σ_{G on n−k vertices} 2^(−iso(G)) is the exact number of unate Krom
functions (a unate function is counted once per valid sign vector, and the valid sign vectors of
a function are 2^(number of variables it does not depend on)).

| n | G(n)/Unate(n) | A109459(n)·2^n n!/G(n) | A109458(n)·n!/G(n) |
|---|---|---|---|
| 5 | 5.82 | 3.35 | 1.62 |
| 6 | 10.58 | 2.56 | 1.46 |
| 7 | 18.21 | 1.95 | 1.32 |
| 8 | 29.07 | 1.58 | 1.22 |
| 9 | 42.37 | 1.36 | 1.14 |
| 10 | 55.70 | 1.22 | 1.10 |

So at n = 10 only 1.8 % of the Krom functions are unate and the ratio is still growing: the
"almost all unate" theorem is a statement about very large n (the dominant small-n structures
are unate functions plus a few "middle" variables, whose relative weight is like n^3/2^n).
The orbit ratios in the last two columns, on the other hand, converge visibly.

## 2. The asymmetry theorem (conjecture, strongly supported)

Burnside's lemma for B_n acting on Krom functions gives the exact identity
```
2^n n! · A109459(n) = G(n) + Σ_{g ≠ 1} Fix(g),        Fix(g) = # g-invariant Krom functions.
```
The measured relative mass Σ_{g≠1} Fix(g)/G(n) is 0.951, 0.577, 0.356, 0.225 for n = 7..10,
shrinking by a factor ≈ 0.63 per step.

**Conjecture A.** Σ_{g≠1} Fix(g) = o(G(n)). Consequently
```
A109459(n) ~ G(n)/(2^n n!) ~ 2^(C(n,2))/n!,      A109458(n) ~ G(n)/n! ~ 2^(C(n+1,2))/n!.
```
**Conjecture B (second order).** Σ_{g≠1} Fix(g)/G(n) = Θ(n^2 2^(−n)). The two dominant classes
are the n single flips, with Fix = G(n−1) − 1 exactly (the flipped variable must be isolated),
and the 2·C(n,2) signed transpositions of order 2, each with Fix ≈ 3·H(n−1) (Lemma 3.2). From
the exact data at n = 10 the flips carry 3.3 % of the mass and a transposition has
Fix ≈ 3.2·G(9), in line with this picture.

## 3. What is proved

Throughout, for an elementary function (skew poset) P and g ∈ B_n, "g-invariant" means the
relation is preserved: R(gl, gm) ⇔ R(l, m). Let g move m variables (a variable is moved if it is
permuted or complemented) in c non-trivial cycles; write j = m − c.

**Lemma 3.1 (flips).** If g complements a variable v and fixes everything else, the g-invariant
Krom functions are exactly those in which v is isolated (unforced, in no block, in no
implication). Hence Fix(g) = G(n−1) − 1 for all Krom functions and Fix_H(g) = H(n−1) for skew
posets.
*Proof.* If v is forced its value flips. If v ≡ ±w for w ≠ v, the flip turns it into v ≡ ∓w,
which cannot also hold. If v < l (or ¬v < l) for a literal l of another variable, the flip gives
¬v < l as well, so l is forced, contradiction; symmetrically for l < v. Conversely an isolated v
is unaffected by the flip. □

**Lemma 3.2 (transposition bound).** Let g swap variables v, w (possibly with a complementation
of both). Then Fix_H(g) ≤ 3 H(n−1).
*Proof.* Delete w. The map P ↦ (P|_{V∖w}, relations inside {v, w}) is injective on g-invariant
skew posets, because the relations between a literal of w and a literal l of a third variable
are the g-images of those between v and l. Inside {v, w} at most three configurations are
g-invariant. For g = (v w): the relation v < w has g-image w < v, and together they make
v ≡ w; so only "none", "v < ¬w" (whose contrapositive w < ¬v is its own g-image) and "¬v < w"
remain, and these two cannot coexist (v → ¬w and ¬v → w say v ≡ ¬w). For the signed swap
v ↦ ¬w, w ↦ ¬v the same argument leaves "none", "v < w", "w < v". □

**Lemma 3.3 (injection for general g).** Fix_H(g) ≤ H(n − j) · 2^(8 c j).
*Proof.* Choose one representative variable per non-trivial cycle; with the fixed variables
this is a set of n − j variables. A g-invariant P is determined by (i) its restriction to the
representatives (a skew poset on n − j points) and (ii) for each cycle representative v and each
of the j non-representative variables u, the relations between the literals of v and the
literals of u (at most 8 bits). Indeed for literals a of g^s(v) and b of g^t(v') (v, v' cycle
representatives), R(a, b) = R(g^(−s) a, g^(−s) b) and g^(−s) b is a literal of g^(t−s)(v'),
which is a representative or a non-representative; in either case the value is in (i) or (ii).
Fixed variables are their own representatives. □

**Proposition 3.4 (small support).** Let U(n) denote the set of g ≠ 1 with j(g) ≤ n/20. Then
Σ_{g ∈ U(n)} Fix_H(g) = o(H(n)).
*Proof.* By Allen/Ilinca–Kahn there is n_0 with H(k) ≤ 2 · 2^(C(k+1,2)) for k ≥ n_0, and
H(k) ≤ H(n_0) for k < n_0; by the unate construction (graphs with no isolated vertex, all
2^k sign vectors) H(k) ≥ 2^(C(k+1,2)) (1 − k 2^(1−k)). With Lemma 3.3 and c ≤ j,
```
Fix_H(g)/H(n) ≤ 4 · 2^(C(n−j+1,2) − C(n+1,2) + 8 j^2) = 4 · 2^(−jn + j(j−1)/2 + 8j^2) ≤ 4 · 2^(−jn + 8.5 j^2).
```
For j ≤ n/20 the exponent is ≤ −0.57 jn. The number of g with a given j is at most
C(n, 2j) (2j)! 2^(2j) ≤ (2n)^(2j) 2^(2j) (choose the moved variables and their images), so
Σ_{g∈U(n)} Fix_H(g)/H(n) ≤ Σ_{j≥1} 4 · 2^(2j log2(4n) − 0.57 jn) → 0. □

**Proposition 3.5 (large support, unate part).** For g with j(g) > n/20, the number of
g-invariant *unate* skew posets is ≤ 2^(C(n+1,2) − n^2/160 + O(n)), and summing over all such
g gives o(H(n)) even after multiplying by |B_n|.
*Proof.* A unate skew poset P is P(s, G) for a sign vector s and a graph G (edges = clauses in
the signed variables; isolated vertices of G are the isolated points of P). Write g = (π, ε)
with π the underlying permutation, moving m' variables, and let f be the number of variables
fixed by π but complemented by ε; m = m' + f > n/20. If P is g-invariant then
P(s, G) = P(g·s, πG); the two sign vectors agree on the variables P depends on, so G = πG. Also
each of the f complemented fixed variables is isolated in P (Lemma 3.1 applied to g^(ord π),
which complements them and moves nothing else), so it is an isolated vertex of G.
Case f ≥ n/40: G is a graph on the other n − f vertices, so there are at most
2^(n + C(n−f,2)) ≤ 2^(C(n+1,2) − n^2/160 + n) such P.
Case m' ≥ n/40: the pairs fixed by π are the C(n−m',2) pairs of fixed vertices and at most
m'/2 swapped pairs, every other pair lies in an orbit of size ≥ 2, so π has at most
(C(n,2) + C(n−m',2))/2 + m'/4 ≤ C(n,2) − m'(n−1)/4 ≤ C(n,2) − n^2/160 + n orbits on pairs, and
there are at most 2^n · 2^(orbits) invariant (s, G).
Either way the count is ≤ 2^(C(n+1,2) − n^2/160 + O(n)), and |B_n| ≤ 2^(n log2 n + n) does not
disturb the 2^(−n^2/160). □

*Remark.* The propositions are stated for skew posets (elementary functions). Conjecture A is
about all Krom functions; the forced variables and equivalence blocks of a g-invariant function
are permuted by g, so the same injections apply layer by layer, but we have not written this
out.

## 4. What is open (the honest gap)

Conjecture A would follow from Propositions 3.4 and 3.5 plus the missing piece:

**Missing Lemma.** For g with j(g) > n/20, the number of g-invariant **non-unate** skew posets
is ≤ 2^(C(n+1,2) − ε n^2) for some ε > 0.

Ilinca–Kahn only give ≤ 2^(C(n+1,2) − Ω(n)) for the number of *all* non-unate skew posets, and
summing that over the 2^n n! group elements loses a factor 2^(n log n). The invariance has to be
used *inside* their entropy/decision-tree argument (Sections 3–5 of arXiv:1005.2863), i.e. one
has to count π-invariant odd-blue-triangle-free graphs, which the trivial bound
3^(orbits on pairs) does not control (it exceeds 2^(C(n+1,2)) as soon as j = Θ(n)). This is a
genuine research problem — a "symmetric" version of the Allen/Ilinca–Kahn theorem, analogous
to the step from Kleitman–Rothschild's enumeration of posets to Prömel's "almost all posets are
rigid". We did not attempt it.

## 5. Consequences for the OEIS entries (all with references)

- A109457: FORMULA `a(n) = (1 + o(1)) * 2^(n(n+1)/2)` [Allen 2007; Ilinca–Kahn 2012]; COMMENT:
  the "elementary" functions of Bollobás–Brightwell–Leader are our skew posets (T(n)).
- A109458: COMMENT (conjecture, supported by n ≤ 10): `a(n) ~ 2^(n(n+1)/2)/n!`.
- A109459: COMMENT (conjecture, supported by n ≤ 10): `a(n) ~ 2^(n(n-1)/2)/n!`, i.e. almost all
  Krom functions have trivial stabiliser in the hyperoctahedral group; Burnside data above.

*Oct 2026:* the OEIS edits were trimmed per the editor's preferences (terms + link only). Of the
lines above only the A109457 FORMULA (the proved asymptotic, with its references) stays in the
entries; the BBL comment and the two conjectures are in the README's "Asymptotics" section.
