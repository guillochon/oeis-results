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
Unate(n) = 1 + 2^n Σ_k C(n,k) Σ_{G on n−k vertices} 2^(−iso(G)) is the exact number of unate Krom
functions. The sum counts the satisfiable ones (a unate function is counted once per valid sign
vector, and the valid sign vectors of a function are 2^(number of variables it does not depend
on)); the 1 is the constant false function, which the sum misses. By enumeration Unate(n) = 4, 14,
96 for n = 1, 2, 3. (An earlier version omitted the 1; the ratios below do not change at this
precision.)

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
are the n single flips, with Fix = G(n−1) exactly (the function must not depend on the flipped variable),
and the 2·C(n,2) signed transpositions of order 2, each with Fix ≈ 3·H(n−1) (Lemma 3.2). From
the exact data at n = 10 the flips carry 3.3 % of the mass and a transposition has
Fix ≈ 3.2·G(9), in line with this picture.

## 3. What is proved

Throughout, for an elementary function (skew poset) P and g ∈ B_n, "g-invariant" means the
relation is preserved: R(gl, gm) ⇔ R(l, m). Write g = (π, ε), π the permutation of the variables.
Let g move m variables (a variable is moved if it is permuted or complemented), let c be the number
of cycles of π of length ≥ 2, and write j = m − c. A variable fixed by π but complemented counts in
m and not in c. Then c ≤ j (a cycle of length L ≥ 2 adds L − 1 ≥ 1 to j, a complemented fixed
variable adds 1), m ≤ 2j, and j = 0 only for g = 1.

**Lemma 3.1 (flips).** If g complements a variable v and fixes everything else, the g-invariant
Krom functions are exactly those that do not depend on v. Hence Fix(g) = G(n−1) for all Krom
functions (both constant functions included) and Fix_H(g) = H(n−1) for skew posets.
*Proof.* Invariance under complementing v is independence from v. Restricting to x_v = 0 is a
bijection from these functions to the Krom functions on the other n − 1 variables (restriction keeps
a 2-CNF a 2-CNF), with inverse "add v as an unused variable". For n = 1, 2, 3 direct enumeration
gives 2, 4, 16 = G(0), G(1), G(2). For a satisfiable function, independence means v is isolated
(unforced, in no block, in no implication): if v is forced its value flips; if v ≡ ±w the flip
gives v ≡ ∓w; if v < l (or ¬v < l) for a literal l of another variable, the flip gives ¬v < l as
well, so l is forced. (An earlier version had G(n−1) − 1: it left out the constant false function,
which every flip fixes.) □

**Lemma 3.1′ (complemented fixed variables).** Let v be fixed by π and complemented by ε, and let
P be a g-invariant skew poset. Then P has no relation between a literal of v and a literal of any
other variable w fixed by π. In particular, if π = 1 every complemented variable is isolated.
*Proof.* Let v < l, l a literal of w (the case l < v is the contrapositive of ¬l < ¬v). If w is
not complemented, the image is ¬v < l, so l is forced; if it is, the image is ¬v < ¬l, i.e.
l < v, so v ≡ l. A skew poset has neither. (Checked by enumeration for all sign flips and all
skew posets with n ≤ 4 in the case π = 1.) □

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
*Proof.* The representatives are one variable per cycle of π of length ≥ 2 and every variable
fixed by π and not complemented: c + (n − m) = n − j variables. The j others are the rest of the
long cycles and the complemented fixed variables. A g-invariant P is determined by (i) its
restriction to the representatives (a skew poset on n − j points) and (ii) for each representative
v of a long cycle and each of the j non-representatives u, the relations between the literals of v
and of u (at most 8 bits). Indeed, for literals a of g^s(v) and b of g^t(v') (v, v' long-cycle
representatives), R(a, b) = R(g^(−s) a, g^(−s) b) and g^(−s) b is a literal of g^(t−s)(v'), a
representative or not: the value is in (i) or (ii). If a is a literal of a variable w fixed by π,
shift by g^(−t) instead: R(a, b) = R(g^(−t) a, g^(−t) b), a pair (w, v'), which is in (i) when w is
not complemented and in (ii) when it is. Two variables fixed by π, at least one complemented, have
no relations at all (Lemma 3.1′). □

**Proposition 3.4 (small support).** Let U(n) denote the set of g ≠ 1 with j(g) ≤ n/20. Then
Σ_{g ∈ U(n)} Fix_H(g) = o(H(n)).
*Proof.* By Allen/Ilinca–Kahn there is n_0 with H(k) ≤ 2 · 2^(C(k+1,2)) for k ≥ n_0, and
H(k) ≤ H(n_0) for k < n_0; by the unate construction (graphs with no isolated vertex, all
2^k sign vectors) H(k) ≥ 2^(C(k+1,2)) (1 − k 2^(1−k)). With Lemma 3.3 and c ≤ j,
```
Fix_H(g)/H(n) ≤ 4 · 2^(C(n−j+1,2) − C(n+1,2) + 8 j^2) = 4 · 2^(−jn + j(j−1)/2 + 8j^2) ≤ 4 · 2^(−jn + 8.5 j^2).
```
For j ≤ n/20 the exponent is ≤ −0.57 jn. Since m ≤ 2j, the number of g with a given j is at most
(2j + 1) C(n, 2j) (2j)! 2^(2j) ≤ (2j + 1)(2n)^(2j) 2^(2j) (choose at most 2j moved variables, their
images and their signs), so
Σ_{g∈U(n)} Fix_H(g)/H(n) ≤ Σ_{j≥1} 4(2j + 1) · 2^(2j log2(4n) − 0.57 jn) → 0. □

**Proposition 3.5 (large support, unate part).** For g with j(g) > n/20, the number of
g-invariant *unate* skew posets is ≤ 2^(C(n+1,2) − n^2/160 + O(n)), and summing over all such
g gives o(H(n)) even after multiplying by |B_n|.
*Proof.* A unate skew poset P is P(s, G) for a sign vector s and a graph G (edges = clauses in
the signed variables; isolated vertices of G are the isolated points of P). Write g = (π, ε)
with π the underlying permutation, moving m' variables, and let f be the number of variables
fixed by π but complemented by ε; m = m' + f > n/20. If P is g-invariant then
P(s, G) = P(g·s, πG); the two sign vectors agree on the variables P depends on, so G = πG. Also
each of the f complemented fixed variables v is an isolated vertex of G. A unate function that
depends on v is monotone in v in one direction only. Composing with g replaces x_v by ¬x_v and
permutes the other coordinates among themselves, so it reverses that direction; invariance then
makes the function monotone both ways in v, i.e. independent of v. (An earlier version applied
Lemma 3.1 to g^(ord π) instead; that fails when ord π is even: g swapping variables 1, 2 and
complementing 3 has g² = 1.)
Case f ≥ n/40: G is a graph on the other n − f vertices, so there are at most
2^(n + C(n−f,2)) ≤ 2^(C(n+1,2) − n^2/160 + n) such P.
Case m' ≥ n/40: the pairs fixed by π are the C(n−m',2) pairs of fixed vertices and at most
m'/2 swapped pairs, every other pair lies in an orbit of size ≥ 2, so π has at most
(C(n,2) + C(n−m',2))/2 + m'/4 ≤ C(n,2) − m'(n−2)/4 ≤ C(n,2) − n^2/160 + n orbits on pairs, and
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

*Nov 2026:* corrected after an external audit: Fix of a single flip is G(n−1), not G(n−1) − 1
(Lemma 3.1); the unate count includes the constant false function; complemented fixed variables
are now explicitly not cycles in the definition of c (so c ≤ j and j = 0 only for g = 1), with
Lemma 3.1′ added for the Lemma 3.3 injection; and the isolation step of Proposition 3.5 now
uses uniqueness of the sign of an essential variable. Conjectures A and B and the table are
unchanged; no finite count was affected.
