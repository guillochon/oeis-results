# A000080: minimal triangle graphs (and A048781)

**A000080(n)** counts, up to isomorphism, the sets T of triples on n points whose graph G(T) (points
joined when they share a triple) is connected and *minimal*: dropping any triple disconnects it. In
hypergraph language these are the minimal connected 3-uniform hypergraphs on n vertices. **A048781(n)**
adds up their triples. Bowen (Math. Comp. 21, 1967) generated them by hand-checked computer search.
Martin Fuller (2015) took both sequences to n = 25 with an isomorph-free generator (nauty, about 24 h for
n ≤ 25).

**Result:** an exact counting formula from a structure theorem plus Pólya / species theory. There is no
enumeration and no isomorphism testing. All terms to **n = 80** come out in about 2 minutes per prime (n = 100
is running). The formula reproduces all 23 known terms of both sequences. New terms start at

    A000080(26) = 4987668376,  A000080(27) = 16341943122,  ...,
    A000080(80) = 2617646724739221521059589519172655683909550

(results/a000080_n80.txt). Pictures: https://guillochon.github.io/minimal-triangle-graphs/

## 1. Structure theorem

Let I be the bipartite incidence graph: one node per point, one node per triple, with an edge when the
point lies on the triple.

1. T is connected (G(T) connected) iff I is connected. T minus a triple e is connected iff I − e is
   connected, because the points stay. So **T is minimal iff every triple node is a cut vertex of I.**
2. A node of degree 3 is a cut vertex iff at least one of its edges is a bridge. If it is a cut vertex, its 3
   edges fall into at least 2 blocks, so some block holds exactly one of them, and that block is a bridge.
   The converse is clear.
3. So in a 2-connected block of I (one that is not a bridge), every triple node has degree exactly 2: it
   has degree at least 2 inside the block and at least one bridge outside. Such a block is the subdivision
   of a **loopless multigraph on the points that is 2-connected with ≥ 2 edges** (a digon of parallel
   edges counts). Its edges are triples, and the third corner of each triple hangs off the block through a
   bridge. Two arms of one block cannot share a point, and an arm cannot land back on the block, because
   either would close a cycle and so belong to the block.
4. Conversely every connected I with these blocks is the incidence graph of a minimal system. Repeated
   triples cannot occur, since they would make a K(2,3) block.

Hence **a minimal triangle graph is a tree of points and units glued at points**. A unit is either a
**lone triple** (species E3) or a **block unit**: a 2-connected loopless multigraph B with each edge carrying its
own pendant arm point, i.e. the two-sort species B(X, Y) with Y := X.

## 2. Counting

Everything is cycle-index (species) algebra evaluated by plethystic substitution into one-variable power
series. Write F(x^k) for the substitution p_k → F(x^k).

* **All multigraphs, edges labelled** (sort X vertices, sort Y edges, no loops): the cycle index is
  explicit, Z_G = Σ_λ p_λ / z_λ · exp(Σ_d o_d(λ) L_d), where o_d(λ) is the number of orbits of length d of a
  permutation of type λ on unordered pairs, and evaluated at q_k = W(x^k) the factor is
  L = Σ_j W(x^j)/j taken at x^d. That gives G̃(S, W), and from the same sum p1∂/∂p1 G and ∂/∂q1 G (weights
  m1 and C(m1,2) + m2).
* **Connected:** C = PLog G. **Blocks**, with B* = E2(X)·Y + B (a single edge or a 2-connected block):
  C• = X·E(B*_X(C•, Y)), C_Y = B*_Y(C•, Y), and C = C• + B*(C•, Y) − C•·B*_X(C•, Y).
  Rather than inverting these as species, choose S with C•(S, W) = R. Then
  B*_X(R, W) = PLog(R/S), B*_Y(R, W) = (∂_{q1}G/G)(S, W), and B*(R, W) = C(S, W) − R + R·B*_X(R, W).
* **Tree of units:** R (point-rooted) = x·PExp(U'(R)) with U = E3 + B(X, X), so
  U' = E2 + B_X + B_Y = B*_X − R·W + B*_Y with W = R. Unrooted, by the dissymmetry theorem:
  **H = R + U(R) − R·U'(R)**, with U(R) = E3(R) + B*(R, R) − E2(R)·R.
* **A048781:** weight every triple by y (W = yR, lone triple y·E3). Work over dual numbers y = 1 + ε,
  where plethysm sends ε → kε. The ε-part of H counts triples.

R and S are found order by order. At degree k one evaluation of G̃ truncated at x^k gives C•, hence
R_k − S_k, hence U' to degree k − 1, hence R_k. The cost is a sum over all partitions λ of m ≤ k
(Σ p(m) nodes, each O((k − m)²)), done for every k. Arithmetic is modulo 62-bit primes. The terms come from
CRT over 3 or 4 primes, and the remaining primes check the result.

## 3. Checks

* All 23 known terms of A000080 and of A048781 (n = 3..25) are reproduced exactly.
* Block counts: check/blocks.py splits all graphs with n ≤ 15 points (lists from Fuller's program) into
  incidence-graph blocks. It finds 1, 1, 1, 2, 2, 4, 6, 12, 20, 42, 84, 187, 424 distinct blocks with 3..15
  points, matching B*(x, x) from the multigraph series.
* species.py (exact rationals, independent code) agrees with the Rust engine for n ≤ 16.
* CRT: every term is reconstructed from 3 primes and agrees with the other 2 (n ≤ 80).
* The random examples on the web page are built from the theorem in reverse. 11,400 of them (n = 3..40)
  were tested to be connected and minimal on exactly n points.
* Growth: a(n+1)/a(n) rises smoothly from 3.25 (n = 25) to 5.03 (n = 79). The blocks make growth slightly
  faster than exponential, because many triples can share one pair of points.

## 4. Files

* `rust/` – engine. `mtg N PRIME THREADS` prints A000080 mod a prime. `mtg2` (dual numbers) prints A000080
  and A048781. It sets below-normal priority and opts out of Windows EcoQoS.
* `run.sh N` – all 5 primes in turn. `crt.py TAG K` – CRT with K primes, the rest as checks.
* `species.py` – exact-rational Python prototype of the same formula (slow, n ≲ 20).
* `check/blocks.py`, `check/a000080_fuller.c` – the block check and Martin Fuller's program (from the
  OEIS entry; the copy in `check/` has its argument parsing patched for nauty 2.8).
* `examples/small.json` – every graph for n ≤ 9 with block labels, for the web page (`viz/A000080`).
* `results/` – per-prime outputs, logs, CRT tables.

## 5. Timing (i5-12600K, 6 threads, one prime)

n = 25: 0.2 s · n = 50: 0.6 s · n = 60: 3.6 s · n = 70: 19 s · n = 80: ~100 s. Each +10 in n costs about 5×.
