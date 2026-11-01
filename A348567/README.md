# A348567: fully symmetric recurrent sandpiles on the n × n grid

a(n) counts the recurrent configurations of the abelian sandpile on the n × n grid (every cell of degree 4, edges
leaving the grid go to the sink) that are invariant under all eight symmetries of the square. This page gives a
closed product formula for a(n), proves the formula a(n)² = A348566(n, n) · 2^A030451(n) conjectured in the entry,
and lists a(0)–a(100) (`b348567.txt`; a(100) has 635 digits). James Guillochon, Nov 1 2026.

Plain-English version: [layman.md](layman.md).

## Result

**Theorem.** For n ≥ 0,

    a(n) = Π_{1 ≤ h ≤ k ≤ n, h and k odd} ( 4 − 2cos(hπ/(n+1)) − 2cos(kπ/(n+1)) ).

**Corollary (conjecture in the entry, Zabolotskii, Nov 2025).** a(n)² = A348566(n, n) · 2^A030451(n).

### Proof

1. *Symmetric recurrents = a cokernel.* Florescu, Morar, Perkinson, Salter, Xu, *Sandpiles and Dominos*
   (EJC 22 (2015) P1.66, arXiv:1406.0100), Prop. 2.10: for any finite group G of graph automorphisms fixing the
   sink, the G-symmetric recurrent configurations form a group isomorphic to Z^O / image(Δ_G), where O is the set
   of orbits of non-sink vertices and Δ_G is the symmetrized reduced Laplacian. So a(n) = |det Δ_G| for
   G = D4 acting on the sandpile grid graph.
2. *Δ_G is Δ restricted to the invariant subspace.* By its definition (2.2), Δ_G is the matrix of the reduced
   Laplacian Δ acting on G-invariant functions, in the basis of orbit indicator vectors. Its determinant is
   basis-free: det(Δ on V^G).
3. *Eigenbasis of V^G.* Δ = T ⊗ I + I ⊗ T with T = tridiag(−1, 2, −1) of size n, eigenvectors
   s_h(i) = sin(hiπ/(n+1)), eigenvalues 2 − 2cos(hπ/(n+1)). Since s_h(n+1−i) = (−1)^(h+1) s_h(i), the
   functions s_h ⊗ s_k invariant under both axis reflections are exactly those with h, k odd; the diagonal
   reflection swaps s_h ⊗ s_k ↔ s_k ⊗ s_h. Hence V^G has basis f_hk = s_h⊗s_k + s_k⊗s_h (h ≤ k odd), each an
   eigenvector with eigenvalue 4 − 2cos(hπ/(n+1)) − 2cos(kπ/(n+1)). Dimension check: m(m+1)/2 with m = ⌈n/2⌉,
   which is the number of D4 orbits of cells. Taking the product gives the theorem.
4. *Corollary.* For Z/2 × Z/2 the same argument gives b(n) = A348566(n, n) = Π over all ordered odd pairs
   (h, k) (this is Florescu et al. Thm 4.2/4.10 in another form). So a(n)² = b(n) · Π_{h odd ≤ n} D_h with
   D_h = 4 − 4cos(hπ/(n+1)) = 8 sin²(hπ/(2N)), N = n+1. Pairing h with 2N − h, the odd h ≤ n together with
   their mirrors (and h = N when N is odd, where sin = 1) are exactly the N odd numbers below 2N, and
   Π_{j=0}^{N−1} sin((2j+1)π/(2N)) = 2^(1−N). So Π_{h odd ≤ n} D_h = 8^⌈n/2⌉ · 2^(−n) = 2^(3⌈n/2⌉−n),
   which is 2^m for n = 2m and 2^(m+2) for n = 2m+1, i.e. 2^A030451(n). ∎

Also A256046(n) | a(n) | b(n) (entry's FORMULA) is consistent: a(n) | b(n) since b/a = Π_{h<k} factors, an
algebraic integer that is rational.

**Asymptotics.** log a(n) ~ (n²/8) · (4G/π) = G n²/(2π), G = Catalan's constant (the spanning-tree entropy of
Z², over the 1/8 of the modes that survive). At n = 100: log a = 1462.4 vs G n²/(2π) = 1457.8.

## Verification

| check | range |
|---|---|
| `brute.py`: enumerate D4-symmetric configurations, Dhar burning test | n ≤ 7, matches DATA |
| `folded.py`: exact integer det (Bareiss) of the folded Laplacian | n ≤ 40, equals the product |
| `product.py`: product at precision ≥ digits + 60, nearest integer within 1e−30 | n ≤ 100 |
| sqrt identity a² = b·2^A030451 | n ≤ 100 |
| b(n) from the product vs A348566 diagonal in DATA | n ≤ 7 |
| a(2m) = 2^m · A065072(m) against the A065072 b-file (independent, domino tilings) | m ≤ 50 (n ≤ 100) |

All 11 terms in DATA reproduce.

**Even n and domino tilings.** a(2m) = 2^m · A065072(m), where 2^m · A065072(m)² = A004003(m) is the number of
domino tilings of the 2m × 2m square (Pachter 1997, Lemma 2). Checked against the A065072 b-file (Heinz,
n = 0..50) for all m ≤ 50, i.e. every even term to a(100), an independent source. With Florescu et al. Thm 4.2
(b(2m) = A004003(m)), the even case of the sqrt conjecture is Pachter's factorization in sandpile form; the
odd case has no such analogue in the literature we found, and the product proof covers both at once.

## Prior art

- L. Florescu, *New connections between the Abelian Sandpile Model and domino tilings*, Reed College bachelor's
  thesis, 2011. Chapter 4, Thm 4.0.29: D4-symmetric
  recurrents on the 2m × 2m grid equal *weighted* domino tilings of Pachter's "Ciucu graph" H_{m−1} (red dominoes
  weight 2), via the folded Laplacian (block tridiagonal, blocks of decreasing size) and the generalized
  Temperley bijection (KPW). Proof is a sketch ("same operations as … the other proofs"); even sizes only; no
  product formula, no numbers, no odd case, no relation like a² = b·2^k. Credit it for the tiling reduction.
- Florescu, Morar, Perkinson, Salter, Xu 2015 (EJC): Prop. 2.10 (used above); §6 points to the thesis.
- Durgin, *Abelian Sandpile Model on Symmetric Graphs*, Harvey Mudd 2009: not read.
- New here: the closed product for all n, the proof of Zabolotskii's conjecture, the relation to A065072, and
  a(11)–a(100).

## Programs and how to reproduce

PARI/GP, all terms (checked against `b348567.txt` for n = 0..100):

    a(n) = localprec(n^2\8 + 40); my(m=(n+1)\2, c=vector(m, i, 2*cos((2*i-1)*Pi/(n+1)))); round(prod(i=1, m, prod(j=i, m, 4-c[i]-c[j])));

The working precision n²/8 + 40 digits exceeds the size of a(n) (about 0.064 n² digits) by a wide margin.

Python (needs mpmath), run from this folder:

    python brute.py 7                    # brute force, n = 1..7 (about a minute)
    python folded.py 40                  # exact integer determinant of the folded Laplacian
    python product.py 100 b348567.txt    # product formula, checks folded.py for n <= 40 and the sqrt
                                         # identity for n <= 100, writes the b-file

- `brute.py`: enumerates the 4^(orbits) symmetric configurations and tests each with Dhar's burning algorithm.
- `folded.py`: builds the D4-folded reduced Laplacian (rows and columns are cell orbits) and takes its
  determinant with fraction-free Bareiss elimination, so no floating point is involved.
- `product.py`: the product at precision (number of digits + 60), rounded only when within 1e-30 of an integer;
  also computes b(n) = A348566(n, n) and asserts a(n)² = b(n) · 2^A030451(n).
- `b348567.txt`: n = 0..100.
