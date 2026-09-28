# A114601 — Symmetric positive definite matrices with 2 on the diagonal and entries in {−1, 0, 1}

**Status (2026-09-27):** analytical track. **Formula complete: all three E constants are computed (E_8 by two independent routes), so every term follows. a(8) = 4925177553 and the new term a(9) = 63545326593 are confirmed by direct enumeration. The proof is written ([`proofs.md`](proofs.md)) and the closed-form EGF is verified. Asymptotics a(n) ~ 0.1473·n!(2e)ⁿn^(−7/8) proved (proofs.md §5). OEIS draft (`oeis_submission.md`, b-file n ≤ 200, PARI checked) and `layman.md` written. Remaining: A086215.**
**Known terms:** a(1)..a(7) = 1, 3, 23, 393, 13089, 737595, 58969079 (Max Alekseyev, 2005/2006). There is no formula in the entry.
**Sibling:** A086215 (the same with non-symmetric M; see "A086215" below).
**Literature (2026-09-27):**
- Greaves–Koolen–Munemasa–Sano–Taniguchi (2013) classify edge-signed graphs with λ_min > −2.
- Abarca–Rivera (2014) characterise positive definite quasi-Cartan matrices of types A and D.
- No enumeration was found, and the connected sequence 1, 2, 16, 304, 11008, … is not in the OEIS.

## The idea: count root bases lattice by lattice
Let M be such a matrix. Then M is the Gram matrix of n linearly independent vectors of norm 2 with pairwise inner products in {−1, 0, 1}. These vectors generate an integral lattice L of rank n spanned by norm-2 vectors. So:

- **L is a root lattice**, a direct sum of copies of A_k, D_k, E_6, E_7 and E_8.
- The n vectors are a **basis of L consisting of roots**.
- Conversely, any basis of an ADE lattice made of roots works. Distinct roots that are not opposite have inner product in {−1, 0, 1}, and the Gram matrix of a basis is positive definite.
- Two ordered bases give the same M iff an isometry of L maps one to the other. Aut(L) acts freely on bases, so
  **#(matrices with lattice L) = #(ordered root bases of L) / |Aut(L)|.**
- M is block diagonal along the connected components of its graph of nonzero entries, and a connected M has an irreducible L. So with c(k) = the number of *connected* matrices of size k:

  **a(n) = n! [xⁿ] exp( Σ_k c(k) xᵏ / k! ),   c(k) = c_A(k) + c_D(k) + c_E(k).**

### The three families
- **A_k** (roots e_i − e_j in Z^(k+1)). A set of k roots is a basis iff the corresponding edges form a spanning tree of K_(k+1), in either orientation. That gives 2ᵏ(k+1)ᵏ⁻¹ sets. With |Aut(A_k)| = 2(k+1)! for k ≥ 2:
  **c_A(k) = 2ᵏ⁻¹(k+1)ᵏ⁻²** (k ≥ 2), and c_A(1) = 1.
- **D_k** (roots ±e_i ± e_j in Zᵏ, k ≥ 4). A set of k roots is a basis iff the corresponding signed multigraph on the k coordinates is connected and unicyclic, and its cycle is *unbalanced* (the roots are independent). The claim is that it then generates all of D_k, not a sublattice of index > 1.
  - Cycles of length ≥ 3: U(k)·4ᵏ/2 sets, where U(k) = A057500(k) counts labelled connected unicyclic graphs.
  - A double edge (both e_i + e_j and e_i − e_j): (k−1)k^(k−2)·4^(k−1) sets.
  - With |Aut(D_k)| = 2ᵏk! for k ≥ 5, and 3 times that for k = 4 (triality):
  **c_D(k) = 2ᵏ⁻¹U(k) + (k−1)kᵏ⁻²2ᵏ⁻²** for k ≥ 5, and c_D(4) = 104. D_3 = A_3 is not counted twice.
- **E_6, E_7, E_8**: one constant each, in ranks 6, 7 and 8 only. **For n ≥ 9 no new E terms appear.** Every term then follows from the formula and the three constants.

## Step 1 results (2026-09-27)

| Check | Result |
|---|---|
| `lattice_check.py`: root bases of A_n counted directly (n = 2..5) | = 2ⁿ(n+1)ⁿ⁻¹ exactly |
| `lattice_check.py`: root bases of D_n counted directly (n = 3..5) | = U(n)4ⁿ/2 + (n−1)nⁿ⁻²4ⁿ⁻¹ exactly. The only other full-rank root sets have det 16: index-2 sublattices, which are not bases |
| `brute.py`: all PD matrices, n ≤ 6 | totals = OEIS a(1..6) |
| `brute.py`: connected matrices classified by determinant | A_n (det n+1), D_n (det 4) and E_6 (det 3) counts match the formulas for every n ≤ 6 |
| **E_6 term** (n = 6, det 3) | **c_E(6) = 355104** |
| `brute.py`, n = 7 (58,969,079 matrices) | total = OEIS a(7); A_7 = 2097152 = 2⁶·8⁵ and D_7 = 7597824 match the formulas |
| **E_7 term** (n = 7, det 2) | **c_E(7) = 43733504**, the same value the OEIS a(7) implies |

n = 6 in detail: 652736 connected = 76832 (A_6 = 2⁵·7⁴) + 220800 (D_6 = 2⁵·3660 + 5·6⁴·2⁴) + 355104 (E_6).

## Step 2 results: the E constants by two independent routes (2026-09-27)
The program is `rust/` (`a114601`), with exact integer arithmetic. Appending a row r to M with d = det M and A = adj M gives det' = 2d − rᵀAr, and the adjugate updates in closed form.

**Route 1, `gram n`:** enumerate every PD matrix of size n and classify the connected ones by determinant.
- It reproduces a(1..7) and every A/D/E count above.
- **n = 8** (68 s): **a(8) = 4925177553**. The connected matrices split into:
  - det 1 (E_8): **4053973888**;
  - det 4 (D_8): 301321216 = 2⁷·1436568 + 7·8⁶·2⁶ ✓;
  - det 9 (A_8): 68024448 = 2⁷·9⁶ ✓.

**Route 2, `roots E6|E7|E8`:** enumerate subsets of the positive roots of E_n (actual vectors: E_8 in its standard coordinates, E_7 and E_6 as the roots orthogonal to one root and to an A_2) and count those with Gram determinant det(E_n), i.e. bases.
- The E term is #bases × 2ⁿ × n! / |Aut(E_n)|.
- Cauchy–Binet checks the whole enumeration: the sum over all n-subsets of det(G_S) equals det(Σ r rᵀ) = hⁿ, with Coxeter number h = 12, 18, 30.
- For E_8, only subsets containing a fixed root are enumerated (W is transitive on roots), then scaled by 120/8. This is validated on E_7, where both ways give the same answer.

| Lattice | Bases (unordered sets of positive roots) | Cauchy–Binet | E term | Agrees with route 1 |
|---|---|---|---|---|
| E_6 | 798984 | 12⁶ ✓ | **355104** (remainder 0) | ✓ |
| E_7 | 196800768 | 18⁷ ✓ | **43733504** (remainder 0) | ✓ |
| E_8 | 273643237440 | 30⁸ ✓ | **4053973888** (remainder 0) | ✓ |

**Result:** c_E(6) = 355104, c_E(7) = 43733504, c_E(8) = 4053973888. The formula gives a(8) = 871203665 + 4053973888 = 4925177553, exactly the direct count.

**New terms from the formula** (`python formula.py N`):
a(8) = 4925177553, a(9) = 63545326593, a(10) = 1772133404723, a(11) = 82461147540119,
a(12) = 5016500407650521, a(13) = 369357572525622689, a(14) = 32049842608755985387, …

These rest on the ADE structure and the A/D basis counts, all proved in [`proofs.md`](proofs.md) (step 3) and verified numerically for n ≤ 8. **a(9) confirmed independently (2026-09-27):** route 1 (`a114601 gram 9 --threads 12 --ckpt logs/gram9.ckpt`, 22256 s on 12 threads, 738 checkpointed chunks) gives a(9) = 63545326593, with connected matrices only in det 10 (2560000000 = c_A(9)) and det 4 (13545271296 = c_D(9)), exactly as the formula predicts (no E-type components exist at size 9). See `logs/gram9.log`.

## Next steps
1. ~~Verify the A and D counts numerically~~ (done above).
2. ~~E_8 term~~ (done, two routes).
3. ~~**Proofs.**~~ Done in [`proofs.md`](proofs.md):
   - reduction to root bases of ADE lattices modulo Aut (Witt), with connected matrices ⇔ irreducible lattices, hence the exp formula;
   - A_k: root bases ⇔ spanning trees;
   - D_k independence ⇔ every component is unicyclic with an unbalanced cycle (signed-graph rank);
   - D_k basis ⇔ additionally connected (an unbalanced cycle gives 2e_v, connectivity gives every root; a disconnected set sits in a sublattice).

   The closed-form EGF in terms of the tree function (§5) is checked against c(k) for k ≤ 30 by `egf_check.py`.

   **Lean** (`lean/OeisLean/A114601.lean`, standard axioms only): `span_eq_D_iff` and `span_eq_A_iff` formalize the D_k and A_k basis characterizations (walk form). The classical inputs (Witt, the Aut orders, tree and unicyclic counts) and the computed E constants are not formalized.
4. **Closed form.** Assemble the EGF. U(k) has the known EGF log(1/(1−T))/2 − T/2 − T²/4, where T is the tree function. Asymptotics then follow from the singularity at 1/e.
5. **A086215.** In M, a zero pair (m_ij, m_ji) comes in 3 ways ((0,0), (1,−1), (−1,1)) and a pair with sum ±1 in 2 ways. So A086215(n) = 3^C(n,2) · Σ over matrices S of (2/3)^(edges of S), which is still multiplicative over components. The edge count of a connected piece is the number of adjacent root pairs: for A_k, the sum over vertices of C(deg, 2) in the tree; for D_k, similar. That gives degree-weighted tree functions.
6. ~~Lean formalization~~ of the A/D basis characterizations: done (see step 3).

## Files
- `brute.py`: brute force for n ≤ 7, with the classification by determinant.
- `lattice_check.py`: root-basis counts in A_n and D_n.
- `formula.py`: a(n) from the formula; also compares against the OEIS data.
- `proofs.md`: proof of the formula, and the closed-form EGF.
- `egf_check.py`: checks the closed-form EGF against c(k).
- `rust/`: `a114601 gram n` (route 1) and `a114601 roots E6|E7|E8 [--fix]` (route 2).
- `logs/`: run logs.
