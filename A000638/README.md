# A000638 — a(21) = 245393739 from A000637(21)

**Result.** A000638(n) counts permutation groups of degree n, that is, conjugacy classes of subgroups of the
symmetric group S_n (equivalently, molecular species of degree n). Its DATA stops at a(20). The companion sequence
A000637 (fixed-point-free permutation groups of degree n) already has a(21) = 140399143, Jeffrey Ketchersid's
computation of the conjugacy classes of subgroups of S_21. Since A000638 is the sequence of partial sums of
A000637,

    A000638(21) = A000638(20) + A000637(21) = 104994596 + 140399143 = 245393739.

No new enumeration was done here. The new term is Ketchersid's work, propagated to the sequence that had not yet
received it, and it is exactly as reliable as A000637(21).

## Why A000638(n) = Σ_{k=0..n} A000637(k)

Let G ≤ S_n and let M ⊆ {1..n} be the set of points that some element of G moves, |M| = k. G fixes every point
outside M and acts on M without a common fixed point, so G is a fixed-point-free subgroup of Sym(M) ≅ S_k.

Conjugation by σ ∈ S_n carries the moved set M of G to σ(M), the moved set of σGσ⁻¹. So every class of
subgroups of S_n with |M| = k contains a subgroup with M = {1..k}, and two such subgroups are conjugate in S_n
exactly when they are conjugate by an element of S_k × S_{n−k}, whose S_{n−k} factor acts trivially on them,
that is, exactly when they are conjugate in S_k. Hence the classes with |M| = k are in bijection with the
classes of fixed-point-free subgroups of S_k, counted by A000637(k), and summing over k = 0..n gives the
identity. (k = 0 is the trivial group, A000637(0) = 1; A000637(1) = 0 since no group on one point is fixed-point-free.)

The identity is in both entries: A000637 has "a(n) = A000638(n) − A000638(n−1)" (Christian G. Bower, 2006) and
A000638 says "Partial sums of A000637".

## Checks

- The partial sums of A000637 reproduce all 21 known terms A000638(0..20) (`derive.py`).
- On the live OEIS (Oct 30 2026), A000637 still ends with a(21) = 140399143 and A000638 still ends at a(20).

## Reproduce

```
python derive.py
```

Pure Python, instant. Its output is in `derive_output.txt`.

## Sources

- Jeffrey Ketchersid, [S1-S21 conjugacy classes and verification code](https://github.com/jeff87654/Lifting),
  GitHub (the source of A000637(19..21), hence of A000638(19..21)).
- D. Holt, *Enumerating subgroups of the symmetric group*, in Computational Group Theory and the Theory of
  Groups II, Contemp. Math. 511 (2010), 33-37 (earlier terms).

The same input also gives three new terms of A005226 (atomic species), see
[../A005226](https://github.com/guillochon/oeis-results/tree/main/A005226).

## Working log

- Oct 2026: found by the transform search of the OEIS-sniper project (`a(n) − a(n−1)` of this sequence matches
  A000637 with one term to spare). Ketchersid's S21 had been added to A000637 on Oct 12 2026 but not to A000638.
