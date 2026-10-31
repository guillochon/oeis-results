# A005226 — a(19)..a(21) from the inverse Euler transform of A000638

**Result.** A005226(n) counts atomic species of degree n (equivalently, connected permutation groups of degree n:
those that are not a direct product of smaller permutation groups). Its DATA stops at a(18). The Euler transform
relation with A000638 (molecular species = permutation groups of degree n) gives

| n | A005226(n) |
|---|---|
| 19 | 6772030 |
| 20 | 76121103 |
| 21 | 110719742 |

No new enumeration was done here. a(19) and a(20) use A000638(19), A000638(20), which are already in the OEIS;
a(21) uses A000638(21) = A000638(20) + A000637(21) = 245393739 (see
[../A000638](https://github.com/guillochon/oeis-results/tree/main/A000638)). All three inputs come from Jeffrey
Ketchersid's computation of the conjugacy classes of subgroups of S_19, S_20 and S_21, and the new terms are
exactly as reliable as those.

## Why A000638 is the Euler transform of A005226

A *molecular* species is one that cannot be written as a sum of two nonzero species; the molecular species of
degree n correspond to the conjugacy classes of subgroups of S_n, so there are A000638(n) of them. An *atomic*
species is a molecular species of positive degree that is not a product of two species of positive degree.
Yeh's unique factorization theorem says that every molecular species is a product of atomic species in exactly
one way, up to the order of the factors (Yeh 1986; Labelle–Yeh 1989). Degrees add under products, so the
molecular species of degree n are exactly the multisets of atomic species whose degrees sum to n, and

    1 + Σ_{n≥1} A000638(n) x^n = Π_{k≥1} (1 − x^k)^(−A005226(k)),

which is the Euler transform. In group terms, a permutation group is uniquely a direct product (on disjoint
point sets) of connected permutation groups. Both entries state the relation (Christian G. Bower, 2006).

The transform is inverted term by term: with c(n) = Σ_{d|n} d·A005226(d), one has n·A000638(n) =
Σ_{k=1..n} c(k)·A000638(n−k), so each c(n), and then each A005226(n), follows from A000638(0..n). Each step is an
exact division, which `derive.py` asserts.

## Checks

- The inversion reproduces all 19 known terms A005226(0..18) from A000638(0..18).
- The Euler transform of the new A005226(0..21) gives back A000638(0..21) exactly.
- All divisions are exact.
- On the live OEIS (Oct 30 2026), A000637 still ends with a(21) = 140399143, A000638 still ends at a(20), and
  A005226 still ends at a(18).

## Reproduce

```
python derive.py
```

Pure Python, instant. Its output is in `derive_output.txt`.

## Sources

- Jeffrey Ketchersid, [S1-S21 conjugacy classes and verification code](https://github.com/jeff87654/Lifting),
  GitHub (the source of A000638(19..21)).
- Y. N. Yeh, *The calculus of virtual species and K-species*, in Combinatoire énumérative, Lecture Notes in
  Math. 1234, Springer (1986), 351-369.
- J. Labelle and Y. N. Yeh, [The relation between Burnside rings and combinatorial
  species](https://doi.org/10.1016/0097-3165(89)90019-8), J. Combin. Theory A 50 (1989), 269-284.

## Working log

- Oct 2026: found by the transform search of the OEIS-sniper project (the inverse Euler transform of A000638
  matches A005226 with terms to spare).
