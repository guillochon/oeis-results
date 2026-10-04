# A227867 — Lipschitz quaternions X with X^2 ≡ 1 (mod n): a closed form

**Status:** **proved and formally verified in Lean 4 + Mathlib** (`lean/OeisLean/A227867.lean`, theorem `A227867.A_closed_form`; axioms: propext, Classical.choice, Quot.sound only). Checked independently by direct count for every n ≤ 256 and against all OEIS data.

**Sequence:** A227867(n) is the number of Lipschitz quaternions X = a + bi + cj + dk (i² = j² = k² = −1) with X² ≡ 1 (mod n). Expanding X², this is the number of solutions mod n of

    a² − b² − c² − d² = 1,   2ab = 2ac = 2ad = 0.

Keyword `mult`; the entry has no formula.

**Related:** A236553 is the same count in the "split" ring i² = j² = +1 (see `../A236553/`). The two sequences agree at every odd n.

## Result
A227867 is multiplicative, with

| | odd p, k ≥ 1 | 2 | 2^k, k ≥ 2 |
|---|---|---|---|
| **A227867(p^k)** | p^(2k−1)·(p+1) + 2 | 8 | 32 |

So for odd n, A227867(n) = A236553(n) = Π_{p^k || n} (p^(2k−1)(p+1) + 2). Writing n = 2^e·m with m odd:

    A227867(n) = A227867(m) · (1, 8, 32, 32, 32, …)[e].

Spot checks: a(3) = 14, a(9) = 110, a(16) = 32, a(48) = 32·14 = 448. All match the OEIS data.

## Proof
Write Q'(b, c, d) = −b² − c² − d². The system is a² + Q'(b, c, d) ≡ 1 with 2ab ≡ 2ac ≡ 2ad ≡ 0.

**Multiplicativity.** This follows from the Chinese remainder theorem, as for A236553.

### Odd p
The case split is the same as for A236553: a unit forces a = ±1 and b = c = d = 0, giving 2 solutions; a non-unit a ≠ 0 is impossible; and a = 0 leaves M' = #{Q'(b, c, d) ≡ 1 mod p^k}.

**Counting mod p: points on the sphere b² + c² + d² ≡ −1.** Every element of F_p is a sum of two squares, so write −1 = s² + t². Then

    (b, c, d) ↦ (d, s·b + t·c, t·b − s·c)

is a bijection of F_p³, with inverse (b', c', d') ↦ (−s c' − t d', −t c' + s d', b'). It turns −b² − c² − d² into c'² + d'² − b'², because (sb + tc)² + (tb − sc)² = (s² + t²)(b² + c²) = −(b² + c²). So the sphere has exactly as many points as the A236553 surface: **p² + p**.

**Lifting.** Exactly as for A236553, with the linear form −2(x_b t_b + x_c t_c + x_d t_d): each solution mod p^j (j ≥ 1) has exactly p² lifts. So M' = p^(2k−1)(p + 1), and **a(p^k) = p^(2k−1)(p+1) + 2**.

### p = 2, k ≥ 3
- **a odd:** a is a unit, so 2b ≡ 2c ≡ 2d ≡ 0 and b, c, d ∈ {0, 2^(k−1)}. Their squares are ≡ 0, so a² ≡ 1, which has exactly 4 solutions mod 2^k. That gives **4 · 2³ = 32** solutions.
- **a even with 2a ≢ 0:** then b, c, d are all even, the left side is ≡ 0 (mod 4), and there are **no** solutions.
- **a ∈ {0, 2^(k−1)}:** then a² ≡ 0, and we would need b² + c² + d² ≡ −1 ≡ 7 (mod 8). But **a sum of three squares is never 7 mod 8**, so there are **no** solutions.

So a(2^k) = 32 for k ≥ 3. Direct counts give a(2) = 8 and a(4) = 32.

**Contrast with A236553:** at odd primes the signs of i² and j² make no difference. At 2, the sum-of-three-squares obstruction wipes out every "a even" solution. That is why A227867(2^k) is stuck at 32, while A236553(2^k) = 2^(2k+2) + 32 keeps growing.

## Verification
- **Direct count.** `verify.py` counts the system directly for every n = 1..256, compares the count with the closed form, and checks the closed form against all 52 OEIS data terms. Everything matches.
- **Formal proof.** `lean/OeisLean/A227867.lean` proves:
  - `A_mul` and `A_odd_prime_pow`, the latter with the sphere count reduced to the A236553 count through the change of variables above (using Mathlib's `ZMod.sq_add_sq`);
  - `A_two`, `A_four` and `A_two_pow`, where the kernel checks "never 7 mod 8" by exhausting ZMod 8;
  - `A_closed_form`.

  Its statement is the count of solutions over `ZMod n`, taken straight from the definition, written as the product of the local factors over the factorization of n. `#print axioms` shows only propext, Classical.choice and Quot.sound.
- **b-file.** `bfile.py` writes n = 1..10000 from the closed form, after checking it against the direct count for n ≤ 256.
- **PARI.** The PARI line below was run in gp 2.15.4 for every n in the b-file; all 10000 terms match. The Python program was also checked against the b-file.

## Programs
Both compute every term, from the closed form.

    (PARI) a(n) = my(f=factor(n)); prod(i=1, #f~, my(p=f[i,1], e=f[i,2]); if(p==2, if(e==1, 8, 32), p^(2*e-1)*(p+1)+2));

    (Python)
    from sympy import factorint
    def A227867(n):
        r = 1
        for p, e in factorint(n).items():
            r *= (8 if e == 1 else 32) if p == 2 else p**(2*e-1)*(p+1)+2
        return r

## How to reproduce
- `python verify.py 256` — direct count of A227867 (and A236553) for n = 1..256, compared with the closed forms and the OEIS data.
- `python bfile.py` — rebuild `b227867.txt`.
- `cd lean && lake exe cache get && lake build` — check the Lean proof.

## Files
| File | What it is |
|---|---|
| `README.md` | This write-up: the result and the proof. |
| `layman.md` | Plain-English explanation. |
| `verify.py` | Direct count of A236553 and A227867 for n ≤ 256, compared with the closed forms and the OEIS data (shared with A236553). |
| `bfile.py`, `b227867.txt` | b-file for n = 1..10000 from the closed form, checked against the direct count for n ≤ 256. |
| `lean/` | Lean 4 + Mathlib formal proof (`OeisLean/A227867.lean`, plus `OeisLean/A236553.lean` and `OeisLean/A236554.lean`, which it uses). Build with `cd lean && lake exe cache get && lake build`. |

## Log
- 2026-09-26: Closed form derived alongside A236553. `verify.py` confirms it for n ≤ 256 and against all OEIS data.
- 2026-09-26: Formalized in Lean (`A_closed_form`); the sphere count reduced to the A236553 count; b-file written.
