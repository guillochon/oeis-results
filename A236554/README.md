# A236554 — Involutions in the quaternion ring over Z/(2^n)Z with i^2 = j^2 = 1

**Status:** **proved and formally verified in Lean 4 + Mathlib** (`lean/OeisLean/A236554.lean`, theorem `A236554.card_A236554`; axioms: propext, Classical.choice, Quot.sound only). Ready for OEIS submission.
**Known terms:** a(1)..a(8) = 8, 64, 288, 1056, 4128, 16416, 65568, 262176 (Grau Ribas, 2014).
**Result:** **a(n) = 2^(2n+2) + 32 for n ≥ 3** (proved). The first new terms are a(9) = 1048608, a(10) = 4194336, a(11) = 16777248.

## What is being counted
The ring has elements a + b·i + c·j + d·k with a, b, c, d ∈ Z/2^n and i² = j² = 1. It is a "split" quaternion ring, not the usual i² = j² = −1 one. An *involution* is an element X with X² = 1. Expanding X² = 1 gives the system from A236553:

    a² − b² + c² + d² ≡ 1,   2ab ≡ 0,   2ac ≡ 0,   2ad ≡ 0   (mod 2^n)

and a(n) is the number of its solutions (a, b, c, d).

**Why it's tagged `hard`:** the only program is Mathematica's `Reduce` over all 2^(4n) tuples. The data stops at n = 8 because that program is slow, not because the problem is hard. `count.py` counts the same system in about 2 seconds for n ≤ 13 by grouping on a and taking a convolution of square counts. It works straight from the definition and uses no idea from the proof below.

## The proof (split by the 2-adic valuation of a)
Let N = 2^n with n ≥ 3.

**Case 1: a odd.** Then 2a·x ≡ 0 means 2x ≡ 0, so b, c, d ∈ {0, N/2}. Their squares are 0 or 2^(2n−2) ≡ 0 (mod 2^n), so the equation reduces to a² ≡ 1 (mod 2^n). That has exactly 4 odd solutions for n ≥ 3 (±1 and ±1 + 2^(n−1)). Each of b, c, d has 2 choices, so this case gives **4 · 8 = 32**, which is the "+32".

**Case 2: a even,** a = 2^v·u with u odd (take v = n if a = 0). The condition 2ax ≡ 0 means v₂(x) ≥ n − 1 − v.
- If v ≤ n − 2, then b, c, d are all even, so the left side is ≡ 0 (mod 4) while the right side is 1. No solutions.
- So v ≥ n − 1, which means a ∈ {0, 2^(n−1)}. Then a² ≡ 0, and 2ax ≡ 0 holds for every x, so b, c, d are unconstrained. This case gives **2·M(n)**, where

      M(n) = #{(b, c, d) mod 2^n : c² + d² − b² ≡ 1}.

**Lemma: M(n+1) = 4·M(n) for n ≥ 3.** Write Q = c² + d² − b².
1. Lifting. If x ≡ x′ (mod 2^n), then Q(x) ≡ Q(x′) (mod 2^(n+1)), because the cross terms carry a factor 2^(n+1) and the square terms 2^(2n). So each solution mod 2^n has 8 lifts mod 2^(n+1), and they all share one value of Q mod 2^(n+1): either 1 or 1 + 2^n. Hence M(n+1) = 8·#{x ∈ S : Q(x) ≡ 1 mod 2^(n+1)}, where S is the solution set mod 2^n.
2. An involution on S that swaps the two classes. Q(x) is odd, so at least one coordinate is odd. Add 2^(n−1) to the first odd coordinate t. Q changes by ±(2^n·t + 2^(2n−2)), which is ≡ 2^n (mod 2^(n+1)) because t is odd and 2n − 2 ≥ n + 1 when n ≥ 3. The change is ≡ 0 (mod 2^n), so the new point is still in S. Coordinate parities are unchanged, so doing it twice returns to x mod 2^n.
3. So exactly half of S has Q ≡ 1 (mod 2^(n+1)), and M(n+1) = 8·M(n)/2 = 4·M(n). ∎

From the computed a(3) = 288 we get M(3) = (288 − 32)/2 = 128, so M(n) = 2^(2n+1) and **a(n) = 32 + 2^(2n+2) for n ≥ 3.** (a(1) = 8 and a(2) = 64 are outside the lemma's range.)

**Equivalent forms.** Since a(n) = 4^(n+1) + 32 for n ≥ 3, with a(1) = 8 and a(2) = 64:
- generating function: Σ a(n)·xⁿ = 8x(1 + 3x − 16x³) / ((1 − x)(1 − 4x));
- recurrence: a(n) = 5·a(n−1) − 4·a(n−2) for n ≥ 5.

Both follow directly from the closed form (the recurrence has characteristic roots 1 and 4), and the generating function was checked against the closed form for n ≤ 29 with sympy.

**Formal verification (2026-09-26):** `lean/OeisLean/A236554.lean` (about 580 lines) proves

    theorem card_A236554 (n : ℕ) (hn : 3 ≤ n) : (T n).card = 2 ^ (2 * n + 2) + 32

where `T n` is the solution set of the system above over `ZMod (2^n)`, taken straight from the definition. It follows the paper proof. M(3) = 128 and a(1) = 8, a(2) = 64 are checked by the Lean kernel (`decide +kernel`, no `native_decide`). Build it with `cd lean && lake build`.

## Checks
- **Direct count.** `count.py` counts the system straight from the definition, with no idea from the proof. It reproduces the published a(1)..a(8) and gives a(9)..a(13); all agree with the formula.
- **Formal proof.** The Lean theorem above, whose statement is the defining system itself; `#print axioms` shows only propext, Classical.choice and Quot.sound.
- **b-file.** `bfile.py` writes n = 1..500 from the formula, after asserting that it matches `count.py` for every n ≤ 13.
- **PARI.** The PARI line below was run in gp 2.15.4 for every n in the b-file; all 500 terms match.

## Programs
Both compute every term, from the closed form.

    (PARI) a(n) = if(n < 3, [8, 64][n], 4^(n+1) + 32);

    (Python)
    def A236554(n): return (8, 64)[n-1] if n < 3 else 4**(n+1) + 32

## How to reproduce
- `python count.py 13` — direct count for n = 1..13, compared with the formula (about 2 seconds).
- `python bfile.py` — rebuild `b236554.txt`.
- `cd lean && lake exe cache get && lake build` — check the Lean proof.

## Related sequences
A236553(n) counts the same system modulo any n. It is multiplicative, and A236554(n) = A236553(2^n). A closed form for all of A236553, and for its Lipschitz-quaternion sibling A227867, is proved in `../A236553/` and `../A227867/`. Those proofs reuse this one for the powers of 2.

## Files
| File | What it is |
|---|---|
| `README.md` | This write-up: the result and the proof. |
| `layman.md` | Plain-English explanation. |
| `count.py` | Direct count from the definition (no proof ideas used). It reproduces a(1..8) and gives a(9..13). |
| `bfile.py`, `b236554.txt` | b-file for n = 1..500 from the formula, checked against `count.py` for n ≤ 13. |
| `lean/` | Lean 4 + Mathlib formal proof (`OeisLean/A236554.lean`, theorem `A236554.card_A236554`). Build with `cd lean && lake exe cache get && lake build`. |

## Log
- 2026-09-26: `count.py` reproduces a(1..8) and gives a(9..13), all matching 2^(2n+2) + 32. Proof drafted (valuation split plus a halving lemma).
- 2026-09-26: Proof formalized in Lean (`card_A236554`); b-file written.
