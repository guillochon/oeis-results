import Mathlib.Order.Lattice.Nat
import Mathlib.Order.ConditionallyCompleteLattice.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Tactic

/-!
# A157615: longest alternating self-avoiding path on an `n × n` board, odd `n`

A157615(n) is the maximal number of squares covered by a self-avoiding path on an `n × n` board
made of alternating horizontal and vertical unit steps. David Wilson conjectured
`a(n) = n^2 - n + 2` for even `n` and `a(n) = n^2 - 2n + 4` for odd `n > 1`.

We prove the odd case (`a_odd`):
* **upper bound** (`le_of_valid`): colour the squares by the parities of their coordinates. Along
  an alternating path the colours cycle with period 4, so every 4 consecutive squares contain an
  (odd, odd) square. There are only `((n-1)/2)^2` of those, so a path has at most
  `4((n-1)/2)^2 + 3 = n^2 - 2n + 4` squares;
* **construction** (`exists_path`): frame induction. Add a width-2 L-shaped frame to the board,
  route the path through it (gaining `4n` squares), and rotate by 180 degrees.
-/

namespace A157615

/-- A horizontal unit step from `u` to `v`. -/
def HStep (u v : ℕ × ℕ) : Prop := v.2 = u.2 ∧ (v.1 = u.1 + 1 ∨ u.1 = v.1 + 1)

/-- A vertical unit step from `u` to `v`. -/
def VStep (u v : ℕ × ℕ) : Prop := v.1 = u.1 ∧ (v.2 = u.2 + 1 ∨ u.2 = v.2 + 1)

/-- `p 0, …, p (L-1)` is a self-avoiding path on the `n × n` board whose unit steps alternate
between horizontal and vertical (step `i` is horizontal iff `i + b` is even, for some fixed `b`). -/
structure Valid (n L : ℕ) (p : ℕ → ℕ × ℕ) : Prop where
  lt : ∀ i < L, (p i).1 < n ∧ (p i).2 < n
  inj : ∀ i < L, ∀ j < L, p i = p j → i = j
  alt : ∃ b, ∀ i, i + 1 < L →
    ((i + b) % 2 = 0 → HStep (p i) (p (i + 1))) ∧ ((i + b) % 2 = 1 → VStep (p i) (p (i + 1)))

/-- **A157615(n)**: the maximal length of such a path. -/
noncomputable def a (n : ℕ) : ℕ := sSup {L | ∃ p, Valid n L p}

/-! ### Upper bound -/

/-- Every window of 4 consecutive squares contains a square with both coordinates odd. -/
theorem window {n L : ℕ} {p : ℕ → ℕ × ℕ} (hv : Valid n L p) (i : ℕ) (hi : i + 3 < L) :
    ∃ j, i ≤ j ∧ j ≤ i + 3 ∧ (p j).1 % 2 = 1 ∧ (p j).2 % 2 = 1 := by
  obtain ⟨b, hb⟩ := hv.alt
  have h0 := hb i (by omega)
  have h1 := hb (i + 1) (by omega)
  have h2 := hb (i + 2) (by omega)
  simp only [show i + 1 + 1 = i + 2 from rfl, show i + 2 + 1 = i + 3 from rfl, HStep,
    VStep] at h0 h1 h2
  have key : ((p i).1 % 2 = 1 ∧ (p i).2 % 2 = 1) ∨ ((p (i + 1)).1 % 2 = 1 ∧ (p (i + 1)).2 % 2 = 1) ∨
      ((p (i + 2)).1 % 2 = 1 ∧ (p (i + 2)).2 % 2 = 1) ∨
      ((p (i + 3)).1 % 2 = 1 ∧ (p (i + 3)).2 % 2 = 1) := by omega
  rcases key with h | h | h | h
  · exact ⟨i, le_refl _, by omega, h⟩
  · exact ⟨i + 1, by omega, by omega, h⟩
  · exact ⟨i + 2, by omega, by omega, h⟩
  · exact ⟨i + 3, by omega, le_refl _, h⟩

/-- There are `k` odd numbers below `2k + 1`. -/
theorem card_odd_lt (k : ℕ) :
    ((Finset.range (2 * k + 1)).filter (fun x => x % 2 = 1)).card = k := by
  have : (Finset.range (2 * k + 1)).filter (fun x => x % 2 = 1) =
      (Finset.range k).image (fun i => 2 * i + 1) := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_image]
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨x / 2, by omega, by omega⟩
    · rintro ⟨i, hi, rfl⟩; omega
  rw [this, Finset.card_image_of_injective _ (fun a b (h : 2 * a + 1 = 2 * b + 1) => by omega),
    Finset.card_range]

/-- **Upper bound**: on the `(2k+1) × (2k+1)` board, a path has at most `4k^2 + 3` squares. -/
theorem le_of_valid {k L : ℕ} {p : ℕ → ℕ × ℕ} (hv : Valid (2 * k + 1) L p) : L ≤ 4 * k ^ 2 + 3 := by
  obtain ⟨K, hK⟩ : ∃ K, K = k ^ 2 := ⟨_, rfl⟩
  rw [← hK]
  by_contra hL
  push Not at hL
  have hw : ∀ w, w < K + 1 → ∃ j, 4 * w ≤ j ∧ j ≤ 4 * w + 3 ∧ (p j).1 % 2 = 1 ∧ (p j).2 % 2 = 1 :=
    fun w hw => window hv (4 * w) (by omega)
  choose! g hg using hw
  set S := (Finset.range (2 * k + 1)).filter (fun x => x % 2 = 1)
  have hcard := Finset.card_le_card_of_injOn (fun w => p (g w)) (s := Finset.range (K + 1))
    (t := S ×ˢ S) ?maps ?inj
  · rw [Finset.card_range, Finset.card_product, card_odd_lt] at hcard
    have : k * k = K := by rw [hK, sq]
    omega
  · intro w hw
    rw [Finset.coe_range, Set.mem_Iio] at hw
    obtain ⟨h1, h2, h3, h4⟩ := hg w hw
    have hb := hv.lt (g w) (by omega)
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, S, Finset.mem_filter,
      Finset.mem_range]
    exact ⟨⟨hb.1, h3⟩, hb.2, h4⟩
  · intro w₁ hw₁ w₂ hw₂ heq
    rw [Finset.coe_range, Set.mem_Iio] at hw₁ hw₂
    have e := hv.inj (g w₁) (by have := hg w₁ hw₁; omega) (g w₂) (by have := hg w₂ hw₂; omega) heq
    have := hg w₁ hw₁
    have := hg w₂ hw₂
    omega

/-! ### Construction by frame induction -/

/-- The route through the width-2 frame (new columns `n, n+1` and rows `n, n+1`), for odd `n`:
`(n, 1)`; a zigzag down the right strip (rows `2..n`); `(n, n+1)`; a zigzag left along the bottom
strip (columns `n-1..0`). It has `4n` squares. -/
def route (n j : ℕ) : ℕ × ℕ :=
  if j = 0 then (n, 1)
  else if j ≤ 2 * n - 2 then (n + ((j - 1) % 2 + ((j - 1) / 2 + 2)) % 2, (j - 1) / 2 + 2)
  else if j = 2 * n - 1 then (n, n + 1)
  else (n - 1 - (j - 2 * n) / 2, n + 1 - ((j - 2 * n) / 2 + (j - 2 * n) % 2) % 2)

section route

variable {n : ℕ} (hn : n % 2 = 1) (h3 : 3 ≤ n)
include hn h3

theorem route_lt {j : ℕ} (hj : j < 4 * n) :
    (route n j).1 < n + 2 ∧ (route n j).2 < n + 2 := by
  unfold route; split_ifs <;> first | contradiction | (simp only; omega)

omit h3 in
theorem route_out {j : ℕ} (hj : j < 4 * n) : n ≤ (route n j).1 ∨ n ≤ (route n j).2 := by
  unfold route; split_ifs <;> first | contradiction | (simp only; omega)

theorem route_inj {i j : ℕ} (hi : i < 4 * n) (hj : j < 4 * n) (h : route n i = route n j) :
    i = j := by
  unfold route at h
  split_ifs at h <;> first | contradiction | (simp only [Prod.mk.injEq] at h; omega)

theorem route_step {j : ℕ} (hj : j + 1 < 4 * n) :
    (j % 2 = 1 → HStep (route n j) (route n (j + 1))) ∧
      (j % 2 = 0 → VStep (route n j) (route n (j + 1))) := by
  unfold route HStep VStep; split_ifs <;> first | contradiction | (simp only; omega)

omit hn h3 in
theorem route_zero : route n 0 = (n, 1) := by simp [route]

theorem route_end1 : route n (4 * n - 2) = (0, n + 1) := by
  unfold route; split_ifs <;> first | contradiction | (simp only [Prod.mk.injEq]; omega)

theorem route_end2 : route n (4 * n - 1) = (0, n) := by
  unfold route; split_ifs <;> first | contradiction | (simp only [Prod.mk.injEq]; omega)

end route

/-- The path followed by the frame route. -/
def ext (n L : ℕ) (p : ℕ → ℕ × ℕ) (i : ℕ) : ℕ × ℕ := if i < L then p i else route n (i - L)

/-- Rotation by 180 degrees of the `N × N` board. -/
def rot (N : ℕ) (u : ℕ × ℕ) : ℕ × ℕ := (N - 1 - u.1, N - 1 - u.2)

theorem valid_rot {N L : ℕ} {p : ℕ → ℕ × ℕ} (hv : Valid N L p) :
    Valid N L (fun i => rot N (p i)) := by
  refine ⟨fun i hi => ?_, fun i hi j hj h => ?_, ?_⟩
  · have := hv.lt i hi; simp only [rot]; omega
  · have := hv.lt i hi; have := hv.lt j hj
    simp only [rot, Prod.mk.injEq] at h
    exact hv.inj i hi j hj (Prod.ext (by omega) (by omega))
  · obtain ⟨b, hb⟩ := hv.alt
    refine ⟨b, fun i hi => ?_⟩
    have h := hb i hi
    have := hv.lt i (by omega); have := hv.lt (i + 1) hi
    simp only [HStep, VStep, rot] at h ⊢
    omega

/-- Appending the frame route to a path that ends `(n-1, 0) → (n-1, 1)` gives a valid path on the
`(n+2) × (n+2)` board with `4n` more squares. -/
theorem valid_ext {n L : ℕ} {p : ℕ → ℕ × ℕ} (hn : n % 2 = 1) (h3 : 3 ≤ n) (hv : Valid n L p)
    (hL : 2 ≤ L) (he1 : p (L - 2) = (n - 1, 0)) (he2 : p (L - 1) = (n - 1, 1)) :
    Valid (n + 2) (L + 4 * n) (ext n L p) := by
  refine ⟨fun i hi => ?_, fun i hi j hj h => ?_, ?_⟩
  · unfold ext; split_ifs with h
    · have := hv.lt i h; omega
    · exact route_lt hn h3 (by omega)
  · unfold ext at h; split_ifs at h with h1 h2 h2
    · exact hv.inj i h1 j h2 h
    · have := hv.lt i h1; have := route_out hn (j := j - L) (by omega); rw [← h] at this; omega
    · have := hv.lt j h2; have := route_out hn (j := i - L) (by omega); rw [h] at this; omega
    · have := route_inj hn h3 (i := i - L) (j := j - L) (by omega) (by omega) h; omega
  · obtain ⟨b, hb⟩ := hv.alt
    -- the last step of `p` is vertical, which fixes the parity of the route's steps
    have hpar : (L - 2 + b) % 2 = 1 := by
      by_contra hc
      have := (hb (L - 2) (by omega)).1 (by omega)
      rw [show L - 2 + 1 = L - 1 by omega, he1, he2] at this
      simp [HStep] at this
    refine ⟨b, fun i hi => ?_⟩
    by_cases h1 : i + 1 < L
    · simp only [ext, if_pos h1, if_pos (show i < L by omega)]
      exact hb i h1
    · by_cases h2 : i < L
      · obtain rfl : i = L - 1 := by omega
        simp only [ext, if_pos h2, if_neg h1, show L - 1 + 1 - L = 0 by omega, route_zero, he2,
          HStep, VStep, true_and]
        omega
      · simp only [ext, if_neg h1, if_neg h2]
        have hs := route_step hn h3 (j := i - L) (by omega)
        rw [show i - L + 1 = i + 1 - L by omega] at hs
        constructor
        · intro h; exact hs.1 (by omega)
        · intro h; exact hs.2 (by omega)

/-- The optimal 5 × 5 path, ending `(4, 0) → (4, 1)`. -/
def P5 (i : ℕ) : ℕ × ℕ :=
  ([(1, 0), (0, 0), (0, 1), (1, 1), (1, 2), (0, 2), (0, 3), (1, 3), (1, 4), (2, 4), (2, 3),
    (3, 3), (3, 2), (2, 2), (2, 1), (3, 1), (3, 0), (4, 0), (4, 1)] : List (ℕ × ℕ)).getD i (0, 0)

set_option maxRecDepth 10000 in
theorem valid_P5 : Valid 5 19 P5 := by
  refine ⟨by decide +kernel, by decide +kernel, ⟨0, fun i hi => ?_⟩⟩
  have : i < 18 := by omega
  interval_cases i <;> (unfold HStep VStep; decide)

/-- An optimal 3 × 3 path (7 squares). -/
def P3 (i : ℕ) : ℕ × ℕ :=
  ([(1, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2), (0, 1)] : List (ℕ × ℕ)).getD i (0, 0)

theorem valid_P3 : Valid 3 7 P3 := by
  refine ⟨by decide +kernel, by decide +kernel, ⟨0, fun i hi => ?_⟩⟩
  have : i < 6 := by omega
  interval_cases i <;> (unfold HStep VStep; decide)

/-- **Construction**: for every `k`, an alternating path of `(2k+4)^2 + 3` squares on the
`(2k+5) × (2k+5)` board, ending `(2k+4, 0) → (2k+4, 1)`. -/
theorem exists_path (k : ℕ) : ∃ p, Valid (2 * k + 5) ((2 * k + 4) ^ 2 + 3) p ∧
    p ((2 * k + 4) ^ 2 + 1) = (2 * k + 4, 0) ∧ p ((2 * k + 4) ^ 2 + 2) = (2 * k + 4, 1) := by
  induction k with
  | zero => exact ⟨P5, valid_P5, by decide, by decide⟩
  | succ k ih =>
    obtain ⟨p, hv, e1, e2⟩ := ih
    obtain ⟨S, hS⟩ : ∃ S, S = (2 * k + 4) ^ 2 := ⟨_, rfl⟩
    have hS' : (2 * (k + 1) + 4) ^ 2 = S + 8 * k + 20 := by rw [hS]; ring
    rw [← hS] at e1 e2 hv
    rw [hS']
    have hn : (2 * k + 5) % 2 = 1 := by omega
    have h3 : 3 ≤ 2 * k + 5 := by omega
    have he1 : p (S + 3 - 2) = (2 * k + 5 - 1, 0) := by
      rw [show S + 3 - 2 = S + 1 by omega, e1]; congr 1
    have he2 : p (S + 3 - 1) = (2 * k + 5 - 1, 1) := by
      rw [show S + 3 - 1 = S + 2 by omega, e2]; congr 1
    have hq := valid_rot (valid_ext hn h3 hv (by omega) he1 he2)
    have hext : ∀ j, ext (2 * k + 5) (S + 3) p (S + 3 + j) = route (2 * k + 5) j :=
      fun j => by simp [ext]
    refine ⟨fun i => rot (2 * k + 5 + 2) (ext (2 * k + 5) (S + 3) p i), ?_, ?_, ?_⟩
    · rw [show 2 * (k + 1) + 5 = 2 * k + 5 + 2 by omega,
        show S + 8 * k + 20 + 3 = S + 3 + 4 * (2 * k + 5) by omega]
      exact hq
    · rw [show S + 8 * k + 20 + 1 = S + 3 + (4 * (2 * k + 5) - 2) by omega]
      simp only [hext, route_end1 hn h3, rot]
      ext <;> simp only <;> omega
    · rw [show S + 8 * k + 20 + 2 = S + 3 + (4 * (2 * k + 5) - 1) by omega]
      simp only [hext, route_end2 hn h3, rot]
      ext <;> simp only <;> omega

/-! ### The odd case of Wilson's conjecture -/

/-- **A157615 for odd `n`**: `a(n) = n^2 - 2n + 4` for every odd `n ≥ 3`. -/
theorem a_odd {n : ℕ} (hn : n % 2 = 1) (h3 : 3 ≤ n) : a n = n ^ 2 - 2 * n + 4 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = 2 * k + 1 := ⟨n / 2, by omega⟩
  have hval : (2 * k + 1) ^ 2 - 2 * (2 * k + 1) + 4 = 4 * k ^ 2 + 3 := by
    have e : (2 * k + 1) ^ 2 = 4 * k ^ 2 + 4 * k + 1 := by ring
    have hk2 : 1 ≤ k ^ 2 := Nat.one_le_pow _ _ (by omega)
    rw [e]; omega
  rw [hval]
  apply IsGreatest.csSup_eq
  refine ⟨?_, fun L ⟨p, hv⟩ => le_of_valid hv⟩
  rcases Nat.lt_or_ge k 2 with hk | hk
  · obtain rfl : k = 1 := by omega
    exact ⟨P3, valid_P3⟩
  · obtain ⟨j, rfl⟩ : ∃ j, k = j + 2 := ⟨k - 2, by omega⟩
    obtain ⟨p, hv, -, -⟩ := exists_path j
    refine ⟨p, ?_⟩
    rw [show 2 * (j + 2) + 1 = 2 * j + 5 by ring,
      show 4 * (j + 2) ^ 2 + 3 = (2 * j + 4) ^ 2 + 3 by ring]
    exact hv

end A157615
