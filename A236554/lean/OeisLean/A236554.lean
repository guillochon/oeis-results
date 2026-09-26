import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

/-!
# A236554: involutions in the quaternion ring over `ZMod (2^n)` with `i^2 = j^2 = 1`

A236554(n) is the number of solutions `(a, b, c, d)` modulo `2^n` of
  `a^2 - b^2 + c^2 + d^2 = 1`, `2ab = 0`, `2ac = 0`, `2ad = 0`.
We prove `A236554(n) = 2^(2n+2) + 32` for `n ≥ 3` (`card_A236554`).

Proof outline (see `attacks/A236554/README.md`):
* the solutions are `{a | a^2 = 1} × {0, h}^3` together with `{0, h} × S n`, where `h = 2^(n-1)`
  and `S n = {(b, c, d) | c^2 + d^2 - b^2 = 1}`;
* `a^2 = 1` has exactly 4 solutions mod `2^n`;
* `|S (n+1)| = 4 |S n|` for `n ≥ 3`: every solution mod `2^n` has 8 lifts, which all have the same
  value of `Q` mod `2^(n+1)`, and an involution shows exactly half of them are solutions.
-/

open Finset

namespace A236554

/-- `2^n = 0` in `ZMod (2^n)`. -/
theorem two_pow_self (n : ℕ) : (2 : ZMod (2 ^ n)) ^ n = 0 := by
  have := ZMod.natCast_self (2 ^ n)
  push_cast at this
  exact this

theorem natCast_ne_zero_of_lt {n k : ℕ} (h0 : 0 < k) (hk : k < 2 ^ n) :
    ((k : ℕ) : ZMod (2 ^ n)) ≠ 0 := by
  rw [Ne, ZMod.natCast_eq_zero_iff]
  exact Nat.not_dvd_of_pos_of_lt h0 hk

/-- `2^k * z = 0` in `ZMod (2^n)` iff `z` is a multiple of `2^(n-k)`. -/
theorem two_pow_mul_eq_zero_iff {n k : ℕ} (hk : k ≤ n) (z : ZMod (2 ^ n)) :
    (2 : ZMod (2 ^ n)) ^ k * z = 0 ↔ ∃ j : ZMod (2 ^ n), z = 2 ^ (n - k) * j := by
  constructor
  · intro h
    have hz := ZMod.natCast_zmod_val z
    rw [← hz] at h
    have h' : (((2 ^ k * z.val : ℕ)) : ZMod (2 ^ n)) = 0 := by push_cast; exact h
    rw [ZMod.natCast_eq_zero_iff] at h'
    have key : 2 ^ k * 2 ^ (n - k) ∣ 2 ^ k * z.val := by
      rwa [← pow_add, Nat.add_sub_cancel' hk]
    obtain ⟨q, hq⟩ := Nat.dvd_of_mul_dvd_mul_left (by positivity) key
    refine ⟨q, ?_⟩
    rw [← hz, hq]
    push_cast
    ring
  · rintro ⟨j, rfl⟩
    rw [← mul_assoc, ← pow_add, Nat.add_sub_cancel' hk, two_pow_self, zero_mul]

/-! ### Parity -/

/-- Reduction mod 2. -/
def ρ (n : ℕ) (hn : n ≠ 0) : ZMod (2 ^ n) →+* ZMod 2 :=
  ZMod.castHom (dvd_pow_self 2 hn) (ZMod 2)

theorem zmod2_cases (x : ZMod 2) : x = 0 ∨ x = 1 := by
  revert x; decide

theorem ρ_natCast {n : ℕ} (hn : n ≠ 0) (k : ℕ) : ρ n hn (k : ZMod (2 ^ n)) = (k : ZMod 2) :=
  map_natCast _ k

theorem eq_two_mul_of_ρ_eq_zero {n : ℕ} (hn : n ≠ 0) {z : ZMod (2 ^ n)} (h : ρ n hn z = 0) :
    ∃ w, z = 2 * w := by
  have hz := ZMod.natCast_zmod_val z
  rw [← hz, ρ_natCast, ZMod.natCast_eq_zero_iff] at h
  obtain ⟨q, hq⟩ := h
  exact ⟨q, by rw [← hz, hq]; push_cast; ring⟩

theorem eq_two_mul_add_one_of_ρ_eq_one {n : ℕ} (hn : n ≠ 0) {z : ZMod (2 ^ n)}
    (h : ρ n hn z = 1) : ∃ w, z = 2 * w + 1 := by
  obtain ⟨w, hw⟩ := eq_two_mul_of_ρ_eq_zero hn (z := z - 1) (by rw [map_sub, h, map_one, sub_self])
  exact ⟨w, by rw [← hw]; ring⟩

theorem isUnit_of_ρ_eq_one {n : ℕ} (hn : n ≠ 0) {z : ZMod (2 ^ n)} (h : ρ n hn z = 1) :
    IsUnit z := by
  obtain ⟨w, rfl⟩ := eq_two_mul_add_one_of_ρ_eq_one hn h
  rw [add_comm]
  apply IsNilpotent.isUnit_one_add
  exact ⟨n, by rw [mul_pow, two_pow_self, zero_mul]⟩

/-- Multiplying by `2^(n-1)`, which squares to 0 times 2, only sees the parity. -/
theorem half_mul_mem {n : ℕ} (hn : n ≠ 0) (j : ZMod (2 ^ n)) :
    (2 : ZMod (2 ^ n)) ^ (n - 1) * j = 0 ∨ (2 : ZMod (2 ^ n)) ^ (n - 1) * j = 2 ^ (n - 1) := by
  have h2 : (2 : ZMod (2 ^ n)) ^ (n - 1) * 2 = 0 := by
    rw [← pow_succ, Nat.sub_add_cancel (Nat.one_le_iff_ne_zero.mpr hn), two_pow_self]
  rcases zmod2_cases (ρ n hn j) with h | h
  · obtain ⟨w, rfl⟩ := eq_two_mul_of_ρ_eq_zero hn h
    left; rw [← mul_assoc, h2, zero_mul]
  · obtain ⟨w, rfl⟩ := eq_two_mul_add_one_of_ρ_eq_one hn h
    right; rw [mul_add, ← mul_assoc, h2, zero_mul, zero_add, mul_one]

/-- `2z = 0` iff `z ∈ {0, 2^(n-1)}`. -/
theorem two_mul_eq_zero_iff {n : ℕ} (hn : n ≠ 0) (z : ZMod (2 ^ n)) :
    2 * z = 0 ↔ z = 0 ∨ z = 2 ^ (n - 1) := by
  constructor
  · intro h
    obtain ⟨j, rfl⟩ := (two_pow_mul_eq_zero_iff (k := 1) (Nat.one_le_iff_ne_zero.mpr hn) z).1
      (by rw [pow_one]; exact h)
    exact half_mul_mem hn j
  · rintro (rfl | rfl)
    · rw [mul_zero]
    · rw [← pow_succ', Nat.sub_add_cancel (Nat.one_le_iff_ne_zero.mpr hn), two_pow_self]

/-! ### Square roots of 1 modulo `2^(m+3)` -/

section roots

variable (m : ℕ)

/-- `h = 2^(n-1)` for `n = m + 3`. -/
abbrev hh : ZMod (2 ^ (m + 3)) := 2 ^ (m + 2)

theorem hh_sq : hh m ^ 2 = 0 := by
  unfold hh; linear_combination (2 : ZMod (2 ^ (m + 3))) ^ (m + 1) * two_pow_self (m + 3)

theorem two_mul_hh : 2 * hh m = 0 := by
  unfold hh; linear_combination two_pow_self (m + 3)

theorem ne_of_sub_natCast {x y : ZMod (2 ^ (m + 3))} {k : ℕ} (h0 : 0 < k) (hk : k < 2 ^ (m + 3))
    (h : x - y = k) : x ≠ y := by
  intro hxy
  apply natCast_ne_zero_of_lt h0 hk
  rw [← h, hxy, sub_self]

theorem ρ_eq_one_of_sq_eq_one {n : ℕ} (hn : n ≠ 0) {a : ZMod (2 ^ n)} (ha : a ^ 2 = 1) :
    ρ n hn a = 1 := by
  rcases zmod2_cases (ρ n hn a) with h | h
  · have := congrArg (ρ n hn) ha
    rw [map_pow, h, map_one] at this
    exact absurd this (by decide)
  · exact h

theorem sq_eq_one_iff (a : ZMod (2 ^ (m + 3))) :
    a ^ 2 = 1 ↔ a = 1 ∨ a = -1 ∨ a = 1 + hh m ∨ a = -1 + hh m := by
  have hn : m + 3 ≠ 0 := by omega
  constructor
  · intro ha
    obtain ⟨k, rfl⟩ := eq_two_mul_add_one_of_ρ_eq_one hn (ρ_eq_one_of_sq_eq_one hn ha)
    have h4 : (2 : ZMod (2 ^ (m + 3))) ^ 2 * (k * (k + 1)) = 0 := by linear_combination ha
    rcases zmod2_cases (ρ _ hn k) with hk | hk
    · -- k even, so k + 1 is a unit and 4k = 0
      have hu : IsUnit (k + 1) := isUnit_of_ρ_eq_one hn (by rw [map_add, hk, map_one, zero_add])
      have h4' : (k + 1) * ((2 : ZMod (2 ^ (m + 3))) ^ 2 * k) = 0 := by linear_combination h4
      obtain ⟨j, hj0⟩ := (two_pow_mul_eq_zero_iff (by omega) k).1 ((hu.mul_right_eq_zero).1 h4')
      rw [show m + 3 - 2 = m + 1 from rfl] at hj0
      subst hj0
      rcases half_mul_mem hn j with hj | hj <;> simp only [show m + 3 - 1 = m + 2 from rfl] at hj
      · left; linear_combination hj
      · right; right; left; linear_combination hj
    · -- k odd, so k is a unit and 4(k + 1) = 0
      have hu : IsUnit k := isUnit_of_ρ_eq_one hn hk
      have h4' : k * ((2 : ZMod (2 ^ (m + 3))) ^ 2 * (k + 1)) = 0 := by linear_combination h4
      obtain ⟨j, hj'⟩ :=
        (two_pow_mul_eq_zero_iff (by omega) (k + 1)).1 ((hu.mul_right_eq_zero).1 h4')
      rw [show m + 3 - 2 = m + 1 from rfl] at hj'
      have hk' : k = 2 ^ (m + 1) * j - 1 := by linear_combination hj'
      subst hk'
      rcases half_mul_mem hn j with hj | hj <;> simp only [show m + 3 - 1 = m + 2 from rfl] at hj
      · right; left; linear_combination hj
      · right; right; right; linear_combination hj
  · rintro (rfl | rfl | rfl | rfl)
    · ring
    · ring
    · linear_combination two_mul_hh m + hh_sq m
    · linear_combination -two_mul_hh m + hh_sq m

theorem card_roots :
    (univ.filter (fun a : ZMod (2 ^ (m + 3)) => a ^ 2 = 1)).card = 4 := by
  have hP : 4 ≤ 2 ^ (m + 2) := by
    calc 4 = 2 ^ 2 := by norm_num
      _ ≤ 2 ^ (m + 2) := Nat.pow_le_pow_right (by norm_num) (by omega)
  have hN : 2 ^ (m + 3) = 2 * 2 ^ (m + 2) := by ring
  have e : univ.filter (fun a : ZMod (2 ^ (m + 3)) => a ^ 2 = 1) =
      {1, -1, 1 + hh m, -1 + hh m} := by
    ext a; simp [sq_eq_one_iff]
  have hh_cast : hh m = ((2 ^ (m + 2) : ℕ) : ZMod (2 ^ (m + 3))) := by push_cast; rfl
  have d1 : (1 : ZMod (2 ^ (m + 3))) ≠ -1 :=
    ne_of_sub_natCast m (k := 2) (by norm_num) (by omega) (by push_cast; ring)
  have d2 : (1 : ZMod (2 ^ (m + 3))) ≠ 1 + hh m :=
    (ne_of_sub_natCast m (k := 2 ^ (m + 2)) (by omega) (by omega)
      (by rw [← hh_cast]; ring)).symm
  have d3 : (1 : ZMod (2 ^ (m + 3))) ≠ -1 + hh m :=
    (ne_of_sub_natCast m (k := 2 ^ (m + 2) - 2) (by omega) (by omega)
      (by rw [Nat.cast_sub (by omega), ← hh_cast]; push_cast; ring)).symm
  have d4 : (-1 : ZMod (2 ^ (m + 3))) ≠ 1 + hh m :=
    (ne_of_sub_natCast m (k := 2 ^ (m + 2) + 2) (by omega) (by omega)
      (by push_cast; ring)).symm
  have d5 : (-1 : ZMod (2 ^ (m + 3))) ≠ -1 + hh m :=
    (ne_of_sub_natCast m (k := 2 ^ (m + 2)) (by omega) (by omega)
      (by rw [← hh_cast]; ring)).symm
  have d6 : (1 + hh m : ZMod (2 ^ (m + 3))) ≠ -1 + hh m :=
    ne_of_sub_natCast m (k := 2) (by norm_num) (by omega) (by push_cast; ring)
  rw [e, card_insert_of_notMem (by simp [d1, d2, d3]), card_insert_of_notMem (by simp [d4, d5]),
    card_insert_of_notMem (by simp [d6]), card_singleton]

end roots

/-! ### The ternary form `Q(b, c, d) = c^2 + d^2 - b^2` and its lifting behaviour -/

/-- `Q (b, c, d) = c^2 + d^2 - b^2`. -/
def Q {α : Type*} [CommRing α] (x : α × α × α) : α := x.2.1 ^ 2 + x.2.2 ^ 2 - x.1 ^ 2

/-- Solutions of `Q = 1` modulo `2^n`. -/
def S (n : ℕ) : Finset (ZMod (2 ^ n) × ZMod (2 ^ n) × ZMod (2 ^ n)) :=
  univ.filter (fun x => Q x = 1)

theorem mem_S {n : ℕ} (x : ZMod (2 ^ n) × ZMod (2 ^ n) × ZMod (2 ^ n)) :
    x ∈ S n ↔ x.2.1 ^ 2 + x.2.2 ^ 2 - x.1 ^ 2 = 1 := by
  unfold S; rw [mem_filter]; simp [Q]

section lift

variable (m : ℕ)

/-- Reduction `ZMod (2^(m+4)) → ZMod (2^(m+3))`. -/
def π : ZMod (2 ^ (m + 4)) →+* ZMod (2 ^ (m + 3)) :=
  ZMod.castHom (pow_dvd_pow 2 (by omega)) _

/-- The reduction on triples, as an additive group hom. -/
def π3 : (ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4))) →+
    (ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3))) :=
  (π m).toAddMonoidHom.prodMap ((π m).toAddMonoidHom.prodMap (π m).toAddMonoidHom)

theorem π3_apply (y : ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4))) :
    π3 m y = (π m y.1, π m y.2.1, π m y.2.2) := rfl

/-- The canonical lift `ZMod (2^(m+3)) → ZMod (2^(m+4))` (not additive). -/
def σ (x : ZMod (2 ^ (m + 3))) : ZMod (2 ^ (m + 4)) := (x.val : ℕ)

def σ3 (x : ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3))) :
    ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4)) :=
  (σ m x.1, σ m x.2.1, σ m x.2.2)

theorem π_σ (x : ZMod (2 ^ (m + 3))) : π m (σ m x) = x := by
  unfold π σ; rw [map_natCast, ZMod.natCast_zmod_val]

theorem π3_σ3 (x) : π3 m (σ3 m x) = x := by
  simp [π3_apply, σ3, π_σ]

theorem Q_map {α β : Type*} [CommRing α] [CommRing β] (f : α →+* β) (x : α × α × α) :
    Q (f x.1, f x.2.1, f x.2.2) = f (Q x) := by
  simp [Q, map_add, map_sub, map_pow]

theorem eq_of_π_eq_zero {z : ZMod (2 ^ (m + 4))} (h : π m z = 0) :
    ∃ j, z = 2 ^ (m + 3) * j := by
  have hz := ZMod.natCast_zmod_val z
  rw [← hz] at h
  unfold π at h
  rw [map_natCast, ZMod.natCast_eq_zero_iff] at h
  obtain ⟨q, hq⟩ := h
  exact ⟨q, by rw [← hz, hq]; push_cast; ring⟩

/-- `Q` modulo `2^(m+4)` only depends on the argument modulo `2^(m+3)`. -/
theorem Q_congr {y y' : ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4)) × ZMod (2 ^ (m + 4))}
    (h : π3 m y = π3 m y') : Q y = Q y' := by
  obtain ⟨b, c, d⟩ := y
  obtain ⟨b', c', d'⟩ := y'
  simp only [π3_apply, Prod.mk.injEq] at h
  obtain ⟨hb, hc, hd⟩ := h
  obtain ⟨jb, hjb⟩ := eq_of_π_eq_zero m (z := b - b') (by rw [map_sub, hb, sub_self])
  obtain ⟨jc, hjc⟩ := eq_of_π_eq_zero m (z := c - c') (by rw [map_sub, hc, sub_self])
  obtain ⟨jd, hjd⟩ := eq_of_π_eq_zero m (z := d - d') (by rw [map_sub, hd, sub_self])
  rw [sub_eq_iff_eq_add] at hjb hjc hjd
  subst hjb hjc hjd
  simp only [Q]
  linear_combination (c' * jc + 2 ^ (m + 2) * jc ^ 2 + d' * jd + 2 ^ (m + 2) * jd ^ 2
    - b' * jb - 2 ^ (m + 2) * jb ^ 2) * two_pow_self (m + 4)

/-- Every fibre of the reduction map on triples has 8 elements. -/
theorem card_fiber (x : ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3))) :
    (univ.filter (fun y => π3 m y = x)).card = 8 := by
  have hsurj : ∀ z, z ∈ Set.range (π3 m) := fun z => ⟨σ3 m z, π3_σ3 m z⟩
  have hall : ∀ z, (univ.filter (fun y => π3 m y = z)).card =
      (univ.filter (fun y => π3 m y = x)).card :=
    fun z => AddMonoidHom.card_fiber_eq_of_mem_range (π3 m) (hsurj z) (hsurj x)
  have hsum := card_eq_sum_card_fiberwise (f := π3 m) (s := univ) (t := univ)
    (fun _ _ => mem_univ _)
  simp only [hall, sum_const, card_univ, Fintype.card_prod, ZMod.card, smul_eq_mul] at hsum
  have hpos : 0 < (2 ^ (m + 3)) * ((2 ^ (m + 3)) * (2 ^ (m + 3))) := by positivity
  have h8 : 2 ^ (m + 4) * (2 ^ (m + 4) * 2 ^ (m + 4)) =
      8 * ((2 ^ (m + 3)) * ((2 ^ (m + 3)) * (2 ^ (m + 3)))) := by ring
  rw [h8] at hsum
  exact (Nat.eq_of_mul_eq_mul_left hpos (by linarith)).symm

/-- Solutions mod `2^(m+3)` whose canonical lift is still a solution mod `2^(m+4)`. -/
def C1 : Finset (ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3))) :=
  univ.filter (fun x => Q (σ3 m x) = 1)

theorem card_S_succ : (S (m + 4)).card = 8 * (C1 m).card := by
  rw [card_eq_sum_card_fiberwise (f := π3 m) (t := univ) (fun _ _ => mem_univ _)]
  have hterm : ∀ x, ((S (m + 4)).filter (fun y => π3 m y = x)).card =
      if Q (σ3 m x) = 1 then 8 else 0 := by
    intro x
    have hQ : ∀ y, π3 m y = x → Q y = Q (σ3 m x) := fun y hy =>
      Q_congr m (by rw [hy, π3_σ3])
    split_ifs with hx
    · rw [← card_fiber m x]
      congr 1
      ext y
      simp only [S, mem_filter, mem_univ, true_and, and_iff_right_iff_imp]
      intro hy; rw [hQ y hy, hx]
    · rw [card_eq_zero, filter_eq_empty_iff]
      intro y hy hy'
      simp only [S, mem_filter, mem_univ, true_and] at hy
      exact hx (by rw [← hQ y hy', hy])
  simp only [hterm]
  rw [C1, card_filter, mul_sum]
  congr 1; ext x; split_ifs <;> simp

/-! #### An involution on `S (m+3)` swapping lifts that are / are not solutions -/

/-- `ρ` on the two levels. -/
abbrev ρ3 := ρ (m + 3) (by omega)
abbrev ρ4 := ρ (m + 4) (by omega)

theorem ρ4_σ (x : ZMod (2 ^ (m + 3))) : ρ4 m (σ m x) = ρ3 m x := by
  conv_rhs => rw [← ZMod.natCast_zmod_val x]
  unfold σ; rw [ρ_natCast, ρ_natCast]

theorem ρ3_hh : ρ3 m (hh m) = 0 := by
  unfold hh; rw [map_pow, map_ofNat, show (2 : ZMod 2) = 0 by decide, zero_pow (by omega)]

/-- Add `h = 2^(m+2)` to the first odd coordinate. -/
def φ (x : ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3))) :
    ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) :=
  if ρ3 m x.1 = 1 then (x.1 + hh m, x.2.1, x.2.2)
  else if ρ3 m x.2.1 = 1 then (x.1, x.2.1 + hh m, x.2.2)
  else (x.1, x.2.1, x.2.2 + hh m)

theorem φ_φ (x) : φ m (φ m x) = x := by
  obtain ⟨b, c, d⟩ := x
  have e : ∀ z : ZMod (2 ^ (m + 3)), z + hh m + hh m = z := fun z => by
    linear_combination two_mul_hh m
  unfold φ
  by_cases hb : ρ3 m b = 1
  · simp [hb, map_add, ρ3_hh, e]
  · by_cases hc : ρ3 m c = 1
    · simp [hb, hc, map_add, ρ3_hh, e]
    · simp [hb, hc, e]

theorem ρ_d_of_mem_S {x} (hx : x ∈ S (m + 3)) (hb : ¬ ρ3 m x.1 = 1) (hc : ¬ ρ3 m x.2.1 = 1) :
    ρ3 m x.2.2 = 1 := by
  have hx' : Q x = 1 := (mem_filter.1 hx).2
  unfold Q at hx'
  have h := congrArg (ρ3 m) hx'
  rw [map_sub, map_add, map_pow, map_pow, map_pow, map_one] at h
  have hb0 : ρ3 m x.1 = 0 := (zmod2_cases _).resolve_right hb
  have hc0 : ρ3 m x.2.1 = 0 := (zmod2_cases _).resolve_right hc
  rw [hb0, hc0] at h
  rcases zmod2_cases (ρ3 m x.2.2) with hd | hd
  · rw [hd] at h; exact absurd h (by decide)
  · exact hd

theorem φ_mem_S {x} (hx : x ∈ S (m + 3)) : φ m x ∈ S (m + 3) := by
  obtain ⟨b, c, d⟩ := x
  simp only [S, mem_filter, mem_univ, true_and] at hx ⊢
  unfold φ
  split_ifs
  · simp only [Q] at hx ⊢; linear_combination hx - b * two_mul_hh m - hh_sq m
  · simp only [Q] at hx ⊢; linear_combination hx + c * two_mul_hh m + hh_sq m
  · simp only [Q] at hx ⊢; linear_combination hx + d * two_mul_hh m + hh_sq m

/-- The lift of `φ x` changes `Q` by exactly `2^(m+3)` modulo `2^(m+4)`. -/
theorem Q_σ3_φ {x} (hx : x ∈ S (m + 3)) :
    Q (σ3 m (φ m x)) = Q (σ3 m x) + 2 ^ (m + 3) := by
  obtain ⟨b, c, d⟩ := x
  have hH : π m ((2 : ZMod (2 ^ (m + 4))) ^ (m + 2)) = hh m := by
    unfold π hh; rw [map_pow, map_ofNat]
  unfold φ
  split_ifs with hb hc
  · -- first coordinate odd
    have hq : Q (σ3 m (b + hh m, c, d)) = Q (σ m b + 2 ^ (m + 2), σ m c, σ m d) :=
      Q_congr m (by simp [π3_apply, σ3, π_σ, map_add, hH])
    rw [hq]
    obtain ⟨w, hw⟩ := eq_two_mul_add_one_of_ρ_eq_one (n := m + 4) (by omega)
      (show ρ4 m (σ m b) = 1 by rw [ρ4_σ]; exact hb)
    simp only [σ3, Q, hw]
    linear_combination (-(w + 1 + 2 ^ m)) * two_pow_self (m + 4)
  · -- second coordinate odd
    have hq : Q (σ3 m (b, c + hh m, d)) = Q (σ m b, σ m c + 2 ^ (m + 2), σ m d) :=
      Q_congr m (by simp [π3_apply, σ3, π_σ, map_add, hH])
    rw [hq]
    obtain ⟨w, hw⟩ := eq_two_mul_add_one_of_ρ_eq_one (n := m + 4) (by omega)
      (show ρ4 m (σ m c) = 1 by rw [ρ4_σ]; exact hc)
    simp only [σ3, Q, hw]
    linear_combination (w + 2 ^ m) * two_pow_self (m + 4)
  · -- third coordinate odd
    have hd := ρ_d_of_mem_S m hx hb hc
    have hq : Q (σ3 m (b, c, d + hh m)) = Q (σ m b, σ m c, σ m d + 2 ^ (m + 2)) :=
      Q_congr m (by simp [π3_apply, σ3, π_σ, map_add, hH])
    rw [hq]
    obtain ⟨w, hw⟩ := eq_two_mul_add_one_of_ρ_eq_one (n := m + 4) (by omega)
      (show ρ4 m (σ m d) = 1 by rw [ρ4_σ]; exact hd)
    simp only [σ3, Q, hw]
    linear_combination (w + 2 ^ m) * two_pow_self (m + 4)

theorem Q_σ3_cases {x} (hx : x ∈ S (m + 3)) :
    Q (σ3 m x) = 1 ∨ Q (σ3 m x) = 1 + 2 ^ (m + 3) := by
  simp only [S, mem_filter, mem_univ, true_and] at hx
  have h1 : π m (Q (σ3 m x) - 1) = 0 := by
    rw [map_sub, map_one, ← Q_map]
    have := π3_σ3 m x
    rw [π3_apply] at this
    rw [this, hx, sub_self]
  obtain ⟨j, hj⟩ := eq_of_π_eq_zero m h1
  rcases half_mul_mem (n := m + 4) (by omega) j with h | h <;>
    simp only [show m + 4 - 1 = m + 3 from rfl] at h
  · left; linear_combination hj + h
  · right; linear_combination hj + h

theorem two_pow_ne_zero' : (2 : ZMod (2 ^ (m + 4))) ^ (m + 3) ≠ 0 := by
  have := natCast_ne_zero_of_lt (n := m + 4) (k := 2 ^ (m + 3)) (by positivity)
    (Nat.pow_lt_pow_right (by norm_num) (by omega))
  push_cast at this; exact this

theorem C1_subset : C1 m ⊆ S (m + 3) := by
  intro x hx
  simp only [C1, mem_filter, mem_univ, true_and] at hx
  simp only [S, mem_filter, mem_univ, true_and]
  have := congrArg (π m) hx
  rw [← Q_map, map_one] at this
  have h' := π3_σ3 m x
  rw [π3_apply] at h'
  rwa [h'] at this

theorem two_mul_card_C1 : 2 * (C1 m).card = (S (m + 3)).card := by
  have hbij : (C1 m).card = (S (m + 3) \ C1 m).card := by
    apply card_bij (fun x _ => φ m x)
    · intro x hx
      have hxS := C1_subset m hx
      simp only [C1, mem_filter, mem_univ, true_and] at hx
      simp only [mem_sdiff, C1, mem_filter, mem_univ, true_and]
      refine ⟨φ_mem_S m hxS, ?_⟩
      rw [Q_σ3_φ m hxS, hx]
      intro h
      exact two_pow_ne_zero' m (by linear_combination h)
    · intro a _ b _ hab
      rw [← φ_φ m a, ← φ_φ m b, hab]
    · intro y hy
      simp only [mem_sdiff, C1, mem_filter, mem_univ, true_and] at hy
      obtain ⟨hyS, hy1⟩ := hy
      refine ⟨φ m y, ?_, φ_φ m y⟩
      simp only [C1, mem_filter, mem_univ, true_and]
      rw [Q_σ3_φ m hyS, (Q_σ3_cases m hyS).resolve_left hy1]
      linear_combination two_pow_self (m + 4)
  have := card_sdiff_add_card_eq_card (C1_subset m)
  omega

theorem card_S_rec : (S (m + 4)).card = 4 * (S (m + 3)).card := by
  rw [card_S_succ, ← two_mul_card_C1]; ring

end lift

set_option maxRecDepth 100000 in
theorem card_S_three : (S 3).card = 128 := by decide +kernel

theorem card_S (m : ℕ) : (S (m + 3)).card = 2 ^ (2 * m + 7) := by
  induction m with
  | zero => exact card_S_three
  | succ k ih =>
    rw [show k + 1 + 3 = k + 4 by omega, card_S_rec, ih]; ring

/-! ### The main count -/

/-- The involutions of the quaternion ring over `ZMod (2^n)` with `i^2 = j^2 = 1`, written as
solutions `(a, b, c, d)` of the system defining A236553 / A236554. -/
def T (n : ℕ) : Finset (ZMod (2 ^ n) × ZMod (2 ^ n) × ZMod (2 ^ n) × ZMod (2 ^ n)) :=
  univ.filter (fun x => x.1 ^ 2 - x.2.1 ^ 2 + x.2.2.1 ^ 2 + x.2.2.2 ^ 2 = 1 ∧
    2 * x.1 * x.2.1 = 0 ∧ 2 * x.1 * x.2.2.1 = 0 ∧ 2 * x.1 * x.2.2.2 = 0)

section main

variable (m : ℕ)

/-- `{0, 2^(n-1)}`, the solutions of `2z = 0`. -/
def H : Finset (ZMod (2 ^ (m + 3))) := {0, hh m}

theorem mem_H (z : ZMod (2 ^ (m + 3))) : z ∈ H m ↔ z = 0 ∨ z = hh m := by simp [H]

theorem sq_of_mem_H {z : ZMod (2 ^ (m + 3))} (hz : z ∈ H m) : z ^ 2 = 0 := by
  rcases (mem_H m z).1 hz with rfl | rfl
  · ring
  · exact hh_sq m

theorem two_mul_of_mem_H {z : ZMod (2 ^ (m + 3))} (hz : z ∈ H m) : 2 * z = 0 := by
  rcases (mem_H m z).1 hz with rfl | rfl
  · ring
  · exact two_mul_hh m

theorem hh_ne_zero : hh m ≠ 0 := by
  have := natCast_ne_zero_of_lt (n := m + 3) (k := 2 ^ (m + 2)) (by positivity)
    (Nat.pow_lt_pow_right (by norm_num) (by omega))
  push_cast at this; exact this

theorem card_H : (H m).card = 2 := card_pair (hh_ne_zero m).symm

theorem T_eq : T (m + 3) = (univ.filter (fun a : ZMod (2 ^ (m + 3)) => a ^ 2 = 1) ×ˢ
    (H m ×ˢ H m ×ˢ H m)) ∪ (H m ×ˢ S (m + 3)) := by
  have hn : m + 3 ≠ 0 := by omega
  ext ⟨a, b, c, d⟩
  simp only [T, mem_filter, mem_univ, true_and, mem_union, mem_product, mem_S]
  constructor
  · rintro ⟨heq, hb, hc, hd⟩
    rcases zmod2_cases (ρ3 m a) with ha | ha
    · -- `a` even
      by_cases h2a : 2 * a = 0
      · right
        have haH : a ∈ H m := (mem_H m a).2 ((two_mul_eq_zero_iff hn a).1 h2a)
        exact ⟨haH, by linear_combination heq - sq_of_mem_H m haH⟩
      · -- then `b, c, d` are even too, and the equation fails modulo 4
        exfalso
        have even_of : ∀ z, 2 * a * z = 0 → ρ3 m z = 0 := by
          intro z hz
          rcases zmod2_cases (ρ3 m z) with h | h
          · exact h
          · exact absurd ((isUnit_of_ρ_eq_one hn h).mul_right_eq_zero.1
              (by linear_combination hz)) h2a
        obtain ⟨b', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hb)
        obtain ⟨c', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hc)
        obtain ⟨d', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hd)
        obtain ⟨a', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn ha
        have h1 : (1 : ZMod (2 ^ (m + 3))) = 2 * (2 * (a' ^ 2 - b' ^ 2 + c' ^ 2 + d' ^ 2)) := by
          linear_combination -heq
        have := congrArg (ρ3 m) h1
        rw [map_one, map_mul, map_ofNat, show (2 : ZMod 2) = 0 by decide, zero_mul] at this
        exact absurd this (by decide)
    · -- `a` odd, hence a unit, so `b, c, d ∈ {0, h}`
      left
      have hu := isUnit_of_ρ_eq_one hn ha
      have hH : ∀ z, 2 * a * z = 0 → z ∈ H m := fun z hz =>
        (mem_H m z).2 ((two_mul_eq_zero_iff hn z).1
          (hu.mul_right_eq_zero.1 (by linear_combination hz)))
      have hbH := hH b hb
      have hcH := hH c hc
      have hdH := hH d hd
      exact ⟨by linear_combination heq + sq_of_mem_H m hbH - sq_of_mem_H m hcH
        - sq_of_mem_H m hdH, hbH, hcH, hdH⟩
  · rintro (⟨ha, hb, hc, hd⟩ | ⟨ha, hQ⟩)
    · exact ⟨by linear_combination ha - sq_of_mem_H m hb + sq_of_mem_H m hc + sq_of_mem_H m hd,
        by linear_combination a * two_mul_of_mem_H m hb,
        by linear_combination a * two_mul_of_mem_H m hc,
        by linear_combination a * two_mul_of_mem_H m hd⟩
    · exact ⟨by linear_combination hQ + sq_of_mem_H m ha,
        by linear_combination b * two_mul_of_mem_H m ha,
        by linear_combination c * two_mul_of_mem_H m ha,
        by linear_combination d * two_mul_of_mem_H m ha⟩

theorem card_T : (T (m + 3)).card = 2 ^ (2 * m + 8) + 32 := by
  have hdisj : Disjoint (univ.filter (fun a : ZMod (2 ^ (m + 3)) => a ^ 2 = 1) ×ˢ
      (H m ×ˢ H m ×ˢ H m)) (H m ×ˢ S (m + 3)) := by
    rw [disjoint_left]
    rintro ⟨a, x⟩ h1 h2
    simp only [mem_product, mem_filter, mem_univ, true_and] at h1 h2
    have h10 : ((1 : ℕ) : ZMod (2 ^ (m + 3))) ≠ 0 :=
      natCast_ne_zero_of_lt (by norm_num) (Nat.one_lt_two_pow (by omega))
    rw [Nat.cast_one] at h10
    exact h10 (by rw [← h1.1, sq_of_mem_H m h2.1])
  rw [T_eq, card_union_of_disjoint hdisj, card_product, card_product, card_product,
    card_product, card_roots, card_H, card_S]
  ring

end main

/-- **A236554**: the number of involutions in the quaternion ring over `ZMod (2^n)` with
`i^2 = j^2 = 1` is `2^(2n+2) + 32` for every `n ≥ 3`. -/
theorem card_A236554 (n : ℕ) (hn : 3 ≤ n) : (T n).card = 2 ^ (2 * n + 2) + 32 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 3 := ⟨n - 3, by omega⟩
  rw [card_T, show 2 * (m + 3) + 2 = 2 * m + 8 by ring]

/-- The two initial terms, which the formula does not cover. -/
theorem card_T_one : (T 1).card = 8 := by decide +kernel

set_option maxRecDepth 100000 in
theorem card_T_two : (T 2).card = 64 := by decide +kernel

end A236554
