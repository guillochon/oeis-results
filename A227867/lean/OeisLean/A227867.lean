import OeisLean.A236553
import Mathlib.FieldTheory.Finite.Basic

/-!
# A227867: Lipschitz quaternions `X` with `X^2 ≡ 1 (mod n)`

A227867(n) is the number of solutions `(a, b, c, d)` modulo `n` of
  `a^2 - b^2 - c^2 - d^2 = 1`, `2ab = 0`, `2ac = 0`, `2ad = 0`
(the involutions of the Lipschitz quaternions `i^2 = j^2 = -1` over `ZMod n`).

Main results:
* `A_mul`: `A` is multiplicative;
* `A_odd_prime_pow`: `A (p^k) = p^(2k-1) (p+1) + 2` for odd primes `p` and `k ≥ 1`
  (the same values as A236553);
* `A_two`, `A_four`, `A_two_pow`: `A 2 = 8`, and `A (2^k) = 32` for `k ≥ 2`
  (a sum of three squares is never `7 mod 8`);
* `A_closed_form`: `A n = ∏_{p^k ∥ n} localFactor p k` for `n ≠ 0`.

The odd-prime base count (points on the sphere `b^2 + c^2 + d^2 = -1` mod `p`) is reduced to the
A236553 count by writing `-1 = s^2 + t^2` and changing variables.
-/

open Finset

namespace A227867

open A236553 (ρ ρ_natCast two_ne_zero_p isUnit_iff_ρ two_isUnit card_roots_p p_pow_self
  πp π3p σp liftp ρ1 ρ2 πp_σp ρ2_σp ρ2_liftp πp_P1 eq_of_πp_eq_zero P1_mul_eq_zero_iff
  card_SQ_prime mem_SQ map4)

/-- The defining system, over any commutative ring. -/
abbrev P {R : Type*} [CommRing R] (x : R × R × R × R) : Prop :=
  x.1 ^ 2 - x.2.1 ^ 2 - x.2.2.1 ^ 2 - x.2.2.2 ^ 2 = 1 ∧
    2 * x.1 * x.2.1 = 0 ∧ 2 * x.1 * x.2.2.1 = 0 ∧ 2 * x.1 * x.2.2.2 = 0

/-- Number of solutions of the system over a ring `R`. -/
noncomputable def cnt (R : Type*) [CommRing R] : ℕ := Nat.card {x : R × R × R × R // P x}

/-- **A227867(n)**. -/
noncomputable def A (n : ℕ) : ℕ := cnt (ZMod n)

/-! ### Multiplicativity -/

section general

variable {R S : Type*} [CommRing R] [CommRing S]

theorem P_map_iff (f : R →+* S) (hf : Function.Injective f) (x : R × R × R × R) :
    P (map4 f x) ↔ P x := by
  obtain ⟨a, b, c, d⟩ := x
  simp only [P, map4]
  constructor
  · rintro ⟨h1, h2, h3, h4⟩
    refine ⟨hf ?_, hf ?_, hf ?_, hf ?_⟩ <;> simpa [map_ofNat]
  · rintro ⟨h1, h2, h3, h4⟩
    refine ⟨?_, ?_, ?_, ?_⟩
    · simpa [map_ofNat] using congrArg f h1
    · simpa [map_ofNat] using congrArg f h2
    · simpa [map_ofNat] using congrArg f h3
    · simpa [map_ofNat] using congrArg f h4

theorem cnt_congr (e : R ≃+* S) : cnt R = cnt S := by
  unfold cnt
  apply Nat.card_congr
  refine Equiv.subtypeEquiv
    (e.toEquiv.prodCongr (e.toEquiv.prodCongr (e.toEquiv.prodCongr e.toEquiv))) ?_
  intro x
  exact (P_map_iff e.toRingHom e.injective x).symm

theorem cnt_prod : cnt (R × S) = cnt R * cnt S := by
  unfold cnt
  rw [← Nat.card_prod]
  apply Nat.card_congr
  exact {
    toFun := fun x => (⟨map4 (RingHom.fst R S) x.1, by
        obtain ⟨⟨a, b, c, d⟩, h⟩ := x
        simp only [P, map4, Prod.ext_iff] at h ⊢
        simp at h ⊢; tauto⟩,
      ⟨map4 (RingHom.snd R S) x.1, by
        obtain ⟨⟨a, b, c, d⟩, h⟩ := x
        simp only [P, map4, Prod.ext_iff] at h ⊢
        simp at h ⊢; tauto⟩)
    invFun := fun y => ⟨((y.1.1.1, y.2.1.1), (y.1.1.2.1, y.2.1.2.1), (y.1.1.2.2.1, y.2.1.2.2.1),
        (y.1.1.2.2.2, y.2.1.2.2.2)), by
        obtain ⟨⟨⟨a, b, c, d⟩, h⟩, ⟨⟨a', b', c', d'⟩, h'⟩⟩ := y
        simp only [P, Prod.ext_iff] at h h' ⊢
        simp; tauto⟩
    left_inv := fun x => by obtain ⟨⟨a, b, c, d⟩, h⟩ := x; rfl
    right_inv := fun y => by obtain ⟨⟨⟨a, b, c, d⟩, h⟩, ⟨⟨a', b', c', d'⟩, h'⟩⟩ := y; rfl }

end general

theorem A_mul {m n : ℕ} (h : m.Coprime n) : A (m * n) = A m * A n := by
  unfold A; rw [cnt_congr (ZMod.chineseRemainder h), cnt_prod]

theorem A_one : A 1 = 1 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide

/-! ### The negative definite form `Q'(b, c, d) = -b^2 - c^2 - d^2` -/

/-- `Q' (b, c, d) = -b^2 - c^2 - d^2`. -/
def Q' {α : Type*} [CommRing α] (x : α × α × α) : α := -x.1 ^ 2 - x.2.1 ^ 2 - x.2.2 ^ 2

theorem Q'_map {α β : Type*} [CommRing α] [CommRing β] (f : α →+* β) (x : α × α × α) :
    Q' (f x.1, f x.2.1, f x.2.2) = f (Q' x) := by
  simp [Q', map_sub, map_pow, map_neg]

/-- Solutions of `Q' = 1` in a finite ring. -/
def SQ' (R : Type*) [CommRing R] [Fintype R] [DecidableEq R] : Finset (R × R × R) :=
  univ.filter (fun y => Q' y = 1)

theorem mem_SQ' {R : Type*} [CommRing R] [Fintype R] [DecidableEq R] (y : R × R × R) :
    y ∈ SQ' R ↔ -y.1 ^ 2 - y.2.1 ^ 2 - y.2.2 ^ 2 = 1 := by
  unfold SQ'; rw [mem_filter]; simp [Q']

theorem card_SQ'_congr {R S : Type*} [CommRing R] [Fintype R] [DecidableEq R] [CommRing S]
    [Fintype S] [DecidableEq S] (e : R ≃+* S) : (SQ' R).card = (SQ' S).card := by
  apply card_nbij' (fun y => (e y.1, e y.2.1, e y.2.2)) (fun z => (e.symm z.1, e.symm z.2.1,
    e.symm z.2.2))
  · intro y hy
    rw [mem_coe, mem_SQ'] at hy ⊢
    have := congrArg e hy
    simpa [map_sub, map_neg, map_pow] using this
  · intro z hz
    rw [mem_coe, mem_SQ'] at hz ⊢
    have := congrArg e.symm hz
    simpa [map_sub, map_neg, map_pow] using this
  · intro y _; simp
  · intro z _; simp

/-! ### Odd prime powers -/

section oddprime

variable (p : ℕ) [hp : Fact p.Prime]

/-- For odd `p`, a solution has either `a = ±1, b = c = d = 0`, or `a = 0` and
`Q' (b,c,d) = 1`. -/
theorem filter_P_eq (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) :
    univ.filter (fun x : ZMod (p ^ e) × ZMod (p ^ e) × ZMod (p ^ e) × ZMod (p ^ e) => P x) =
      (univ.filter (fun a : ZMod (p ^ e) => a ^ 2 = 1)) ×ˢ {(0, 0, 0)} ∪
        {0} ×ˢ SQ' (ZMod (p ^ e)) := by
  ext ⟨a, b, c, d⟩
  simp only [mem_filter, mem_univ, true_and, mem_union, mem_product, mem_singleton, mem_SQ',
    Prod.mk.injEq]
  constructor
  · rintro ⟨heq, hb, hc, hd⟩
    by_cases ha : ρ p e he a = 0
    · by_cases ha0 : a = 0
      · subst ha0; right; exact ⟨rfl, by linear_combination heq⟩
      · exfalso
        have hnu : ∀ z, 2 * a * z = 0 → ρ p e he z = 0 := by
          intro z hz
          by_contra hz'
          have h2a : 2 * a = 0 := ((isUnit_iff_ρ p he z).2 hz').mul_left_eq_zero.1 hz
          exact ha0 ((two_isUnit p hp2 he).mul_right_eq_zero.1 h2a)
        have := congrArg (ρ p e he) heq
        simp [ha, hnu b hb, hnu c hc, hnu d hd] at this
    · have h2a : IsUnit (2 * a) := (two_isUnit p hp2 he).mul ((isUnit_iff_ρ p he a).2 ha)
      have hb0 := h2a.mul_right_eq_zero.1 hb
      have hc0 := h2a.mul_right_eq_zero.1 hc
      have hd0 := h2a.mul_right_eq_zero.1 hd
      subst hb0 hc0 hd0
      left; exact ⟨by linear_combination heq, rfl, rfl, rfl⟩
  · rintro (⟨ha, rfl, rfl, rfl⟩ | ⟨rfl, hQ⟩)
    · exact ⟨by linear_combination ha, by ring, by ring, by ring⟩
    · exact ⟨by linear_combination hQ, by ring, by ring, by ring⟩

theorem cnt_eq (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) :
    cnt (ZMod (p ^ e)) = 2 + (SQ' (ZMod (p ^ e))).card := by
  unfold cnt
  rw [Nat.card_eq_fintype_card, Fintype.card_subtype, filter_P_eq p hp2 he,
    card_union_of_disjoint, card_product, card_product, card_roots_p p hp2 he, card_singleton,
    card_singleton]
  · ring
  · rw [disjoint_left]
    rintro ⟨a, x⟩ h1 h2
    simp only [mem_product, mem_filter, mem_univ, true_and, mem_singleton] at h1 h2
    have h := h1.1
    rw [h2.1] at h
    have := congrArg (ρ p e he) h
    simp at this

end oddprime

section hensel

variable (p : ℕ) [hp : Fact p.Prime] (f : ℕ)

omit hp in
theorem Q'_add (u j : ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2))) :
    Q' (u.1 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.1,
      u.2.1 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.2.1,
      u.2.2 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.2.2) =
    Q' u + (p : ZMod (p ^ (f + 2))) ^ (f + 1) *
      (-2 * (u.1 * j.1 + u.2.1 * j.2.1 + u.2.2 * j.2.2)) := by
  simp only [Q']
  linear_combination -(j.1 ^ 2 + j.2.1 ^ 2 + j.2.2 ^ 2) * (p : ZMod (p ^ (f + 2))) ^ f *
    p_pow_self p (f + 2)

/-- Every solution of `Q' = 1` modulo `p^(f+1)` lifts to exactly `p^2` solutions modulo
`p^(f+2)`. -/
theorem card_fiber_SQ' (hp2 : p ≠ 2) {x} (hx : x ∈ SQ' (ZMod (p ^ (f + 1)))) :
    ((SQ' (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card = p ^ 2 := by
  set P1 : ZMod (p ^ (f + 2)) := (p : ZMod (p ^ (f + 2))) ^ (f + 1) with hP1
  set u : ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) :=
    (σp p f x.1, σp p f x.2.1, σp p f x.2.2) with hu
  have hx' : Q' x = 1 := by rw [mem_SQ'] at hx; simpa [Q'] using hx
  have hQu : πp p f (Q' u - 1) = 0 := by
    rw [map_sub, map_one, ← Q'_map]
    simp only [hu, πp_σp]
    rw [hx', sub_self]
  obtain ⟨e0, he0⟩ := eq_of_πp_eq_zero p f hQu
  set xb := ρ1 p f x.1
  set xc := ρ1 p f x.2.1
  set xd := ρ1 p f x.2.2
  have hQbar : -xb ^ 2 - xc ^ 2 - xd ^ 2 = 1 := by
    have := congrArg (ρ1 p f) hx'
    simpa [Q', map_sub, map_neg, map_pow] using this
  let ℓ : (ZMod p × ZMod p × ZMod p) →+ ZMod p :=
    AddMonoidHom.mk' (fun t => -2 * (xb * t.1 + xc * t.2.1 + xd * t.2.2))
      (by intro a b; simp only [Prod.fst_add, Prod.snd_add]; ring)
  have hℓ : ∀ t, ℓ t = -2 * (xb * t.1 + xc * t.2.1 + xd * t.2.2) := fun t => rfl
  set β := -ρ2 p f e0
  let τ : ZMod p × ZMod p × ZMod p →
      ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) :=
    fun t => (u.1 + P1 * liftp p f t.1, u.2.1 + P1 * liftp p f t.2.1,
      u.2.2 + P1 * liftp p f t.2.2)
  have key : ∀ t, Q' (τ t) = 1 ↔ ℓ t = β := by
    intro t
    have hexp := Q'_add p f u (liftp p f t.1, liftp p f t.2.1, liftp p f t.2.2)
    simp only at hexp
    have h1 : Q' (τ t) = 1 ↔ P1 * (e0 + -2 * (u.1 * liftp p f t.1 +
        u.2.1 * liftp p f t.2.1 + u.2.2 * liftp p f t.2.2)) = 0 := by
      change Q' (u.1 + P1 * liftp p f t.1, u.2.1 + P1 * liftp p f t.2.1,
        u.2.2 + P1 * liftp p f t.2.2) = 1 ↔ _
      rw [hexp]
      constructor
      · intro h; linear_combination h - he0
      · intro h; linear_combination h + he0
    rw [h1, P1_mul_eq_zero_iff, hℓ]
    simp only [map_add, map_mul, map_neg, map_ofNat, hu, ρ2_σp, ρ2_liftp]
    constructor
    · intro h; linear_combination h
    · intro h; linear_combination h
  have hbij : (univ.filter (fun t => ℓ t = β)).card =
      ((SQ' (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card := by
    apply card_bij (fun t _ => τ t)
    · intro t ht
      simp only [mem_filter, mem_univ, true_and] at ht
      rw [mem_filter, mem_SQ']
      refine ⟨?_, ?_⟩
      · have := (key t).2 ht; simpa [Q'] using this
      · simp only [π3p, τ, hu, map_add, map_mul, πp_σp, hP1, πp_P1, zero_mul, add_zero]
    · intro t _ t' _ htt
      simp only [τ, Prod.mk.injEq] at htt
      obtain ⟨h1, h2, h3⟩ := htt
      have e1 := (P1_mul_eq_zero_iff p f (liftp p f t.1 - liftp p f t'.1)).1
        (by linear_combination h1)
      have e2 := (P1_mul_eq_zero_iff p f (liftp p f t.2.1 - liftp p f t'.2.1)).1
        (by linear_combination h2)
      have e3 := (P1_mul_eq_zero_iff p f (liftp p f t.2.2 - liftp p f t'.2.2)).1
        (by linear_combination h3)
      rw [map_sub, ρ2_liftp, ρ2_liftp, sub_eq_zero] at e1 e2 e3
      exact Prod.ext e1 (Prod.ext e2 e3)
    · intro y hy
      rw [mem_filter, mem_SQ'] at hy
      obtain ⟨hyQ, hyx⟩ := hy
      simp only [π3p, Prod.ext_iff] at hyx
      obtain ⟨hb, hc, hd⟩ := hyx
      obtain ⟨jb, hjb⟩ := eq_of_πp_eq_zero p f (z := y.1 - u.1)
        (by rw [map_sub, hb, hu, πp_σp, sub_self])
      obtain ⟨jc, hjc⟩ := eq_of_πp_eq_zero p f (z := y.2.1 - u.2.1)
        (by rw [map_sub, hc, hu, πp_σp, sub_self])
      obtain ⟨jd, hjd⟩ := eq_of_πp_eq_zero p f (z := y.2.2 - u.2.2)
        (by rw [map_sub, hd, hu, πp_σp, sub_self])
      have hfix : ∀ j : ZMod (p ^ (f + 2)), P1 * liftp p f (ρ2 p f j) = P1 * j := by
        intro j
        have := (P1_mul_eq_zero_iff p f (liftp p f (ρ2 p f j) - j)).2
          (by rw [map_sub, ρ2_liftp, sub_self])
        linear_combination this
      have hτ : τ (ρ2 p f jb, ρ2 p f jc, ρ2 p f jd) = y := by
        simp only [τ, hfix]
        ext
        · linear_combination -hjb
        · linear_combination -hjc
        · linear_combination -hjd
      refine ⟨(ρ2 p f jb, ρ2 p f jc, ρ2 p f jd), ?_, hτ⟩
      simp only [mem_filter, mem_univ, true_and]
      rw [← key, hτ]
      simpa [Q'] using hyQ
  rw [← hbij]
  have h2 := two_ne_zero_p p hp2
  have hsurj : ∀ r, r ∈ Set.range ℓ := by
    intro r
    by_cases hb : xb = 0
    · by_cases hc : xc = 0
      · have hd : xd ≠ 0 := by
          intro hd; rw [hb, hc, hd] at hQbar; simp at hQbar
        refine ⟨(0, 0, -r / (2 * xd)), ?_⟩
        rw [hℓ]; field_simp; ring
      · refine ⟨(0, -r / (2 * xc), 0), ?_⟩
        rw [hℓ]; field_simp; ring
    · refine ⟨(-r / (2 * xb), 0, 0), ?_⟩
      rw [hℓ]; field_simp; ring
  have hall : ∀ r, (univ.filter (fun t => ℓ t = r)).card =
      (univ.filter (fun t => ℓ t = β)).card :=
    fun r => AddMonoidHom.card_fiber_eq_of_mem_range ℓ (hsurj r) (hsurj β)
  have hsum := card_eq_sum_card_fiberwise (f := ℓ) (s := univ) (t := univ)
    (fun _ _ => mem_univ _)
  simp only [hall, sum_const, card_univ, Fintype.card_prod, ZMod.card, smul_eq_mul] at hsum
  have hpos : 0 < p := hp.out.pos
  have : p * (p * p) = p * p ^ 2 := by ring
  rw [this] at hsum
  exact (Nat.eq_of_mul_eq_mul_left hpos (by linarith)).symm

theorem card_SQ'_succ (hp2 : p ≠ 2) :
    (SQ' (ZMod (p ^ (f + 2)))).card = p ^ 2 * (SQ' (ZMod (p ^ (f + 1)))).card := by
  rw [card_eq_sum_card_fiberwise (f := π3p p f) (t := univ) (fun _ _ => mem_univ _)]
  have hterm : ∀ x, ((SQ' (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card =
      if x ∈ SQ' (ZMod (p ^ (f + 1))) then p ^ 2 else 0 := by
    intro x
    split_ifs with hx
    · exact card_fiber_SQ' p f hp2 hx
    · rw [card_eq_zero, filter_eq_empty_iff]
      intro y hy hyx
      apply hx
      rw [mem_SQ'] at hy ⊢
      have := congrArg (πp p f) hy
      rw [← hyx]
      simpa [π3p, map_sub, map_neg, map_pow] using this
  simp only [hterm]
  rw [sum_ite_mem, univ_inter, sum_const, smul_eq_mul, mul_comm]

end hensel

section base

variable (p : ℕ) [hp : Fact p.Prime]

/-- Points on the sphere `b^2 + c^2 + d^2 = -1` mod `p`: write `-1 = s^2 + t^2`; then
`(b, c, d) ↦ (d, s b + t c, t b - s c)` turns `-b^2 - c^2 - d^2` into `c^2 + d^2 - b^2`. -/
theorem card_SQ'_prime (hp2 : p ≠ 2) : (SQ' (ZMod p)).card = p ^ 2 + p := by
  obtain ⟨s, t, hst⟩ := ZMod.sq_add_sq p (-1)
  rw [← card_SQ_prime p hp2]
  apply card_nbij' (fun y => (y.2.2, s * y.1 + t * y.2.1, t * y.1 - s * y.2.1))
    (fun z => (-s * z.2.1 - t * z.2.2, -t * z.2.1 + s * z.2.2, z.1))
  · intro y hy
    rw [mem_coe, mem_SQ'] at hy
    rw [mem_coe, mem_SQ]
    linear_combination hy + (y.1 ^ 2 + y.2.1 ^ 2) * hst
  · intro z hz
    rw [mem_coe, mem_SQ] at hz
    rw [mem_coe, mem_SQ']
    linear_combination hz - (z.2.1 ^ 2 + z.2.2 ^ 2) * hst
  · intro y _
    ext
    · simp only; linear_combination (-y.1) * hst
    · simp only; linear_combination (-y.2.1) * hst
    · rfl
  · intro z _
    ext
    · rfl
    · simp only; linear_combination (-z.2.1) * hst
    · simp only; linear_combination (-z.2.2) * hst

theorem card_SQ'_odd_prime_pow (hp2 : p ≠ 2) (f : ℕ) :
    (SQ' (ZMod (p ^ (f + 1)))).card = p ^ (2 * f) * (p ^ 2 + p) := by
  induction f with
  | zero =>
    rw [card_SQ'_congr (ZMod.ringEquivCongr (pow_one p)), card_SQ'_prime p hp2]; ring
  | succ k ih =>
    rw [show k + 1 + 1 = k + 2 by omega, card_SQ'_succ p k hp2, ih]; ring

/-- **A227867 at odd prime powers.** -/
theorem A_odd_prime_pow (hp2 : p ≠ 2) {k : ℕ} (hk : k ≠ 0) :
    A (p ^ k) = p ^ (2 * k - 1) * (p + 1) + 2 := by
  obtain ⟨f, rfl⟩ : ∃ f, k = f + 1 := ⟨k - 1, by omega⟩
  unfold A
  rw [cnt_eq p hp2 (by omega), card_SQ'_odd_prime_pow p hp2,
    show 2 * (f + 1) - 1 = 2 * f + 1 by omega]
  ring

end base

/-! ### Powers of 2 -/

theorem A_two : A 2 = 8 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide +kernel

set_option maxRecDepth 100000 in
theorem A_four : A 4 = 32 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide +kernel

/-- A sum of three squares is never `7 mod 8`. -/
theorem neg_sum_three_sq_ne_one : ∀ x y z : ZMod 8, -x ^ 2 - y ^ 2 - z ^ 2 ≠ 1 := by decide

section twopow

open A236554 (hh H mem_H sq_of_mem_H two_mul_of_mem_H card_H card_roots zmod2_cases
  isUnit_of_ρ_eq_one eq_two_mul_of_ρ_eq_zero two_mul_eq_zero_iff)

variable (m : ℕ)

/-- Modulo `2^(m+3)`, the solutions are exactly `{a | a^2 = 1} × {0, h}^3`. -/
theorem filter_P_eq_two :
    univ.filter (fun x : ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) × ZMod (2 ^ (m + 3)) ×
      ZMod (2 ^ (m + 3)) => P x) =
    univ.filter (fun a : ZMod (2 ^ (m + 3)) => a ^ 2 = 1) ×ˢ (H m ×ˢ H m ×ˢ H m) := by
  have hn : m + 3 ≠ 0 := by omega
  ext ⟨a, b, c, d⟩
  simp only [mem_filter, mem_univ, true_and, mem_product]
  constructor
  · rintro ⟨heq, hb, hc, hd⟩
    rcases zmod2_cases (A236554.ρ (m + 3) hn a) with ha | ha
    · -- `a` even: impossible
      exfalso
      by_cases h2a : 2 * a = 0
      · -- `a ∈ {0, h}`, so `-(b^2 + c^2 + d^2) = 1`, impossible modulo 8
        have haH : a ∈ H m := (mem_H m a).2 ((two_mul_eq_zero_iff hn a).1 h2a)
        have hQ : -b ^ 2 - c ^ 2 - d ^ 2 = 1 := by linear_combination heq - sq_of_mem_H m haH
        let φ8 : ZMod (2 ^ (m + 3)) →+* ZMod 8 := ZMod.castHom (Dvd.intro (2 ^ m) (by ring)) _
        have := congrArg φ8 hQ
        simp only [map_sub, map_neg, map_pow, map_one] at this
        exact neg_sum_three_sq_ne_one _ _ _ this
      · -- then `b, c, d` are even too, and the equation fails modulo 2
        have even_of : ∀ z, 2 * a * z = 0 → A236554.ρ (m + 3) hn z = 0 := by
          intro z hz
          rcases zmod2_cases (A236554.ρ (m + 3) hn z) with h | h
          · exact h
          · exact absurd ((isUnit_of_ρ_eq_one hn h).mul_right_eq_zero.1
              (by linear_combination hz)) h2a
        obtain ⟨b', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hb)
        obtain ⟨c', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hc)
        obtain ⟨d', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn (even_of _ hd)
        obtain ⟨a', rfl⟩ := eq_two_mul_of_ρ_eq_zero hn ha
        have h1 : (1 : ZMod (2 ^ (m + 3))) = 2 * (2 * (a' ^ 2 - b' ^ 2 - c' ^ 2 - d' ^ 2)) := by
          linear_combination -heq
        have := congrArg (A236554.ρ (m + 3) hn) h1
        rw [map_one, map_mul, map_ofNat, show (2 : ZMod 2) = 0 by decide, zero_mul] at this
        exact absurd this (by decide)
    · -- `a` odd, hence a unit, so `b, c, d ∈ {0, h}`
      have hu := isUnit_of_ρ_eq_one hn ha
      have hH : ∀ z, 2 * a * z = 0 → z ∈ H m := fun z hz =>
        (mem_H m z).2 ((two_mul_eq_zero_iff hn z).1
          (hu.mul_right_eq_zero.1 (by linear_combination hz)))
      have hbH := hH b hb
      have hcH := hH c hc
      have hdH := hH d hd
      exact ⟨by linear_combination heq + sq_of_mem_H m hbH + sq_of_mem_H m hcH
        + sq_of_mem_H m hdH, hbH, hcH, hdH⟩
  · rintro ⟨ha, hb, hc, hd⟩
    exact ⟨by linear_combination ha - sq_of_mem_H m hb - sq_of_mem_H m hc - sq_of_mem_H m hd,
      by linear_combination a * two_mul_of_mem_H m hb,
      by linear_combination a * two_mul_of_mem_H m hc,
      by linear_combination a * two_mul_of_mem_H m hd⟩

theorem A_two_pow_add_three : A (2 ^ (m + 3)) = 32 := by
  unfold A cnt
  rw [Nat.card_eq_fintype_card, Fintype.card_subtype, filter_P_eq_two, card_product,
    card_product, card_product, card_roots, card_H]

end twopow

theorem A_two_pow {k : ℕ} (hk : 2 ≤ k) : A (2 ^ k) = 32 := by
  rcases Nat.lt_or_ge k 3 with h | h
  · obtain rfl : k = 2 := by omega
    exact A_four
  · obtain ⟨m, rfl⟩ : ∃ m, k = m + 3 := ⟨k - 3, by omega⟩
    exact A_two_pow_add_three m

/-! ### The closed form -/

/-- The value of A227867 at the prime power `p^k` (`k ≥ 1`). -/
def localFactor (p k : ℕ) : ℕ :=
  if p = 2 then (if k = 1 then 8 else 32) else p ^ (2 * k - 1) * (p + 1) + 2

theorem A_prime_pow {p k : ℕ} (hp : p.Prime) (hk : k ≠ 0) : A (p ^ k) = localFactor p k := by
  unfold localFactor
  split_ifs with h2 h1
  · subst h2 h1; exact A_two
  · subst h2; exact A_two_pow (by omega)
  · have := Fact.mk hp; exact A_odd_prime_pow p h2 hk

/-- **A227867, closed form**: for `n ≠ 0`, `A n = ∏_{p^k ∥ n} localFactor p k`, where
`localFactor p k = p^(2k-1) (p+1) + 2` for odd `p`, `8` for `2`, and `32` for `2^k`, `k ≥ 2`. -/
theorem A_closed_form {n : ℕ} (hn : n ≠ 0) :
    A n = n.factorization.prod (fun p k => localFactor p k) := by
  rw [Nat.multiplicative_factorization A (fun x y h => A_mul h) A_one hn]
  apply Finsupp.prod_congr
  intro p hp
  rw [Nat.support_factorization] at hp
  exact A_prime_pow (Nat.prime_of_mem_primeFactors hp) (Finsupp.mem_support_iff.1
    (by rwa [Nat.support_factorization]))

end A227867
