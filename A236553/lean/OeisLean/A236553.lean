import OeisLean.A236554
import Mathlib.Data.ZMod.Units
import Mathlib.Data.Nat.Factorization.Induction

/-!
# A236553: involutions in the quaternion ring over `ZMod n` with `i^2 = j^2 = 1`

A236553(n) is the number of solutions `(a, b, c, d)` modulo `n` of
  `a^2 - b^2 + c^2 + d^2 = 1`, `2ab = 0`, `2ac = 0`, `2ad = 0`.

Main results:
* `A_mul`: `A` is multiplicative (Chinese remainder theorem);
* `A_odd_prime_pow`: `A (p^k) = p^(2k-1) (p+1) + 2` for odd primes `p` and `k ≥ 1`;
* `A_two_pow`: `A 2 = 8`, `A 4 = 64`, `A (2^k) = 2^(2k+2) + 32` for `k ≥ 3` (from A236554);
* `A_closed_form`: `A n = ∏_{p^k ∥ n} localFactor p k` for `n ≠ 0`.
-/

open Finset

namespace A236553

open A236554 (Q Q_map)

/-- The defining system, over any commutative ring. -/
abbrev P {R : Type*} [CommRing R] (x : R × R × R × R) : Prop :=
  x.1 ^ 2 - x.2.1 ^ 2 + x.2.2.1 ^ 2 + x.2.2.2 ^ 2 = 1 ∧
    2 * x.1 * x.2.1 = 0 ∧ 2 * x.1 * x.2.2.1 = 0 ∧ 2 * x.1 * x.2.2.2 = 0

/-- Number of solutions of the system over a ring `R`. -/
noncomputable def cnt (R : Type*) [CommRing R] : ℕ := Nat.card {x : R × R × R × R // P x}

/-- **A236553(n)**. -/
noncomputable def A (n : ℕ) : ℕ := cnt (ZMod n)

section general

variable {R S : Type*} [CommRing R] [CommRing S]

def map4 (f : R →+* S) (x : R × R × R × R) : S × S × S × S :=
  (f x.1, f x.2.1, f x.2.2.1, f x.2.2.2)

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

/-- `A` is multiplicative. -/
theorem A_mul {m n : ℕ} (h : m.Coprime n) : A (m * n) = A m * A n := by
  unfold A; rw [cnt_congr (ZMod.chineseRemainder h), cnt_prod]

theorem A_one : A 1 = 1 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide

/-! ### Odd prime powers: the case split -/

/-- Solutions of `Q = 1` in a finite ring. -/
def SQ (R : Type*) [CommRing R] [Fintype R] [DecidableEq R] : Finset (R × R × R) :=
  univ.filter (fun y => Q y = 1)

theorem mem_SQ {R : Type*} [CommRing R] [Fintype R] [DecidableEq R] (y : R × R × R) :
    y ∈ SQ R ↔ y.2.1 ^ 2 + y.2.2 ^ 2 - y.1 ^ 2 = 1 := by
  unfold SQ; rw [mem_filter]; simp [Q]

section oddprime

variable (p : ℕ) [hp : Fact p.Prime]

/-- Reduction `ZMod (p^e) → ZMod p`. -/
def ρ (e : ℕ) (he : e ≠ 0) : ZMod (p ^ e) →+* ZMod p := ZMod.castHom (dvd_pow_self p he) _

theorem ρ_natCast {e : ℕ} (he : e ≠ 0) (k : ℕ) : ρ p e he (k : ZMod (p ^ e)) = (k : ZMod p) :=
  map_natCast _ k

theorem two_ne_zero_p (hp2 : p ≠ 2) : (2 : ZMod p) ≠ 0 := by
  intro h
  have h' : ((2 : ℕ) : ZMod p) = 0 := by exact_mod_cast h
  rw [ZMod.natCast_eq_zero_iff] at h'
  exact hp2 ((Nat.prime_dvd_prime_iff_eq hp.out Nat.prime_two).1 h')

theorem isUnit_iff_ρ {e : ℕ} (he : e ≠ 0) (z : ZMod (p ^ e)) : IsUnit z ↔ ρ p e he z ≠ 0 := by
  constructor
  · intro hu h0
    have := hu.map (ρ p e he)
    rw [h0] at this
    exact not_isUnit_zero this
  · intro h
    rw [← ZMod.natCast_zmod_val z] at h ⊢
    rw [ρ_natCast, Ne, ZMod.natCast_eq_zero_iff] at h
    rw [ZMod.isUnit_iff_coprime]
    exact Nat.Coprime.pow_right e ((Nat.Prime.coprime_iff_not_dvd hp.out).2 h).symm

theorem two_isUnit (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) : IsUnit (2 : ZMod (p ^ e)) :=
  (isUnit_iff_ρ p he 2).2 (by rw [map_ofNat]; exact two_ne_zero_p p hp2)

theorem sq_eq_one_iff_p (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) (a : ZMod (p ^ e)) :
    a ^ 2 = 1 ↔ a = 1 ∨ a = -1 := by
  constructor
  · intro h
    have hprod : (a - 1) * (a + 1) = 0 := by linear_combination h
    by_cases hu : ρ p e he (a - 1) = 0
    · have hu' : IsUnit (a + 1) := (isUnit_iff_ρ p he _).2 (by
        rw [show a + 1 = (a - 1) + 2 by ring, map_add, hu, map_ofNat, zero_add]
        exact two_ne_zero_p p hp2)
      left; linear_combination hu'.mul_left_eq_zero.1 hprod
    · right; linear_combination ((isUnit_iff_ρ p he _).2 hu).mul_right_eq_zero.1 hprod
  · rintro (rfl | rfl) <;> ring

theorem one_ne_neg_one_p (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) : (1 : ZMod (p ^ e)) ≠ -1 := by
  intro h
  have h2 : (2 : ZMod (p ^ e)) = 0 := by linear_combination h
  have := congrArg (ρ p e he) h2
  rw [map_ofNat, map_zero] at this
  exact two_ne_zero_p p hp2 this

theorem card_roots_p (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) :
    (univ.filter (fun a : ZMod (p ^ e) => a ^ 2 = 1)).card = 2 := by
  have : univ.filter (fun a : ZMod (p ^ e) => a ^ 2 = 1) = {1, -1} := by
    ext a; simp [sq_eq_one_iff_p p hp2 he]
  rw [this, card_pair (one_ne_neg_one_p p hp2 he)]

/-- For odd `p`, a solution has either `a = ±1, b = c = d = 0`, or `a = 0` and `Q (b,c,d) = 1`. -/
theorem filter_P_eq (hp2 : p ≠ 2) {e : ℕ} (he : e ≠ 0) :
    univ.filter (fun x : ZMod (p ^ e) × ZMod (p ^ e) × ZMod (p ^ e) × ZMod (p ^ e) => P x) =
      (univ.filter (fun a : ZMod (p ^ e) => a ^ 2 = 1)) ×ˢ {(0, 0, 0)} ∪
        {0} ×ˢ SQ (ZMod (p ^ e)) := by
  ext ⟨a, b, c, d⟩
  simp only [mem_filter, mem_univ, true_and, mem_union, mem_product, mem_singleton, mem_SQ,
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
        rw [map_add, map_add, map_sub, map_pow, map_pow, map_pow, map_pow, ha, hnu b hb,
          hnu c hc, hnu d hd, map_one] at this
        simp at this
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
    cnt (ZMod (p ^ e)) = 2 + (SQ (ZMod (p ^ e))).card := by
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

/-! ### Odd prime powers: Hensel lifting, counted -/

section hensel

variable (p : ℕ) [hp : Fact p.Prime]

omit hp in
theorem p_pow_self (e : ℕ) : (p : ZMod (p ^ e)) ^ e = 0 := by
  have := ZMod.natCast_self (p ^ e)
  push_cast at this
  exact this

theorem eq_p_mul_of_ρ_eq_zero {e : ℕ} (he : e ≠ 0) {z : ZMod (p ^ e)} (h : ρ p e he z = 0) :
    ∃ w, z = p * w := by
  have hz := ZMod.natCast_zmod_val z
  rw [← hz, ρ_natCast, ZMod.natCast_eq_zero_iff] at h
  obtain ⟨q, hq⟩ := h
  exact ⟨q, by rw [← hz, hq]; push_cast; ring⟩

variable (f : ℕ)

/-- Reduction `ZMod (p^(f+2)) → ZMod (p^(f+1))`. -/
def πp : ZMod (p ^ (f + 2)) →+* ZMod (p ^ (f + 1)) := ZMod.castHom (pow_dvd_pow p (by omega)) _

def π3p (y : ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2))) :
    ZMod (p ^ (f + 1)) × ZMod (p ^ (f + 1)) × ZMod (p ^ (f + 1)) :=
  (πp p f y.1, πp p f y.2.1, πp p f y.2.2)

/-- Canonical lift `ZMod (p^(f+1)) → ZMod (p^(f+2))`. -/
def σp (x : ZMod (p ^ (f + 1))) : ZMod (p ^ (f + 2)) := (x.val : ℕ)

/-- Canonical lift `ZMod p → ZMod (p^(f+2))`. -/
def liftp (t : ZMod p) : ZMod (p ^ (f + 2)) := (t.val : ℕ)

abbrev ρ1 := ρ p (f + 1) (by omega)
abbrev ρ2 := ρ p (f + 2) (by omega)

theorem πp_σp (x : ZMod (p ^ (f + 1))) : πp p f (σp p f x) = x := by
  unfold πp σp; rw [map_natCast, ZMod.natCast_zmod_val]

theorem ρ2_σp (x : ZMod (p ^ (f + 1))) : ρ2 p f (σp p f x) = ρ1 p f x := by
  conv_rhs => rw [← ZMod.natCast_zmod_val x]
  unfold σp; rw [ρ_natCast, ρ_natCast]

theorem ρ2_liftp (t : ZMod p) : ρ2 p f (liftp p f t) = t := by
  unfold liftp; rw [ρ_natCast, ZMod.natCast_zmod_val]

omit hp in
theorem πp_P1 : πp p f ((p : ZMod (p ^ (f + 2))) ^ (f + 1)) = 0 := by
  unfold πp; rw [map_pow, map_natCast, p_pow_self]

theorem eq_of_πp_eq_zero {z : ZMod (p ^ (f + 2))} (h : πp p f z = 0) :
    ∃ j, z = (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j := by
  have hz := ZMod.natCast_zmod_val z
  rw [← hz] at h
  unfold πp at h
  rw [map_natCast, ZMod.natCast_eq_zero_iff] at h
  obtain ⟨q, hq⟩ := h
  exact ⟨q, by rw [← hz, hq]; push_cast; ring⟩

/-- `p^(f+1) z = 0` modulo `p^(f+2)` iff `z ≡ 0 (mod p)`. -/
theorem P1_mul_eq_zero_iff (z : ZMod (p ^ (f + 2))) :
    (p : ZMod (p ^ (f + 2))) ^ (f + 1) * z = 0 ↔ ρ2 p f z = 0 := by
  constructor
  · intro h
    have hz := ZMod.natCast_zmod_val z
    rw [← hz] at h
    have h' : ((p ^ (f + 1) * z.val : ℕ) : ZMod (p ^ (f + 2))) = 0 := by push_cast; exact h
    rw [ZMod.natCast_eq_zero_iff] at h'
    have key : p ^ (f + 1) * p ∣ p ^ (f + 1) * z.val := by rwa [← pow_succ]
    have hdvd := Nat.dvd_of_mul_dvd_mul_left (by have := hp.out.pos; positivity) key
    rw [← hz, ρ_natCast, ZMod.natCast_eq_zero_iff]
    exact hdvd
  · intro h
    obtain ⟨w, rfl⟩ := eq_p_mul_of_ρ_eq_zero p (by omega) h
    rw [← mul_assoc, ← pow_succ, p_pow_self, zero_mul]

omit hp in
/-- The quadratic form on a coset of the kernel is affine. -/
theorem Q_add (u j : ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2))) :
    Q (u.1 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.1,
      u.2.1 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.2.1,
      u.2.2 + (p : ZMod (p ^ (f + 2))) ^ (f + 1) * j.2.2) =
    Q u + (p : ZMod (p ^ (f + 2))) ^ (f + 1) *
      (2 * (u.2.1 * j.2.1 + u.2.2 * j.2.2 - u.1 * j.1)) := by
  simp only [Q]
  linear_combination (j.2.1 ^ 2 + j.2.2 ^ 2 - j.1 ^ 2) * (p : ZMod (p ^ (f + 2))) ^ f *
    p_pow_self p (f + 2)

/-- Every solution modulo `p^(f+1)` lifts to exactly `p^2` solutions modulo `p^(f+2)`. -/
theorem card_fiber_SQ (hp2 : p ≠ 2) {x} (hx : x ∈ SQ (ZMod (p ^ (f + 1)))) :
    ((SQ (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card = p ^ 2 := by
  set P1 : ZMod (p ^ (f + 2)) := (p : ZMod (p ^ (f + 2))) ^ (f + 1) with hP1
  set u : ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) :=
    (σp p f x.1, σp p f x.2.1, σp p f x.2.2) with hu
  have hx' : Q x = 1 := by rw [mem_SQ] at hx; simpa [Q] using hx
  have hQu : πp p f (Q u - 1) = 0 := by
    rw [map_sub, map_one, ← Q_map]
    simp only [hu, πp_σp]
    rw [hx', sub_self]
  obtain ⟨e0, he0⟩ := eq_of_πp_eq_zero p f hQu
  set xb := ρ1 p f x.1
  set xc := ρ1 p f x.2.1
  set xd := ρ1 p f x.2.2
  have hQbar : xc ^ 2 + xd ^ 2 - xb ^ 2 = 1 := by
    have := congrArg (ρ1 p f) hx'
    simpa [Q, map_sub, map_add, map_pow] using this
  let ℓ : (ZMod p × ZMod p × ZMod p) →+ ZMod p :=
    AddMonoidHom.mk' (fun t => 2 * (xc * t.2.1 + xd * t.2.2 - xb * t.1))
      (by intro a b; simp only [Prod.fst_add, Prod.snd_add]; ring)
  have hℓ : ∀ t, ℓ t = 2 * (xc * t.2.1 + xd * t.2.2 - xb * t.1) := fun t => rfl
  set β := -ρ2 p f e0
  let τ : ZMod p × ZMod p × ZMod p →
      ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) × ZMod (p ^ (f + 2)) :=
    fun t => (u.1 + P1 * liftp p f t.1, u.2.1 + P1 * liftp p f t.2.1,
      u.2.2 + P1 * liftp p f t.2.2)
  have key : ∀ t, Q (τ t) = 1 ↔ ℓ t = β := by
    intro t
    have hexp := Q_add p f u (liftp p f t.1, liftp p f t.2.1, liftp p f t.2.2)
    simp only at hexp
    have h1 : Q (τ t) = 1 ↔ P1 * (e0 + 2 * (u.2.1 * liftp p f t.2.1 +
        u.2.2 * liftp p f t.2.2 - u.1 * liftp p f t.1)) = 0 := by
      change Q (u.1 + P1 * liftp p f t.1, u.2.1 + P1 * liftp p f t.2.1,
        u.2.2 + P1 * liftp p f t.2.2) = 1 ↔ _
      rw [hexp]
      constructor
      · intro h; linear_combination h - he0
      · intro h; linear_combination h + he0
    rw [h1, P1_mul_eq_zero_iff, hℓ]
    simp only [map_add, map_mul, map_sub, map_ofNat, hu, ρ2_σp, ρ2_liftp]
    constructor
    · intro h; linear_combination h
    · intro h; linear_combination h
  -- the fibre is in bijection with the solutions of the affine equation `ℓ t = β`
  have hbij : (univ.filter (fun t => ℓ t = β)).card =
      ((SQ (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card := by
    apply card_bij (fun t _ => τ t)
    · intro t ht
      simp only [mem_filter, mem_univ, true_and] at ht
      rw [mem_filter, mem_SQ]
      refine ⟨?_, ?_⟩
      · have := (key t).2 ht; simpa [Q] using this
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
      rw [mem_filter, mem_SQ] at hy
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
      simpa [Q] using hyQ
  rw [← hbij]
  -- `ℓ` is surjective, so every fibre has `p^3 / p = p^2` elements
  have h2 := two_ne_zero_p p hp2
  have hsurj : ∀ r, r ∈ Set.range ℓ := by
    intro r
    by_cases hc : xc = 0
    · by_cases hd : xd = 0
      · have hb : xb ≠ 0 := by
          intro hb; rw [hb, hc, hd] at hQbar; simp at hQbar
        refine ⟨(-r / (2 * xb), 0, 0), ?_⟩
        rw [hℓ]; field_simp; ring
      · refine ⟨(0, 0, r / (2 * xd)), ?_⟩
        rw [hℓ]; field_simp; ring
    · refine ⟨(0, r / (2 * xc), 0), ?_⟩
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

theorem card_SQ_succ (hp2 : p ≠ 2) :
    (SQ (ZMod (p ^ (f + 2)))).card = p ^ 2 * (SQ (ZMod (p ^ (f + 1)))).card := by
  rw [card_eq_sum_card_fiberwise (f := π3p p f) (t := univ) (fun _ _ => mem_univ _)]
  have hterm : ∀ x, ((SQ (ZMod (p ^ (f + 2)))).filter (fun y => π3p p f y = x)).card =
      if x ∈ SQ (ZMod (p ^ (f + 1))) then p ^ 2 else 0 := by
    intro x
    split_ifs with hx
    · exact card_fiber_SQ p f hp2 hx
    · rw [card_eq_zero, filter_eq_empty_iff]
      intro y hy hyx
      apply hx
      rw [mem_SQ] at hy ⊢
      have := congrArg (πp p f) hy
      rw [← hyx]
      simpa [π3p, map_sub, map_add, map_pow] using this
  simp only [hterm]
  rw [sum_ite_mem, univ_inter, sum_const, smul_eq_mul, mul_comm]

end hensel

/-! ### Odd primes: the base count modulo `p` -/

section base

variable (p : ℕ) [hp : Fact p.Prime]

theorem card_SQ_congr {R S : Type*} [CommRing R] [Fintype R] [DecidableEq R] [CommRing S]
    [Fintype S] [DecidableEq S] (e : R ≃+* S) : (SQ R).card = (SQ S).card := by
  apply card_nbij' (fun y => (e y.1, e y.2.1, e y.2.2)) (fun z => (e.symm z.1, e.symm z.2.1,
    e.symm z.2.2))
  · intro y hy
    rw [mem_coe, mem_SQ] at hy ⊢
    have := congrArg e hy
    simpa [map_sub, map_add, map_pow] using this
  · intro z hz
    rw [mem_coe, mem_SQ] at hz ⊢
    have := congrArg e.symm hz
    simpa [map_sub, map_add, map_pow] using this
  · intro y _; simp
  · intro z _; simp

/-- Number of solutions of `u w = r` in `ZMod p`. -/
def N (r : ZMod p) : ℕ := (univ.filter (fun z : ZMod p × ZMod p => z.1 * z.2 = r)).card

theorem N_of_ne_zero {r : ZMod p} (hr : r ≠ 0) : N p r = p - 1 := by
  unfold N
  have hcard : (univ.filter (fun u : ZMod p => u ≠ 0)).card = p - 1 := by
    rw [filter_ne', card_erase_of_mem (mem_univ _), card_univ, ZMod.card]
  rw [← hcard]
  apply card_nbij' (fun z => z.1) (fun u => (u, r / u))
  · intro z hz
    rw [mem_coe, mem_filter] at hz ⊢
    refine ⟨mem_univ _, ?_⟩
    intro h0; simp only at h0; rw [h0, zero_mul] at hz; exact hr hz.2.symm
  · intro u hu
    rw [mem_coe, mem_filter] at hu ⊢
    exact ⟨mem_univ _, mul_div_cancel₀ r hu.2⟩
  · intro z hz
    rw [mem_coe, mem_filter] at hz
    have h1 : z.1 ≠ 0 := by intro h0; rw [h0, zero_mul] at hz; exact hr hz.2.symm
    ext
    · rfl
    · change r / z.1 = z.2
      rw [div_eq_iff h1]; linear_combination -hz.2
  · intro u _; rfl

theorem N_zero : N p 0 = p - 1 + p := by
  have hsum := card_eq_sum_card_fiberwise (f := fun z : ZMod p × ZMod p => z.1 * z.2)
    (s := univ) (t := univ) (fun _ _ => mem_univ _)
  rw [← add_sum_erase _ _ (mem_univ (0 : ZMod p))] at hsum
  have hrest : ∑ r ∈ univ.erase (0 : ZMod p),
      (univ.filter (fun z : ZMod p × ZMod p => z.1 * z.2 = r)).card = (p - 1) * (p - 1) := by
    have hc : ∀ r ∈ univ.erase (0 : ZMod p),
        (univ.filter (fun z : ZMod p × ZMod p => z.1 * z.2 = r)).card = p - 1 :=
      fun r hr => N_of_ne_zero p (ne_of_mem_erase hr)
    rw [sum_congr rfl hc, sum_const,
      card_erase_of_mem (mem_univ _), card_univ, ZMod.card, smul_eq_mul]
  rw [hrest, card_univ, Fintype.card_prod, ZMod.card] at hsum
  have h1 : 1 ≤ p := hp.out.one_lt.le
  unfold N
  zify [h1] at hsum ⊢
  linear_combination -hsum

theorem card_sq_eq_one_prime (hp2 : p ≠ 2) :
    (univ.filter (fun c : ZMod p => c ^ 2 = 1)).card = 2 := by
  have hne : (1 : ZMod p) ≠ -1 := by
    intro h
    exact two_ne_zero_p p hp2 (by linear_combination h)
  have : univ.filter (fun c : ZMod p => c ^ 2 = 1) = {1, -1} := by
    ext c; simp [pow_two, mul_self_eq_one_iff]
  rw [this, card_pair hne]

theorem card_SQ_prime (hp2 : p ≠ 2) : (SQ (ZMod p)).card = p ^ 2 + p := by
  have h2 := two_ne_zero_p p hp2
  set i2 : ZMod p := (2 : ZMod p)⁻¹
  have hi2 : 2 * i2 = 1 := mul_inv_cancel₀ h2
  -- change variables `u = d - b`, `w = d + b`
  have step1 : (SQ (ZMod p)).card =
      (univ.filter (fun z : ZMod p × ZMod p × ZMod p => z.1 ^ 2 + z.2.1 * z.2.2 = 1)).card := by
    apply card_nbij' (fun y => (y.2.1, y.2.2 - y.1, y.2.2 + y.1))
      (fun z => ((z.2.2 - z.2.1) * i2, z.1, (z.2.1 + z.2.2) * i2))
    · intro y hy
      rw [mem_coe, mem_SQ] at hy
      rw [mem_coe, mem_filter]
      exact ⟨mem_univ _, by linear_combination hy⟩
    · intro z hz
      rw [mem_coe, mem_filter] at hz
      rw [mem_coe, mem_SQ]
      linear_combination hz.2 + z.2.1 * z.2.2 * (2 * i2 + 1) * hi2
    · intro y _
      ext
      · simp only; linear_combination y.1 * hi2
      · rfl
      · simp only; linear_combination y.2.2 * hi2
    · intro z _
      ext
      · rfl
      · simp only; linear_combination z.2.1 * hi2
      · simp only; linear_combination z.2.2 * hi2
  -- sum over `c` of the number of `(u, w)` with `u w = 1 - c^2`
  have step2 : (univ.filter (fun z : ZMod p × ZMod p × ZMod p => z.1 ^ 2 + z.2.1 * z.2.2 = 1)).card
      = ∑ c : ZMod p, N p (1 - c ^ 2) := by
    rw [card_filter, Fintype.sum_prod_type]
    apply sum_congr rfl
    intro c _
    unfold N
    rw [card_filter]
    apply sum_congr rfl
    intro z _
    simp only [eq_sub_iff_add_eq']
  have hN : ∀ c : ZMod p, N p (1 - c ^ 2) = (p - 1) + if c ^ 2 = 1 then p else 0 := by
    intro c
    split_ifs with hc
    · rw [hc, sub_self, N_zero]
    · rw [N_of_ne_zero p (by intro h; apply hc; linear_combination -h), add_zero]
  rw [step1, step2, sum_congr rfl (fun c _ => hN c), sum_add_distrib, ← sum_filter, sum_const,
    sum_const, card_sq_eq_one_prime p hp2, card_univ, ZMod.card, smul_eq_mul, smul_eq_mul]
  have h1 : 1 ≤ p := hp.out.one_lt.le
  zify [h1]
  ring

theorem card_SQ_odd_prime_pow (hp2 : p ≠ 2) (f : ℕ) :
    (SQ (ZMod (p ^ (f + 1)))).card = p ^ (2 * f) * (p ^ 2 + p) := by
  induction f with
  | zero =>
    rw [card_SQ_congr (ZMod.ringEquivCongr (pow_one p)), card_SQ_prime p hp2]; ring
  | succ k ih =>
    rw [show k + 1 + 1 = k + 2 by omega, card_SQ_succ p k hp2, ih]; ring

/-- **A236553 at odd prime powers.** -/
theorem A_odd_prime_pow (hp2 : p ≠ 2) {k : ℕ} (hk : k ≠ 0) :
    A (p ^ k) = p ^ (2 * k - 1) * (p + 1) + 2 := by
  obtain ⟨f, rfl⟩ : ∃ f, k = f + 1 := ⟨k - 1, by omega⟩
  unfold A
  rw [cnt_eq p hp2 (by omega), card_SQ_odd_prime_pow p hp2,
    show 2 * (f + 1) - 1 = 2 * f + 1 by omega]
  ring

end base

/-! ### Powers of 2, from A236554 -/

theorem cnt_eq_T (n : ℕ) : cnt (ZMod (2 ^ n)) = (A236554.T n).card := by
  unfold cnt A236554.T
  rw [Nat.card_eq_fintype_card, Fintype.card_subtype]

theorem A_two : A 2 = 8 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide +kernel

set_option maxRecDepth 100000 in
theorem A_four : A 4 = 64 := by
  unfold A cnt; rw [Nat.card_eq_fintype_card]; decide +kernel

theorem A_two_pow {k : ℕ} (hk : 3 ≤ k) : A (2 ^ k) = 2 ^ (2 * k + 2) + 32 := by
  unfold A; rw [cnt_eq_T, A236554.card_A236554 k hk]

/-! ### The closed form -/

/-- The value of A236553 at the prime power `p^k` (`k ≥ 1`). -/
def localFactor (p k : ℕ) : ℕ :=
  if p = 2 then (if k = 1 then 8 else if k = 2 then 64 else 2 ^ (2 * k + 2) + 32)
  else p ^ (2 * k - 1) * (p + 1) + 2

theorem A_prime_pow {p k : ℕ} (hp : p.Prime) (hk : k ≠ 0) : A (p ^ k) = localFactor p k := by
  unfold localFactor
  split_ifs with h2 h1 h22
  · subst h2 h1; exact A_two
  · subst h2 h22; exact A_four
  · subst h2; exact A_two_pow (by omega)
  · have := Fact.mk hp; exact A_odd_prime_pow p h2 hk

/-- **A236553, closed form**: for `n ≠ 0`, `A n = ∏_{p^k ∥ n} localFactor p k`, where
`localFactor p k = p^(2k-1) (p+1) + 2` for odd `p`, and `8, 64, 2^(2k+2) + 32` for
`2, 4, 2^k (k ≥ 3)`. -/
theorem A_closed_form {n : ℕ} (hn : n ≠ 0) :
    A n = n.factorization.prod (fun p k => localFactor p k) := by
  rw [Nat.multiplicative_factorization A (fun x y h => A_mul h) A_one hn]
  apply Finsupp.prod_congr
  intro p hp
  rw [Nat.support_factorization] at hp
  exact A_prime_pow (Nat.prime_of_mem_primeFactors hp) (Finsupp.mem_support_iff.1
    (by rwa [Nat.support_factorization]))

end A236553
