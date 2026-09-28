import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Fintype.Pi
import Mathlib.Data.Fintype.Sigma
import Mathlib.Tactic

/-!
# A006455: the split identity for naturally labeled posets

A *naturally labeled poset* on `Fin N` is a strict partial order contained in the usual order
(`x R y → x < y`). A006455(N) counts them. We store the strict order as a Boolean relation.

Cut `Fin (n + k)` into the bottom part `Fin n` (via `Fin.castAdd`) and the top part `Fin k`
(via `Fin.natAdd`). A natural order `r` on `Fin (n + k)` restricts to natural orders `P` on the
bottom and `Q` on the top. It has no relation from top to bottom, and it records, for each top
element `y`, the set `f y` of bottom elements below it. Each `f y` is a down-set of `P`, and `f`
is order-preserving from `Q` into down-sets. Conversely, any such triple glues back uniquely.
This is the case of natural labelings of the "generalized vertical sums" of a Campo and Erné
(arXiv:1802.01419, Theorem 2.1).

Main results:
* `equiv`: natural orders on `Fin (n + k)` ≃ triples `(P, Q, f)` as above;
* `a_add`: A006455(n + k) = Σ_P Σ_Q #{order-preserving f from Q to the down-sets of P}.

This is the identity behind the program that computed a(13)–a(16). Standard axioms only.
-/

namespace A006455

open Finset

/-- `r` is a natural order on `Fin N`: contained in `<`, and transitive. -/
def IsNat {N : ℕ} (r : Fin N → Fin N → Bool) : Prop :=
  (∀ a b, r a b = true → a < b) ∧ ∀ a b c, r a b = true → r b c = true → r a c = true

instance {N : ℕ} (r : Fin N → Fin N → Bool) : Decidable (IsNat r) := by
  unfold IsNat; infer_instance

/-- A006455(N): the number of naturally labeled posets on `Fin N`. -/
def a (N : ℕ) : ℕ := Fintype.card {r : Fin N → Fin N → Bool // IsNat r}

/-- `D` is a down-set of the order `P`. -/
def IsDown {n : ℕ} (P : Fin n → Fin n → Bool) (D : Fin n → Bool) : Prop :=
  ∀ x x', P x x' = true → D x' = true → D x = true

/-- Gluing data for bottom order `P` and top order `Q`: every `f y` is a down-set of `P`, and
`f` is order-preserving from `Q` (pointwise inclusion). -/
def Glue {n k : ℕ} (P : Fin n → Fin n → Bool) (Q : Fin k → Fin k → Bool)
    (f : Fin k → Fin n → Bool) : Prop :=
  (∀ y, IsDown P (f y)) ∧ ∀ y y', Q y y' = true → ∀ x, f y x = true → f y' x = true

instance {n k : ℕ} (P : Fin n → Fin n → Bool) (Q : Fin k → Fin k → Bool)
    (f : Fin k → Fin n → Bool) : Decidable (Glue P Q f) := by
  unfold Glue IsDown; infer_instance

/-- A triple (bottom order, top order, gluing map). -/
abbrev Triple (n k : ℕ) :=
  (Fin n → Fin n → Bool) × (Fin k → Fin k → Bool) × (Fin k → Fin n → Bool)

/-- The triples that come from natural orders. -/
def Good {n k : ℕ} (t : Triple n k) : Prop := IsNat t.1 ∧ IsNat t.2.1 ∧ Glue t.1 t.2.1 t.2.2

/-- Restrict a relation on `Fin (n + k)` to the bottom, the top, and bottom-below-top. -/
def split {n k : ℕ} (r : Fin (n + k) → Fin (n + k) → Bool) : Triple n k :=
  (fun i i' => r (Fin.castAdd k i) (Fin.castAdd k i'),
   fun j j' => r (Fin.natAdd n j) (Fin.natAdd n j'),
   fun j i => r (Fin.castAdd k i) (Fin.natAdd n j))

/-- Glue a triple into a relation on `Fin (n + k)` (nothing from top to bottom). -/
def glue {n k : ℕ} (t : Triple n k) : Fin (n + k) → Fin (n + k) → Bool := fun u v =>
  Fin.addCases (motive := fun _ => Bool)
    (fun i => Fin.addCases (motive := fun _ => Bool) (fun i' => t.1 i i') (fun j => t.2.2 j i) v)
    (fun j => Fin.addCases (motive := fun _ => Bool) (fun _ => false) (fun j' => t.2.1 j j') v) u

section glue_simp
variable {n k : ℕ} (t : Triple n k) (i i' : Fin n) (j j' : Fin k)

@[simp] lemma glue_cc : glue t (Fin.castAdd k i) (Fin.castAdd k i') = t.1 i i' := by simp [glue]
@[simp] lemma glue_cn : glue t (Fin.castAdd k i) (Fin.natAdd n j) = t.2.2 j i := by simp [glue]
@[simp] lemma glue_nc : glue t (Fin.natAdd n j) (Fin.castAdd k i) = false := by simp [glue]
@[simp] lemma glue_nn : glue t (Fin.natAdd n j) (Fin.natAdd n j') = t.2.1 j j' := by simp [glue]

end glue_simp

section lt_lemmas
variable {n k : ℕ}

lemma castAdd_lt_natAdd (i : Fin n) (j : Fin k) : Fin.castAdd k i < Fin.natAdd n j := by
  rw [Fin.lt_def]; simp; omega

lemma castAdd_lt_iff (i i' : Fin n) : Fin.castAdd k i < Fin.castAdd k i' ↔ i < i' := by
  rw [Fin.lt_def, Fin.lt_def, Fin.val_castAdd, Fin.val_castAdd]

lemma natAdd_lt_iff (j j' : Fin k) : Fin.natAdd n j < Fin.natAdd n j' ↔ j < j' := by
  rw [Fin.lt_def, Fin.lt_def, Fin.val_natAdd, Fin.val_natAdd]; omega

lemma not_natAdd_lt_castAdd (i : Fin n) (j : Fin k) : ¬ Fin.natAdd n j < Fin.castAdd k i :=
  not_lt.2 (castAdd_lt_natAdd i j).le

end lt_lemmas

variable {n k : ℕ}

theorem split_good {r : Fin (n + k) → Fin (n + k) → Bool} (hr : IsNat r) : Good (split r) := by
  obtain ⟨hlt, htr⟩ := hr
  refine ⟨⟨fun i i' h => (castAdd_lt_iff i i').1 (hlt _ _ h), fun _ _ _ h1 h2 => htr _ _ _ h1 h2⟩,
    ⟨fun j j' h => (natAdd_lt_iff j j').1 (hlt _ _ h), fun _ _ _ h1 h2 => htr _ _ _ h1 h2⟩,
    fun _ _ _ h1 h2 => htr _ _ _ h1 h2, fun _ _ h x hx => htr _ _ _ hx h⟩

theorem glue_nat {t : Triple n k} (ht : Good t) : IsNat (glue t) := by
  obtain ⟨⟨hPlt, hPtr⟩, ⟨hQlt, hQtr⟩, hdown, hmono⟩ := ht
  constructor
  · intro u v h
    induction u using Fin.addCases with
    | left i =>
      induction v using Fin.addCases with
      | left i' => exact (castAdd_lt_iff i i').2 (hPlt _ _ (by simpa using h))
      | right j => exact castAdd_lt_natAdd i j
    | right j =>
      induction v using Fin.addCases with
      | left i' => simp at h
      | right j' => exact (natAdd_lt_iff j j').2 (hQlt _ _ (by simpa using h))
  · intro u v w h1 h2
    induction u using Fin.addCases with
    | left i =>
      induction v using Fin.addCases with
      | left i' =>
        induction w using Fin.addCases with
        | left i'' => simp only [glue_cc] at *; exact hPtr _ _ _ h1 h2
        | right j => simp only [glue_cc, glue_cn] at *; exact hdown j i i' h1 h2
      | right j =>
        induction w using Fin.addCases with
        | left i'' => simp at h2
        | right j' => simp only [glue_cn, glue_nn] at *; exact hmono j j' h2 i h1
    | right j =>
      induction v using Fin.addCases with
      | left i' => simp at h1
      | right j' =>
        induction w using Fin.addCases with
        | left i'' => simp at h2
        | right j'' => simp only [glue_nn] at *; exact hQtr _ _ _ h1 h2

theorem glue_split {r : Fin (n + k) → Fin (n + k) → Bool} (hr : IsNat r) :
    glue (split r) = r := by
  funext u v
  induction u using Fin.addCases with
  | left i =>
    induction v using Fin.addCases with
    | left i' => simp [split]
    | right j => simp [split]
  | right j =>
    induction v using Fin.addCases with
    | left i' =>
      rw [glue_nc]
      by_contra h
      exact not_natAdd_lt_castAdd i' j (hr.1 _ _ (by simpa using Ne.symm h))
    | right j' => simp [split]

theorem split_glue (t : Triple n k) : split (glue t) = t := by
  obtain ⟨P, Q, f⟩ := t
  simp [split]

/-- **The split bijection.** Natural orders on `Fin (n + k)` correspond to triples
(natural order on the bottom, natural order on the top, order-preserving map from the top into
the down-sets of the bottom). -/
def equiv (n k : ℕ) :
    {r : Fin (n + k) → Fin (n + k) → Bool // IsNat r} ≃ {t : Triple n k // Good t} where
  toFun r := ⟨split r.1, split_good r.2⟩
  invFun t := ⟨glue t.1, glue_nat t.2⟩
  left_inv r := Subtype.ext (glue_split r.2)
  right_inv t := Subtype.ext (split_glue t.1)

/-- Regroup good triples as a dependent sum over the bottom and top orders. -/
def goodEquivSigma (n k : ℕ) :
    {t : Triple n k // Good t} ≃
      Σ P : {P : Fin n → Fin n → Bool // IsNat P}, Σ Q : {Q : Fin k → Fin k → Bool // IsNat Q},
        {f : Fin k → Fin n → Bool // Glue P.1 Q.1 f} where
  toFun t := ⟨⟨t.1.1, t.2.1⟩, ⟨t.1.2.1, t.2.2.1⟩, ⟨t.1.2.2, t.2.2.2⟩⟩
  invFun s := ⟨(s.1.1, s.2.1.1, s.2.2.1), s.1.2, s.2.1.2, s.2.2.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

/-- **Split identity.** A006455(n + k) is the sum, over natural orders `P` on `Fin n` and `Q` on
`Fin k`, of the number of order-preserving maps from `Q` to the down-sets of `P`. -/
theorem a_add (n k : ℕ) :
    a (n + k) = ∑ P : {P : Fin n → Fin n → Bool // IsNat P},
      ∑ Q : {Q : Fin k → Fin k → Bool // IsNat Q},
        Fintype.card {f : Fin k → Fin n → Bool // Glue P.1 Q.1 f} := by
  rw [a, Fintype.card_congr ((equiv n k).trans (goodEquivSigma n k)), Fintype.card_sigma]
  simp only [Fintype.card_sigma]

/-! ## Weights of unlabeled posets

The program sums over one representative `s` per isomorphism class, with weight
`e(s) / |Aut s|`. The next theorem shows that this weight is the number of natural orders on
`Fin n` isomorphic to `s`: `e(s) = |Aut s| · #{natural orders isomorphic to s}`.
-/

section weights

variable {α : Type*} [Fintype α] [DecidableEq α] {m : ℕ}

/-- The relation `s` transported to `Fin m` along `σ`. -/
def push (σ : α ≃ Fin m) (s : α → α → Bool) : Fin m → Fin m → Bool :=
  fun u v => s (σ.symm u) (σ.symm v)

/-- Linear extensions of `s`, as bijections onto `Fin m` that respect `s`. -/
def LinExt (s : α → α → Bool) (m : ℕ) := {σ : α ≃ Fin m // ∀ x y, s x y = true → σ x < σ y}

/-- Automorphisms of `s`. -/
def Aut (s : α → α → Bool) := {τ : α ≃ α // ∀ x y, s (τ x) (τ y) = s x y}

/-- Natural orders on `Fin m` isomorphic to `s`. -/
def Labelings (s : α → α → Bool) (m : ℕ) :=
  {r : Fin m → Fin m → Bool // IsNat r ∧ ∃ σ : α ≃ Fin m, push σ s = r}

instance (s : α → α → Bool) : Fintype (LinExt s m) := by unfold LinExt; infer_instance
instance (s : α → α → Bool) : Fintype (Aut s) := by unfold Aut; infer_instance
instance (s : α → α → Bool) : Fintype (Labelings s m) := by unfold Labelings; infer_instance

omit [Fintype α] [DecidableEq α] in
@[simp] lemma push_apply (σ : α ≃ Fin m) (s : α → α → Bool) (u v : Fin m) :
    push σ s u v = s (σ.symm u) (σ.symm v) := rfl

omit [Fintype α] [DecidableEq α] in
lemma push_isNat {s : α → α → Bool} (hs : ∀ x y z, s x y = true → s y z = true → s x z = true)
    (σ : LinExt s m) : IsNat (push σ.1 s) := by
  refine ⟨fun u v h => ?_, fun u v w h1 h2 => hs _ _ _ h1 h2⟩
  simpa using σ.2 _ _ h

omit [Fintype α] [DecidableEq α] in
/-- A bijection whose transported order is natural is a linear extension. -/
lemma linExt_of_push {s : α → α → Bool} {σ : α ≃ Fin m} (h : IsNat (push σ s)) :
    ∀ x y, s x y = true → σ x < σ y := by
  intro x y hxy
  exact h.1 _ _ (by simpa using hxy)

/-- **Weight lemma.** For a strict order `s` (transitivity is all that is needed),
`e(s) = |Aut s| · #{natural orders on Fin m isomorphic to s}`. -/
theorem card_linExt {s : α → α → Bool}
    (hs : ∀ x y z, s x y = true → s y z = true → s x z = true) :
    Fintype.card (LinExt s m) = Fintype.card (Aut s) * Fintype.card (Labelings s m) := by
  classical
  let f : LinExt s m → Labelings s m := fun σ => ⟨push σ.1 s, push_isNat hs σ, σ.1, rfl⟩
  rw [← Fintype.card_congr (Equiv.sigmaFiberEquiv f), Fintype.card_sigma]
  have hfib : ∀ L : Labelings s m, Fintype.card {σ // f σ = L} = Fintype.card (Aut s) := by
    rintro ⟨r, hr, σ0, h0⟩
    subst h0
    refine Fintype.card_congr
      { toFun := fun σ => ⟨σ.1.1.trans σ0.symm, fun x y => ?_⟩
        invFun := fun τ => ⟨⟨τ.1.trans σ0, fun x y hxy => ?_⟩, ?_⟩
        left_inv := fun σ => ?_
        right_inv := fun τ => ?_ }
    · have := congrFun (congrFun (congrArg Subtype.val σ.2) (σ.1.1 x)) (σ.1.1 y)
      simpa [f] using this.symm
    · exact linExt_of_push hr _ _ (by simpa [τ.2] using hxy)
    · apply Subtype.ext
      funext u v
      simp only [f]
      have := τ.2 (τ.1.symm (σ0.symm u)) (τ.1.symm (σ0.symm v))
      simpa using this.symm
    · apply Subtype.ext; apply Subtype.ext; ext x; simp
    · apply Subtype.ext; ext x; simp
  simp [hfib, Finset.sum_const, Finset.card_univ, mul_comm]

end weights

/-! Small values, checked by the kernel. -/

set_option maxRecDepth 100000 in
theorem a_three : a 3 = 7 := by decide +kernel

-- The split identity instantiated at n = 1, k = 2 (the sum itself is not evaluated).
theorem a_add_1_2 : (∑ P : {P : Fin 1 → Fin 1 → Bool // IsNat P},
    ∑ Q : {Q : Fin 2 → Fin 2 → Bool // IsNat Q},
      Fintype.card {f : Fin 2 → Fin 1 → Bool // Glue P.1 Q.1 f}) = 7 := by
  rw [← a_three, show (3 : ℕ) = 1 + 2 from rfl, a_add]

-- Expect only the standard axioms propext, Classical.choice, Quot.sound.
#print axioms a_add
#print axioms a_three
#print axioms card_linExt

end A006455
