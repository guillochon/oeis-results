import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Finset.Card
import Mathlib.SetTheory.Cardinal.Finite
import Mathlib.Tactic

/-!
# A109459 / A109458 / A109457: the structure theorem for Krom (2-SAT) functions

A *Krom function* on a set of variables `V` is a set `S` of assignments `V → Bool` that is closed
under bitwise majority; equivalently (Krom, Schaefer) the solution set of a 2-CNF. The counts
A109457 (labelled), A109458 (up to permuting variables) and A109459 (up to permuting and
complementing variables) were extended to `n = 10` by enumerating "skew posets" up to
isomorphism. This file proves the mathematics that justifies that reduction:

* `sol_valid_eq`: a **nonempty majority-closed** `S` equals the solution set of the implications
  `l → m` between literals that hold on `S` (the 2-clauses valid on `S`);
* `valid_sol_iff`: for a **skew preorder** `R` on literals (reflexive, transitive, reversed by
  negation, and with no literal implying its own negation) the implications valid on `Sol R` are
  exactly `R`, and `exists_sat` shows that no literal is forced;
* `kromEquivSkew`: hence the Krom functions with no forced variable are in bijection with the
  skew preorders on the literals. The equivalence classes of a skew preorder are the blocks of
  equivalent literals of the attack write-up, and its quotient is the skew poset that the
  program enumerates.

Everything is over an arbitrary finite type of variables. Standard axioms only.
-/

namespace A109459

variable {V : Type*}

/-- A literal: a variable together with the value it asserts. -/
abbrev Lit (V : Type*) := V × Bool

/-- Negation of a literal. -/
def neg (l : Lit V) : Lit V := (l.1, !l.2)

@[simp] lemma neg_neg (l : Lit V) : neg (neg l) = l := by simp [neg]
@[simp] lemma neg_fst (l : Lit V) : (neg l).1 = l.1 := rfl
@[simp] lemma neg_snd (l : Lit V) : (neg l).2 = !l.2 := rfl
@[simp] lemma neg_mk (v : V) (b : Bool) : neg (v, b) = (v, !b) := rfl

/-- `Sat x l`: the assignment `x` makes the literal `l` true. -/
def Sat (x : V → Bool) (l : Lit V) : Prop := x l.1 = l.2

lemma sat_neg_iff (x : V → Bool) (l : Lit V) : Sat x (neg l) ↔ ¬ Sat x l := by
  obtain ⟨v, b⟩ := l
  simp only [Sat, neg]
  cases b <;> cases x v <;> simp

/-- Majority of three booleans. -/
def maj3 (a b c : Bool) : Bool := (a && b) || (b && c) || (a && c)

/-- Bitwise majority of three assignments. -/
def maj (x y z : V → Bool) : V → Bool := fun v => maj3 (x v) (y v) (z v)

lemma maj3_of_two {a b c d : Bool}
    (h : (a = d ∧ b = d) ∨ (b = d ∧ c = d) ∨ (a = d ∧ c = d)) : maj3 a b c = d := by
  cases a <;> cases b <;> cases c <;> cases d <;> simp_all [maj3]

lemma two_of_maj3 {a b c d : Bool} (h : maj3 a b c = d) :
    (a = d ∧ b = d) ∨ (b = d ∧ c = d) ∨ (a = d ∧ c = d) := by
  cases a <;> cases b <;> cases c <;> cases d <;> simp_all [maj3]

/-- `S` is closed under bitwise majority: the defining property of a Krom function. -/
def MajClosed (S : Set (V → Bool)) : Prop := ∀ x ∈ S, ∀ y ∈ S, ∀ z ∈ S, maj x y z ∈ S

/-- The implication `l → m` holds at every point of `S`, i.e. the 2-clause `¬l ∨ m` is valid
on `S`. Unit clauses are the case `m = neg l`. -/
def Valid (S : Set (V → Bool)) (l m : Lit V) : Prop := ∀ x ∈ S, Sat x l → Sat x m

/-- The solution set of a set of implications between literals. -/
def Sol (R : Lit V → Lit V → Prop) : Set (V → Bool) := {x | ∀ l m, R l m → Sat x l → Sat x m}

lemma mem_sol {R : Lit V → Lit V → Prop} {x : V → Bool} :
    x ∈ Sol R ↔ ∀ l m, R l m → Sat x l → Sat x m := Iff.rfl

/-- Solution sets of 2-clauses are majority-closed (the easy direction of Krom's theorem). -/
theorem majClosed_sol (R : Lit V → Lit V → Prop) : MajClosed (Sol R) := by
  intro x hx y hy z hz l m hlm hl
  have hl' := two_of_maj3 hl
  apply maj3_of_two
  rcases hl' with ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exact Or.inl ⟨hx l m hlm h1, hy l m hlm h2⟩
  · exact Or.inr (Or.inl ⟨hy l m hlm h1, hz l m hlm h2⟩)
  · exact Or.inr (Or.inr ⟨hx l m hlm h1, hz l m hlm h2⟩)

/-! ### Krom's theorem: majority-closed sets are 2-CNF definable -/

/-- If `x` satisfies every 2-clause valid on `S`, then for any two coordinates some point of `S`
agrees with `x` on both. -/
theorem pair_agree {S : Set (V → Bool)} {x : V → Bool} (hx : x ∈ Sol (Valid S)) (a b : V) :
    ∃ y ∈ S, y a = x a ∧ y b = x b := by
  by_contra h
  push Not at h
  have hv : Valid S (a, x a) (b, !x b) := by
    intro y hy hya
    simp only [Sat] at hya ⊢
    have := h y hy hya
    cases hxb : x b <;> cases hyb : y b <;> simp_all
  have := hx _ _ hv rfl
  simp [Sat] at this

/-- The key induction: a point of `S` agrees with `x` on any finite set of coordinates. -/
theorem agree_on {S : Set (V → Bool)} (hM : MajClosed S) (hS : S.Nonempty) {x : V → Bool}
    (hx : x ∈ Sol (Valid S)) [DecidableEq V] (A : Finset V) :
    ∃ y ∈ S, ∀ v ∈ A, y v = x v := by
  induction A using Finset.strongInduction with
  | H A ih =>
    rcases A.eq_empty_or_nonempty with rfl | ⟨a, ha⟩
    · obtain ⟨y, hy⟩ := hS
      exact ⟨y, hy, by simp⟩
    rcases (A.erase a).eq_empty_or_nonempty with he | ⟨b, hb⟩
    · have hA : ∀ v ∈ A, v = a := by
        intro v hv
        by_contra hva
        have : v ∈ A.erase a := Finset.mem_erase.2 ⟨hva, hv⟩
        simp [he] at this
      obtain ⟨y, hy, hya, -⟩ := pair_agree hx a a
      exact ⟨y, hy, fun v hv => by rw [hA v hv]; exact hya⟩
    · have hba : b ≠ a := (Finset.mem_erase.1 hb).1
      have hbA : b ∈ A := (Finset.mem_erase.1 hb).2
      obtain ⟨y, hy, hyA⟩ := ih (A.erase a) (Finset.erase_ssubset ha)
      obtain ⟨z, hz, hzA⟩ := ih (A.erase b) (Finset.erase_ssubset hbA)
      obtain ⟨w, hw, hwa, hwb⟩ := pair_agree hx a b
      refine ⟨maj y z w, hM y hy z hz w hw, ?_⟩
      intro v hv
      by_cases hva : v = a
      · subst hva
        exact maj3_of_two (Or.inr (Or.inl ⟨hzA v (Finset.mem_erase.2 ⟨hba.symm, ha⟩), hwa⟩))
      by_cases hvb : v = b
      · subst hvb
        exact maj3_of_two (Or.inr (Or.inr ⟨hyA v (Finset.mem_erase.2 ⟨hba, hbA⟩), hwb⟩))
      exact maj3_of_two (Or.inl ⟨hyA v (Finset.mem_erase.2 ⟨hva, hv⟩),
        hzA v (Finset.mem_erase.2 ⟨hvb, hv⟩)⟩)

/-- **Krom's theorem.** A nonempty majority-closed set is the solution set of the implications
(2-clauses) valid on it. -/
theorem sol_valid_eq [Fintype V] [DecidableEq V] {S : Set (V → Bool)} (hM : MajClosed S)
    (hS : S.Nonempty) : Sol (Valid S) = S := by
  ext x
  constructor
  · intro hx
    obtain ⟨y, hy, hyx⟩ := agree_on hM hS hx Finset.univ
    have : y = x := funext fun v => hyx v (Finset.mem_univ v)
    exact this ▸ hy
  · intro hx l m hlm hl
    exact hlm x hx hl

/-! ### Skew preorders and their solution sets -/

/-- A skew preorder on the literals: reflexive, transitive, reversed by negation, and no literal
implies its own negation (no forced literal). It is the closed implication structure of a 2-CNF
with no forced variable. -/
structure IsSkew (R : Lit V → Lit V → Prop) : Prop where
  refl : ∀ l, R l l
  trans : ∀ {l m k}, R l m → R m k → R l k
  contra : ∀ {l m}, R l m → R (neg m) (neg l)
  cons : ∀ l, ¬ R l (neg l)

/-- `U` is closed upwards under `R`. -/
def Closed (R : Lit V → Lit V → Prop) (U : Set (Lit V)) : Prop := ∀ l ∈ U, ∀ m, R l m → m ∈ U

/-- `U` contains no complementary pair. -/
def Consistent (U : Set (Lit V)) : Prop := ∀ l ∈ U, neg l ∉ U

/-- `U` decides the variable `v`. -/
def Assigned (U : Set (Lit V)) (v : V) : Prop := (v, true) ∈ U ∨ (v, false) ∈ U

/-- The assignment read off from a set of literals that decides every variable. -/
lemma sat_decide_iff {U : Set (Lit V)} (hcons : Consistent U) (hass : ∀ v, Assigned U v)
    [DecidablePred (· ∈ U)] (l : Lit V) :
    Sat (fun v => decide ((v, true) ∈ U)) l ↔ l ∈ U := by
  obtain ⟨v, b⟩ := l
  simp only [Sat]
  cases b
  · simp only [decide_eq_false_iff_not]
    constructor
    · intro h
      rcases hass v with h' | h'
      · exact absurd h' h
      · exact h'
    · intro h h'
      exact hcons _ h' (by simpa using h)
  · simp

/-- **Extension lemma.** A closed, consistent set of literals extends to a full solution of `R`.
The induction is over the set `D` of variables still allowed to be undecided. -/
theorem extend {R : Lit V → Lit V → Prop} (hR : IsSkew R) [DecidableEq V] (D : Finset V) :
    ∀ U : Set (Lit V), Closed R U → Consistent U → (∀ v, v ∉ D → Assigned U v) →
      ∃ x ∈ Sol R, ∀ l ∈ U, Sat x l := by
  classical
  induction D using Finset.induction_on with
  | empty =>
    intro U hc hcons hass
    have hass' : ∀ v, Assigned U v := fun v => hass v (by simp)
    refine ⟨fun v => decide ((v, true) ∈ U), ?_, ?_⟩
    · intro l m hlm hl
      exact (sat_decide_iff hcons hass' m).2 (hc l ((sat_decide_iff hcons hass' l).1 hl) m hlm)
    · intro l hl
      exact (sat_decide_iff hcons hass' l).2 hl
  | insert v D hv ih =>
    intro U hc hcons hass
    by_cases hA : Assigned U v
    · refine ih U hc hcons ?_
      intro w hw
      by_cases hwv : w = v
      · exact hwv ▸ hA
      · exact hass w (by simp [hwv, hw])
    · -- make `v` true and close up
      let U' : Set (Lit V) := U ∪ {m | R (v, true) m}
      have hU : ∀ l, l ∈ U' ↔ l ∈ U ∨ R (v, true) l := fun l => Iff.rfl
      have hc' : Closed R U' := by
        intro l hl m hlm
        rcases (hU l).1 hl with h | h
        · exact (hU m).2 (Or.inl (hc l h m hlm))
        · exact (hU m).2 (Or.inr (hR.trans h hlm))
      have hnf : (v, false) ∉ U := fun h => hA (Or.inr h)
      have hcons' : Consistent U' := by
        intro l hl hnl
        rcases (hU l).1 hl with h1 | h1 <;> rcases (hU (neg l)).1 hnl with h2 | h2
        · exact hcons l h1 h2
        · -- l ∈ U, (v,true) ≤ neg l  ⟹  l ≤ (v,false) ∈ U
          have := hR.contra h2
          rw [neg_neg] at this
          exact hnf (hc l h1 _ this)
        · -- (v,true) ≤ l, neg l ∈ U  ⟹  neg l ≤ (v,false) ∈ U
          have := hR.contra h1
          exact hnf (hc _ h2 _ this)
        · -- (v,true) ≤ l and (v,true) ≤ neg l  ⟹  (v,true) ≤ (v,false)
          have := hR.contra h2
          rw [neg_neg] at this
          exact hR.cons (v, true) (hR.trans h1 this)
      have hass' : ∀ w, w ∉ D → Assigned U' w := by
        intro w hw
        by_cases hwv : w = v
        · subst hwv
          exact Or.inl ((hU _).2 (Or.inr (hR.refl _)))
        · rcases hass w (by simp [hwv, hw]) with h | h
          · exact Or.inl ((hU _).2 (Or.inl h))
          · exact Or.inr ((hU _).2 (Or.inl h))
      obtain ⟨x, hx, hxU⟩ := ih U' hc' hcons' hass'
      exact ⟨x, hx, fun l hl => hxU l ((hU l).2 (Or.inl hl))⟩

/-- The up-set of a literal is closed. -/
lemma closed_upset {R : Lit V → Lit V → Prop} (hR : IsSkew R) (l : Lit V) :
    Closed R {k | R l k} := fun _ hk _ hkm => hR.trans hk hkm

/-- The up-set of a literal is consistent. -/
lemma consistent_upset {R : Lit V → Lit V → Prop} (hR : IsSkew R) (l : Lit V) :
    Consistent {k | R l k} := by
  intro k hk hnk
  have := hR.contra hnk
  rw [neg_neg] at this
  exact hR.cons l (hR.trans hk this)

/-- No literal is forced on `Sol R`. -/
theorem exists_sat {R : Lit V → Lit V → Prop} (hR : IsSkew R) [Fintype V] [DecidableEq V]
    (l : Lit V) : ∃ x ∈ Sol R, Sat x l := by
  obtain ⟨x, hx, hxU⟩ := extend hR Finset.univ {k | R l k} (closed_upset hR l)
    (consistent_upset hR l) (fun v hv => absurd (Finset.mem_univ v) hv)
  exact ⟨x, hx, hxU l (hR.refl l)⟩

/-- `Sol R` is nonempty (also when `V` is empty). -/
theorem sol_nonempty {R : Lit V → Lit V → Prop} (hR : IsSkew R) [Fintype V] [DecidableEq V] :
    (Sol R).Nonempty := by
  obtain ⟨x, hx, -⟩ := extend hR Finset.univ ∅ (fun _ h => absurd h (Set.notMem_empty _))
    (fun _ h => absurd h (Set.notMem_empty _)) (fun v hv => absurd (Finset.mem_univ v) hv)
  exact ⟨x, hx⟩

/-- **The implications valid on `Sol R` are exactly `R`.** -/
theorem valid_sol_iff {R : Lit V → Lit V → Prop} (hR : IsSkew R) [Fintype V] [DecidableEq V]
    (l m : Lit V) : Valid (Sol R) l m ↔ R l m := by
  constructor
  · intro hval
    by_contra hnot
    let U : Set (Lit V) := {k | R l k} ∪ {k | R (neg m) k}
    have hU : ∀ k, k ∈ U ↔ R l k ∨ R (neg m) k := fun k => Iff.rfl
    have hc : Closed R U := by
      intro k hk j hkj
      rcases (hU k).1 hk with h | h
      · exact (hU j).2 (Or.inl (hR.trans h hkj))
      · exact (hU j).2 (Or.inr (hR.trans h hkj))
    have hcons : Consistent U := by
      intro k hk hnk
      rcases (hU k).1 hk with h1 | h1 <;> rcases (hU (neg k)).1 hnk with h2 | h2
      · exact consistent_upset hR l k h1 h2
      · -- l ≤ k and neg m ≤ neg k  ⟹  k ≤ m  ⟹  l ≤ m
        have := hR.contra h2
        rw [neg_neg, neg_neg] at this
        exact hnot (hR.trans h1 this)
      · -- neg m ≤ k and l ≤ neg k  ⟹  k ≤ neg l  ⟹  neg m ≤ neg l  ⟹  l ≤ m
        have := hR.contra h2
        rw [neg_neg] at this
        have := hR.contra (hR.trans h1 this)
        rw [neg_neg, neg_neg] at this
        exact hnot this
      · exact consistent_upset hR (neg m) k h1 h2
    obtain ⟨x, hx, hxU⟩ := extend hR Finset.univ U hc hcons
      (fun v hv => absurd (Finset.mem_univ v) hv)
    have hl : Sat x l := hxU l ((hU l).2 (Or.inl (hR.refl l)))
    have hm : Sat x (neg m) := hxU (neg m) ((hU _).2 (Or.inr (hR.refl _)))
    exact (sat_neg_iff x m).1 hm (hval x hx hl)
  · intro h x hx hl
    exact hx l m h hl

/-! ### The bijection -/

/-- The implication relation of a nonforced Krom set is a skew preorder. -/
theorem isSkew_valid {S : Set (V → Bool)} (hnf : ∀ l, ∃ x ∈ S, Sat x l) : IsSkew (Valid S) where
  refl := fun _ _ _ h => h
  trans := fun h1 h2 x hx hl => h2 x hx (h1 x hx hl)
  contra := fun {l m} h x hx hm => by
    rw [sat_neg_iff] at hm ⊢
    exact fun hl => hm (h x hx hl)
  cons := fun l h => by
    obtain ⟨x, hx, hl⟩ := hnf l
    exact (sat_neg_iff x l).1 (h x hx hl) hl

/-- Krom functions on `V` with no forced variable. -/
abbrev KromNF (V : Type*) : Type _ :=
  {S : Set (V → Bool) // MajClosed S ∧ S.Nonempty ∧ ∀ l, ∃ x ∈ S, Sat x l}

/-- Skew preorders on the literals of `V`. -/
abbrev Skew (V : Type*) : Type _ := {R : Lit V → Lit V → Prop // IsSkew R}

/-- **Structure theorem.** Krom functions with no forced variable correspond bijectively to
skew preorders on the literals, via `S ↦ Valid S` and `R ↦ Sol R`. -/
def kromEquivSkew [Fintype V] [DecidableEq V] : KromNF V ≃ Skew V where
  toFun S := ⟨Valid S.1, isSkew_valid S.2.2.2⟩
  invFun R := ⟨Sol R.1, majClosed_sol _, sol_nonempty R.2, exists_sat R.2⟩
  left_inv S := Subtype.ext (sol_valid_eq S.2.1 S.2.2.1)
  right_inv R := Subtype.ext (funext fun l => funext fun m => propext (valid_sol_iff R.2 l m))

/-- The two sides have the same number of elements. -/
theorem card_kromNF_eq [Fintype V] [DecidableEq V] : Nat.card (KromNF V) = Nat.card (Skew V) :=
  Nat.card_congr kromEquivSkew

/-! ### The forced-variable layer: all Krom functions -/

section Forced

/-- `ok f x`: the assignment `x` takes the forced values prescribed by `f`. -/
def ok (f : V → Option Bool) (x : V → Bool) : Prop := ∀ v b, f v = some b → x v = b

/-- Restriction of an assignment to the unforced variables. -/
def restr (f : V → Option Bool) (x : V → Bool) : {v // f v = none} → Bool := fun v => x v.1

/-- Extension of an assignment of the unforced variables by the forced values. -/
def ext (f : V → Option Bool) (x' : {v // f v = none} → Bool) : V → Bool :=
  fun v => if h : f v = none then x' ⟨v, h⟩ else (f v).getD false

lemma restr_ext (f : V → Option Bool) (x' : {v // f v = none} → Bool) :
    restr f (ext f x') = x' := by
  funext ⟨v, hv⟩
  simp [restr, ext, hv]

lemma ok_ext (f : V → Option Bool) (x' : {v // f v = none} → Bool) : ok f (ext f x') := by
  intro v b hb
  simp [ext, hb]

lemma ext_restr {f : V → Option Bool} {x : V → Bool} (hx : ok f x) : ext f (restr f x) = x := by
  funext v
  by_cases h : f v = none
  · simp [ext, restr, h]
  · obtain ⟨b, hb⟩ := Option.ne_none_iff_exists'.1 h
    simp [ext, hb, hx v b hb]

lemma restr_maj (f : V → Option Bool) (x y z : V → Bool) :
    restr f (maj x y z) = maj (restr f x) (restr f y) (restr f z) := rfl

lemma ok_maj {f : V → Option Bool} {x y z : V → Bool} (hx : ok f x) (hy : ok f y)
    (hz : ok f z) : ok f (maj x y z) := by
  intro v b hb
  simp only [maj]
  rw [hx v b hb, hy v b hb, hz v b hb]
  cases b <;> rfl

open Classical in
/-- The forced value of `v` on `S`, if any. -/
noncomputable def forced (S : Set (V → Bool)) (v : V) : Option Bool :=
  if ∀ x ∈ S, x v = true then some true else if ∀ x ∈ S, x v = false then some false else none

lemma forced_eq_some_iff {S : Set (V → Bool)} (hS : S.Nonempty) (v : V) (b : Bool) :
    forced S v = some b ↔ ∀ x ∈ S, x v = b := by
  unfold forced
  obtain ⟨y, hy⟩ := hS
  by_cases h1 : ∀ x ∈ S, x v = true
  · rw [if_pos h1]
    constructor
    · intro h
      cases b
      · simp at h
      · exact h1
    · intro h
      cases b
      · have := h y hy
        have := h1 y hy
        simp_all
      · rfl
  · rw [if_neg h1]
    by_cases h2 : ∀ x ∈ S, x v = false
    · rw [if_pos h2]
      constructor
      · intro h
        cases b
        · exact h2
        · simp at h
      · intro h
        cases b
        · rfl
        · exact absurd h h1
    · rw [if_neg h2]
      constructor
      · intro h
        simp at h
      · intro h
        cases b
        · exact absurd h h2
        · exact absurd h h1

lemma forced_eq_none_iff {S : Set (V → Bool)} (hS : S.Nonempty) (v : V) :
    forced S v = none ↔ ∀ b, ∃ x ∈ S, x v = b := by
  constructor
  · intro h b
    by_contra hne
    push Not at hne
    have : ∀ x ∈ S, x v = !b := fun x hx => by
      have := hne x hx
      cases b <;> cases hxv : x v <;> simp_all
    rw [← forced_eq_some_iff hS] at this
    rw [h] at this
    simp at this
  · intro h
    unfold forced
    obtain ⟨x, hx, hxv⟩ := h false
    obtain ⟨y, hy, hyv⟩ := h true
    have h1 : ¬ ∀ z ∈ S, z v = true := fun h1 => by simpa [hxv] using h1 x hx
    have h2 : ¬ ∀ z ∈ S, z v = false := fun h2 => by simpa [hyv] using h2 y hy
    simp [h1, h2]

/-- Nonempty Krom functions. -/
abbrev KromNE (V : Type*) : Type _ := {S : Set (V → Bool) // MajClosed S ∧ S.Nonempty}

/-- The set of assignments with the forced values `f` whose restriction lies in `S'`. -/
def build (f : V → Option Bool) (S' : Set ({v // f v = none} → Bool)) : Set (V → Bool) :=
  {x | ok f x ∧ restr f x ∈ S'}

/-- **Forced-variable decomposition.** For a fixed pattern `f` of forced values, the nonempty
Krom functions whose forced variables are exactly those of `f` correspond to the forced-free Krom
functions on the remaining variables. -/
def fiberEquiv (f : V → Option Bool) :
    {S : KromNE V // forced S.1 = f} ≃ KromNF {v // f v = none} where
  toFun S := ⟨restr f '' S.1.1, by
    obtain ⟨⟨S, hM, hS⟩, hf⟩ := S
    refine ⟨?_, hS.image _, ?_⟩
    · rintro _ ⟨x, hx, rfl⟩ _ ⟨y, hy, rfl⟩ _ ⟨z, hz, rfl⟩
      exact ⟨maj x y z, hM x hx y hy z hz, rfl⟩
    · rintro ⟨⟨v, hv⟩, b⟩
      have hv' : forced S v = none := hf ▸ hv
      obtain ⟨x, hx, hxv⟩ := (forced_eq_none_iff hS v).1 hv' b
      exact ⟨restr f x, ⟨x, hx, rfl⟩, hxv⟩⟩
  invFun S' := ⟨⟨build f S'.1, by
    obtain ⟨S', hM, hS, hnf⟩ := S'
    refine ⟨?_, ?_⟩
    · rintro x ⟨hx, hx'⟩ y ⟨hy, hy'⟩ z ⟨hz, hz'⟩
      exact ⟨ok_maj hx hy hz, by rw [restr_maj]; exact hM _ hx' _ hy' _ hz'⟩
    · obtain ⟨x', hx'⟩ := hS
      exact ⟨ext f x', ok_ext f x', by rw [restr_ext]; exact hx'⟩⟩, by
    obtain ⟨S', hM, hS, hnf⟩ := S'
    have hne : (build f S').Nonempty := by
      obtain ⟨x', hx'⟩ := hS
      exact ⟨ext f x', ok_ext f x', by rw [restr_ext]; exact hx'⟩
    funext v
    cases hfv : f v with
    | some b =>
      rw [forced_eq_some_iff hne]
      intro x hx
      exact hx.1 v b hfv
    | none =>
      rw [forced_eq_none_iff hne]
      intro b
      obtain ⟨x', hx', hxv⟩ := hnf (⟨v, hfv⟩, b)
      refine ⟨ext f x', ⟨ok_ext f x', by rw [restr_ext]; exact hx'⟩, ?_⟩
      simpa [ext, hfv, Sat] using hxv⟩
  left_inv S := by
    obtain ⟨⟨S, hM, hS⟩, hf⟩ := S
    apply Subtype.ext
    apply Subtype.ext
    ext x
    constructor
    · rintro ⟨hx, y, hy, hxy⟩
      have : x = y := by
        funext v
        cases hfv : f v with
        | some b =>
          have hyv : y v = b :=
            ((forced_eq_some_iff hS v b).1 (hf ▸ hfv : forced S v = some b)) y hy
          rw [hx v b hfv, hyv]
        | none =>
          have := congrFun hxy ⟨v, hfv⟩
          exact this.symm
      exact this ▸ hy
    · intro hx
      refine ⟨?_, x, hx, rfl⟩
      intro v b hb
      exact ((forced_eq_some_iff hS v b).1 (hf ▸ hb : forced S v = some b)) x hx
  right_inv S' := by
    obtain ⟨S', hM, hS, hnf⟩ := S'
    apply Subtype.ext
    ext x'
    constructor
    · rintro ⟨x, ⟨-, hx'⟩, rfl⟩
      exact hx'
    · intro hx'
      exact ⟨ext f x', ⟨ok_ext f x', by rw [restr_ext]; exact hx'⟩, restr_ext f x'⟩

/-- All Krom functions, including the empty one. -/
abbrev Krom (V : Type*) : Type _ := {S : Set (V → Bool) // MajClosed S}

open Classical in
/-- A Krom function is empty or a nonempty Krom function. -/
noncomputable def kromEquivSum : Krom V ≃ KromNE V ⊕ Unit where
  toFun S := if h : S.1.Nonempty then Sum.inl ⟨S.1, S.2, h⟩ else Sum.inr ()
  invFun := fun
    | Sum.inl S => ⟨S.1, S.2.1⟩
    | Sum.inr _ => ⟨∅, fun x hx => absurd hx (Set.notMem_empty x)⟩
  left_inv S := by
    by_cases h : S.1.Nonempty
    · simp only [dif_pos h]
    · simp only [dif_neg h]
      exact Subtype.ext (Set.not_nonempty_iff_eq_empty.1 h).symm
  right_inv := fun
    | Sum.inl S => by simp only [dif_pos S.2.2]
    | Sum.inr () => by simp only [dif_neg Set.not_nonempty_empty]

instance [Finite V] : Finite (KromNE V) := Subtype.finite
instance [Finite V] : Finite (Krom V) := Subtype.finite
instance [Finite V] : Finite (KromNF V) := Subtype.finite
instance [Finite V] : Finite (Skew V) := Subtype.finite

/-- **Counting Krom functions.** The number of Krom functions on `V` is one (the empty function)
plus, for every pattern `f` of forced values, the number of skew preorders on the literals of
the unforced variables. Summing over the patterns with `k` forced variables gives the
`binomial(n,k) 2^k` of the A109457 formula. -/
theorem card_krom [Fintype V] [DecidableEq V] :
    Nat.card (Krom V) = 1 + ∑ f : V → Option Bool, Nat.card (Skew {v // f v = none}) := by
  rw [Nat.card_congr kromEquivSum, Nat.card_sum, Nat.card_unique (α := Unit), add_comm]
  congr 1
  rw [← Nat.card_congr (Equiv.sigmaFiberEquiv (fun S : KromNE V => forced S.1)), Nat.card_sigma]
  refine Finset.sum_congr rfl fun f _ => ?_
  rw [Nat.card_congr (fiberEquiv f), card_kromNF_eq]

end Forced

end A109459

