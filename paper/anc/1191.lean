/-
Copyright 2026 Lucas Parobek.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

# Erdős Problem #1191 — Lean 4 Formal Verification Suite
Paper: "Unconditional Diophantine Haar Closure, Exact Harmonic-Difference Conservation,
        and Multi-Octave Dyadic Budgets for Erdős Problem #1191" (Parobek, 2026)
Repo:  https://github.com/doublegate/Erdos_1191

This self-contained Lean 4 module (zero `sorry`s, zero custom axioms) machine-checks the
discrete arithmetic lemmas, dyadic interval bounds, powers-of-two recursive sum identity and
bound, KKT integer quadratic inequality, binary sub-block variance identity, discrete
Legendre–Fenchel remainder identities, Order-2 Cesàro–Fejér recurrence step, exponent/power
conservation identities, and finite Sidon set witness checks used in [Pa26].
-/

namespace Erdos1191

local notation "ℕ" => Nat

/-! ## 1. Theorem 2.1(a): Exact Harmonic-Integral Conservation Algebraic Kernel -/

/--
**Arithmetic Lemma for Theorem 2.1(a) (Harmonic-Difference Numerator Identity):**
For any natural numbers `1 ≤ a_i < a_j` with difference `d = a_j - a_i`,
the numerator identity `(d + 1) + (a_i - 1) = a_j` holds on `ℕ` (used to derive the
partial-fraction decomposition `1/a_j + (a_i - 1)/((a_j - a_i + 1) a_j) = 1/(a_j - a_i + 1)`
in Theorem 2.1(a) of [Pa26]).
-/
theorem erdos_1191.variants.harmonic_conservation_numerator (a_i a_j d : ℕ)
    (h_pos : 1 ≤ a_i) (h_lt : a_i < a_j) (h_diff : d = a_j - a_i) :
    (d + 1) + (a_i - 1) = a_j := by
  omega

/-! ## 2. Theorem 2.2(a)–(b), Remark 2.3 & Lemma 2.4: Cross-Shell Localization, Powers-of-Two Verification & KKT Convexity -/

/--
**Arithmetic Lemma for Theorem 2.2(a) (Two-Shell Difference Localization Bound):**
For natural numbers `x, y, pow_km2` satisfying `1 ≤ x ≤ pow_km2` and
`2 · pow_km2 < y ≤ 4 · pow_km2`, the difference `y - x` satisfies
`pow_km2 < y - x < 4 · pow_km2` (used to localize cross-shell differences
into `D_k-1 ∪ D_k` when `pow_km2 = 2^k-2` in Theorem 2.2(a) of [Pa26]).
-/
theorem erdos_1191.variants.two_shell_difference_localization (x y pow_km2 : ℕ)
    (hx_pos : 1 ≤ x)
    (hx_hi : x ≤ pow_km2)
    (hy_lo : 2 * pow_km2 < y)
    (hy_hi : y ≤ 4 * pow_km2) :
    pow_km2 < y - x ∧ y - x < 4 * pow_km2 := by
  omega

/--
Recursive integer sum `S(0) = 0`, `S(M + 1) = 2 S(M) + 2(M + 1)`, representing the
`2^M`-scaled sum `sum_{k=1..M} k · 2^{M - k + 1}` evaluated on the powers-of-two
sequence `A = {2^n : n ≥ 0}` (`m_k = 1, A(2^k-1) = k`) in Remark 2.3 of [Pa26].
-/
def powersOfTwoScaledBilinear : ℕ → ℕ
  | 0 => 0
  | M + 1 => 2 * powersOfTwoScaledBilinear M + 2 * (M + 1)

/--
**Arithmetic Lemma for Remark 2.3 (Closed-Form Identity for `powersOfTwoScaledBilinear`):**
For every `M ≥ 0`, the recursively defined integer sum `powersOfTwoScaledBilinear M`
satisfies the exact closed-form identity
`powersOfTwoScaledBilinear(M) + 2M + 4 = 2^{M + 2}`
(used to evaluate the powers-of-two bilinear sum in Remark 2.3 of [Pa26]).
-/
theorem erdos_1191.variants.powers_of_two_bilinear_identity (M : ℕ) :
    powersOfTwoScaledBilinear M + 2 * M + 4 = 2 ^ (M + 2) := by
  induction M with
  | zero => rfl
  | succ m ih =>
    have h_def : powersOfTwoScaledBilinear (m + 1) = 2 * powersOfTwoScaledBilinear m + 2 * (m + 1) := rfl
    have h_pow : 2 ^ ((m + 1) + 2) = 2 * 2 ^ (m + 2) := by
      change 2 ^ ((m + 2) + 1) = 2 * 2 ^ (m + 2)
      rw [Nat.pow_succ, Nat.mul_comm]
    omega

/--
**Arithmetic Lemma for Theorem 2.2(b) & Remark 2.3 (Upper Bound on `powersOfTwoScaledBilinear`):**
For every `M ≥ 1`, the recursively defined sum `powersOfTwoScaledBilinear M` satisfies
`powersOfTwoScaledBilinear(M) ≤ (M + 1) 2^M - 2`
(certifying `22 ≤ 30` at `M = 3` for the `2^M`-scaled powers-of-two sum in Remark 2.3 of [Pa26]).
-/
theorem erdos_1191.variants.powers_of_two_satisfies_corrected_bilinear_bound (M : ℕ) (hM : 1 ≤ M) :
    powersOfTwoScaledBilinear M ≤ (M + 1) * 2 ^ M - 2 := by
  have h_id := erdos_1191.variants.powers_of_two_bilinear_identity M
  have h_pow_split : 2 ^ (M + 2) = 4 * 2 ^ M := by
    rw [Nat.pow_add, Nat.mul_comm]
  by_cases h1 : M = 1
  · subst h1; decide
  · by_cases h2 : M = 2
    · subst h2; decide
    · have h3 : 3 ≤ M := by omega
      have h_ge4 : 4 * 2 ^ M ≤ (M + 1) * 2 ^ M :=
        Nat.mul_le_mul_right (2 ^ M) (by omega)
      omega

/--
**Arithmetic Lemma for Lemma 2.4 (Quadratic Polynomial Inequality for `k ≥ 4`):**
For every natural number `k ≥ 4`, `16 k^2 + 2(k^2 - 9) = 18(k^2 - 1)` and
`16 k^2 < 18(k^2 - 1)`. In Lemma 2.4 of [Pa26], this integer inequality is used to verify
that `1/(sqrt(k-1)) + 1/(sqrt(k+1)) < 3/sqrt(2k)` for all `k ≥ 4`
(alongside `9 > 8` at `k = 3`) when establishing positivity of the KKT multipliers
`λ_k = 2^1-k(2ρ_k - ρ_k+1)`.
-/
theorem erdos_1191.variants.kkt_increment_convexity_k_ge_4 (k : ℕ) (hk : 4 ≤ k) :
    16 * (k * k) + 2 * (k * k - 9) = 18 * (k * k - 1) ∧
    16 * (k * k) < 18 * (k * k - 1) := by
  have hk2 : 16 ≤ k * k := by
    have h1 : 4 * 4 ≤ k * k := Nat.mul_le_mul hk hk
    exact h1
  omega

/-! ## 3. Theorem 2.5(b): Binary Sub-Block Variance Refinement -/

/--
**Arithmetic Lemma for Theorem 2.5(b) & Theorem 2.11(c) (Binary Sub-Block Variance Refinement Identity):**
For any natural numbers `S_1, S_2 ∈ ℕ`, `2(S_1^2 + S_2^2) = (S_1 + S_2)^2 + |S_1 - S_2|^2`
(encoded in Lean via truncated subtraction on `ℕ` as `(S_1 - S_2)^2 + (S_2 - S_1)^2`),
and hence `(S_1 + S_2)^2 ≤ 2(S_1^2 + S_2^2)`. This provides the arithmetic step used in the
sub-block refinement bound `L(S_1, B) + L(S_2, B) ≥ L(S_1 + S_2, 2B)` of [Pa26].
-/
theorem erdos_1191.variants.binary_subblock_variance_refinement (S1 S2 : ℕ) :
    2 * (S1 * S1 + S2 * S2) =
      (S1 + S2) * (S1 + S2) + (S1 - S2) * (S1 - S2) + (S2 - S1) * (S2 - S1) ∧
    (S1 + S2) * (S1 + S2) ≤ 2 * (S1 * S1 + S2 * S2) := by
  have h_le : S1 ≤ S2 ∨ S2 ≤ S1 := Nat.le_total S1 S2
  rcases h_le with h12 | h21
  · have h_z : S1 - S2 = 0 := Nat.sub_eq_zero_of_le h12
    rw [h_z]
    have h_id : (S1 + S2) * (S1 + S2) + (S2 - S1) * (S2 - S1) = 2 * (S1 * S1 + S2 * S2) := by
      have h_sub : S1 + (S2 - S1) = S2 := Nat.add_sub_of_le h12
      calc
        (S1 + S2) * (S1 + S2) + (S2 - S1) * (S2 - S1)
          = (2 * S1 + (S2 - S1)) * (2 * S1 + (S2 - S1)) + (S2 - S1) * (S2 - S1) := by
            have h_s : S1 + S2 = 2 * S1 + (S2 - S1) := by omega
            rw [h_s]
        _ = 4 * (S1 * S1) + 4 * (S1 * (S2 - S1)) + 2 * ((S2 - S1) * (S2 - S1)) := by
            have h_exp : (2 * S1 + (S2 - S1)) * (2 * S1 + (S2 - S1)) =
                4 * (S1 * S1) + 4 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1) := by
              calc
                (2 * S1 + (S2 - S1)) * (2 * S1 + (S2 - S1))
                  = (2 * S1) * (2 * S1) + (2 * S1) * (S2 - S1) + ((S2 - S1) * (2 * S1) + (S2 - S1) * (S2 - S1)) := by
                    rw [Nat.add_mul, Nat.mul_add, Nat.mul_add]
                _ = 4 * (S1 * S1) + 2 * (S1 * (S2 - S1)) + 2 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1) := by
                    have ha : (2 * S1) * (2 * S1) = (2 * 2) * (S1 * S1) := by ac_rfl
                    have hb : (2 * S1) * (S2 - S1) = 2 * (S1 * (S2 - S1)) := by ac_rfl
                    have hc : (S2 - S1) * (2 * S1) = 2 * (S1 * (S2 - S1)) := by ac_rfl
                    rw [ha, hb, hc]
                    omega
                _ = 4 * (S1 * S1) + 4 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1) := by omega
            omega
        _ = 2 * (S1 * S1 + (S1 * S1 + 2 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1))) := by omega
        _ = 2 * (S1 * S1 + S2 * S2) := by
            have h_s2 : S2 * S2 = S1 * S1 + 2 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1) := by
              calc
                S2 * S2 = (S1 + (S2 - S1)) * (S1 + (S2 - S1)) := by rw [h_sub]
                _ = S1 * S1 + S1 * (S2 - S1) + ((S2 - S1) * S1 + (S2 - S1) * (S2 - S1)) := by
                    rw [Nat.add_mul, Nat.mul_add, Nat.mul_add]
                _ = S1 * S1 + 2 * (S1 * (S2 - S1)) + (S2 - S1) * (S2 - S1) := by
                    rw [Nat.mul_comm (S2 - S1) S1]
                    omega
            rw [h_s2]
    omega
  · have h_z : S2 - S1 = 0 := Nat.sub_eq_zero_of_le h21
    rw [h_z]
    have h_id : (S1 + S2) * (S1 + S2) + (S1 - S2) * (S1 - S2) = 2 * (S1 * S1 + S2 * S2) := by
      have h_sub : S2 + (S1 - S2) = S1 := Nat.add_sub_of_le h21
      calc
        (S1 + S2) * (S1 + S2) + (S1 - S2) * (S1 - S2)
          = (2 * S2 + (S1 - S2)) * (2 * S2 + (S1 - S2)) + (S1 - S2) * (S1 - S2) := by
            have h_s : S1 + S2 = 2 * S2 + (S1 - S2) := by omega
            rw [h_s]
        _ = 4 * (S2 * S2) + 4 * (S2 * (S1 - S2)) + 2 * ((S1 - S2) * (S1 - S2)) := by
            have h_exp : (2 * S2 + (S1 - S2)) * (2 * S2 + (S1 - S2)) =
                4 * (S2 * S2) + 4 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2) := by
              calc
                (2 * S2 + (S1 - S2)) * (2 * S2 + (S1 - S2))
                  = (2 * S2) * (2 * S2) + (2 * S2) * (S1 - S2) + ((S1 - S2) * (2 * S2) + (S1 - S2) * (S1 - S2)) := by
                    rw [Nat.add_mul, Nat.mul_add, Nat.mul_add]
                _ = 4 * (S2 * S2) + 2 * (S2 * (S1 - S2)) + 2 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2) := by
                    have ha : (2 * S2) * (2 * S2) = (2 * 2) * (S2 * S2) := by ac_rfl
                    have hb : (2 * S2) * (S1 - S2) = 2 * (S2 * (S1 - S2)) := by ac_rfl
                    have hc : (S1 - S2) * (2 * S2) = 2 * (S2 * (S1 - S2)) := by ac_rfl
                    rw [ha, hb, hc]
                    omega
                _ = 4 * (S2 * S2) + 4 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2) := by omega
            omega
        _ = 2 * ((S2 * S2 + 2 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2)) + S2 * S2) := by omega
        _ = 2 * (S1 * S1 + S2 * S2) := by
            have h_s1 : S1 * S1 = S2 * S2 + 2 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2) := by
              calc
                S1 * S1 = (S2 + (S1 - S2)) * (S2 + (S1 - S2)) := by rw [h_sub]
                _ = S2 * S2 + S2 * (S1 - S2) + ((S1 - S2) * S2 + (S1 - S2) * (S1 - S2)) := by
                    rw [Nat.add_mul, Nat.mul_add, Nat.mul_add]
                _ = S2 * S2 + 2 * (S2 * (S1 - S2)) + (S1 - S2) * (S1 - S2) := by
                    rw [Nat.mul_comm (S1 - S2) S2]
                    omega
            rw [h_s1]
    omega

/-! ## 4. Theorem 2.6(a)–(d): Discrete Legendre–Fenchel Sawtooth Integrality Excess & Scaling -/

/--
**Theorem 2.6(a) (Discrete Legendre–Fenchel Remainder Decomposition):**
For any bucket count `B ≥ 1` and total occupancy `S = B q + r` with integer quotient `q`
and remainder `r ≤ B`, the minimum ordered-pair sum over `B` integer buckets summing to `S`
satisfies the exact algebraic decomposition
`(B - r) q(q - 1) + r(q + 1) q = B q(q - 1) + 2 r q`.
-/
theorem erdos_1191.variants.discrete_legendre_fenchel_excess (B q r : ℕ)
    (hr : r ≤ B) :
    (B - r) * (q * (q - 1)) + r * ((q + 1) * q) =
      B * (q * (q - 1)) + 2 * (r * q) := by
  have h1 : (q + 1) * q = q * (q - 1) + 2 * q := by
    cases q with
    | zero => rfl
    | succ m =>
      change (m + 2) * (m + 1) = (m + 1) * m + 2 * (m + 1)
      rw [Nat.add_mul, Nat.mul_comm m (m + 1)]
  rw [h1, Nat.mul_add, ← Nat.add_assoc, ← Nat.add_mul, Nat.sub_add_cancel hr, Nat.mul_left_comm r 2 q]

/--
**Theorem 2.6(b) (Squared-Total Discrete Legendre–Fenchel Integrality Identity):**
Multiplying the minimum ordered-pair count by `B` yields the exact additive identity
`B(B q^2 + 2 r q) + B r = (B q + r)^2 + r(B - r)`,
separating the smooth quadratic term `(B q + r)^2` from the strictly non-negative
discrete integrality excess `r(B - r) = B^2 θ(1 - θ)`.
-/
theorem erdos_1191.variants.legendre_fenchel_squared_identity (B q r : ℕ)
    (hr : r ≤ B) :
    B * (B * (q * q) + 2 * (r * q)) + B * r =
      (B * q + r) * (B * q + r) + r * (B - r) := by
  have h_r : r * r + r * (B - r) = r * B := by
    rw [← Nat.mul_add, Nat.add_sub_cancel' hr]
  have h_expand : (B * q + r) * (B * q + r) = B * (B * (q * q) + 2 * (r * q)) + r * r := by
    calc
      (B * q + r) * (B * q + r)
        = (B * q) * (B * q) + (B * q) * r + (r * (B * q) + r * r) := by
          rw [Nat.add_mul, Nat.mul_add, Nat.mul_add]
      _ = B * (B * (q * q)) + B * (r * q) + B * (r * q) + r * r := by
          ac_rfl
      _ = B * (B * (q * q) + 2 * (r * q)) + r * r := by
          rw [Nat.two_mul, Nat.mul_add, Nat.mul_add, ← Nat.add_assoc]
  rw [h_expand, Nat.add_assoc, h_r, Nat.mul_comm r B]

/--
**Theorem 2.6(c) (Karamata Integer Convexity / Robin-Hood Smoothing Step):**
Whenever two buckets have occupancies `x_1 ≥ x_2 + 2`, transferring one element from the
heavier bucket (`x_1 - 1`) to the lighter bucket (`x_2 + 1`) preserves the sum and strictly
reduces the ordered-pair count `x_1(x_1 - 1) + x_2(x_2 - 1)` by `2(x_1 - x_2 - 1) ≥ 2`.
-/
theorem erdos_1191.variants.robin_hood_pair_smoothing (x1 x2 : ℕ) (h_gap : x2 + 2 ≤ x1) :
    (x1 - 1) * (x1 - 2) + (x2 + 1) * x2 + 2 * (x1 - x2 - 1) =
      x1 * (x1 - 1) + x2 * (x2 - 1) := by
  have h_step1 : x1 * (x1 - 1) = (x1 - 1) * (x1 - 2) + 2 * (x1 - 1) := by
    cases x1 with
    | zero => omega
    | succ m =>
      cases m with
      | zero => omega
      | succ k =>
        change (k + 2) * (k + 1) = (k + 1) * k + 2 * (k + 1)
        rw [Nat.add_mul, Nat.mul_comm k (k + 1)]
  have h_step2 : (x2 + 1) * x2 = x2 * (x2 - 1) + 2 * x2 := by
    cases x2 with
    | zero => rfl
    | succ m =>
      change (m + 2) * (m + 1) = (m + 1) * m + 2 * (m + 1)
      rw [Nat.add_mul, Nat.mul_comm m (m + 1)]
  rw [h_step1, h_step2]
  omega

/--
**Arithmetic Lemma for Theorem 2.6(d) & Theorem 2.11(c) (Quarter-Square Identity & Bound):**
For any natural numbers `B, r ∈ ℕ` with `r ≤ B`, the scaled product `4 r(B - r)`
plus the single square `(2r - B)^2` (encoded in Lean via truncated subtraction on `ℕ`
as `(2r - B)^2 + (B - 2r)^2`) equals `B^2`, and hence `4 r(B - r) ≤ B^2`.
This provides the arithmetic bound `r(B - r) ≤ B^2 / 4` used to bound the relative
integrality excess `r(B - r)/S^2 ≤ 1/(4μ^2)` (`μ = S/B`) in [Pa26].
-/
theorem erdos_1191.variants.normalized_integrality_excess_bound (B r : ℕ) (hr : r ≤ B) :
    4 * (r * (B - r)) + (2 * r - B) * (2 * r - B) + (B - 2 * r) * (B - 2 * r) = B * B ∧
    4 * (r * (B - r)) ≤ B * B := by
  have h_split : 2 * r ≤ B ∨ B ≤ 2 * r := Nat.le_total (2 * r) B
  rcases h_split with h_lo | h_hi
  · have h_zero : 2 * r - B = 0 := Nat.sub_eq_zero_of_le h_lo
    rw [h_zero]
    have h_sum : (B - 2 * r) + 2 * r = B := Nat.sub_add_cancel h_lo
    have h_alg : 4 * (r * (B - r)) + (B - 2 * r) * (B - 2 * r) = B * B := by
      calc
        4 * (r * (B - r)) + (B - 2 * r) * (B - 2 * r)
          = (2 * r) * (2 * (B - r)) + (B - 2 * r) * (B - 2 * r) := by
            have h4 : 4 * (r * (B - r)) = (2 * r) * (2 * (B - r)) := by
              calc
                4 * (r * (B - r)) = (2 * 2) * (r * (B - r)) := rfl
                _ = (2 * r) * (2 * (B - r)) := by ac_rfl
            rw [h4]
        _ = (2 * r) * ((B - 2 * r) + B) + (B - 2 * r) * (B - 2 * r) := by
            have h2 : 2 * (B - r) = (B - 2 * r) + B := by omega
            rw [h2]
        _ = (2 * r) * (B - 2 * r) + (2 * r) * B + (B - 2 * r) * (B - 2 * r) := by
            rw [Nat.mul_add]
        _ = ((B - 2 * r) + 2 * r) * (B - 2 * r) + (2 * r) * B := by
            rw [Nat.add_mul]
            omega
        _ = B * (B - 2 * r) + B * (2 * r) := by
            rw [h_sum, Nat.mul_comm (2 * r) B]
        _ = B * ((B - 2 * r) + 2 * r) := by
            rw [← Nat.mul_add]
        _ = B * B := by
            rw [h_sum]
    omega
  · have h_zero : B - 2 * r = 0 := Nat.sub_eq_zero_of_le h_hi
    rw [h_zero]
    have h_2r : (2 * r - B) + B = 2 * r := Nat.sub_add_cancel h_hi
    have h_alg : 4 * (r * (B - r)) + (2 * r - B) * (2 * r - B) = B * B := by
      calc
        4 * (r * (B - r)) + (2 * r - B) * (2 * r - B)
          = (2 * r) * (2 * (B - r)) + (2 * r - B) * (2 * r - B) := by
            have h4 : 4 * (r * (B - r)) = (2 * r) * (2 * (B - r)) := by
              calc
                4 * (r * (B - r)) = (2 * 2) * (r * (B - r)) := rfl
                _ = (2 * r) * (2 * (B - r)) := by ac_rfl
            rw [h4]
        _ = ((2 * r - B) + B) * (2 * (B - r)) + (2 * r - B) * (2 * r - B) := by
            rw [h_2r]
        _ = (2 * r - B) * (2 * (B - r)) + B * (2 * (B - r)) + (2 * r - B) * (2 * r - B) := by
            rw [Nat.add_mul]
        _ = (2 * r - B) * (2 * (B - r) + (2 * r - B)) + B * (2 * (B - r)) := by
            rw [Nat.mul_add]
            omega
        _ = (2 * r - B) * B + B * (2 * (B - r)) := by
            have h_in : 2 * (B - r) + (2 * r - B) = B := by omega
            rw [h_in]
        _ = B * ((2 * r - B) + 2 * (B - r)) := by
            rw [Nat.mul_comm (2 * r - B) B, ← Nat.mul_add]
        _ = B * B := by
            have h_in : (2 * r - B) + 2 * (B - r) = B := by omega
            rw [h_in]
    omega

/-! ## 5. Theorem 2.7(a): Order-2 Cesàro–Fejér Triangular Kernel Recurrence -/

/--
**Arithmetic Lemma for Theorem 2.7(a) (Order-2 Cesàro–Fejér Hockey-Stick Recurrence):**
For every `m ∈ ℕ`, `m(m - 1) + 2 m = (m + 1) m`, giving the arithmetic recurrence step
for the unweighted Order-2 Cesàro–Fejér sum `sum_{w=1..W} max(w - d, 0) = C(W - d + 1, 2)`
(with `d = |x - y| < W`) in Theorem 2.7(a) of [Pa26].
-/
theorem erdos_1191.variants.cesaro_fejer_order2_step (m : ℕ) :
    m * (m - 1) + 2 * m = (m + 1) * m := by
  cases m with
  | zero => rfl
  | succ k =>
    change (k + 1) * k + 2 * (k + 1) = (k + 2) * (k + 1)
    rw [Nat.mul_comm (k + 1) k, ← Nat.add_mul]

/-! ## 6. Theorem 2.9(a)–(b): Scale-Invariant Block-Window Conservation -/

/--
**Theorem 2.9(a) (Scale-Invariant Block-Window Exponent Identity):**
For any shell index `k`, scale depth `s`, and octave offset `r` with `1 ≤ s + r ≤ k`,
the sliding-window width exponent `k - (s + r)` and sub-block count exponent `(s + r) - 1`
sum to `k - 1` independently of `s` and `r`.
-/
theorem erdos_1191.variants.scale_invariant_block_window_exponents (k s r : ℕ)
    (hsr_pos : 1 ≤ s + r) (hsr_le : s + r ≤ k) :
    (k - (s + r)) + ((s + r) - 1) = k - 1 := by
  omega

/--
**Theorem 2.9(b) (Scale-Invariant Block-Window Product Identity):**
For any shell index `k`, scale depth `s`, and octave offset `r` with `1 ≤ s + r ≤ k`,
the sliding-window width `w = 2^{k - s - r}` and sub-block count `B = 2^{s + r - 1}`
satisfy `w · B = 2^{k - 1} = |I_k|` identically across all scales `s` and octaves `r`.
-/
theorem erdos_1191.variants.scale_invariant_block_window_product (k s r : ℕ)
    (hsr_pos : 1 ≤ s + r) (hsr_le : s + r ≤ k) :
    2 ^ (k - (s + r)) * 2 ^ ((s + r) - 1) = 2 ^ (k - 1) := by
  rw [← Nat.pow_add, erdos_1191.variants.scale_invariant_block_window_exponents k s r hsr_pos hsr_le]

/-! ## 7. Theorem 2.11(b): Hyperfactorial Sawtooth Abel Step -/

/--
**Theorem 2.11(b) (Hyperfactorial Sawtooth Abel Coefficient Identity):**
For every integer period `m ≥ 1`, `m(m + 1) - (m - 1)m = 2m`.
This exact identity telescopes the piecewise logarithmic integrals of `{u}(1 - {u})/u`
into the logarithm of the hyperfactorial `H(M - 1) = prod_{m=1..M-1} m^m`, yielding the
Glaisher–Kinkelin sawtooth constant `C_GK = 2 ln A_GK = 1/6 - 2 ζ'(-1)`.
-/
theorem erdos_1191.variants.hyperfactorial_sawtooth_abel_step (m : ℕ) (hm : 1 ≤ m) :
    m * (m + 1) - (m - 1) * m = 2 * m := by
  cases m with
  | zero => omega
  | succ k =>
    simp [Nat.add_mul, Nat.mul_add]
    omega

/-! ## 8. Section 4: Loop-Free O(1) 64-Bit Bit-Shift Sidon Legality & Certified Witnesses -/

/--
**Section 4 (Bit-Reversed Shift Difference Identity):**
In a 64-bit register holding `rev_mask = sum_{a ∈ A} 2^{63 - a}`, right-shifting by `63 - x`
places the bit for `a < x ≤ 63` at bit position `(63 - a) - (63 - x) = x - a`.
-/
theorem erdos_1191.variants.bitshift_sidon_difference_index (a x : ℕ) (h_lt : a < x) (h_hi : x ≤ 63) :
    (63 - a) - (63 - x) = x - a := by
  omega

/-- Boolean Sidon legality checker over a finite list of strictly increasing naturals. -/
def allPositiveDiffs : List ℕ → List ℕ
  | [] => []
  | a :: rest => (rest.map (fun b => b - a)) ++ allPositiveDiffs rest

def isStrictlyIncreasing : List ℕ → Bool
  | [] => true
  | [_] => true
  | a :: b :: rest => (a < b) && isStrictlyIncreasing (b :: rest)

def hasNoDuplicates : List ℕ → Bool
  | [] => true
  | d :: rest => !(rest.contains d) && hasNoDuplicates rest

def isCertifiedSidon (A : List ℕ) : Bool :=
  isStrictlyIncreasing A && hasNoDuplicates (allPositiveDiffs A)

/--
**Table 1 (Mian–Chowla Greedy Prefix Verification):**
Machine-checked proof that the first 10 terms of the Mian–Chowla sequence (OEIS `A005282`)
`[1, 2, 4, 8, 13, 21, 31, 45, 66, 81]` form a valid Sidon set with `45` distinct differences.
-/
def mianChowlaPrefix10 : List ℕ := [1, 2, 4, 8, 13, 21, 31, 45, 66, 81]

theorem mianChowlaPrefix10_is_sidon :
    isCertifiedSidon mianChowlaPrefix10 = true ∧
    (allPositiveDiffs mianChowlaPrefix10).length = 45 := by
  decide

/--
**Table 2 (Exact Maximin Lower-Envelope Witness Verification at `N = 32`):**
Machine-checked proof that the Maximin optimal witness `A^* = [1, 2, 5, 7, 14, 22, 32]`
of cardinality `|A^*| = 7` is a valid Sidon set with all `21` positive differences distinct.
-/
def maximinWitnessN32 : List ℕ := [1, 2, 5, 7, 14, 22, 32]

theorem maximinWitnessN32_is_sidon :
    isCertifiedSidon maximinWitnessN32 = true ∧
    (allPositiveDiffs maximinWitnessN32).length = 21 := by
  decide

end Erdos1191
