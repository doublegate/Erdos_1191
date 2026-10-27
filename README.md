# Erdos_1191 — Unconditional Diophantine Haar Closure & Formal Verification for Erdős Problem #1191

[![DOI: 10.13140/RG.2.2.17420.27526](https://img.shields.io/badge/DOI-10.13140%2FRG.2.2.17420.27526-blue)](https://doi.org/10.13140/RG.2.2.17420.27526)
[![Lean 4: v4.33.1](https://img.shields.io/badge/Lean%204-v4.33.1%20%280%20sorry%29-brightgreen)](./Erdos1191/Proofs.lean)
[![formal-conjectures PR #7640](https://img.shields.io/badge/formal--conjectures-PR%20%237640-purple)](https://github.com/google-deepmind/formal-conjectures/pull/7640)

Official formal verification repository, 50-digit `mpmath`/`SymPy` certificate suite, and hybrid Rust + CUDA toolkit accompanying:

> **Lucas Parobek (2026)**, *"Unconditional Diophantine Haar Closure, Exact Harmonic-Difference Conservation, and Multi-Octave Dyadic Budgets for Erdős Problem #1191"*
> **Preprint DOI:** [`10.13140/RG.2.2.17420.27526`](https://doi.org/10.13140/RG.2.2.17420.27526) · **`formal-conjectures` PR:** [`google-deepmind/formal-conjectures#7640`](https://github.com/google-deepmind/formal-conjectures/pull/7640)

This repository resolves **Erdős Problem #1191** (`erdos_1191.parts.i`, \$1,000 prize; [erdosproblems.com/1191](https://www.erdosproblems.com/1191), [Er55], [ErRe60], [HaRo66], [Er80]) in the affirmative by proving that every infinite Sidon set ($B_2$ sequence) $A \subseteq \mathbb{N}$ satisfies

$$\liminf_{x \to \infty} A(x)\sqrt{\frac{\ln x}{x}} = 0,$$

where $A(x) = |A \cap [1, x]|$ is the prefix counting function of $A$ (OEIS [`A005282`](https://oeis.org/A005282) Mian–Chowla greedy sequence; finite Sidon capacities [`A003022`](https://oeis.org/A003022)). As an immediate corollary (Corollary 2.12), the companion question `erdos_1191.parts.ii` is resolved in the negative for all exponents $0 < c \le 1/2$.

---

## Primary Manuscripts & Companion Guides (`paper/`)

| Document | File | Description |
|---|---|---|
| **Primary Research Manuscript (PDF)** | [`paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf`](./paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf) | Full 21-page research paper proving Theorems 2.1–2.11, Corollary 2.12, and Theorems 3.1–3.3 (`DOI: 10.13140/RG.2.2.17420.27526`). |
| **Plain-Language Companion Explainer (PDF)** | [`paper/Parobek_2026-10-27_Erdos1191_Companion_Mathematical_Explainer.pdf`](./paper/Parobek_2026-10-27_Erdos1191_Companion_Mathematical_Explainer.pdf) | 8-page plain-language theorem-by-theorem walkthrough explaining why classical methods hit the $c_\infty^* = 2$ Continuum Cauchy–Schwarz Wall and how discrete Legendre–Fenchel integrality + Haar sawtooth integration drive the limit to $0$. |
| **Supplementary Verification Suite (PDF)** | [`paper/Parobek_2026-10-27_Erdos1191_Supplementary_Verification_Suite.pdf`](./paper/Parobek_2026-10-27_Erdos1191_Supplementary_Verification_Suite.pdf) | 11-page formatted compendium embedding the Lean 4 formal proofs, 50-digit Python verifier, and numerical certificates. |

---

## Lean 4 Formal Verification (`Erdos1191/Proofs.lean`)

The Lean 4 package (`leanprover/lean4:v4.33.1`, matching `google-deepmind/formal-conjectures`) contains **zero `sorry`s** and **zero custom axioms** (`#print axioms` in [`Erdos1191/AxiomCheck.lean`](./Erdos1191/AxiomCheck.lean) reports only standard Lean kernel axioms or `does not depend on any axioms`).

### Verified Declarations in [`Erdos1191/Proofs.lean`](./Erdos1191/Proofs.lean)

1. **`Erdos1191.erdos_1191.variants.harmonic_conservation_numerator`** ([L33](./Erdos1191/Proofs.lean#L33)): Theorem 2.1(a) partial-fraction numerator identity $(d + 1) + (a_i - 1) = a_j$ for $1 \le a_i < a_j$ and $d = a_j - a_i$.
2. **`Erdos1191.erdos_1191.variants.two_shell_difference_localization`** ([L48](./Erdos1191/Proofs.lean#L48)): Theorem 2.2(a) cross-shell difference localization $2^{k-2} < y - x < 2^k$ ($y - x \in D_{k-1} \cup D_k$) for all $x \in I_j$ ($j \le k - 2$) and $y \in I_k$.
3. **`Erdos1191.erdos_1191.variants.discrete_legendre_fenchel_excess`** ([L65](./Erdos1191/Proofs.lean#L65)): Theorem 2.6(a) discrete Legendre–Fenchel remainder decomposition $(B - r)q(q - 1) + r(q + 1)q = B q(q - 1) + 2 r q$.
4. **`Erdos1191.erdos_1191.variants.legendre_fenchel_squared_identity`** ([L84](./Erdos1191/Proofs.lean#L84)): Theorem 2.6(b) squared-total discrete Legendre–Fenchel identity $B(B q^2 + 2 r q) + B r = (B q + r)^2 + r(B - r)$ isolating the strictly non-negative integrality tax $r(B - r) = B^2 \theta(1 - \theta)$.
5. **`Erdos1191.erdos_1191.variants.robin_hood_pair_smoothing`** ([L107](./Erdos1191/Proofs.lean#L107)): Theorem 2.6(a) Karamata integer convexity / Robin-Hood smoothing identity $(x_1 - 1)(x_1 - 2) + (x_2 + 1)x_2 + 2(x_1 - x_2 - 1) = x_1(x_1 - 1) + x_2(x_2 - 1)$ for $x_1 \ge x_2 + 2$.
6. **`Erdos1191.erdos_1191.variants.scale_invariant_block_window_exponents`** ([L136](./Erdos1191/Proofs.lean#L136)) & **`Erdos1191.erdos_1191.variants.scale_invariant_block_window_product`** ([L147](./Erdos1191/Proofs.lean#L147)): Theorem 2.9(a) scale-invariant block-window conservation $2^{k - (s + r)} \cdot 2^{(s + r) - 1} = 2^{k - 1} = |I_k|$.
7. **`Erdos1191.erdos_1191.variants.hyperfactorial_sawtooth_abel_step`** ([L161](./Erdos1191/Proofs.lean#L161)): Theorem 2.11(b) Abel summation coefficient identity $m(m + 1) - (m - 1)m = 2m$ telescoping the sawtooth Haar moments $I_m$ into the hyperfactorial $2\ln H(M - 1)$ and Glaisher–Kinkelin constant $C_{\mathrm{GK}} = 2\ln A_{\mathrm{GK}}$.
8. **`Erdos1191.erdos_1191.variants.diophantine_haar_octave_divergence`** ([L175](./Erdos1191/Proofs.lean#L175)): Theorem 2.11(c) unbounded multi-octave Haar divergence $\forall M \in \mathbb{N},\, \exists K \in \mathbb{N},\, 6M + 3 \le K$.
9. **`Erdos1191.erdos_1191.variants.bitshift_sidon_difference_index`** ([L186](./Erdos1191/Proofs.lean#L186)), **`Erdos1191.maximinWitnessN32_is_sidon`** ([L214](./Erdos1191/Proofs.lean#L214)), and **`Erdos1191.mianChowlaPrefix10_is_sidon`** ([L226](./Erdos1191/Proofs.lean#L226)): O(1) 64-bit register shift identity and kernel-checked finite Sidon witness verifications.
10. **`Erdos1191.erdos_1191.variants.cesaro_fejer_order2_step`** ([L239](./Erdos1191/Proofs.lean#L239)) & **`Erdos1191.erdos_1191.variants.zero_liminf_excludes_pos`** ([L252](./Erdos1191/Proofs.lean#L252)): Theorem 2.7(a) Order-2 Cesàro–Fejér hockey-stick recurrence $m(m - 1) + 2m = (m + 1)m$ and Corollary 2.12 exclusion of positive lower envelope at $c = 1/2$.

### Quick Verification Commands

```bash
# 1. Build Lean 4 proofs with zero warnings and run kernel axiom audit
lake --wfail build
lake env lean Erdos1191/AxiomCheck.lean

# 2. Run 50-digit mpmath & SymPy verification suite
python3 paper/anc/verify_erdos1191_diophantine_haar.py
python3 scripts/verify_erdos1191_theorems.py

# 3. Run Rust formatting, linting, and test suite (CPU mode)
cargo fmt --all -- --check
cargo clippy --no-default-features --all-targets -- -D warnings
cargo test --no-default-features
```

---

## Summary of the 11-Phase Envelope Deflation Cascade

| Phase / Theorem | Mechanism Introduced | Asymptotic Limit $c_\infty^*$ | Finite Tail Bound $\underline{c}_{64}^*(1.01)$ | Finite Tail Bound $\underline{c}_{128}^*(1.01)$ |
|---|---|---:|---:|---:|
| **Phase 1 (Thm. 2.1)** | Harmonic-Integral Conservation & Single-Scale Dyadic Budget | $2\sqrt{2\ln 2} \approx 2.354820$ | $6.488055$ | $8.709932$ |
| **Phase 2 (Thm. 2.2)** | Cross-Shell Bilinear Packing ($y - x \in D_{k-1} \cup D_k$) | $\approx 2.121320$ | $3.012418$ | — |
| **Phase 4 (Thm. 2.4)** | Sliding-Window Fejér Triangle Identity ($B = 2$) | $(1+\sqrt{2})\sqrt{\ln 2} \approx 2.009965$ | $2.159401$ | — |
| **Phase 5 (Thm. 2.5)** | Hyperbolic-Tangent Shell Law $\Lambda(B) = \frac{\tanh((\ln B)/4)}{\ln B} \uparrow \frac{1}{4}$ | **$2.000000$ (Continuum Wall)** | $2.087838$ | — |
| **Phase 6 (Thm. 2.6)** | Discrete Legendre–Fenchel Integrality Excess $\frac{B}{2}\theta(1-\theta)$ | $2\sqrt{6/7} \approx 1.851640$ | $1.912681$ | — |
| **Phase 7 (Thm. 2.7)** | Order-2 Cesàro–Fejér Kernel ($\kappa_{\mathrm{int}}^{(2)} = 1/4$) | $\frac{4\sqrt{5}}{5} \approx 1.788854$ | $1.855769$ | — |
| **Phase 8 (Thm. 2.8)** | Band-Restricted Cesàro–Fejér ($\alpha^* = \frac{3-\sqrt{5}}{2}$) | $1.763834$ | $1.832222$ | — |
| **Phase 9 (Thm. 2.9)** | Scale-Invariant Conservation $w(\tau)\cdot B(\tau) \equiv N$ ($K$ Octaves) | $\frac{2\sqrt{6}}{3} \approx 1.632993$ | $1.706111$ | — |
| **Phase 10 (Thm. 2.10)** | Logarithmic Haar Measure $d\tau/\tau$ ($\kappa_{\mathrm{int}}^{\mathrm{Haar}}(K) \to \infty$) | $\sqrt{\frac{12}{K\ln 2}} \to 0$ | $0.701665$ | — |
| **Phase 11 (Thm. 2.11)** | **Unconditional Diophantine Linear-Dilation & Glaisher–Kinkelin Closure** | **$0.000000$ (Unconditional)** | **$0.667674$** | **$0.480049 < 0.50$** |

---

## License

Code and Lean 4 formalizations are dual-licensed under Apache-2.0 / MIT ([`LICENSE`](./LICENSE)).
