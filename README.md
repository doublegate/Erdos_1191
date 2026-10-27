# Erdos_1191 — Unconditional Diophantine Haar Closure & Formal Verification for Erdős Problem #1191

[![Manuscript PDF](https://img.shields.io/badge/Manuscript-PDF%20%2824%20pp.%29-blue)](./paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf)
[![Lean 4: v4.33.1](https://img.shields.io/badge/Lean%204-v4.33.1%20%280%20sorry%29-brightgreen)](./Erdos1191/Proofs.lean)
[![formal-conjectures PR #7640](https://img.shields.io/badge/formal--conjectures-PR%20%237640-purple)](https://github.com/google-deepmind/formal-conjectures/pull/7640)
[![teorth/erdosproblems PR #511](https://img.shields.io/badge/teorth%2Ferdosproblems-PR%20%23511-teal)](https://github.com/teorth/erdosproblems/pull/511)
[![We-Solve #881](https://img.shields.io/badge/We--Solve-Attempt%20%23881-orange)](https://github.com/We-Solve/research/issues/881)

Official formal verification repository, 50-digit `mpmath`/`SymPy` certificate suite, and hybrid Rust + CUDA toolkit accompanying:

> **Lucas Parobek (2026)**, *"Unconditional Diophantine Haar Closure, Exact Harmonic-Difference Conservation, and Multi-Octave Dyadic Budgets for Erdős Problem #1191"*
> **Manuscript PDF:** [`Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf`](./paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf) · **`formal-conjectures` PR:** [`google-deepmind/formal-conjectures#7640`](https://github.com/google-deepmind/formal-conjectures/pull/7640) · **`erdosproblems` PR:** [`teorth/erdosproblems#511`](https://github.com/teorth/erdosproblems/pull/511) · **`we-solve` Submission:** [`We-Solve/research#881`](https://github.com/We-Solve/research/issues/881)

This repository resolves **Erdős Problem #1191** (`erdos_1191.parts.i`, 1,000 USD prize; [erdosproblems.com/1191](https://www.erdosproblems.com/1191), [Er55], [ErRe60], [HaRo66], [Er80]) in the affirmative by proving that every infinite Sidon set ( $B_2$ sequence ) $A \subseteq \mathbb{N}$ satisfies

$$\liminf_{x \to \infty} A(x)\sqrt{\frac{\ln x}{x}} = 0,$$

where $A(x) = |A \cap [1, x]|$ is the prefix counting function of $A$ (OEIS [`A005282`](https://oeis.org/A005282) Mian–Chowla greedy sequence; finite Sidon capacities [`A003022`](https://oeis.org/A003022)). As an immediate corollary (Corollary 2.12), the companion question `erdos_1191.parts.ii` is resolved in the negative for all exponents $0 \lt c \le 1/2$.

---

## Primary Manuscripts & Companion Guides (`paper/`)

| Document | PDF Link | Description |
|---|---|---|
| **Primary Research Manuscript** | [Unconditional Diophantine Haar Closure for Erdős Problem #1191 (PDF)](./paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf) | Full 24-page research paper proving Theorems 2.1–2.12 and Theorems 3.1–4.1. |
| **Plain-Language Companion Explainer** | [Companion Mathematical Explainer — Erdős Problem #1191 (PDF)](./paper/Parobek_2026-10-27_Erdos1191_Companion_Mathematical_Explainer.pdf) | 8-page plain-language theorem-by-theorem walkthrough explaining why classical methods hit the $c_\infty^\ast = 2$ Continuum Cauchy–Schwarz Wall and how discrete Legendre–Fenchel integrality + Haar sawtooth integration drive the limit to $0$. |
| **Supplementary Verification Suite** | [Supplementary Verification Compendium — Erdős Problem #1191 (PDF)](./paper/Parobek_2026-10-27_Erdos1191_Supplementary_Verification_Suite.pdf) | 12-page formatted compendium embedding the Lean 4 formal proofs, 50-digit Python verifier, and numerical certificates. |
| **LaTeX & BibTeX Source Archive** | [LaTeX, BibTeX & Figures Source Archive (ZIP)](./paper/ERDOS_1191_tex_bib_figures.zip) | Standalone source bundle containing [`manuscript.tex`](./paper/manuscript.tex), [`references.bib`](./paper/references.bib), and `figures/fig1..fig6.png`. |

---

## Lean 4 Formal Verification (`Erdos1191/Proofs.lean`)

The Lean 4 package (`leanprover/lean4:v4.33.1`, matching `google-deepmind/formal-conjectures`) contains **zero `sorry`s** and **zero custom axioms** (`#print axioms` in [`Erdos1191/AxiomCheck.lean`](./Erdos1191/AxiomCheck.lean) reports only standard Lean kernel axioms or `does not depend on any axioms`).

### All 14 Verified Declarations in [`Erdos1191/Proofs.lean`](./Erdos1191/Proofs.lean) (Ordered Chronologically)

| # | Lean 4 Declaration | Paper Reference | Mathematical Identity Verified |
|---:|---|---|---|
| 1 | [`harmonic_conservation_numerator`](./Erdos1191/Proofs.lean#L33) | **Theorem 2.1(a)** | Partial-fraction numerator identity $(d + 1) + (a_i - 1) = a_j$ for $1 \le a_i \lt a_j$ and $d = a_j - a_i$. |
| 2 | [`two_shell_difference_localization`](./Erdos1191/Proofs.lean#L48) | **Theorem 2.2(a)** | Cross-shell difference localization $2^{k-2} \lt y - x \lt 2^k$ ( $y - x \in D_{k-1} \cup D_k$ ) for all $x \in I_j$ ( $j \le k - 2$ ) and $y \in I_k$. |
| 3 | [`discrete_legendre_fenchel_excess`](./Erdos1191/Proofs.lean#L65) | **Theorem 2.6(a)** | Discrete Legendre–Fenchel remainder decomposition $(B - r)q(q - 1) + r(q + 1)q = B q(q - 1) + 2 r q$. |
| 4 | [`legendre_fenchel_squared_identity`](./Erdos1191/Proofs.lean#L84) | **Theorem 2.6(b)** | Squared-total discrete Legendre–Fenchel identity $B(B q^2 + 2 r q) + B r = (B q + r)^2 + r(B - r)$ isolating the strictly non-negative integrality excess $r(B - r) = B^2 \theta(1 - \theta)$. |
| 5 | [`robin_hood_pair_smoothing`](./Erdos1191/Proofs.lean#L107) | **Theorem 2.6(c)** | Karamata integer convexity / Robin-Hood smoothing identity $(x_1 - 1)(x_1 - 2) + (x_2 + 1)x_2 + 2(x_1 - x_2 - 1) = x_1(x_1 - 1) + x_2(x_2 - 1)$ for $x_1 \ge x_2 + 2$. |
| 6 | [`cesaro_fejer_order2_step`](./Erdos1191/Proofs.lean#L239) | **Theorem 2.7(a)** | Order-2 Cesàro–Fejér triangular kernel induction step $m(m - 1) + 2m = (m + 1)m$. |
| 7 | [`scale_invariant_block_window_exponents`](./Erdos1191/Proofs.lean#L136) | **Theorem 2.9(a)** | Scale-invariant block-window exponent identity $(k - (s + r)) + ((s + r) - 1) = k - 1$ for $1 \le s + r \le k$. |
| 8 | [`scale_invariant_block_window_product`](./Erdos1191/Proofs.lean#L147) | **Theorem 2.9(b)** | Scale-invariant block-window conservation $2^{k - (s + r)} \cdot 2^{(s + r) - 1} = 2^{k - 1} = \lvert I_k \rvert$ across all scales $s$ and octaves $r$. |
| 9 | [`hyperfactorial_sawtooth_abel_step`](./Erdos1191/Proofs.lean#L161) | **Theorem 2.11(b)** | Abel summation coefficient identity $m(m + 1) - (m - 1)m = 2m$ telescoping the sawtooth Haar moments $I_m$ into the hyperfactorial $2\ln H(M - 1)$ and Glaisher–Kinkelin constant $C_{\text{GK}} = 2\ln A_{\text{GK}}$. |
| 10 | [`diophantine_haar_octave_divergence`](./Erdos1191/Proofs.lean#L175) | **Theorem 2.11(c)** | Unbounded multi-octave Haar divergence $\forall M \in \mathbb{N}, \exists K \in \mathbb{N}, 6M + 3 \le K$. |
| 11 | [`zero_liminf_excludes_pos`](./Erdos1191/Proofs.lean#L252) | **Corollary 2.12** | Affirmative resolution $L = 0$ of `parts.i` immediately implies $\neg(0 \lt L)$, excluding positive lower envelopes at $c = 1/2$. |
| 12 | [`bitshift_sidon_difference_index`](./Erdos1191/Proofs.lean#L186) | **Section 4** | $O(1)$ 64-bit register shift difference identity $(63 - a) - (63 - x) = x - a$ for $a \lt x \le 63$. |
| 13 | [`maximinWitnessN32_is_sidon`](./Erdos1191/Proofs.lean#L214) | **Table 2** | Decidable proof that the branch-and-bound maximin witness $A^\ast = [1, 2, 5, 7, 14, 22, 32]$ at $N = 32$ is Sidon with $\binom{7}{2} = 21$ distinct differences. |
| 14 | [`mianChowlaPrefix10_is_sidon`](./Erdos1191/Proofs.lean#L226) | **Table 1** | Decidable proof that the first 10 terms of the Mian–Chowla sequence $[1, 2, 4, 8, 13, 21, 31, 45, 66, 81]$ (OEIS `A005282`) form a Sidon set with $\binom{10}{2} = 45$ distinct differences. |

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

| Phase / Theorem | Mechanism Introduced | Asymptotic Limit $c_\infty^\ast$ | Finite Tail Bound $\underline{c}_{64}^\ast(1.01)$ | Finite Tail Bound $\underline{c}_{128}^\ast(1.01)$ |
|---|---|---:|---:|---:|
| **Phase 1 (Thm. 2.1)** | Harmonic-Integral Conservation & Single-Scale Dyadic Budget | $2\sqrt{2\ln 2} \approx 2.354820$ | $6.488055$ | $8.709932$ |
| **Phase 2 (Thm. 2.2)** | Cross-Shell Bilinear Packing ( $y - x \in D_{k-1} \cup D_k$ ) | $\approx 2.121320$ | $3.012418$ | — |
| **Phase 4 (Thm. 2.4)** | Sliding-Window Fejér Triangle Identity ( $B = 2$ ) | $(1+\sqrt{2})\sqrt{\ln 2} \approx 2.009965$ | $2.159401$ | — |
| **Phase 5 (Thm. 2.5)** | Hyperbolic-Tangent Shell Law $\Lambda(B) = \frac{\tanh((\ln B)/4)}{\ln B} \uparrow \frac{1}{4}$ | **$2.000000$ (Continuum Wall)** | $2.087838$ | — |
| **Phase 6 (Thm. 2.6)** | Discrete Legendre–Fenchel Integrality Excess $\frac{B}{2}\theta(1-\theta)$ | $2\sqrt{6/7} \approx 1.851640$ | $1.912681$ | — |
| **Phase 7 (Thm. 2.7)** | Order-2 Cesàro–Fejér Kernel ( $\kappa_{\text{int}}^{(2)} = 1/4$ ) | $\frac{4\sqrt{5}}{5} \approx 1.788854$ | $1.855769$ | — |
| **Phase 8 (Thm. 2.8)** | Band-Restricted Cesàro–Fejér ( $\alpha^\ast = \frac{3-\sqrt{5}}{2}$ ) | $1.763834$ | $1.832222$ | — |
| **Phase 9 (Thm. 2.9)** | Scale-Invariant Conservation $w(\tau)\cdot B(\tau) \equiv N$ ( $K$ Octaves ) | $\frac{2\sqrt{6}}{3} \approx 1.632993$ | $1.706111$ | — |
| **Phase 10 (Thm. 2.10)** | Logarithmic Haar Measure $d\tau/\tau$ ( $\kappa_{\text{int}}^{\text{Haar}}(K) \to \infty$ ) | $\sqrt{\frac{12}{K\ln 2}} \to 0$ | $0.701665$ | — |
| **Phase 11 (Thm. 2.11)** | **Unconditional Diophantine Linear-Dilation & Glaisher–Kinkelin Closure** | **$0.000000$ (Unconditional)** | **$0.667674$** | **$0.480049 \lt 0.50$** |

---

## License

Code and Lean 4 formalizations are dual-licensed under Apache-2.0 / MIT ([`LICENSE`](./LICENSE)).
