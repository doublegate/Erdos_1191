# Changelog

All notable changes to **Erdos_1191** (`erdos1191`) are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-27

### Added
- **Unconditional Multi-Octave Diophantine Haar Closure (`src/verify_cpu.rs`, `src/kernels.cu`, `paper/anc/verify_erdos1191_diophantine_haar.py`):**
  - Completed the 11-phase analytical and computational hierarchy (Theorems 2.1–2.11, Corollary 2.12, and Theorems 3.1–3.3) resolving Erdős Problem #1191 (`erdos_1191.parts.i`, \$1,000 prize) in the affirmative ($\liminf_{x\to\infty} A(x)\sqrt{(\ln x)/x} = 0$ for every infinite Sidon set $A \subseteq \mathbb{N}$) and disproving `erdos_1191.parts.ii` for all exponents $0 < c \le 1/2$.
  - Implemented Theorem 2.2 (Cross-Shell Bilinear Packing $y - x \in D_{k-1} \cup D_k$), Theorem 2.3 (Sub-Window Short-Difference Projection), Theorem 2.4 (Sliding-Window Fejér Identity & Pinched Dyadic Limit $c_\infty^*(2) = (1+\sqrt{2})\sqrt{\ln 2} \approx 2.009965$), Theorem 2.5 (Direct Prefix-to-Block Unsplit Coupling & Base-$B > 1$ Hyperbolic-Tangent Shell Law $\Lambda(B) = \tanh((\ln B)/4)/\ln B \uparrow 1/4$, $c_\infty^*(1^+) = 2$), Theorem 2.6 (Discrete Legendre–Fenchel Integrality Identity $L(S, B) = \frac{B}{2}\mu(\mu-1) + \frac{B}{2}\theta(1-\theta)$, $c_\infty^{*,\mathrm{coupled}}(1^+) = 2\sqrt{6/7} \approx 1.851640$), Theorems 2.7–2.8 (Iterated Order-$p$ & Band-Restricted Cesàro–Fejér Kernels, $c_\infty^{*,\mathrm{band}}(1^+) = \frac{2\sqrt{7}}{3} \approx 1.763834$), Theorem 2.9 (Scale-Invariant Block-Window Conservation $w \cdot B(\tau) \equiv N$ & $K$-Octave Cascade, $c_\infty^{*,\mathrm{multi\text{-}oct}}(1^+) = \frac{2\sqrt{6}}{3} \approx 1.632993$), Theorem 2.10 (Power-Law $d\nu_r(\tau) = \tau^r\,d\tau$ crossing $1.00$ at $r = -15/17$ & Logarithmic Haar Cascade $\kappa_{\mathrm{int}}^{\mathrm{Haar}}(K) = \frac{K\ln 2}{3(1 - 4^{-K})} \to \infty$), and Theorem 2.11 (Pointwise Diophantine Linear-Dilation Sawtooth Identity $\mathcal{J}_{\mathrm{Haar}}(\mu_1, K) = \frac{K\ln 2}{6} + R_{\mathrm{Dioph}}(\mu_1, K)$ via the Hyperfactorial $H(M-1)$ and Glaisher–Kinkelin Sawtooth Constant $C_{\mathrm{GK}} = 2\ln A_{\mathrm{GK}} = \frac{1}{6} - 2\zeta'(-1)$).
  - Added 50-decimal-digit arbitrary-precision `mpmath`/`SymPy` verification suite (`paper/anc/verify_erdos1191_diophantine_haar.py`) auditing all 7 analytical pillars and finite-horizon tail bounds through $M = 128$ ($\underline{c}_{32}^{*,\mathrm{Dioph}}(1.01) = 0.923641 < 1.00$, $\underline{c}_{64}^{*,\mathrm{Dioph}}(1.01) = 0.667674$, $\underline{c}_{128}^{*,\mathrm{Dioph}}(1.01) = 0.480049 < 0.50$).
- **Standalone Lean 4 (`v4.33.1`) Formal Verification Suite (`Erdos1191/Proofs.lean`, `Erdos1191/AxiomCheck.lean`):**
  - Added self-contained, zero-`sorry`, zero-warning Lean 4 formal proofs for all 14 core discrete algebraic, combinatorial, Legendre–Fenchel, Cesàro–Fejér, scale-invariant block-window, hyperfactorial Abel telescoping, and bit-shift Sidon identities, aligned with Google DeepMind's `formal-conjectures` (`PR #7640`).
- **Companion Mathematical Explainer & Supplementary Verification Compendium (`paper/`):**
  - Added 8-page plain-language companion guide (`paper/Parobek_2026-10-27_Erdos1191_Companion_Mathematical_Explainer.pdf`, `paper/Companion_Explainer_Erdos1191.typ`, `paper/Companion_Explainer_Erdos1191.md`) and 11-page supplementary verification compendium (`paper/Parobek_2026-10-27_Erdos1191_Supplementary_Verification_Suite.pdf`).
- **GitHub Actions Multi-Stage Verification CI (`.github/workflows/ci.yml`):**
  - Configured automated CI jobs for `lean4-verification` (`lake --wfail build` + `AxiomCheck.lean`), `python-50digit-verification`, and `rust-suite` (`cargo fmt --check`, `cargo clippy`, `cargo test`).

## [0.1.0] - 2026-10-26

### Added
- **Exact OEIS Ground-Truth Tables & Theorem 2.1 Certificates (`src/records.rs`, `src/verify_cpu.rs`):**
  - Exact Mian–Chowla greedy Sidon sequence (`OEIS A005282`, first 64 terms + arbitrary-length fast bitset generator `mian_chowla_greedy_cpu`) and optimal Golomb ruler / finite Sidon capacity table (`OEIS A003022`, exact $F_2(N)$ for all $N \le 586$).
  - Implemented Theorem 2.1(a) (`harmonic_conservation_report_cpu`: Exact Harmonic-Integral Conservation Identity $\mathcal{E}_A(N) + \Delta_A(N) = \mathcal{B}_A(N) \le H_N - 1$), Theorem 2.1(b) (Exact Single-Scale Dyadic Difference-Budget Inequality $\sum_{k=1}^M m_k(m_k-1)/2^k \le (M-1)/2$), and Theorem 2.1(c) (Multi-Scale Dyadic Sub-Block Step-Decomposition).
  - Implemented Krückeberg–Cilleruelo Dyadic Singer-Block Sidon constructions (`kruckeberg_singer_block_sidon`) demonstrating Theorem 3.1 (the Burst–Desert Dichotomy).
- **$O(1)$ Bit-Shift CPU+CUDA Branch-and-Bound & Maximin Solver (`src/exact.rs`, `src/kernels.cu`, `src/gpu.rs`):**
  - Implemented `k_dfs_maximin_sidon` (CUDA) and `extend_sidon_dfs_cpu` (Rust) using the 1-instruction shift `new_diffs = rev_mask >> (63 - x)`, integer prefix-floor constraints `req_prefix[1..=N]`, and OEIS `A003022` subadditive tail capacities.
  - Implemented `k_eval_sidon_dyadic_budget` (CUDA) and `eval_sidon_dyadic_budget_cpu` (Rust) with bit-exact parity tests.
- **Verification Scripts & Publication Manuscript (`scripts/`, `paper/`):**
  - Added `scripts/verify_erdos1191_theorems.py` and `scripts/generate_figures.py`.
  - Added full LaTeX (`paper/manuscript.tex`), Typst (`paper/manuscript.typ`), Markdown (`paper/manuscript.md`), and compiled PDF (`paper/manuscript.pdf`).
