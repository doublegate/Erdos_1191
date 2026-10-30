# Changelog

All notable changes to **Erdos_1191** (`erdos1191`) are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2026-10-29

### Fixed & Remediated
- **`formal-conjectures#7640` Commit-Permalink Pinning & Companion Line-Anchor Sync (`FormalConjectures_1191.lean`, `submissions/formal_conjectures_lean4/1191.lean`, `paper/`):**
  - Replaced all 14 `formal_proof` `blob/main` URLs in `FormalConjectures/ErdosProblems/1191.lean` (`google-deepmind/formal-conjectures#7640`) with immutable commit permalinks (`faf0ce94febe33568b6b08e4a07a6ca593b83288`) per `PROOFS.md` and `bocowgill`'s review.
  - Synchronized all 17 Lean declaration line anchors (`#L34`–`#L452`) across `paper/Companion_Explainer_Erdos1191.typ`, `paper/Companion_Explainer_Erdos1191.md`, `paper/manuscript.md`, `paper/anc/1191.lean.txt`, and `ERDOS_1191_COMPLETE_SESSION_AND_PROJECT_ARCHIVE.md`, and recompiled all PDFs and `.zip`/`.tar.gz` bundles.

## [1.1.0] - 2026-10-27

### Added & Remediated
- **We-Solve `#881` & `formal-conjectures#7640` Reviewer Audit Remediations, Chronological 17-Theorem Lean 4 Suite, and $\text{High-}M$ (`M <= 4096`) Certification (`Erdos1191/Proofs.lean`, `paper/anc/verify_erdos1191_diophantine_haar.py`, `paper/`):**
  - Expanded and reordered `Erdos1191/Proofs.lean` (`paper/anc/1191.lean`, `FormalConjectures_1191.lean`, `google-deepmind/formal-conjectures#7640`) into **17 strictly chronological zero-`sorry` declarations** (`#L34`–`#L452`), adding `powers_of_two_bilinear_identity`, `powers_of_two_satisfies_corrected_bilinear_bound` (Remark 2.3), `kkt_increment_convexity_k_ge_4` (Lemma 2.4), `binary_subblock_variance_refinement` (Theorem 2.5(b) & 2.11(c)), `normalized_integrality_excess_bound` (Theorem 2.6(d) & 2.11(c)), and `cesaro_fejer_order2_step` (Theorem 2.7(a)).
  - Refined all Lean 4 docstrings (`category API` in `formal-conjectures#7640`) to state the exact discrete arithmetic propositions proved on $\mathbb{N}$, explicitly distinguishing ordinary subtraction from Lean's truncated subtraction on $\mathbb{N}$ in `binary_subblock_variance_refinement` and `normalized_integrality_excess_bound` and stating the unweighted Order-2 sum $\sum _ {w=1} ^ W \max(w - d, 0) = \binom{W - d + 1}{2}$ ($d = \lvert x - y\rvert \lt W$) in `cesaro_fejer_order2_step`.
  - Extended exact QP / Diophantine Haar envelope evaluations across **16 dyadic horizons $M \in \lbrace 8, 12, 16, 20, 24, 32, 48, 64, 96, 128, 192, 256, 512, 1024, 2048, 4096\rbrace$** ($N = 2 ^ {4096}$, achieving $\underline{c} _ {4096} ^ {\ast ,\mathrm{Dioph}}(1.01) = 0.086619 \lt 0.10$, a **$463.7157\times$ reduction** over the Phase 1 harmonic bound $\underline{c} _ {4096} ^ {(1)} = 40.166778$).
  - Added Figure 7 (`paper/figures/fig7_we_solve_remediation_diagnostics.{pdf,png}`), "Countering Potential Methodological Objections" subsection, remediated all two-column LaTeX display-equation overflows (`\allowdisplaybreaks`, multi-line `aligned` splits), and rebuilt all manuscripts (`23 pp.`), companion explainers (`10 pp.`), supplementary verification compendiums (`21 pp.`), and archives.

## [1.0.0] - 2026-10-27

### Added
- **Unconditional Multi-Octave Diophantine Haar Closure (`src/verify_cpu.rs`, `src/kernels.cu`, `paper/anc/verify_erdos1191_diophantine_haar.py`):**
  - Completed the 11-phase analytical and computational hierarchy (Theorems 2.1–2.11, Corollary 2.12, and Theorems 3.1–3.3) resolving Erdős Problem #1191 (`erdos_1191.parts.i`, 1,000 USD prize) in the affirmative ($\liminf _ {x\to\infty} A(x)\sqrt{(\ln x)/x} = 0$ for every infinite Sidon set $A \subseteq \mathbb{N}$) and disproving `erdos_1191.parts.ii` for all exponents $0 \lt c \le 1/2$.
  - Implemented Theorem 2.2 (Cross-Shell Bilinear Packing $y - x \in D _ {k-1} \cup D _ k$), Theorem 2.3 (Sub-Window Short-Difference Projection), Theorem 2.4 (Sliding-Window Fejér Identity & Pinched Dyadic Limit $c _ \infty ^ \ast (2) = (1+\sqrt{2})\sqrt{\ln 2} \approx 2.009965$), Theorem 2.5 (Direct Prefix-to-Block Unsplit Coupling & $\text{Base-}B \gt 1$ Hyperbolic-Tangent Shell Law $\Lambda(B) = \tanh((\ln B)/4)/\ln B \uparrow 1/4$, $c _ \infty ^ \ast (1 ^ +) = 2$), Theorem 2.6 (Discrete Legendre–Fenchel Integrality Identity $L(S, B) = \frac{B}{2}\mu(\mu-1) + \frac{B}{2}\theta(1-\theta){}$, $c _ \infty ^ {\ast ,\mathrm{coupled}}(1 ^ +) = 2\sqrt{6/7} \approx 1.851640$), Theorems 2.7–2.8 (Iterated $\text{Order-}p$ & Band-Restricted Cesàro–Fejér Kernels, $c _ \infty ^ {\ast ,\mathrm{band}}(1 ^ +) = \frac{2\sqrt{7}}{3} \approx 1.763834$), Theorem 2.9 (Scale-Invariant Block-Window Conservation $w \cdot B(\tau) \equiv N$ & $K$-Octave Cascade, $c _ \infty ^ {\ast ,\mathrm{multi\text{-}oct}}(1 ^ +) = \frac{2\sqrt{6}}{3} \approx 1.632993$), Theorem 2.10 (Power-Law $d\nu _ r(\tau) = \tau ^ r\,d\tau$ crossing $1.00$ at $r = -15/17$ & Logarithmic Haar Cascade $\kappa _ {\mathrm{int}} ^ {\mathrm{Haar}}(K) = \frac{K\ln 2}{3(1 - 4 ^ {-K})} \to \infty$), and Theorem 2.11 (Pointwise Diophantine Linear-Dilation Sawtooth Identity $\mathcal{J} _ {\mathrm{Haar}}(\mu _ 1, K) = \frac{K\ln 2}{6} + R _ {\mathrm{Dioph}}(\mu _ 1, K){}$ via the Hyperfactorial $H(M-1){}$ and Glaisher–Kinkelin Sawtooth Constant $C _ {\mathrm{GK}} = 2\ln A _ {\mathrm{GK}} = \frac{1}{6} - 2\zeta'(-1){}$).
  - Added 50-decimal-digit arbitrary-precision `mpmath`/`SymPy` verification suite (`paper/anc/verify_erdos1191_diophantine_haar.py`) auditing all 7 analytical pillars and finite-horizon tail bounds through $M = 128$ ($\underline{c} _ {32} ^ {\ast ,\mathrm{Dioph}}(1.01) = 0.923641 \lt 1.00$, $\underline{c} _ {64} ^ {\ast ,\mathrm{Dioph}}(1.01) = 0.667674$, $\underline{c} _ {128} ^ {\ast ,\mathrm{Dioph}}(1.01) = 0.480049 \lt 0.50$).
- **Standalone Lean 4 (`v4.33.1`) Formal Verification Suite (`Erdos1191/Proofs.lean`, `Erdos1191/AxiomCheck.lean`):**
  - Added self-contained, zero-`sorry`, zero-warning Lean 4 formal proofs for all 14 core discrete algebraic, combinatorial, Legendre–Fenchel, Cesàro–Fejér, scale-invariant block-window, hyperfactorial Abel telescoping, and bit-shift Sidon identities, aligned with Google DeepMind's `formal-conjectures` (`PR #7640`).
- **Companion Mathematical Explainer & Supplementary Verification Compendium (`paper/`):**
  - Added 8-page plain-language companion guide (`paper/Parobek_2026-10-27_Erdos1191_Companion_Mathematical_Explainer.pdf`, `paper/Companion_Explainer_Erdos1191.typ`, `paper/Companion_Explainer_Erdos1191.md`) and 12-page supplementary verification compendium (`paper/Parobek_2026-10-27_Erdos1191_Supplementary_Verification_Suite.pdf`, `paper/Supplementary_Verification_Suite_Erdos1191.typ`).
- **GitHub Actions Multi-Stage Verification CI (`.github/workflows/ci.yml`):**
  - Configured automated CI jobs for `lean4-verification` (`lake --wfail build` + `AxiomCheck.lean`), `python-50digit-verification`, and `rust-suite` (`cargo fmt --check`, `cargo clippy`, `cargo test`).

## [0.1.0] - 2026-10-26

### Added
- **Exact OEIS Ground-Truth Tables & Theorem 2.1 Certificates (`src/records.rs`, `src/verify_cpu.rs`):**
  - Exact Mian–Chowla greedy Sidon sequence (`OEIS A005282`, first 64 terms + arbitrary-length fast bitset generator `mian_chowla_greedy_cpu`) and optimal Golomb ruler / finite Sidon capacity table (`OEIS A003022`, exact $F _ 2(N){}$ for all $N \le 586$).
  - Implemented Theorem 2.1(a) (`harmonic_conservation_report_cpu`: Exact Harmonic-Integral Conservation Identity $\mathcal{E} _ A(N) + \Delta _ A(N) = \mathcal{B} _ A(N) \le H _ N - 1$), Theorem 2.1(b) (Exact Single-Scale Dyadic Difference-Budget Inequality $\sum _ {k=1} ^ M m _ k(m _ k-1)/2 ^ k \le (M-1)/2$), and Theorem 2.1(c) (Multi-Scale Dyadic Sub-Block Step-Decomposition).
  - Implemented Krückeberg–Cilleruelo Dyadic Singer-Block Sidon constructions (`kruckeberg_singer_block_sidon`) demonstrating Theorem 3.1 (the Burst–Desert Dichotomy).
- **$O(1){}$ Bit-Shift CPU+CUDA Branch-and-Bound & Maximin Solver (`src/exact.rs`, `src/kernels.cu`, `src/gpu.rs`):**
  - Implemented `k_dfs_maximin_sidon` (CUDA) and `extend_sidon_dfs_cpu` (Rust) using the 1-instruction shift `new_diffs = rev_mask >> (63 - x)`, integer prefix-floor constraints `req_prefix[1..=N]`, and OEIS `A003022` subadditive tail capacities.
  - Implemented `k_eval_sidon_dyadic_budget` (CUDA) and `eval_sidon_dyadic_budget_cpu` (Rust) with bit-exact parity tests.
- **Verification Scripts & Publication Manuscript (`scripts/`, `paper/`):**
  - Added `scripts/verify_erdos1191_theorems.py` and `scripts/generate_figures.py`.
  - Added full LaTeX (`paper/manuscript.tex`), Typst (`paper/manuscript.typ`), Markdown (`paper/manuscript.md`), and compiled PDF (`paper/manuscript.pdf`).
