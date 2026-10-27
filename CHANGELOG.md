# Changelog

All notable changes to **erdos1191** are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
  - Added full LaTeX (`paper/manuscript.tex`), Typst (`paper/manuscript.typ`), Markdown (`paper/manuscript.md`), compiled PDF (`paper/manuscript.pdf`), and arXiv bundle (`paper/arxiv_submission.tar.gz`).
