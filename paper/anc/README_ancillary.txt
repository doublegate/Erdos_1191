================================================================================
Ancillary Verification Suite — Erdős Problem #1191 ($1,000 Prize)
Paper: "Unconditional Diophantine Haar Closure, Exact Harmonic-Difference
        Conservation, and Multi-Octave Dyadic Budgets for Erdős Problem #1191"
Author: Lucas Parobek (parobek@gmail.com)
Paper PDF: https://github.com/doublegate/Erdos_1191/blob/main/paper/Parobek_2026-10-27_Erdos1191_Unconditional_Diophantine_Haar_Closure_Sidon_Sets.pdf
Repository: https://github.com/doublegate/Erdos_1191
formal-conjectures PR: https://github.com/google-deepmind/formal-conjectures/pull/7640
teorth/erdosproblems PR: https://github.com/teorth/erdosproblems/pull/511
We-Solve/research Issue: https://github.com/We-Solve/research/issues/881

Contents of anc/:
  1. 1191.lean (and 1191.lean.txt for plain-text viewers)
     Self-contained Lean 4 (v4.33.1) formal verification module (17 verified
     declarations in strictly chronological manuscript order: 15 `erdos_1191.variants.*`
     theorems + 2 decidable Sidon witness theorems `mianChowlaPrefix10_is_sidon`
     and `maximinWitnessN32_is_sidon`), zero `sorry`s, zero warnings, zero
     custom axioms. Matches `Erdos1191/Proofs.lean` and `google-deepmind/formal-conjectures#7640`.
     Run via:
       lake --wfail build && lake env lean Erdos1191/AxiomCheck.lean
       (or standalone: lean paper/anc/1191.lean)

  2. verify_erdos1191_diophantine_haar.py
     Standalone zero-dependency Python 3 verification script that checks:
       - Theorem 2.1(a): Exact Harmonic-Integral Conservation Identity
         E_A(N) + Delta_A(N) == B_A(N) <= H_N - 1 on the Mian-Chowla greedy
         Sidon sequence (OEIS A005282).
       - Theorem 2.2(b), Remark 2.3 & Lemma 2.4: Corrected Cross-Shell Bilinear
         Inequality verified on the powers-of-two Sidon set A = {2^n : n >= 0}
         (sum = 4 - (2M+4)/2^M <= M + 1 - 2^{-(M-1)}, 2.75 <= 3.75 at M=3) and
         KKT convex-cone increment positivity 16*k^2 < 18*(k^2 - 1) for k >= 4.
       - Theorem 2.5(b) & Theorem 2.6(a)-(d): Binary sub-block variance
         refinement, Discrete Legendre-Fenchel Integrality Excess
         L(S, B) = (B/2)*mu*(mu-1) + (B/2)*theta*(1-theta), and quarter-square
         sawtooth ceiling 4*r*(B-r) + |2r-B|^2 = B^2 (relative scaling 1/(4*mu^2)).
       - Theorem 2.7(a) & Theorem 2.11(a)-(c): Order-2 Cesaro-Fejer kernel
         recurrence and Pointwise Diophantine Linear-Dilation Sawtooth Identity
         F_Haar(M) = M^2/2 - M(M-1)ln(M) + 2*ln(H(M-1))
         = (1/6)*ln(M) + 2*ln(A_GK) + 1/(360*M^2) - 1/(2520*M^4) + O(M^{-6})
         in 80-digit arbitrary-precision arithmetic, plus the M=8..128
         unit-occupancy collision-scale certificates (< 1.00 at M=32: 0.923641;
         < 0.50 at M=128: 0.480049) and fixed-K asymptotic continuum limit 2.00.
     Run via:
       python3 paper/anc/verify_erdos1191_diophantine_haar.py

  3. diophantine_haar_certificates.json
     Machine-readable JSON certificate output containing the Glaisher-Kinkelin
     sawtooth constants C_GK = +0.497508954055 and C_LF = -0.002491045945,
     80-digit Euler-Maclaurin verification residuals, and the M=8..128
     multi-octave Diophantine Haar bounds.
================================================================================
