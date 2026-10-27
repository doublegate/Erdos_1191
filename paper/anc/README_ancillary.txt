================================================================================
arXiv Ancillary Verification Suite — Erdős Problem #1191 ($1,000 Prize)
Paper: "Unconditional Diophantine Haar Closure, Exact Harmonic-Difference
        Conservation, and Multi-Octave Dyadic Budgets for Erdős Problem #1191"
Author: Lucas Parobek (parobek@gmail.com)
Repository: https://github.com/doublegate/Erdos/tree/main/erdos1191

Contents of anc/:
  1. verify_erdos1191_diophantine_haar.py
     Standalone zero-dependency Python 3 verification script that checks:
       - Theorem 2.1(a): Exact Harmonic-Integral Conservation Identity
         E_A(N) + Delta_A(N) == B_A(N) <= H_N - 1 on the Mian-Chowla greedy
         Sidon sequence (OEIS A005282).
       - Theorem 2.6(a): Discrete Legendre-Fenchel Integrality Excess
         L(S, B) = (B/2)*mu*(mu-1) + (B/2)*theta*(1-theta).
       - Theorem 2.11(a)-(d): Unconditional Diophantine Linear-Dilation
         Sawtooth Identity F_Haar(M) = M^2/2 - M(M-1)ln(M) + 2*ln(H(M-1))
         = (1/6)*ln(M) + 2*ln(A_GK) + 1/(360*M^2) - 1/(2520*M^4) + O(M^{-6})
         in 80-digit arbitrary-precision arithmetic, plus the M=8..128
         finite-horizon envelope certificates (< 1.00 at M=32: 0.923641;
         < 0.50 at M=128: 0.480049).
     Run via:
       python3 verify_erdos1191_diophantine_haar.py

  2. diophantine_haar_certificates.json
     Machine-readable JSON certificate output containing the Glaisher-Kinkelin
     sawtooth constants C_GK = +0.497508954055 and C_LF = -0.002491045945,
     80-digit Euler-Maclaurin verification residuals, and the M=8..128
     multi-octave Diophantine Haar bounds.

  3. 1191.lean
     Lean 4 / Mathlib formalization blueprint for google-deepmind/formal-conjectures.
================================================================================
