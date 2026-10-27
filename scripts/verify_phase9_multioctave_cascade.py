#!/usr/bin/env python3
"""
Phase 9 Symbolic & Numerical Verification Suite for Erdős Problem #1191:
Theorem 2.9 — Scale-Invariant Block-Window Conservation w * B(tau) == N,
K-Octave Integrality Cascade kappa_int^{band}(2^{-K}) = 1/(2*(1 + 2^{-K} + 4^{-K})) -> 1/2,
and Multi-Octave Deflation of c_infty^* to 2*sqrt(6)/3 ≈ 1.632993 & c_16^* < 2.00.

Verifies:
  1. Theorem 2.9(a) Symbolic Proof (Scale-Invariant Block-Window Conservation Law):
     For any macroscopic sub-window width w = tau * W (tau in (0, 1]) tiling a domain of
     length N into B(tau) = N / (tau * W) sub-blocks of bin size w = tau * W, the exact
     Discrete Legendre–Fenchel integrality excess per window is (w / 2)*theta(tau)*(1 - theta(tau)),
     so the total integrality excess across all B(tau) blocks satisfies the Scale-Invariant
     Block-Window Product Identity:
       E_int(w) = B(tau) * (w / 2) * theta(tau)*(1 - theta(tau))
                = (N / (tau * W)) * ((tau * W) / 2) * theta(tau)*(1 - theta(tau))
               == (N / 2) * theta(tau)*(1 - theta(tau)),
     where tau * tau^{-1} == 1 cancels identically across every dyadic octave tau in (0, 1].
  2. Theorem 2.9(b) Symbolic & Numerical Proof (K-Octave Integrality Cascade &
     Sub-1.64 Continuum Limit):
     Summing F_w(A) over the top K dyadic octaves w in [W_0 + 1, W] (W_0 = floor(2^{-K} W),
     alpha = 2^{-K}) harvests the unit-cycle integrality average int_0^1 theta*(1 - theta) d_theta = 1/6
     uniformly across tau in [alpha, 1] while normalizing by the truncated Sidon ceiling
     binom(W + 1, 3) - binom(W_0 + 1, 3) ~ ((1 - alpha^3) / 6) * W^3:
       kappa_int^{band}(2^{-K}) = (3 / (1 - 2^{-3K})) * int_{2^{-K}}^1 (1/6) d_tau
                                = (1 - 2^{-K}) / (2 * (1 - 2^{-3K}))
                               == 1 / (2 * (1 + 2^{-K} + 4^{-K}))
                                = 2^{2K - 1} / (4^K + 2^K + 1).
     Exact rational cascade:
       * K = 1 (alpha = 1/2, Phase 8): kappa_int = 2/7 ≈ 0.28571429,
         mu_scale(1) = 18/7, c_infty^{*,(1)}(1^+) = 2*sqrt(7)/3 ≈ 1.76383421.
       * K = 2 (alpha = 1/4, 2-Octave): kappa_int = 8/21 ≈ 0.38095238,
         mu_scale(2) = 58/21, c_infty^{*,(2)}(1^+) = 2*sqrt(21/29) = 2*sqrt(609)/29 ≈ 1.70186400.
       * K = 3 (alpha = 1/8, 3-Octave): kappa_int = 32/73 ≈ 0.43835616,
         mu_scale(3) = 210/73, c_infty^{*,(3)}(1^+) = 2*sqrt(73/105) ≈ 1.66761878.
       * K = 4 (alpha = 1/16, 4-Octave): kappa_int = 128/273 ≈ 0.46886447,
         mu_scale(4) = 802/273, c_infty^{*,(4)}(1^+) = 2*sqrt(273/401) ≈ 1.65027006.
       * K -> infty (alpha -> 0^+, Full Multi-Octave Continuum Limit):
         lim_{K->infty} kappa_int^{band}(2^{-K}) = 1/2 = 0.50000000,
         mu_scale(infty) = 3, c_infty^{*,(infty)}(1^+) = 2*sqrt(2/3) = 2*sqrt(6)/3 ≈ 1.63299316.
  3. Theorem 2.9 Multi-Octave Envelope Bounds:
     Because the backward-lag integral beta_{s,band}(alpha)^2 = (3/(1-alpha^3)) int_alpha^1 tau^2 beta_s(tau)^2 d_tau
     weights sub-windows by tau^2 (keeping >= 87.5% of its weight in the upper octave [1/2, 1])
     whereas the scale-invariant integrality excess weights sub-windows uniformly by 1 d_tau,
     taking the multi-octave cascade strictly increases the per-scale weight gamma_s^{(K)} for
     every K >= 1 and every s >= 1, pushing the finite-horizon tail bound c_M^{*,multi-oct}(1.01)
     strictly below 2.000000 already at M = 16 (down from M = 20 in Phase 8, M = 24 in Phase 7,
     and M = 32 in Phase 6) and down to ~1.71 at M = 64.
"""

import json
import math

import sympy as sp


def psi_lag(delta: float) -> float:
    if delta <= 1e-12:
        return 1.0 - 0.25 * delta
    if delta >= 1.0:
        return 2.0 / 3.0
    return (2.0 / (3.0 * delta)) * (1.0 - (1.0 - delta) ** 1.5)


def beta_s_tau_sq(s: int, tau: float) -> float:
    if s == 0:
        return (4.0 / 3.0) ** 2
    base_inc = 1.0 - 1.0 / math.sqrt(2.0)
    d1 = tau * (2.0 ** (-1 - s))
    d0 = tau * (2.0 ** (-s))
    lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
    return (lag_inc / base_inc) ** 2


def beta_s_band_sq(s: int, alpha: float = 0.0, n_simp: int = 64) -> float:
    """
    Compute the Band-Restricted / Multi-Octave Order-2 Cesàro–Fejér lag-amplification factor
      beta_{s,band}(alpha)^2 := (3 / (1 - alpha^3)) * int_alpha^1 tau^2 * beta_s(tau)^2 d_tau
    via composite Simpson quadrature on [alpha, 1] for alpha in [0, 1).
    At alpha = 0 (K -> infty multi-octave continuum limit), 3 / (1 - alpha^3) = 3.
    """
    if s == 0:
        return (4.0 / 3.0) ** 2
    a, b = alpha, 1.0
    h = (b - a) / n_simp
    total = (a * a) * beta_s_tau_sq(s, a) + (b * b) * beta_s_tau_sq(s, b)
    for i in range(1, n_simp):
        tau = a + i * h
        weight = 4.0 if (i % 2 == 1) else 2.0
        total += weight * (tau * tau) * beta_s_tau_sq(s, tau)
    integral = (h / 3.0) * total
    return (3.0 / (1.0 - alpha**3)) * integral


def l_pigeon(s_val: int, bins: int) -> int:
    if bins <= 0 or s_val <= bins:
        return 0
    q, t = divmod(s_val, bins)
    return bins * (q * (q - 1) // 2) + t * q


def verify_symbolic_theorem_2_9():
    N, W, tau, theta, alpha, _K = sp.symbols("N W tau theta alpha K", positive=True)

    # 1. Theorem 2.9(a) Symbolic Proof of Scale-Invariant Block-Window Conservation:
    #    w(tau) = tau * W, B(tau) = N / (tau * W) => w(tau) * B(tau) == N,
    #    and E_int(w) = B(tau) * (w(tau) / 2) * theta * (1 - theta) == (N / 2) * theta * (1 - theta).
    w_tau = tau * W
    B_tau = N / (tau * W)
    product_w_B = sp.simplify(w_tau * B_tau)
    assert product_w_B == N, (
        f"Scale-invariant product w(tau)*B(tau) != N: {product_w_B}"
    )

    E_int_tau = sp.simplify(B_tau * (w_tau / 2) * theta * (1 - theta))
    expected_E_int = (N / 2) * theta * (1 - theta)
    assert sp.simplify(E_int_tau - expected_E_int) == 0, (
        "E_int(w) scale invariance failed!"
    )

    # 2. Theorem 2.9(b) Symbolic Proof of the K-Octave Integrality Cascade:
    inner_avg = sp.integrate(theta * (1 - theta), (theta, 0, 1))
    assert inner_avg == sp.Rational(1, 6)

    kappa_alpha = sp.simplify(
        (sp.Rational(3, 1) / (1 - alpha**3)) * sp.integrate(inner_avg, (tau, alpha, 1))
    )
    expected_kappa_alpha = sp.Rational(1, 2) / (1 + alpha + alpha**2)
    assert sp.simplify(kappa_alpha - expected_kappa_alpha) == 0

    # Substitute alpha = 1 / q (with q = 2^K) and verify rational form q^2 / (2*(q^2 + q + 1)) = 2^{2K-1} / (4^K + 2^K + 1)
    q = sp.Symbol("q", positive=True)
    kappa_q = sp.simplify(kappa_alpha.subs(alpha, 1 / q))
    expected_kappa_q = (q**2) / (2 * (q**2 + q + 1))
    assert sp.simplify(kappa_q - expected_kappa_q) == 0

    # Verify exact limit as K -> infty (alpha -> 0^+):
    lim_kappa = sp.limit(kappa_alpha, alpha, 0, dir="+")
    assert lim_kappa == sp.Rational(1, 2), f"Expected limit 1/2, got {lim_kappa}"

    lim_mu_scale = 2 * (1 + lim_kappa)
    assert lim_mu_scale == sp.Rational(3, 1)

    c_inf_multi = sp.simplify(2 / sp.sqrt(1 + lim_kappa))
    assert sp.simplify(c_inf_multi - 2 * sp.sqrt(6) / 3) == 0

    # Verify exact rational cascade for K = 1, 2, 3, 4
    expected_kappas = {
        1: sp.Rational(2, 7),
        2: sp.Rational(8, 21),
        3: sp.Rational(32, 73),
        4: sp.Rational(128, 273),
    }
    cascade_table = []
    prev_k_val = 0.0
    prev_c_val = float("inf")
    for k_oct in range(1, 7):
        alpha_k = sp.Rational(1, 2**k_oct)
        kap_k = sp.simplify(kappa_alpha.subs(alpha, alpha_k))
        if k_oct in expected_kappas:
            assert kap_k == expected_kappas[k_oct], f"Mismatch at K={k_oct}: {kap_k}"
        mu_k = sp.simplify(2 * (1 + kap_k))
        c_inf_k = sp.simplify(2 / sp.sqrt(1 + kap_k))
        kap_f = float(kap_k.evalf())
        c_f = float(c_inf_k.evalf())
        assert kap_f > prev_k_val and kap_f < 0.5
        assert c_f < prev_c_val and c_f > float(c_inf_multi.evalf())
        prev_k_val = kap_f
        prev_c_val = c_f
        cascade_table.append(
            {
                "K_octaves": k_oct,
                "alpha": f"1/{2**k_oct}",
                "kappa_int_exact": str(kap_k),
                "kappa_int_float": round(kap_f, 8),
                "mu_scale_exact": str(mu_k),
                "mu_scale_float": round(float(mu_k.evalf()), 8),
                "c_infty_exact": str(c_inf_k),
                "c_infty_float": round(c_f, 8),
            }
        )

    # 3. Compare per-scale weights gamma_s across Phase 6 (Order-1), Phase 7 (CF2),
    #    Phase 8 (K=1 Band), and Phase 9 (K=2, K=3, and K->infty Multi-Octave)
    scale_weight_rows = []
    for s_idx in range(1, 6):
        b1_sq = beta_s_tau_sq(s_idx, 1.0)
        b_k1_sq = beta_s_band_sq(s_idx, 0.5)
        b_k2_sq = beta_s_band_sq(s_idx, 0.25)
        b_k3_sq = beta_s_band_sq(s_idx, 0.125)
        b_inf_sq = beta_s_band_sq(s_idx, 0.0)

        g_p6 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s_idx))
        g_p7_raw = b_k1_sq + 0.25 * (1.0 - 2.0 ** (-s_idx))
        g_p8_k1_raw = b_k1_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s_idx))
        g_p9_k2_raw = b_k2_sq + (8.0 / 21.0) * (1.0 - 2.0 ** (-s_idx))
        g_p9_k3_raw = b_k3_sq + (32.0 / 73.0) * (1.0 - 2.0 ** (-s_idx))
        g_p9_inf = b_inf_sq + 0.5 * (1.0 - 2.0 ** (-s_idx))

        assert g_p7_raw < g_p8_k1_raw < g_p9_k2_raw < g_p9_k3_raw < g_p9_inf, (
            f"Strict K-octave monotonicity of Order-2 weight failed at s={s_idx}!"
        )
        assert g_p9_inf > g_p9_k3_raw > max(g_p6, g_p8_k1_raw), (
            f"Multi-octave K>=3 failed to strictly dominate Order-1 and Phase 8 at s={s_idx}!"
        )
        scale_weight_rows.append(
            {
                "s": s_idx,
                "beta_s_k1_sq": round(b_k1_sq, 8),
                "beta_s_inf_sq": round(b_inf_sq, 8),
                "gamma_s_p6_order1": round(g_p6, 8),
                "gamma_s_p7_cf2_raw": round(g_p7_raw, 8),
                "gamma_s_p8_k1_raw": round(g_p8_k1_raw, 8),
                "gamma_s_p9_k2_raw": round(g_p9_k2_raw, 8),
                "gamma_s_p9_k3_raw": round(g_p9_k3_raw, 8),
                "gamma_s_p9_multi_inf": round(g_p9_inf, 8),
            }
        )

    return {
        "scale_invariant_product_w_times_B_eq_N": True,
        "kappa_int_K_formula": "2^(2K-1) / (4^K + 2^K + 1)",
        "kappa_int_limit_exact": "1/2",
        "mu_scale_limit_exact": "3",
        "c_infty_multi_oct_exact": "2*sqrt(6)/3",
        "c_infty_multi_oct_float": round(float(c_inf_multi.evalf()), 10),
        "k_octave_cascade_table": cascade_table,
        "per_scale_weights_s1_to_s5": scale_weight_rows,
    }


def generate_mian_chowla(terms: int):
    seq = []
    diffs = set()
    cand = 1
    while len(seq) < terms:
        ok = True
        for a in seq:
            if (cand - a) in diffs:
                ok = False
                break
        if ok:
            for a in seq:
                diffs.add(cand - a)
            seq.append(cand)
        cand += 1
    return seq


def verify_multioctave_cascade_on_mian_chowla():
    """
    Verify Theorem 2.9 exact integer K-Octave Band-Restricted Cesàro–Fejér identities &
    sub-block scale-invariant pigeonhole lower-bound chains on the Mian–Chowla sequence
    across dyadic scales r in {3, ..., 7} for K = 1 (W_0 = W/2), K = 2 (W_0 = W/4),
    and K = 3 (W_0 = W/8).
    """
    mc = generate_mian_chowla(64)
    m_scales = 10
    n_horizon = 1 << m_scales
    in_h = [x for x in mc if x <= n_horizon]
    diffs = [
        in_h[j] - in_h[i] for i in range(len(in_h)) for j in range(i + 1, len(in_h))
    ]

    pref = [0] * (n_horizon + 1)
    idx = 0
    for x in range(1, n_horizon + 1):
        while idx < len(in_h) and in_h[idx] <= x:
            idx += 1
        pref[x] = idx

    def win_c(u: int, w: int) -> int:
        hi = min(u, n_horizon)
        lo = max(0, u - w)
        return pref[hi] - pref[lo]

    def refined_pigeon_lb_for_w(w: int, r_scale: int, ref_levels: int = 2) -> int:
        t_r = min(1 << r_scale, n_horizon)
        blocks = [(1, t_r + w - 1)]
        for k in range(r_scale + 1, m_scales + 1):
            t_prev = 1 << (k - 1)
            t_curr = min(1 << k, n_horizon)
            if t_prev >= n_horizon:
                break
            blocks.append((t_prev + w, t_curr + w - 1))
        lb_tot = 0
        n_sub = 1 << ref_levels
        for a, b in blocks:
            total_len = b - a + 1
            for j in range(n_sub):
                sub_lo = a + (j * total_len) // n_sub
                sub_hi = a + ((j + 1) * total_len) // n_sub - 1
                if sub_hi < sub_lo:
                    continue
                bins = sub_hi - sub_lo + 1
                s_val = sum(win_c(u, w) for u in range(sub_lo, sub_hi + 1))
                lb_tot += l_pigeon(s_val, bins)
        return lb_tot

    rows = []
    for r in range(3, 8):
        w_cap = 1 << r
        row_k = {"r": r, "W": w_cap}
        prev_lb = -1
        prev_rhs = -1
        for k_oct in [1, 2, 3]:
            w0 = w_cap >> k_oct
            lb_k = 0
            win_k = 0
            for w in range(w0 + 1, w_cap + 1):
                lb_w = refined_pigeon_lb_for_w(w, r, ref_levels=2)
                fw_exact = sum(
                    win_c(u, w) * (win_c(u, w) - 1) // 2
                    for u in range(1, n_horizon + w)
                )
                assert lb_w <= fw_exact
                lb_k += lb_w
                win_k += fw_exact

            rhs_k = 0
            for d in diffs:
                if d < w_cap:
                    rem_w = w_cap - d
                    term_w = rem_w * (rem_w + 1) // 2
                    rem_0 = max(w0 - d, 0)
                    term_0 = rem_0 * (rem_0 + 1) // 2
                    rhs_k += term_w - term_0

            cap_w = (w_cap - 1) * w_cap * (w_cap + 1) // 6
            cap_0 = (w0 - 1) * w0 * (w0 + 1) // 6 if w0 >= 2 else 0
            cap_k = cap_w - cap_0

            assert win_k == rhs_k, f"Multi-octave identity failed at r={r}, K={k_oct}"
            assert lb_k <= rhs_k <= cap_k
            assert lb_k >= prev_lb and rhs_k >= prev_rhs
            prev_lb = lb_k
            prev_rhs = rhs_k

            row_k[f"K{k_oct}_w0"] = w0
            row_k[f"K{k_oct}_pigeonhole_lb"] = lb_k
            row_k[f"K{k_oct}_exact_rhs"] = rhs_k
            row_k[f"K{k_oct}_binom_cap"] = cap_k
        rows.append(row_k)
    return rows


def verify_multioctave_envelope_bounds():
    kappa_pref = 2.5 - 1.5 * math.sqrt(2.0) + 0.5 * math.log(1.0 + math.sqrt(2.0))
    lam_2 = math.tanh(math.log(2.0) / 4.0) / math.log(2.0)
    lam_101 = math.tanh(math.log(1.01) / 4.0) / math.log(1.01)
    scale_ratio = lam_101 / lam_2
    beta0_sq = (4.0 / 3.0) ** 2

    c_inf_p8 = 2.0 * math.sqrt(7.0 / 9.0)
    c_inf_p9_k2 = 2.0 * math.sqrt(21.0 / 29.0)
    c_inf_p9_k3 = 2.0 * math.sqrt(73.0 / 105.0)
    c_inf_p9_multi = 2.0 * math.sqrt(2.0 / 3.0)

    rows = []
    prev_p9_tail = float("inf")
    for m_val in [8, 12, 16, 20, 24, 32, 48, 64]:
        r0 = m_val // 2
        q_p8_tail = 0.0
        q_p9_k2_tail = 0.0
        q_p9_k3_tail = 0.0
        q_p9_multi_tail = 0.0

        for k in range(r0 + 1, m_val + 1):
            alpha_k = 1.0 / math.sqrt(k * math.log(2.0)) - 1.0 / math.sqrt(
                2.0 * (k - 1) * math.log(2.0)
            )
            s_k = max(0, k - 2) // 2
            f8 = beta0_sq
            f9_k2 = beta0_sq
            f9_k3 = beta0_sq
            f9_multi = beta0_sq
            for s in range(1, s_k + 1):
                b1_sq = beta_s_tau_sq(s, 1.0)
                b_k1_sq = beta_s_band_sq(s, 0.5)
                b_k2_sq = beta_s_band_sq(s, 0.25)
                b_k3_sq = beta_s_band_sq(s, 0.125)
                b_multi_sq = beta_s_band_sq(s, 0.0)
                g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s))

                f8 += 2.0 * max(g1, b_k1_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s)))
                f9_k2 += 2.0 * max(g1, b_k2_sq + (8.0 / 21.0) * (1.0 - 2.0 ** (-s)))
                f9_k3 += 2.0 * max(g1, b_k3_sq + (32.0 / 73.0) * (1.0 - 2.0 ** (-s)))
                f9_multi += 2.0 * max(g1, b_multi_sq + 0.5 * (1.0 - 2.0 ** (-s)))

            q_p8_tail += scale_ratio * f8 * alpha_k * alpha_k
            q_p9_k2_tail += scale_ratio * f9_k2 * alpha_k * alpha_k
            q_p9_k3_tail += scale_ratio * f9_k3 * alpha_k * alpha_k
            q_p9_multi_tail += scale_ratio * f9_multi * alpha_k * alpha_k

        rhs_tail = 0.5 * (m_val - r0)
        c_p8_tail = math.sqrt(rhs_tail / q_p8_tail)
        c_p9_k2_tail = math.sqrt(rhs_tail / q_p9_k2_tail)
        c_p9_k3_tail = math.sqrt(rhs_tail / q_p9_k3_tail)
        c_p9_multi_tail = math.sqrt(rhs_tail / q_p9_multi_tail)

        assert c_p9_multi_tail < c_p9_k3_tail < c_p9_k2_tail < c_p8_tail
        assert c_p9_multi_tail < prev_p9_tail
        prev_p9_tail = c_p9_multi_tail

        q_p8_full = 0.0
        q_p9_multi_full = 0.0
        prev_req = 0.0
        for k in range(3, m_val + 1):
            req_k = math.sqrt((2.0**k) / (k * math.log(2.0)))
            rho_k = req_k - prev_req
            prev_req = req_k
            rho_sq_2k = (rho_k * rho_k) / (2.0**k)
            s_k = max(0, k - 2) // 2

            f8 = beta0_sq
            f9_multi = beta0_sq
            for s in range(1, s_k + 1):
                b1_sq = beta_s_tau_sq(s, 1.0)
                b_k1_sq = beta_s_band_sq(s, 0.5)
                b_multi_sq = beta_s_band_sq(s, 0.0)
                g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s))
                f8 += 2.0 * max(g1, b_k1_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s)))
                f9_multi += 2.0 * max(g1, b_multi_sq + 0.5 * (1.0 - 2.0 ** (-s)))

            w_bil = 1.0 if k == m_val else 0.5
            q_p8_full += scale_ratio * f8 * rho_sq_2k + (w_bil + 0.5 * kappa_pref) / (
                k * math.log(2.0)
            )
            q_p9_multi_full += scale_ratio * f9_multi * rho_sq_2k + (
                w_bil + 0.5 * kappa_pref
            ) / (k * math.log(2.0))

        rhs_full = 0.5 * (m_val + 2)
        c_p8_full = math.sqrt(rhs_full / q_p8_full)
        c_p9_multi_full = math.sqrt(rhs_full / q_p9_multi_full)
        assert c_p9_multi_full < c_p8_full
        assert c_p9_multi_full < c_inf_p9_multi < c_p9_multi_tail

        rows.append(
            {
                "M": m_val,
                "c_p8_band_tail_b101": round(c_p8_tail, 6),
                "c_p9_k2_tail_b101": round(c_p9_k2_tail, 6),
                "c_p9_k3_tail_b101": round(c_p9_k3_tail, 6),
                "c_p9_multi_tail_b101": round(c_p9_multi_tail, 6),
                "c_p8_band_full_b101": round(c_p8_full, 6),
                "c_p9_multi_full_b101": round(c_p9_multi_full, 6),
                "c_inf_p8_limit": round(c_inf_p8, 6),
                "c_inf_p9_k2_limit": round(c_inf_p9_k2, 6),
                "c_inf_p9_k3_limit": round(c_inf_p9_k3, 6),
                "c_inf_p9_multi_limit": round(c_inf_p9_multi, 6),
            }
        )

    # Verify sub-2.000000 breakthrough already at M >= 16 (down from M = 20 in Phase 8)
    for r in rows:
        if r["M"] >= 16:
            assert r["c_p9_multi_tail_b101"] < 2.0, (
                f"Expected < 2.0 at M={r['M']}, got {r['c_p9_multi_tail_b101']}"
            )

    return rows


def main():
    sym_res = verify_symbolic_theorem_2_9()
    mc_rows = verify_multioctave_cascade_on_mian_chowla()
    env_rows = verify_multioctave_envelope_bounds()
    out = {
        "theorem_2_9_symbolic": sym_res,
        "theorem_2_9_mian_chowla_multioctave": mc_rows,
        "theorem_2_9_multioctave_envelope_bounds": env_rows,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
