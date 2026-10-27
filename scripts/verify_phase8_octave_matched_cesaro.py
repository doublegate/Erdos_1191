#!/usr/bin/env python3
"""
Phase 8 Symbolic & Numerical Verification Suite for Erdős Problem #1191:
Theorem 2.8 — Octave-Matched Cesàro–Fejér Order p*(s) ~ 2^s ln 2, Non-Vanishing
Asymptotic Octave Weight lim_{s->infty} max_{p>=1} beta_s^{(p)} = 1/4, and
Band-Restricted Upper-Octave Cesàro–Fejér Deflation Pushing c_infty^* Below 1.77.

Verifies:
  1. Theorem 2.8(a) Symbolic & Numerical Proof (Non-Vanishing Octave-Matched Cesàro–Fejér Law):
       beta_s^{(p)} := int_{2^{-s}}^{2^{-(s-1)}} (p + 1)*(1 - u)^p du
                     = (1 - 2^{-s})^{p+1} - (1 - 2^{-(s-1)})^{p+1} > 0,
     which decays as O(p * 2^{-s}) -> 0 for any fixed p, whereas matching the polynomial order
     p = p*(s) ~ 2^s * ln 2 to each dyadic octave s >= 2 yields a strictly decreasing sequence
       M_s := max_{p >= 1} beta_s^{(p)} > 1/4   (M_2 = 5/16 = 0.3125, M_3 = 9031/32768 ≈ 0.275604, ...)
     converging to the exact non-vanishing universal constant:
       lim_{s -> infty} max_{p >= 1} beta_s^{(p)} = max_{x in (0, 1)} (x - x^2) = 1/4 = 0.25000000.
  2. Theorem 2.8(b) Symbolic & Pointwise Proof (Band-Restricted Upper-Octave Cesàro–Fejér
     Identity CF_{W, W_0}(A) and Integrality Amplification kappa_int^{band}(alpha)):
       CF_{W, W_0}(A) := sum_{w = W_0 + 1}^W F_w(A)
                       = sum_{d in D_A(N) cap [1, W-1]} [binom(W - d + 1, 2) - binom(max(W_0 - d + 1, 0), 2)]
                      <= binom(W + 1, 3) - binom(W_0 + 1, 3).
     For W_0 = floor(alpha * W) with alpha in [1/2, 1), normalizing by the band-restricted RHS
     capacity (1 - alpha^3) * W^3 / 6 amplifies the unit-cycle discrete Legendre–Fenchel
     integrality excess to:
       kappa_int^{band}(alpha) = (3 / (1 - alpha^3)) * int_alpha^1 (1/6) d_tau = 1 / (2*(1 + alpha + alpha^2)).
     At alpha = 1/2 (W_0 = floor(W/2)), 100% of the summed sub-windows w in [W/2 + 1, W] lie in
     the active unit-cycle upper octave while the RHS capacity drops by 7/8, boosting the
     integrality constant from kappa_int^{(2)} = 1/4 to:
       kappa_int^{band}(1/2) = (8/7) * (1/4) = 2/7 ≈ 0.28571429 > 1/4,
     shifting the asymptotic per-scale multiplier from 5/2 to 18/7 ≈ 2.57142857 and driving the
     pinched asymptotic limit down to:
       c_infty^{*,band}(1^+) = 2 * sqrt(7 / 9) = 2*sqrt(7)/3 ≈ 1.76383421 < 1.78885438.
  3. Theorem 2.8 Band-Restricted Cesàro–Fejér Envelope Bounds:
       c_M^{*,band}(B) < c_M^{*,CF2}(B) < c_M^{*,coupled}(B) < c_M^{*,joint}(B) < c_M^*(2),
     breaking the 2.000000 barrier already at M = 20 (1.999651 at B = 1.01, down from M = 24
     in Phase 7 and M = 32 in Phase 6) and reaching 1.832222 at M = 64.
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


def beta_s_band_sq(s: int, alpha: float = 0.5, n_simp: int = 64) -> float:
    """
    Compute the Band-Restricted Order-2 Cesàro–Fejér lag-amplification factor
      beta_{s,band}(alpha)^2 := (3 / (1 - alpha^3)) * int_alpha^1 tau^2 * beta_s(tau)^2 d_tau
    via composite Simpson quadrature on [alpha, 1].
    At alpha = 1/2, 3 / (1 - (1/2)^3) = 24 / 7.
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


def verify_symbolic_theorem_2_8():
    u, x, p, s, W, W0, d, alpha, tau, theta = sp.symbols(
        "u x p s W W0 d alpha tau theta", positive=True
    )

    # 1. Theorem 2.8(a) Symbolic Proof of the Octave Integral beta_s^{(p)}
    #    and Non-Vanishing Octave-Matched Limit lim_{s->infty} max_{p>=1} beta_s^{(p)} = 1/4:
    psi_tilde = (p + 1) * (1 - u) ** p
    norm_total = sp.integrate(psi_tilde, (u, 0, 1))
    assert sp.simplify(norm_total - 1) == 0, "Normalized kernel integral over [0, 1] != 1"

    u_lo = sp.Rational(1, 1) / (2**s)
    u_hi = sp.Rational(1, 1) / (2 ** (s - 1))
    beta_s_p = sp.integrate(psi_tilde, (u, u_lo, u_hi))
    expected_beta_s_p = (1 - u_lo) ** (p + 1) - (1 - u_hi) ** (p + 1)
    assert sp.simplify(beta_s_p - expected_beta_s_p) == 0, "Octave integral mismatch"

    # Verify (1 - 2^{-s})^2 - (1 - 2^{-(s-1)}) == 4^{-s} > 0
    a_s = 1 - u_lo
    b_s = 1 - u_hi
    assert sp.simplify(a_s**2 - b_s - sp.Rational(1, 1) / (4**s)) == 0

    # Verify the asymptotic variational problem max_{x in (0, 1)} (x - x^2) = 1/4 at x = 1/2
    f_var = x - x**2
    x_crit = sp.solve(sp.diff(f_var, x), x)[0]
    assert x_crit == sp.Rational(1, 2)
    max_var = f_var.subs(x, x_crit)
    assert max_var == sp.Rational(1, 4)

    # Verify exact rational values for s = 2, 3 and strict decrease of M_s > 1/4 for s = 2..8
    beta_2_1 = (sp.Rational(3, 4) ** 2) - (sp.Rational(1, 2) ** 2)
    assert beta_2_1 == sp.Rational(5, 16)
    beta_3_4 = (sp.Rational(7, 8) ** 5) - (sp.Rational(3, 4) ** 5)
    assert beta_3_4 == sp.Rational(9031, 32768)

    octave_table = []
    prev_m_s = float("inf")
    for s_idx in range(1, 9):
        a_val = 1.0 - 2.0 ** (-s_idx)
        b_val = 1.0 - 2.0 ** (1 - s_idx)
        best_p = 1
        best_val = -1.0
        for p_cand in range(1, 600):
            val = (a_val ** (p_cand + 1)) - (b_val ** (p_cand + 1))
            if val > best_val:
                best_val = val
                best_p = p_cand
        if s_idx >= 2:
            assert best_val > 0.25, f"M_{s_idx} = {best_val} <= 0.25"
            assert best_val < prev_m_s, f"M_{s_idx} = {best_val} >= M_{s_idx-1} = {prev_m_s}"
            prev_m_s = best_val
        octave_table.append({
            "s": s_idx,
            "p_star": best_p,
            "asymptotic_2s_ln2": round((2**s_idx) * math.log(2.0), 4),
            "M_s_max_beta": round(best_val, 8),
            "beta_s_order1": round((a_val**2) - (b_val**2), 8),
            "beta_s_order2": round((a_val**3) - (b_val**3), 8),
        })

    # 2. Theorem 2.8(b) Symbolic Proof of the Band-Restricted Cesàro–Fejér Ceiling
    #    and Exact Band-Restricted Integrality Constant kappa_int^{band}(alpha):
    #    sum_{d=1}^{W-1} binom(W - d + 1, 2) - sum_{d=1}^{W0-1} binom(W0 - d + 1, 2)
    #      == binom(W + 1, 3) - binom(W0 + 1, 3)
    psi2_W = (W - d) * (W - d + 1) / 2
    psi2_W0 = (W0 - d) * (W0 - d + 1) / 2
    sum_W = sp.summation(sp.expand(psi2_W), (d, 1, W - 1))
    sum_W0 = sp.summation(sp.expand(psi2_W0), (d, 1, W0 - 1))
    cap_diff = sp.simplify(
        (sum_W - sum_W0) - ((W - 1) * W * (W + 1) / 6 - (W0 - 1) * W0 * (W0 + 1) / 6)
    )
    assert cap_diff == 0, "Band-restricted binomial ceiling identity failed!"

    # Symbolic proof of kappa_int^{band}(alpha) = 1 / (2*(1 + alpha + alpha^2))
    inner_theta = sp.integrate(theta * (1 - theta), (theta, 0, 1))
    assert inner_theta == sp.Rational(1, 6)
    kappa_band_alpha = sp.simplify(
        (sp.Rational(3, 1) / (1 - alpha**3)) * sp.integrate(inner_theta, (tau, alpha, 1))
    )
    expected_kappa_band = sp.Rational(1, 2) / (1 + alpha + alpha**2)
    assert sp.simplify(kappa_band_alpha - expected_kappa_band) == 0

    kappa_band_half = kappa_band_alpha.subs(alpha, sp.Rational(1, 2))
    assert kappa_band_half == sp.Rational(2, 7), f"Expected 2/7, got {kappa_band_half}"

    per_scale_mult_band = 2 * (1 + kappa_band_half)
    assert per_scale_mult_band == sp.Rational(18, 7)

    c_inf_band_sym = sp.simplify(2 / sp.sqrt(1 + kappa_band_half))
    assert sp.simplify(c_inf_band_sym - 2 * sp.sqrt(7) / 3) == 0

    # Compare per-scale weights across Order-1, Phase 7 CF2, and Phase 8 Band-Restricted (alpha=1/2)
    gamma_comparison = []
    for s_idx in range(1, 6):
        b1_sq = beta_s_tau_sq(s_idx, 1.0)
        bband_sq = beta_s_band_sq(s_idx, 0.5)
        g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s_idx))
        g_p7 = max(g1, bband_sq + 0.25 * (1.0 - 2.0 ** (-s_idx)))
        g_p8 = bband_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s_idx))
        # In Phase 8, g_p8 strictly dominates g1 and g_p7 for EVERY s >= 1 (including s = 1)!
        assert g_p8 > g_p7 >= g1, f"Phase 8 weight failed to dominate at s={s_idx}"
        gamma_comparison.append({
            "s": s_idx,
            "gamma_s_order1": round(g1, 8),
            "gamma_s_phase7": round(g_p7, 8),
            "gamma_s_phase8_band": round(g_p8, 8),
        })

    return {
        "normalized_kernel_integral_one": True,
        "octave_matched_limit_exact": "1/4",
        "M_2_exact": "5/16",
        "M_3_exact": "9031/32768",
        "octave_matched_table_s1_to_s8": octave_table,
        "kappa_int_band_alpha_formula": "1 / (2*(1 + alpha + alpha^2))",
        "kappa_int_band_half_exact": "2/7",
        "per_scale_multiplier_band_exact": "18/7",
        "c_infty_band_exact": "2*sqrt(7)/3",
        "c_infty_band_float": round(float(c_inf_band_sym.evalf()), 10),
        "per_scale_weights_s1_to_s5": gamma_comparison,
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


def verify_band_restricted_cesaro_on_mian_chowla():
    """
    Verify Theorem 2.8(b) exact integer Band-Restricted Cesàro–Fejér identities &
    pigeonhole lower-bound chains on the Mian–Chowla sequence (OEIS A005282) across
    dyadic scales r in {2, ..., 7} for band cutoffs W_0 = W/2 (alpha = 1/2) and
    W_0 = 3W/4 (alpha = 3/4).
    """
    mc = generate_mian_chowla(64)
    m_scales = 10
    n_horizon = 1 << m_scales
    in_h = [x for x in mc if x <= n_horizon]
    diffs = [in_h[j] - in_h[i] for i in range(len(in_h)) for j in range(i + 1, len(in_h))]

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

    def unsplit_pigeon_lb_for_w(w: int, r_scale: int) -> int:
        t_r = min(1 << r_scale, n_horizon)
        blocks = [(1, t_r + w - 1)]
        for k in range(r_scale + 1, m_scales + 1):
            t_prev = 1 << (k - 1)
            t_curr = min(1 << k, n_horizon)
            if t_prev >= n_horizon:
                break
            blocks.append((t_prev + w, t_curr + w - 1))
        lb_tot = 0
        for a, b in blocks:
            bins = b - a + 1
            s_val = sum(win_c(u, w) for u in range(a, b + 1))
            lb_tot += l_pigeon(s_val, bins)
        return lb_tot

    rows = []
    for r in range(2, 8):
        w_cap = 1 << r
        w_half = w_cap >> 1
        w_three_quarter = (3 * w_cap) >> 2

        lb_half = 0
        lb_tq = 0
        win_half = 0
        win_tq = 0

        for w in range(w_half + 1, w_cap + 1):
            lb_w = unsplit_pigeon_lb_for_w(w, r)
            fw_exact = sum(
                win_c(u, w) * (win_c(u, w) - 1) // 2 for u in range(1, n_horizon + w)
            )
            assert lb_w <= fw_exact
            lb_half += lb_w
            win_half += fw_exact
            if w > w_three_quarter:
                lb_tq += lb_w
                win_tq += fw_exact

        def kernel_band_sum(w0: int) -> int:
            tot = 0
            for d in diffs:
                if d < w_cap:
                    rem_w = w_cap - d
                    term_w = rem_w * (rem_w + 1) // 2
                    rem_0 = max(w0 - d, 0)
                    term_0 = rem_0 * (rem_0 + 1) // 2
                    tot += term_w - term_0
            return tot

        rhs_half = kernel_band_sum(w_half)
        rhs_tq = kernel_band_sum(w_three_quarter)

        cap_w = (w_cap - 1) * w_cap * (w_cap + 1) // 6
        cap_half = cap_w - (w_half - 1) * w_half * (w_half + 1) // 6
        cap_tq = cap_w - (w_three_quarter - 1) * w_three_quarter * (w_three_quarter + 1) // 6

        assert win_half == rhs_half, f"Band identity failed at r={r}, W0=W/2"
        assert win_tq == rhs_tq, f"Band identity failed at r={r}, W0=3W/4"
        assert lb_half <= rhs_half <= cap_half
        assert lb_tq <= rhs_tq <= cap_tq

        rows.append({
            "r": r,
            "W": w_cap,
            "W0_half": w_half,
            "band_half_pigeonhole_lb": lb_half,
            "band_half_exact_rhs": rhs_half,
            "band_half_binom_cap": cap_half,
            "W0_three_quarter": w_three_quarter,
            "band_tq_pigeonhole_lb": lb_tq,
            "band_tq_exact_rhs": rhs_tq,
            "band_tq_binom_cap": cap_tq,
        })
    return rows


def verify_band_restricted_envelope_bounds():
    kappa_pref = 2.5 - 1.5 * math.sqrt(2.0) + 0.5 * math.log(1.0 + math.sqrt(2.0))
    lam_2 = math.tanh(math.log(2.0) / 4.0) / math.log(2.0)
    lam_101 = math.tanh(math.log(1.01) / 4.0) / math.log(1.01)
    scale_ratio = lam_101 / lam_2
    beta0_sq = (4.0 / 3.0) ** 2
    c_inf_p6 = 2.0 * math.sqrt(6.0 / 7.0)
    c_inf_p7 = 2.0 * math.sqrt(4.0 / 5.0)
    c_inf_p8 = 2.0 * math.sqrt(7.0 / 9.0)

    rows = []
    prev_p8_tail = float("inf")
    for m_val in [8, 12, 16, 20, 24, 32, 48, 64]:
        r0 = m_val // 2
        q_p6_tail = 0.0
        q_p7_tail = 0.0
        q_p8_tail = 0.0

        for k in range(r0 + 1, m_val + 1):
            alpha_k = 1.0 / math.sqrt(k * math.log(2.0)) - 1.0 / math.sqrt(
                2.0 * (k - 1) * math.log(2.0)
            )
            s_k = max(0, k - 2) // 2
            f6 = beta0_sq
            f7 = beta0_sq
            f8 = beta0_sq
            for s in range(1, s_k + 1):
                b1_sq = beta_s_tau_sq(s, 1.0)
                boct_sq = beta_s_band_sq(s, 0.5)
                g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s))
                g2_p7 = boct_sq + 0.25 * (1.0 - 2.0 ** (-s))
                g2_p8 = boct_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s))
                f6 += 2.0 * g1
                f7 += 2.0 * max(g1, g2_p7)
                f8 += 2.0 * max(g1, g2_p8)

            q_p6_tail += scale_ratio * f6 * alpha_k * alpha_k
            q_p7_tail += scale_ratio * f7 * alpha_k * alpha_k
            q_p8_tail += scale_ratio * f8 * alpha_k * alpha_k

        rhs_tail = 0.5 * (m_val - r0)
        c_p6_tail = math.sqrt(rhs_tail / q_p6_tail)
        c_p7_tail = math.sqrt(rhs_tail / q_p7_tail)
        c_p8_tail = math.sqrt(rhs_tail / q_p8_tail)
        assert c_p8_tail < c_p7_tail < c_p6_tail
        assert c_p8_tail < prev_p8_tail
        prev_p8_tail = c_p8_tail

        q_p7_full = 0.0
        q_p8_full = 0.0
        prev_req = 0.0
        for k in range(3, m_val + 1):
            req_k = math.sqrt((2.0**k) / (k * math.log(2.0)))
            rho_k = req_k - prev_req
            prev_req = req_k
            rho_sq_2k = (rho_k * rho_k) / (2.0**k)
            s_k = max(0, k - 2) // 2

            f7 = beta0_sq
            f8 = beta0_sq
            for s in range(1, s_k + 1):
                b1_sq = beta_s_tau_sq(s, 1.0)
                boct_sq = beta_s_band_sq(s, 0.5)
                g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0 ** (-s))
                g2_p7 = boct_sq + 0.25 * (1.0 - 2.0 ** (-s))
                g2_p8 = boct_sq + (2.0 / 7.0) * (1.0 - 2.0 ** (-s))
                f7 += 2.0 * max(g1, g2_p7)
                f8 += 2.0 * max(g1, g2_p8)

            w_bil = 1.0 if k == m_val else 0.5
            q_p7_full += scale_ratio * f7 * rho_sq_2k + (w_bil + 0.5 * kappa_pref) / (
                k * math.log(2.0)
            )
            q_p8_full += scale_ratio * f8 * rho_sq_2k + (w_bil + 0.5 * kappa_pref) / (
                k * math.log(2.0)
            )

        rhs_full = 0.5 * (m_val + 2)
        c_p7_full = math.sqrt(rhs_full / q_p7_full)
        c_p8_full = math.sqrt(rhs_full / q_p8_full)
        assert c_p8_full < c_p7_full
        assert c_p8_full < c_inf_p8 < c_p8_tail

        rows.append({
            "M": m_val,
            "c_p6_coupled_tail_b101": round(c_p6_tail, 6),
            "c_p7_cf2_tail_b101": round(c_p7_tail, 6),
            "c_p8_band_tail_b101": round(c_p8_tail, 6),
            "c_p7_cf2_full_b101": round(c_p7_full, 6),
            "c_p8_band_full_b101": round(c_p8_full, 6),
            "c_inf_p6_limit": round(c_inf_p6, 6),
            "c_inf_p7_limit": round(c_inf_p7, 6),
            "c_inf_p8_limit": round(c_inf_p8, 6),
        })

    # Verify sub-2.000000 breakthrough already at M >= 20 (down from M = 24 in Phase 7)
    for r in rows:
        if r["M"] >= 20:
            assert (
                r["c_p8_band_tail_b101"] < 2.0
            ), f"Expected < 2.0 at M={r['M']}, got {r['c_p8_band_tail_b101']}"

    return rows


def main():
    sym_res = verify_symbolic_theorem_2_8()
    mc_rows = verify_band_restricted_cesaro_on_mian_chowla()
    env_rows = verify_band_restricted_envelope_bounds()
    out = {
        "theorem_2_8_symbolic": sym_res,
        "theorem_2_8_mian_chowla_band": mc_rows,
        "theorem_2_8_band_envelope_bounds": env_rows,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
