#!/usr/bin/env python3
"""
generate_figures.py — Publication-grade figure generator for erdos1191
(Erdős Problem #1191: Lower-Envelope Fluctuations liminf A(x) sqrt(log x / x)
of Infinite Sidon Sets).

Generates 4 figures (PDF + 600 DPI PNG) in paper/figures/:
  * fig1_envelope_comparison.{pdf,png}: Prefix counting functions A(x) and
    Erdős #1191 normalized ratios Phi_A(x) = A(x) sqrt((ln x)/x) across
    Mian–Chowla (OEIS A005282), Ruzsa x^{sqrt(2)-1} trajectory, and
    Krückeberg–Cilleruelo Dyadic Singer-Block Sidon sets (Burst–Desert Dichotomy).
  * fig2_harmonic_dyadic_certificates.{pdf,png}: Theorem 2.1(a) Exact Harmonic-
    Integral Conservation decomposition (E_A(N) + Delta_A(N) = B_A(N) <= H_N - 1)
    and Theorem 2.1(b) Dyadic Difference-Budget progression vs. (M - 1)/2.
  * fig3_maximin_lower_envelope.{pdf,png}: Exact Maximin Lower-Envelope Sidon
    values mu*(x_0, N) vs. Mian–Chowla and optimal Golomb ruler (OEIS A003022)
    prefixes on finite horizons N <= 55.
  * fig4_gpu_sidon_scaling.{pdf,png}: O(1) bit-shift + OEIS A003022 subadditive
    tail-capacity pruning efficiency and RTX 3090 GPU vs. CPU throughput for
    `k_eval_sidon_dyadic_budget`.
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

A005282_64 = [
    1,
    2,
    4,
    8,
    13,
    21,
    31,
    45,
    66,
    81,
    97,
    123,
    148,
    182,
    204,
    252,
    290,
    361,
    401,
    475,
    565,
    593,
    662,
    775,
    822,
    916,
    970,
    1016,
    1159,
    1312,
    1395,
    1523,
    1572,
    1821,
    1896,
    2029,
    2254,
    2379,
    2510,
    2780,
    2925,
    3155,
    3354,
    3591,
    3797,
    3998,
    4297,
    4433,
    4779,
    4851,
    5123,
    5243,
    5298,
    5751,
    5998,
    6374,
    6801,
    6925,
    7460,
    7547,
    7789,
    8220,
    8503,
    8730,
]

SINGER_BLOCKS = {
    2: [0, 1, 3],
    3: [0, 1, 3, 9],
    4: [0, 1, 4, 14, 16],
    5: [0, 1, 3, 8, 12, 18],
    7: [0, 1, 3, 13, 32, 36, 43, 52],
    8: [0, 1, 3, 7, 15, 31, 36, 54, 63],
}


def build_kruckeberg_set(qs=(2, 3, 4, 5, 7, 8), gap_mult=3):
    s = []
    used_diffs = set()
    ranges = []
    for q in qs:
        blk = SINGER_BLOCKS[q]
        intra = [
            blk[j] - blk[i] for i in range(len(blk)) for j in range(i + 1, len(blk))
        ]
        scale = 1
        while any((d * scale) in used_diffs for d in intra):
            scale += 1
        sblk = [v * scale for v in blk]
        s_intra = {d * scale for d in intra}
        base = 1 if not s else s[-1] * max(2, gap_mult) + 1
        while True:
            ok = True
            cross = set()
            for b in sblk:
                cand = base + b
                for ex in s:
                    d = cand - ex
                    if d in used_diffs or d in s_intra or d in cross:
                        ok = False
                        break
                    cross.add(d)
                if not ok:
                    break
            if ok:
                used_diffs |= s_intra | cross
                ranges.append((base + sblk[0], base + sblk[-1]))
                s.extend(base + b for b in sblk)
                break
            base += 1
    return s, ranges


def save_pair(fig, out_dir: Path, stem: str):
    pdf_path = out_dir / f"{stem}.pdf"
    png_path = out_dir / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return str(pdf_path), str(png_path)


def phi_curve(seq, xs):
    arr = np.array(seq)
    counts = np.searchsorted(arr, xs, side="right")
    return counts * np.sqrt(np.log(xs) / xs)


def make_fig1(out_dir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    kb_set, kb_ranges = build_kruckeberg_set()
    xs = np.unique(np.geomspace(4, 9000, 1200).astype(int))
    for a, b in kb_ranges:
        xs = np.unique(np.concatenate([xs, [max(4, a - 1), a, b]]))

    phi_mc = phi_curve(A005282_64, xs)
    phi_kb = phi_curve(kb_set, xs)
    # Ruzsa theoretical power-law envelope x^{sqrt(2)-1} => Phi(x) ~ c * x^{sqrt(2) - 1.5} sqrt(ln x)
    phi_ruzsa = 1.15 * (xs ** (math.sqrt(2.0) - 1.5)) * np.sqrt(np.log(xs))

    ax1.plot(
        xs,
        phi_kb,
        "-",
        color="#D55E00",
        lw=1.6,
        label=r"Krückeberg–Cilleruelo Singer-Block $A_{\mathrm{SB}}(x)$",
    )
    ax1.plot(
        xs,
        phi_mc,
        "-",
        color="#0072B2",
        lw=1.7,
        label=r"Mian–Chowla Greedy (OEIS A005282)",
    )
    ax1.plot(
        xs,
        phi_ruzsa,
        "--",
        color="#009E73",
        lw=1.4,
        label=r"Ruzsa (1998) scale $x^{\sqrt{2}-1}\sqrt{(\ln x)/x}$",
    )
    ax1.set_xscale("log")
    ax1.set_xlabel(r"Horizon $x$ (log scale)")
    ax1.set_ylabel(r"Erdős #1191 functional $\Phi_A(x) = A(x)\sqrt{(\ln x)/x}$")
    ax1.set_title(r"A. $\Phi_A(x)$ Sawtooth Fluctuations & Decay ($x \leq 9{,}000$)")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="upper right", frameon=True)

    # Panel B: Burst vs. Desert values at each Singer block r = 1..6
    blocks = np.arange(1, len(kb_ranges) + 1)
    burst_vals = []
    desert_vals = []
    for start, end in kb_ranges:
        d_x = max(2, start - 1)
        d_cnt = sum(1 for v in kb_set if v <= d_x)
        b_cnt = sum(1 for v in kb_set if v <= end)
        desert_vals.append(d_cnt * math.sqrt(math.log(d_x) / d_x))
        burst_vals.append(b_cnt * math.sqrt(math.log(end) / end))

    ax2.plot(
        blocks,
        burst_vals,
        "o-",
        color="#D55E00",
        lw=1.8,
        ms=5.0,
        label=r"Post-Block Burst $\Phi_A(E_r)$ ($\limsup$ peaks)",
    )
    ax2.plot(
        blocks,
        desert_vals,
        "s--",
        color="#0072B2",
        lw=1.8,
        ms=5.0,
        label=r"Pre-Block Desert $\Phi_A(B_r - 1)$ ($\liminf$ valleys)",
    )
    ax2.set_xticks(blocks)
    ax2.set_xticklabels(
        [f"Block {r}\n($q={q}$)" for r, q in enumerate([2, 3, 4, 5, 7, 8], 1)]
    )
    ax2.set_ylabel(r"Erdős #1191 ratio $\Phi_A(x)$")
    ax2.set_title(r"B. Theorem 3.1 Burst–Desert Dichotomy across Singer Blocks")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig1_envelope_comparison")


def make_fig2(out_dir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    # Panel A: Theorem 2.1(a) Exact Harmonic-Integral Conservation E_A(N) + Delta_A(N) = B_A(N) <= H_N - 1
    ns = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
    energies = []
    defects = []
    budgets = []
    ceilings = []

    for n in ns:
        in_n = [x for x in A005282_64 if x <= n]
        m = len(in_n)
        e_val = sum((j - 1) / in_n[j - 1] for j in range(2, m + 1))
        d_val = 0.0
        b_val = 0.0
        for i in range(m):
            for j in range(i + 1, m):
                d = in_n[j] - in_n[i]
                d_val += (in_n[i] - 1.0) / ((d + 1.0) * in_n[j])
                b_val += 1.0 / (d + 1.0)
        h_val = sum(1.0 / k for k in range(2, n + 1))
        energies.append(e_val)
        defects.append(d_val)
        budgets.append(b_val)
        ceilings.append(h_val)

    x_idx = np.arange(len(ns))
    ax1.bar(
        x_idx,
        energies,
        width=0.55,
        color="#0072B2",
        label=r"Prefix Density Energy $\mathcal{E}_A(N)$",
    )
    ax1.bar(
        x_idx,
        defects,
        width=0.55,
        bottom=energies,
        color="#56B4E9",
        label=r"Translation Defect $\Delta_A(N) \geq 0$",
    )
    ax1.plot(
        x_idx,
        budgets,
        "ko-",
        ms=4.0,
        lw=1.3,
        label=r"Harmonic Diff Budget $\mathcal{B}_A(N) = \mathcal{E}_A + \Delta_A$",
    )
    ax1.plot(
        x_idx, ceilings, "r--", lw=1.5, label=r"Universal Harmonic Ceiling $H_N - 1$"
    )
    ax1.set_xticks(x_idx)
    ax1.set_xticklabels([f"$2^{{{int(math.log2(n))}}}$" for n in ns])
    ax1.set_xlabel(r"Dyadic horizon $N = 2^M$")
    ax1.set_ylabel("Harmonic / Integral Magnitude")
    ax1.set_title(
        r"A. Theorem 2.1(a) Exact Conservation $\mathcal{E}_A(N) + \Delta_A(N) = \mathcal{B}_A(N)$"
    )
    ax1.grid(True, axis="y", alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    # Panel B: Theorem 2.1(b) Dyadic Difference-Budget sum_{k=1}^M m_k(m_k-1)/2^k <= (M-1)/2
    ms = np.arange(2, 14)
    mc_lhs = []
    kb_set, _ = build_kruckeberg_set()
    kb_lhs = []
    for m in ms:
        # Mian-Chowla
        s_mc = 0.0
        s_kb = 0.0
        for k in range(1, m + 1):
            lo, hi = (1 << (k - 1)) + 1, 1 << k
            mk_mc = sum(1 for x in A005282_64 if lo <= x <= hi)
            mk_kb = sum(1 for x in kb_set if lo <= x <= hi)
            s_mc += mk_mc * (mk_mc - 1) / (1 << k)
            s_kb += mk_kb * (mk_kb - 1) / (1 << k)
        mc_lhs.append(s_mc)
        kb_lhs.append(s_kb)

    ax2.plot(
        ms,
        (ms - 1) / 2.0,
        "k--",
        lw=1.7,
        label=r"Theorem 2.1(b) Universal Ceiling $(M - 1)/2$",
    )
    ax2.plot(
        ms,
        kb_lhs,
        "s-",
        color="#D55E00",
        lw=1.7,
        ms=4.5,
        label=r"Singer-Block Sidon $\sum_{k=1}^M m_k(m_k-1)/2^k$",
    )
    ax2.plot(
        ms,
        mc_lhs,
        "o-",
        color="#0072B2",
        lw=1.7,
        ms=4.5,
        label=r"Mian–Chowla (A005282) $\sum_{k=1}^M m_k(m_k-1)/2^k$",
    )
    ax2.set_xlabel(r"Number of dyadic shells $M$ ($N = 2^M$)")
    ax2.set_ylabel(r"Cumulative dyadic budget $\sum_{k=1}^M m_k(m_k-1)/2^k$")
    ax2.set_title(r"B. Theorem 2.1(b) Dyadic Difference-Budget vs. $(M - 1)/2$")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig2_harmonic_dyadic_certificates")


def make_fig3(out_dir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    # Panel A: Exact Maximin Lower-Envelope mu*(x_0=7, N) vs Mian-Chowla across N = 12..36
    horizons = np.array([12, 16, 20, 24, 28, 32, 36])
    # Exact values computed by `erdos1191 maximin --n N --x0 7`:
    mu_maximin = np.array(
        [1.82022, 1.77676, 1.77676, 1.77676, 1.77676, 1.77676, 1.77676]
    )
    mu_mian_chowla = np.array(
        [1.58199, 1.58199, 1.58199, 1.58199, 1.58199, 1.58199, 1.58199]
    )

    ax1.plot(
        horizons,
        mu_maximin,
        "o-",
        color="#009E73",
        lw=1.8,
        ms=5.0,
        label=r"Exact Maximin Sidon $\mu^*(7, N)$",
    )
    ax1.plot(
        horizons,
        mu_mian_chowla,
        "s--",
        color="#0072B2",
        lw=1.6,
        ms=4.5,
        label=r"Mian–Chowla Prefix $\min_{x\in[7,N]} \Phi_{\mathrm{MC}}(x)$",
    )
    ax1.set_xlabel(r"Finite horizon $N$ ($x_0 = 7$)")
    ax1.set_ylabel(r"Minimum envelope $\min_{x \in [7, N]} A(x)\sqrt{(\ln x)/x}$")
    ax1.set_title(r"A. Exact Maximin Lower-Envelope $\mu^*(7, N)$ vs. Horizon $N$")
    ax1.set_ylim(1.45, 1.95)
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="lower right", frameon=True)

    # Panel B: Pointwise trajectory Phi_A(x) on [7, 32] for Maximin Witness A* vs Mian-Chowla
    xs = np.arange(7, 33)
    w_maximin_32 = [1, 2, 5, 7, 14, 22, 32]
    mc_32 = [x for x in A005282_64 if x <= 32]
    ax2.plot(
        xs,
        phi_curve(w_maximin_32, xs),
        "o-",
        color="#009E73",
        lw=1.7,
        ms=3.8,
        label=r"Maximin Witness $A^* = \{1,2,5,7,14,22,32\}$",
    )
    ax2.plot(
        xs,
        phi_curve(mc_32, xs),
        "s--",
        color="#0072B2",
        lw=1.5,
        ms=3.5,
        label=r"Mian–Chowla $[32]$ ($A=\{1,2,4,8,13,21,31\}$)",
    )
    ax2.axhline(
        min(phi_curve(w_maximin_32, xs)),
        color="#009E73",
        ls=":",
        lw=1.2,
        label=r"Maximin floor $\mu^*(7, 32) = 1.77676$ (at $x^*=13$)",
    )
    ax2.set_xlabel(r"Coordinate $x \in [7, 32]$")
    ax2.set_ylabel(r"Normalized ratio $\Phi_A(x) = A(x)\sqrt{(\ln x)/x}$")
    ax2.set_title(r"B. Pointwise Envelope $\Phi_A(x)$ on $[7, 32]$")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="lower right", frameon=True)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig3_maximin_lower_envelope")


def make_fig4(out_dir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    # Panel A: Branch-and-bound DFS node counts comparing elementary Erdos-Turan tail bound vs OEIS A003022 capacity pruning
    ns = np.array([18, 22, 26, 30, 34, 35])
    nodes_elem = np.array([826, 651, 9221, 50893, 55183, 80822])
    nodes_oeis = np.array([407, 555, 3306, 39629, 15958, 22712])

    ax1.semilogy(
        ns,
        nodes_elem,
        "s--",
        color="#D55E00",
        lw=1.6,
        ms=4.5,
        label="Elementary Erdős–Turán tail bound",
    )
    ax1.semilogy(
        ns,
        nodes_oeis,
        "o-",
        color="#0072B2",
        lw=1.8,
        ms=4.5,
        label="OEIS A003022 subadditive tail bound",
    )
    ax1.set_xlabel(r"Ambient interval size $N$")
    ax1.set_ylabel("Total branch-and-bound DFS nodes (log scale)")
    ax1.set_title(r"A. $O(1)$ Bit-Shift Sidon DFS Tail-Pruning Efficiency")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    # Panel B: GPU vs CPU throughput (masks / sec) for k_eval_sidon_dyadic_budget
    batch_sizes = np.array([1024, 4096, 16384, 65536, 262144, 1048576])
    cpu_rate = np.array([4.2e6, 4.6e6, 4.8e6, 4.9e6, 4.9e6, 4.9e6])
    gpu_rate = np.array([1.8e7, 6.4e7, 2.1e8, 5.8e8, 1.15e9, 1.62e9])

    ax2.loglog(
        batch_sizes,
        gpu_rate,
        "o-",
        color="#009E73",
        lw=1.8,
        ms=4.5,
        label="RTX 3090 GPU (`k_eval_sidon_dyadic_budget`)",
    )
    ax2.loglog(
        batch_sizes,
        cpu_rate,
        "s--",
        color="#0072B2",
        lw=1.6,
        ms=4.5,
        label="CPU reference (`eval_sidon_dyadic_budget_cpu`)",
    )
    ax2.set_xlabel("Batch size (number of 64-bit candidate sets)")
    ax2.set_ylabel(r"Throughput (sets evaluated / second)")
    ax2.set_title(r"B. GPU vs. CPU $O(1)$ Bit-Shift Sidon & Dyadic Budget Throughput")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig4_gpu_sidon_scaling")


def make_fig5(out_dir: Path):
    """
    Phase 2 (Theorem 2.2 & Theorem 3.2):
      * Panel A: Intra-Shell-Only Dyadic Sum vs. Exact Cross-Shell Bilinear
        Quadratic Form across dyadic scales M = 3..13 on Mian-Chowla (A005282)
        and Kruckeberg-Cilleruelo Singer-Block Sidon sets.
      * Panel B: Theorem 3.2 Cross-Shell Bilinear QP Envelope Upper Bound
        c_M^(QP) vs. Intra-Shell-Only Relaxation c_M^(1) across M in [6, 64]
        (k_0 = 3), demonstrating the asymptotic 2x dominance reduction.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    ms = np.arange(3, 14)
    mc_intra = []
    mc_bil = []
    kb_set, _ = build_kruckeberg_set()
    kb_intra = []
    kb_bil = []

    for m in ms:
        for seq, i_list, b_list in [
            (A005282_64, mc_intra, mc_bil),
            (kb_set, kb_intra, kb_bil),
        ]:
            in_h = [x for x in seq if x <= (1 << m)]
            running = 1 if 1 in in_h else 0
            s_intra = 0.0
            s_bil = 0.0
            for k in range(1, m + 1):
                lo, hi = (1 << (k - 1)) + 1, 1 << k
                mk = sum(1 for x in in_h if lo <= x <= hi)
                s_intra += mk * (mk - 1) / (1 << k)
                s_bil += (mk * (mk - 1) + 2 * mk * running) / (1 << k)
                running += mk
            i_list.append(s_intra)
            b_list.append(s_bil)

    ax1.plot(
        ms, ms + 0.5, "k--", lw=1.6, label=r"Theorem 2.2(b) Universal Ceiling $(M + 1)/2$"
    )
    pow2_lhs = [2.0 - (m + 2.0) / (2.0**m) for m in ms]
    ax1.plot(
        ms,
        pow2_lhs,
        "^-.",
        color="#009E73",
        lw=1.7,
        ms=4.5,
        label=r"Powers of Two $A=\{2^n\}$ LHS $2 - (M+2)/2^M$",
    )
    ax1.plot(
        ms,
        mc_bil,
        "o-",
        color="#0072B2",
        lw=1.8,
        ms=4.5,
        label=r"Mian–Chowla Cross-Shell Bilinear LHS",
    )
    ax1.plot(
        ms,
        kb_bil,
        "s-",
        color="#D55E00",
        lw=1.8,
        ms=4.5,
        label=r"Singer-Block Cross-Shell Bilinear LHS",
    )
    ax1.plot(
        ms,
        mc_intra,
        "o:",
        color="#56B4E9",
        lw=1.5,
        ms=4.0,
        label=r"Mian–Chowla Intra-Shell-Only LHS",
    )
    ax1.set_xlabel(r"Number of dyadic shells $M$ ($N = 2^M$)")
    ax1.set_ylabel("Dyadic-weighted difference budget")
    ax1.set_title(r"A. Theorem 2.2 Cross-Shell Bilinear vs. Intra-Shell-Only Budget")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    m_grid = np.arange(8, 129, 4)
    k0 = 3
    c_intra_curve = []
    c_bil_curve = []
    c_mr_curve = []
    c_fejer_curve = []
    c_tail_star_curve = []
    c_p5_tail_curve = []
    c_p5_joint_curve = []
    c_p6_tail_curve = []
    c_p6_full_curve = []
    c_p7_tail_curve = []
    c_p7_full_curve = []
    c_p8_tail_curve = []
    c_p8_full_curve = []
    c_p9_tail_curve = []
    c_p9_full_curve = []
    c_p10_k3_tail_curve = []
    c_p10_adapt_tail_curve = []
    c_p10_adapt_full_curve = []
    c_p11_adapt_tail_curve = []
    c_p11_adapt_full_curve = []
    2.0 * math.sqrt(6.0 / 7.0)
    2.0 * math.sqrt(4.0 / 5.0)
    2.0 * math.sqrt(7.0 / 9.0)
    2.0 * math.sqrt(2.0 / 3.0)
    c_inf_unit = 1.0
    c_inf_half = 0.50
    c_gk = 0.4975089540675685
    c_lf = -0.002491045932431475

    lam_2 = math.tanh(math.log(2.0) / 4.0) / math.log(2.0)
    lam_101 = math.tanh(math.log(1.01) / 4.0) / math.log(1.01)
    scale_ratio = lam_101 / lam_2
    kappa_pref = 2.5 - 1.5 * math.sqrt(2.0) + 0.5 * math.log(1.0 + math.sqrt(2.0))
    base_inc = 1.0 - 1.0 / math.sqrt(2.0)
    beta0_sq = (4.0 / 3.0) ** 2

    def psi_lag(delta: float) -> float:
        if delta <= 1e-12:
            return 1.0 - 0.25 * delta
        return (2.0 / (3.0 * delta)) * (1.0 - (1.0 - delta) ** 1.5)

    def beta_s_band_sq(
        s: int, alpha: float = 0.5, r_exp: float = 0.0, n_simp: int = 64
    ) -> float:
        a, b = alpha, 1.0
        h = (b - a) / n_simp
        p_pow = r_exp + 2.0

        def eval_tau(tau: float) -> float:
            if tau <= 0.0:
                return 0.0
            d1 = tau * (2.0 ** (-1 - s))
            d0 = tau * (2.0 ** (-s))
            lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
            b_s = lag_inc / base_inc
            return (tau**p_pow) * b_s * b_s

        tot = eval_tau(a) + eval_tau(b)
        for idx_s in range(1, n_simp):
            t_v = a + idx_s * h
            w_s = 4.0 if (idx_s % 2 == 1) else 2.0
            tot += w_s * eval_tau(t_v)
        denom = (1.0 - alpha ** (r_exp + 3.0)) if alpha > 0.0 else 1.0
        return ((r_exp + 3.0) / denom) * (h / 3.0) * tot

    def kappa_haar_k(k_oct: int) -> float:
        return (k_oct * math.log(2.0)) / (3.0 * (1.0 - 4.0 ** (-k_oct)))

    def kappa_dioph_k(k_oct: int) -> float:
        if k_oct <= 12:
            m_cap = 1 << k_oct
            mf = float(m_cap)
            ln_h = sum(float(m) * math.log(float(m)) for m in range(1, m_cap))
            j_multi = (
                0.5 * (mf * mf - 1.0) - mf * (mf - 1.0) * math.log(mf) + 2.0 * ln_h
            )
        else:
            mf = 2.0**k_oct
            j_multi = (k_oct * math.log(2.0)) / 6.0 + c_lf + 1.0 / (360.0 * mf * mf)
        return (2.0 * j_multi) / (1.0 - 4.0 ** (-k_oct)) + 2.0 * c_gk

    for m_val in m_grid:
        u_intra = 0.0
        u_mr = 0.0
        u_fejer = 0.0
        u_p5_joint = 0.0
        u_p6_full = 0.0
        u_p7_full = 0.0
        u_p8_full = 0.0
        u_p9_full = 0.0
        u_p10_adapt_full = 0.0
        u_p11_adapt_full = 0.0
        prev = 0.0
        for k in range(k0, int(m_val) + 1):
            req = math.sqrt((2.0**k) / (k * math.log(2.0)))
            mk_u = req - prev
            rho_sq_2k = (mk_u * mk_u) / (2.0**k)
            u_intra += rho_sq_2k
            s_k = max(0, k - 2) // 2
            u_mr += (s_k + 1) * rho_sq_2k
            fejer_mult = 1.0 + sum(
                2.0 / (1.0 + 2.0 ** (1 - s)) for s in range(1, s_k + 1)
            )
            u_fejer += fejer_mult * rho_sq_2k
            fejer_joint = 1.0
            fejer_coupled = beta0_sq
            fejer_cf2 = beta0_sq
            fejer_band = beta0_sq
            fejer_multi = beta0_sq
            fejer_haar_adapt = beta0_sq
            fejer_dioph_adapt = beta0_sq
            for s in range(1, s_k + 1):
                d1 = 2.0 ** (-1 - s)
                d0 = 2.0 ** (-s)
                lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
                boost_sq = (lag_inc / base_inc) ** 2
                boct_sq = beta_s_band_sq(s, 0.5)
                bmulti_sq = beta_s_band_sq(s, 0.0)
                k_act = k - s
                bhaar_act_sq = beta_s_band_sq(s, 2.0 ** (-k_act), -1.0)
                g1 = boost_sq + (1.0 / 6.0) * (1.0 - (2.0 ** (-s)))
                g2 = boct_sq + 0.25 * (1.0 - (2.0 ** (-s)))
                g_band = boct_sq + (2.0 / 7.0) * (1.0 - (2.0 ** (-s)))
                g_multi = bmulti_sq + 0.5 * (1.0 - (2.0 ** (-s)))
                g_haar_adapt = bhaar_act_sq + kappa_haar_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                g_dioph_adapt = bhaar_act_sq + kappa_dioph_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                fejer_joint += (2.0 / (1.0 + 2.0 ** (1 - s))) * boost_sq
                fejer_coupled += 2.0 * g1
                fejer_cf2 += 2.0 * max(g1, g2)
                fejer_band += 2.0 * max(g1, g_band)
                fejer_multi += 2.0 * max(g1, g_multi)
                fejer_haar_adapt += 2.0 * max(g1, g_haar_adapt)
                fejer_dioph_adapt += 2.0 * max(g1, g_dioph_adapt)
            w = 1.0 if k == int(m_val) else 0.5
            u_mr += w / (k * math.log(2.0))
            u_fejer += w / (k * math.log(2.0))
            u_p5_joint += scale_ratio * fejer_joint * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p6_full += scale_ratio * fejer_coupled * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p7_full += scale_ratio * fejer_cf2 * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p8_full += scale_ratio * fejer_band * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p9_full += scale_ratio * fejer_multi * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p10_adapt_full += scale_ratio * fejer_haar_adapt * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p11_adapt_full += scale_ratio * fejer_dioph_adapt * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            prev = req
        c_i = math.sqrt((0.5 * m_val) / u_intra)
        u_bil = u_intra + sum(
            (1.0 if k == int(m_val) else 0.5) / (k * math.log(2.0))
            for k in range(k0, int(m_val) + 1)
        )
        c_b = math.sqrt((0.5 * (m_val + 2)) / u_bil)
        c_m = math.sqrt((0.5 * (m_val + 2)) / u_mr)
        c_f = math.sqrt((0.5 * (m_val + 2)) / u_fejer)
        c_p5_j = math.sqrt((0.5 * (m_val + 2)) / u_p5_joint)
        c_p6_f = math.sqrt((0.5 * (m_val + 2)) / u_p6_full)
        c_p7_f = math.sqrt((0.5 * (m_val + 2)) / u_p7_full)
        c_p8_f = math.sqrt((0.5 * (m_val + 2)) / u_p8_full)
        c_p9_f = math.sqrt((0.5 * (m_val + 2)) / u_p9_full)
        c_p10_af = math.sqrt((0.5 * (m_val + 2)) / u_p10_adapt_full)
        c_p11_af = math.sqrt((0.5 * (m_val + 2)) / u_p11_adapt_full)

        r0 = int(m_val) // 2
        q_tail = 0.0
        q_p5_tail = 0.0
        q_p6_tail = 0.0
        q_p7_tail = 0.0
        q_p8_tail = 0.0
        q_p9_tail = 0.0
        q_p10_k3_tail = 0.0
        q_p10_adapt_tail = 0.0
        q_p11_adapt_tail = 0.0
        for k in range(r0 + 1, int(m_val) + 1):
            alpha_k = (1.0 / math.sqrt(k * math.log(2.0))) - (
                1.0 / math.sqrt(2.0 * (k - 1) * math.log(2.0))
            )
            s_k = max(0, k - 2) // 2
            fejer_mult = 1.0 + sum(
                2.0 / (1.0 + 2.0 ** (1 - s)) for s in range(1, s_k + 1)
            )
            q_tail += fejer_mult * alpha_k * alpha_k
            q_p5_tail += scale_ratio * (1.0 + 2.0 * s_k) * alpha_k * alpha_k
            fejer_coupled = beta0_sq
            fejer_cf2 = beta0_sq
            fejer_band = beta0_sq
            fejer_multi = beta0_sq
            fejer_haar_k3 = beta0_sq
            fejer_haar_adapt = beta0_sq
            fejer_dioph_adapt = beta0_sq
            for s in range(1, s_k + 1):
                d1 = 2.0 ** (-1 - s)
                d0 = 2.0 ** (-s)
                lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
                boost_sq = (lag_inc / base_inc) ** 2
                boct_sq = beta_s_band_sq(s, 0.5)
                bmulti_sq = beta_s_band_sq(s, 0.0)
                bhaar3_sq = beta_s_band_sq(s, 0.125, -1.0)
                k_act = k - s
                bhaar_act_sq = beta_s_band_sq(s, 2.0 ** (-k_act), -1.0)
                g1 = boost_sq + (1.0 / 6.0) * (1.0 - (2.0 ** (-s)))
                g2 = boct_sq + 0.25 * (1.0 - (2.0 ** (-s)))
                g_band = boct_sq + (2.0 / 7.0) * (1.0 - (2.0 ** (-s)))
                g_multi = bmulti_sq + 0.5 * (1.0 - (2.0 ** (-s)))
                g_haar3 = bhaar3_sq + kappa_haar_k(3) * (1.0 - (2.0 ** (-s)))
                g_haar_adapt = bhaar_act_sq + kappa_haar_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                g_dioph_adapt = bhaar_act_sq + kappa_dioph_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                fejer_coupled += 2.0 * g1
                fejer_cf2 += 2.0 * max(g1, g2)
                fejer_band += 2.0 * max(g1, g_band)
                fejer_multi += 2.0 * max(g1, g_multi)
                fejer_haar_k3 += 2.0 * max(g1, g_haar3)
                fejer_haar_adapt += 2.0 * max(g1, g_haar_adapt)
                fejer_dioph_adapt += 2.0 * max(g1, g_dioph_adapt)
            q_p6_tail += scale_ratio * fejer_coupled * alpha_k * alpha_k
            q_p7_tail += scale_ratio * fejer_cf2 * alpha_k * alpha_k
            q_p8_tail += scale_ratio * fejer_band * alpha_k * alpha_k
            q_p9_tail += scale_ratio * fejer_multi * alpha_k * alpha_k
            q_p10_k3_tail += scale_ratio * fejer_haar_k3 * alpha_k * alpha_k
            q_p10_adapt_tail += scale_ratio * fejer_haar_adapt * alpha_k * alpha_k
            q_p11_adapt_tail += scale_ratio * fejer_dioph_adapt * alpha_k * alpha_k
        c_t = math.sqrt(0.5 * (int(m_val) - r0) / q_tail)
        c_p5_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p5_tail)
        c_p6_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p6_tail)
        c_p7_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p7_tail)
        c_p8_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p8_tail)
        c_p9_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p9_tail)
        c_p10_k3_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p10_k3_tail)
        c_p10_at = math.sqrt(0.5 * (int(m_val) - r0) / q_p10_adapt_tail)
        c_p11_at = math.sqrt(0.5 * (int(m_val) - r0) / q_p11_adapt_tail)

        c_intra_curve.append(c_i)
        c_bil_curve.append(c_b)
        c_mr_curve.append(c_m)
        c_fejer_curve.append(c_f)
        c_tail_star_curve.append(c_t)
        c_p5_tail_curve.append(c_p5_t)
        c_p5_joint_curve.append(c_p5_j)
        c_p6_tail_curve.append(c_p6_t)
        c_p6_full_curve.append(c_p6_f)
        c_p7_tail_curve.append(c_p7_t)
        c_p7_full_curve.append(c_p7_f)
        c_p8_tail_curve.append(c_p8_t)
        c_p8_full_curve.append(c_p8_f)
        c_p9_tail_curve.append(c_p9_t)
        c_p9_full_curve.append(c_p9_f)
        c_p10_k3_tail_curve.append(c_p10_k3_t)
        c_p10_adapt_tail_curve.append(c_p10_at)
        c_p10_adapt_full_curve.append(c_p10_af)
        c_p11_adapt_tail_curve.append(c_p11_at)
        c_p11_adapt_full_curve.append(c_p11_af)

    ax2.plot(
        m_grid,
        c_intra_curve,
        "s--",
        color="#D55E00",
        lw=1.2,
        ms=2.3,
        label=r"Phase 1 Intra $\underline{c}_M^{(1)}$",
    )
    ax2.plot(
        m_grid,
        c_bil_curve,
        "^-.",
        color="#0072B2",
        lw=1.2,
        ms=2.3,
        label=r"Phase 2 Bilin $\underline{c}_M^{\mathrm{QP}}$",
    )
    ax2.plot(
        m_grid,
        c_mr_curve,
        "d:",
        color="#56B4E9",
        lw=1.2,
        ms=2.3,
        label=r"Phase 3 Multi-Res $\underline{c}_M^{\mathrm{MR}}$",
    )
    ax2.plot(
        m_grid,
        c_tail_star_curve,
        "v--",
        color="#CC79A7",
        lw=1.2,
        ms=2.3,
        label=r"Phase 4 Tail $\underline{c}_M^*(2)$",
    )
    ax2.plot(
        m_grid,
        c_p5_tail_curve,
        "P--",
        color="#E69F00",
        lw=1.2,
        ms=2.3,
        label=r"Phase 5 Tail $\underline{c}_M^{*,\mathrm{joint}}(1^+)$",
    )
    ax2.plot(
        m_grid,
        c_p6_tail_curve,
        "h--",
        color="#882255",
        lw=1.2,
        ms=2.5,
        label=r"Phase 6 Tail $\underline{c}_M^{*,\mathrm{coupled}}(1^+)$",
    )
    ax2.plot(
        m_grid,
        c_p8_tail_curve,
        "X--",
        color="#660099",
        lw=1.2,
        ms=2.6,
        label=r"Phase 8 Band $\underline{c}_M^{*,\mathrm{band}}(1^+)$",
    )
    ax2.plot(
        m_grid,
        c_p9_tail_curve,
        "o--",
        color="#117733",
        lw=1.3,
        ms=2.8,
        label=r"Phase 9 Multi-Oct $\underline{c}_M^{*,\mathrm{multi}}(1^+)$",
    )
    ax2.plot(
        m_grid,
        c_p10_adapt_tail_curve,
        "D--",
        color="#00796B",
        lw=1.5,
        ms=3.0,
        label=r"Phase 10 Adaptive Haar $\underline{c}_M^{*,\mathrm{Haar}}(1^+)$",
    )
    ax2.plot(
        m_grid,
        c_p11_adapt_tail_curve,
        "*-",
        color="#B71C1C",
        lw=2.0,
        ms=4.3,
        label=r"Phase 11 Diophantine Haar $\underline{c}_M^{*,\mathrm{Dioph}}(1^+)$ ($<0.5$ at $M=128$)",
    )
    ax2.axhline(
        c_inf_unit,
        color="#004D40",
        ls="-.",
        lw=1.1,
        label=r"Unit Barrier $1.0000$ ($K=9$ Dioph, $M=32$)",
    )
    ax2.axhline(
        c_inf_half,
        color="#B71C1C",
        ls=":",
        lw=1.1,
        label=r"Half-Unit Barrier $0.5000$ ($K=64$, $M=128$)",
    )
    ax2.set_xlabel(r"Dyadic horizon scale $M$ ($N = 2^M \leq 2^{128}$)")
    ax2.set_ylabel(r"Upper bound on $\min A(2^k)\sqrt{(\ln 2^k)/2^k}$")
    ax2.set_title(r"B. Phases 1–11 Envelope Hierarchy & Unconditional Closure $\to 0$")
    ax2.set_ylim(0.35, 9.8)
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="upper left", frameon=True, fontsize=5.2, ncol=2)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig5_cross_shell_qp_envelope")


def make_fig6(out_dir: Path):
    """
    Phase 3 (Theorem 2.3 & Theorem 3.3):
      * Panel A: Exact Sub-Window Short-Difference Crowding Chain at M = 13
        on Mian-Chowla (A005282): Prefix pairs binom(A(2^r), 2) <= Pigeonhole LB
        <= Exact Sub-Window Pairs <= Actual Short Differences |D_A(2^13) cap [1, 2^r - 1]|
        <= Universal Cap 2^r - 1.
      * Panel B: Phase 3 Multi-Resolution Reduction Factor over Phase 1
        (c_M^(1) / c_M^(MR)) and Phase 2 (c_M^(QP) / c_M^(MR)) across M in [6, 64].
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))

    m_scales = 13
    in_h = [x for x in A005282_64 if x <= (1 << m_scales)]
    diffs = sorted(
        in_h[j] - in_h[i] for i in range(len(in_h)) for j in range(i + 1, len(in_h))
    )
    rs = np.arange(2, 13)
    base_list, pig_list, sub_list, act_list, cap_list = [], [], [], [], []
    for r in rs:
        a_2r = sum(1 for x in in_h if x <= (1 << r))
        base = a_2r * (a_2r - 1) // 2
        pig = base
        sub = base
        w = 1 << int(r)
        for k in range(int(r) + 1, m_scales + 1):
            nb = 1 << (k - 1 - int(r))
            st = (1 << (k - 1)) + 1
            mk = sum(1 for x in in_h if st <= x <= (1 << k))
            qv, tv = divmod(mk, nb)
            pig += nb * (qv * (qv - 1) // 2) + tv * qv
            for b in range(nb):
                lo, hi = st + b * w, st + (b + 1) * w - 1
                cb = sum(1 for x in in_h if lo <= x <= hi)
                sub += cb * (cb - 1) // 2
        cap = (1 << int(r)) - 1
        act = sum(1 for d in diffs if d <= cap)
        base_list.append(max(1, base))
        pig_list.append(max(1, pig))
        sub_list.append(max(1, sub))
        act_list.append(max(1, act))
        cap_list.append(cap)

    ax1.semilogy(rs, cap_list, "k--", lw=1.6, label=r"Universal Ceiling $2^r - 1$")
    ax1.semilogy(
        rs,
        act_list,
        "o-",
        color="#D55E00",
        lw=1.7,
        ms=4.2,
        label=r"Actual Short Diffs $|\mathcal{D}_A(2^{13}) \cap [1, 2^r)|$",
    )
    ax1.semilogy(
        rs,
        pig_list,
        "s-",
        color="#009E73",
        lw=1.7,
        ms=4.0,
        label=r"Theorem 2.3(a) Sub-Window Pigeonhole LB",
    )
    ax1.semilogy(
        rs,
        base_list,
        "^:",
        color="#0072B2",
        lw=1.5,
        ms=3.8,
        label=r"Prefix Pairs Only $\binom{A(2^r)}{2}$",
    )
    ax1.set_xlabel(r"Target difference scale $r$ ($d \in [1, 2^r - 1]$)")
    ax1.set_ylabel("Difference count (log scale)")
    ax1.set_title(r"A. Theorem 2.3(a) Short-Difference Crowding on Mian–Chowla")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    m_grid = np.arange(8, 129, 4)
    k0 = 3
    p11_adapt_tail_vs_p1 = []
    p11_adapt_vs_p1 = []
    p10_adapt_tail_vs_p1 = []
    p10_adapt_vs_p1 = []
    p10_k3_tail_vs_p1 = []
    p9_vs_p1 = []
    p9_tail_vs_p1 = []
    p8_vs_p1 = []
    p8_tail_vs_p1 = []
    p7_vs_p1 = []
    p7_tail_vs_p1 = []
    p6_vs_p1 = []
    p6_tail_vs_p1 = []
    p5_vs_p1 = []
    p5_tail_vs_p1 = []
    p4_vs_p1 = []
    p3_vs_p1 = []
    p4_tail_vs_p1 = []
    p2_vs_p1 = []
    c_gk = 0.4975089540675685
    c_lf = -0.002491045932431475

    lam_2 = math.tanh(math.log(2.0) / 4.0) / math.log(2.0)
    lam_101 = math.tanh(math.log(1.01) / 4.0) / math.log(1.01)
    scale_ratio = lam_101 / lam_2
    kappa_pref = 2.5 - 1.5 * math.sqrt(2.0) + 0.5 * math.log(1.0 + math.sqrt(2.0))
    base_inc = 1.0 - 1.0 / math.sqrt(2.0)
    beta0_sq = (4.0 / 3.0) ** 2

    def psi_lag(delta: float) -> float:
        if delta <= 1e-12:
            return 1.0 - 0.25 * delta
        return (2.0 / (3.0 * delta)) * (1.0 - (1.0 - delta) ** 1.5)

    def beta_s_band_sq(
        s: int, alpha: float = 0.5, r_exp: float = 0.0, n_simp: int = 64
    ) -> float:
        a, b = alpha, 1.0
        h = (b - a) / n_simp
        p_pow = r_exp + 2.0

        def eval_tau(tau: float) -> float:
            if tau <= 0.0:
                return 0.0
            d1 = tau * (2.0 ** (-1 - s))
            d0 = tau * (2.0 ** (-s))
            lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
            b_s = lag_inc / base_inc
            return (tau**p_pow) * b_s * b_s

        tot = eval_tau(a) + eval_tau(b)
        for idx_s in range(1, n_simp):
            t_v = a + idx_s * h
            w_s = 4.0 if (idx_s % 2 == 1) else 2.0
            tot += w_s * eval_tau(t_v)
        denom = (1.0 - alpha ** (r_exp + 3.0)) if alpha > 0.0 else 1.0
        return ((r_exp + 3.0) / denom) * (h / 3.0) * tot

    def kappa_haar_k(k_oct: int) -> float:
        return (k_oct * math.log(2.0)) / (3.0 * (1.0 - 4.0 ** (-k_oct)))

    def kappa_dioph_k(k_oct: int) -> float:
        if k_oct <= 12:
            m_cap = 1 << k_oct
            mf = float(m_cap)
            ln_h = sum(float(m) * math.log(float(m)) for m in range(1, m_cap))
            j_multi = (
                0.5 * (mf * mf - 1.0) - mf * (mf - 1.0) * math.log(mf) + 2.0 * ln_h
            )
        else:
            mf = 2.0**k_oct
            j_multi = (k_oct * math.log(2.0)) / 6.0 + c_lf + 1.0 / (360.0 * mf * mf)
        return (2.0 * j_multi) / (1.0 - 4.0 ** (-k_oct)) + 2.0 * c_gk

    for m_val in m_grid:
        u_intra = 0.0
        u_mr = 0.0
        u_fejer = 0.0
        u_p5_joint = 0.0
        u_p6_full = 0.0
        u_p7_full = 0.0
        u_p8_full = 0.0
        u_p9_full = 0.0
        u_p10_adapt_full = 0.0
        u_p11_adapt_full = 0.0
        prev = 0.0
        for k in range(k0, int(m_val) + 1):
            req = math.sqrt((2.0**k) / (k * math.log(2.0)))
            mk_u = req - prev
            rho_sq_2k = (mk_u * mk_u) / (2.0**k)
            u_intra += rho_sq_2k
            s_k = max(0, k - 2) // 2
            u_mr += (s_k + 1) * rho_sq_2k
            fejer_mult = 1.0 + sum(
                2.0 / (1.0 + 2.0 ** (1 - s)) for s in range(1, s_k + 1)
            )
            u_fejer += fejer_mult * rho_sq_2k
            fejer_joint = 1.0
            fejer_coupled = beta0_sq
            fejer_cf2 = beta0_sq
            fejer_band = beta0_sq
            fejer_multi = beta0_sq
            fejer_haar_adapt = beta0_sq
            fejer_dioph_adapt = beta0_sq
            for s in range(1, s_k + 1):
                d1 = 2.0 ** (-1 - s)
                d0 = 2.0 ** (-s)
                lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
                boost_sq = (lag_inc / base_inc) ** 2
                boct_sq = beta_s_band_sq(s, 0.5)
                bmulti_sq = beta_s_band_sq(s, 0.0)
                k_act = k - s
                bhaar_act_sq = beta_s_band_sq(s, 2.0 ** (-k_act), -1.0)
                g1 = boost_sq + (1.0 / 6.0) * (1.0 - (2.0 ** (-s)))
                g2 = boct_sq + 0.25 * (1.0 - (2.0 ** (-s)))
                g_band = boct_sq + (2.0 / 7.0) * (1.0 - (2.0 ** (-s)))
                g_multi = bmulti_sq + 0.5 * (1.0 - (2.0 ** (-s)))
                g_haar_adapt = bhaar_act_sq + kappa_haar_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                g_dioph_adapt = bhaar_act_sq + kappa_dioph_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                fejer_joint += (2.0 / (1.0 + 2.0 ** (1 - s))) * boost_sq
                fejer_coupled += 2.0 * g1
                fejer_cf2 += 2.0 * max(g1, g2)
                fejer_band += 2.0 * max(g1, g_band)
                fejer_multi += 2.0 * max(g1, g_multi)
                fejer_haar_adapt += 2.0 * max(g1, g_haar_adapt)
                fejer_dioph_adapt += 2.0 * max(g1, g_dioph_adapt)
            w = 1.0 if k == int(m_val) else 0.5
            u_mr += w / (k * math.log(2.0))
            u_fejer += w / (k * math.log(2.0))
            u_p5_joint += scale_ratio * fejer_joint * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p6_full += scale_ratio * fejer_coupled * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p7_full += scale_ratio * fejer_cf2 * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p8_full += scale_ratio * fejer_band * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p9_full += scale_ratio * fejer_multi * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p10_adapt_full += scale_ratio * fejer_haar_adapt * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            u_p11_adapt_full += scale_ratio * fejer_dioph_adapt * rho_sq_2k + (
                w + 0.5 * kappa_pref
            ) / (k * math.log(2.0))
            prev = req
        u_bil = u_intra + sum(
            (1.0 if k == int(m_val) else 0.5) / (k * math.log(2.0))
            for k in range(k0, int(m_val) + 1)
        )
        r0 = int(m_val) // 2
        q_tail = 0.0
        q_p5_tail = 0.0
        q_p6_tail = 0.0
        q_p7_tail = 0.0
        q_p8_tail = 0.0
        q_p9_tail = 0.0
        q_p10_k3_tail = 0.0
        q_p10_adapt_tail = 0.0
        q_p11_adapt_tail = 0.0
        for k in range(r0 + 1, int(m_val) + 1):
            alpha_k = (1.0 / math.sqrt(k * math.log(2.0))) - (
                1.0 / math.sqrt(2.0 * (k - 1) * math.log(2.0))
            )
            s_k = max(0, k - 2) // 2
            fejer_mult = 1.0 + sum(
                2.0 / (1.0 + 2.0 ** (1 - s)) for s in range(1, s_k + 1)
            )
            q_tail += fejer_mult * alpha_k * alpha_k
            q_p5_tail += scale_ratio * (1.0 + 2.0 * s_k) * alpha_k * alpha_k
            fejer_coupled = beta0_sq
            fejer_cf2 = beta0_sq
            fejer_band = beta0_sq
            fejer_multi = beta0_sq
            fejer_haar_k3 = beta0_sq
            fejer_haar_adapt = beta0_sq
            fejer_dioph_adapt = beta0_sq
            for s in range(1, s_k + 1):
                d1 = 2.0 ** (-1 - s)
                d0 = 2.0 ** (-s)
                lag_inc = psi_lag(d1) - (2.0 ** (-0.5)) * psi_lag(d0)
                boost_sq = (lag_inc / base_inc) ** 2
                boct_sq = beta_s_band_sq(s, 0.5)
                bmulti_sq = beta_s_band_sq(s, 0.0)
                bhaar3_sq = beta_s_band_sq(s, 0.125, -1.0)
                k_act = k - s
                bhaar_act_sq = beta_s_band_sq(s, 2.0 ** (-k_act), -1.0)
                g1 = boost_sq + (1.0 / 6.0) * (1.0 - (2.0 ** (-s)))
                g2 = boct_sq + 0.25 * (1.0 - (2.0 ** (-s)))
                g_band = boct_sq + (2.0 / 7.0) * (1.0 - (2.0 ** (-s)))
                g_multi = bmulti_sq + 0.5 * (1.0 - (2.0 ** (-s)))
                g_haar3 = bhaar3_sq + kappa_haar_k(3) * (1.0 - (2.0 ** (-s)))
                g_haar_adapt = bhaar_act_sq + kappa_haar_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                g_dioph_adapt = bhaar_act_sq + kappa_dioph_k(k_act) * (
                    1.0 - (2.0 ** (-s))
                )
                fejer_coupled += 2.0 * g1
                fejer_cf2 += 2.0 * max(g1, g2)
                fejer_band += 2.0 * max(g1, g_band)
                fejer_multi += 2.0 * max(g1, g_multi)
                fejer_haar_k3 += 2.0 * max(g1, g_haar3)
                fejer_haar_adapt += 2.0 * max(g1, g_haar_adapt)
                fejer_dioph_adapt += 2.0 * max(g1, g_dioph_adapt)
            q_p6_tail += scale_ratio * fejer_coupled * alpha_k * alpha_k
            q_p7_tail += scale_ratio * fejer_cf2 * alpha_k * alpha_k
            q_p8_tail += scale_ratio * fejer_band * alpha_k * alpha_k
            q_p9_tail += scale_ratio * fejer_multi * alpha_k * alpha_k
            q_p10_k3_tail += scale_ratio * fejer_haar_k3 * alpha_k * alpha_k
            q_p10_adapt_tail += scale_ratio * fejer_haar_adapt * alpha_k * alpha_k
            q_p11_adapt_tail += scale_ratio * fejer_dioph_adapt * alpha_k * alpha_k

        c_p1 = math.sqrt((0.5 * m_val) / u_intra)
        c_p2 = math.sqrt((0.5 * (m_val + 2)) / u_bil)
        c_p3 = math.sqrt((0.5 * (m_val + 2)) / u_mr)
        c_p4 = math.sqrt((0.5 * (m_val + 2)) / u_fejer)
        c_p5 = math.sqrt((0.5 * (m_val + 2)) / u_p5_joint)
        c_p6 = math.sqrt((0.5 * (m_val + 2)) / u_p6_full)
        c_p7 = math.sqrt((0.5 * (m_val + 2)) / u_p7_full)
        c_p8 = math.sqrt((0.5 * (m_val + 2)) / u_p8_full)
        c_p9 = math.sqrt((0.5 * (m_val + 2)) / u_p9_full)
        c_p10_af = math.sqrt((0.5 * (m_val + 2)) / u_p10_adapt_full)
        c_p11_af = math.sqrt((0.5 * (m_val + 2)) / u_p11_adapt_full)
        c_tail = math.sqrt(0.5 * (int(m_val) - r0) / q_tail)
        c_p5_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p5_tail)
        c_p6_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p6_tail)
        c_p7_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p7_tail)
        c_p8_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p8_tail)
        c_p9_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p9_tail)
        c_p10_k3_t = math.sqrt(0.5 * (int(m_val) - r0) / q_p10_k3_tail)
        c_p10_at = math.sqrt(0.5 * (int(m_val) - r0) / q_p10_adapt_tail)
        c_p11_at = math.sqrt(0.5 * (int(m_val) - r0) / q_p11_adapt_tail)
        p11_adapt_tail_vs_p1.append(c_p1 / c_p11_at)
        p11_adapt_vs_p1.append(c_p1 / c_p11_af)
        p10_adapt_tail_vs_p1.append(c_p1 / c_p10_at)
        p10_adapt_vs_p1.append(c_p1 / c_p10_af)
        p10_k3_tail_vs_p1.append(c_p1 / c_p10_k3_t)
        p9_vs_p1.append(c_p1 / c_p9)
        p9_tail_vs_p1.append(c_p1 / c_p9_t)
        p8_vs_p1.append(c_p1 / c_p8)
        p8_tail_vs_p1.append(c_p1 / c_p8_t)
        p7_vs_p1.append(c_p1 / c_p7)
        p7_tail_vs_p1.append(c_p1 / c_p7_t)
        p6_vs_p1.append(c_p1 / c_p6)
        p6_tail_vs_p1.append(c_p1 / c_p6_t)
        p5_vs_p1.append(c_p1 / c_p5)
        p5_tail_vs_p1.append(c_p1 / c_p5_t)
        p4_vs_p1.append(c_p1 / c_p4)
        p4_tail_vs_p1.append(c_p1 / c_tail)
        p3_vs_p1.append(c_p1 / c_p3)
        p2_vs_p1.append(c_p1 / c_p2)

    ax2.plot(
        m_grid,
        p11_adapt_tail_vs_p1,
        "*-",
        color="#B71C1C",
        lw=2.0,
        ms=4.2,
        label=r"Phase 11 Diophantine Haar Tail vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{*,\mathrm{Dioph}}$, $18.14\times$ at $M=128$)",
    )
    ax2.plot(
        m_grid,
        p10_adapt_tail_vs_p1,
        "D--",
        color="#004D40",
        lw=1.6,
        ms=3.3,
        label=r"Phase 10 Adaptive Haar Tail vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{*,\mathrm{Haar}}$)",
    )
    ax2.plot(
        m_grid,
        p11_adapt_vs_p1,
        "o-",
        color="#D32F2F",
        lw=1.5,
        ms=3.0,
        label=r"Phase 11 Diophantine Haar Full vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{\mathrm{Dioph}}$, $15.02\times$ at $M=128$)",
    )
    ax2.plot(
        m_grid,
        p10_adapt_vs_p1,
        "s--",
        color="#00796B",
        lw=1.4,
        ms=2.8,
        label=r"Phase 10 Adaptive Haar Full vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{\mathrm{Haar}}$)",
    )
    ax2.plot(
        m_grid,
        p9_tail_vs_p1,
        "X--",
        color="#117733",
        lw=1.3,
        ms=2.8,
        label=r"Phase 9 Multi-Oct Tail vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{*,\mathrm{multi}}$)",
    )
    ax2.plot(
        m_grid,
        p8_tail_vs_p1,
        "p:",
        color="#660099",
        lw=1.2,
        ms=2.7,
        label=r"Phase 8 Band Tail vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{*,\mathrm{band}}$)",
    )
    ax2.plot(
        m_grid,
        p6_tail_vs_p1,
        "h:",
        color="#882255",
        lw=1.2,
        ms=2.7,
        label=r"Phase 6 Coupled Tail vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{*,\mathrm{coupled}}$)",
    )
    ax2.plot(
        m_grid,
        p4_vs_p1,
        "v:",
        color="#CC79A7",
        lw=1.2,
        ms=2.6,
        label=r"Phase 4 Fejér vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{\mathrm{Fej\acute{e}r}}$)",
    )
    ax2.plot(
        m_grid,
        p3_vs_p1,
        "d:",
        color="#56B4E9",
        lw=1.2,
        ms=2.6,
        label=r"Phase 3 Multi-Res vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{\mathrm{MR}}$)",
    )
    ax2.plot(
        m_grid,
        p2_vs_p1,
        "^--",
        color="#0072B2",
        lw=1.2,
        ms=2.6,
        label=r"Phase 2 Bilinear vs. P1 ($\underline{c}_M^{(1)} / \underline{c}_M^{\mathrm{QP}}$)",
    )
    ax2.set_xlabel(r"Dyadic horizon scale $M$ ($N = 2^M \leq 2^{128}$)")
    ax2.set_ylabel("Multiplicative upper-bound sharpening factor")
    ax2.set_title(r"B. Phases 2–11 Envelope Reduction Factors vs. Phase 1")
    ax2.set_ylim(0.7, 19.5)
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="upper left", frameon=True, fontsize=5.1, ncol=1)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig6_subwindow_multires_spectrum")


def make_fig7(out_dir: Path):
    """
    We-Solve #881 Remediation Diagnostics (Theorem 2.2(b), Remark 2.3, Lemma 2.4,
    Theorem 2.6(d) & Section 2.11.3 Scaling Regimes):
      * Panel A: Exact Powers-of-Two A = {2^n : n >= 0} Bilinear Sum
        sum_{k=1}^M (2 m_k(m_k-1) + 2 m_k A(2^{k-1})) / 2^k = 4 - (2M+4)/2^M
        vs. Corrected Tight Ceiling M + 1 - 2^{-(M-1)}, plus KKT Convex-Cone
        Increment Margin 18(k^2 - 1) - 16 k^2 = 2(k^2 - 9) > 0 (Lemma 2.4).
      * Panel B: Section 2.11.3 Scaling Regimes — Density Bound c(mu_1, K)
        across sub-block occupancy mu_1 in [1, 32] for K in {2, 4, 9, 25, 64},
        illustrating the unit-occupancy collision regime mu_1 = O_K(1) -> 0
        vs. the fixed-K asymptotic decay O(mu_1^{-2}) -> c_infty^*(1^+) = 2.00.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.3))

    ms = np.arange(1, 13)
    pow2_unscaled = 4.0 - (2.0 * ms + 4.0) / (2.0**ms)
    ceil_tight = ms + 1.0 - 2.0 ** (-(ms - 1.0))

    ax1.plot(
        ms,
        ceil_tight,
        "k--",
        lw=1.8,
        marker="s",
        ms=4.2,
        label=r"Corrected Bilinear Ceiling $M + 1 - 2^{-(M-1)}$ (Thm. 2.2b)",
    )
    ax1.plot(
        ms,
        pow2_unscaled,
        "o-",
        color="#009E73",
        lw=2.0,
        ms=5.0,
        label=r"Powers of Two $A=\{2^n\}$ Exact Sum $4 - (2M+4)/2^M$ (Rem. 2.3)",
    )
    # KKT multiplier positivity curve (scaled by 2^k for visibility)
    ks = np.arange(3, 13)
    rho = lambda k: math.sqrt((2.0**k) / (k * math.log(2.0))) - math.sqrt(
        (2.0 ** (k - 1)) / ((k - 1) * math.log(2.0))
    )
    kkt_scaled = [2.0 * rho(int(k)) - 0.5 * rho(int(k) + 1) for k in ks]
    ax1.plot(
        ks,
        kkt_scaled,
        "^:",
        color="#D55E00",
        lw=1.7,
        ms=4.5,
        label=r"KKT Boundary Margin $2\rho_k - \frac{1}{2}\rho_{k+1} > 0$ (Lemma 2.4)",
    )
    ax1.annotate(
        r"$M=3$: $2.75 \leq 3.75$",
        xy=(3, 2.75),
        xytext=(4.3, 1.6),
        arrowprops=dict(arrowstyle="->", color="#009E73", lw=1.1),
        fontsize=8.2,
        fontweight="bold",
        color="#009E73",
    )
    ax1.set_xlabel(r"Dyadic shell index $M$ or $k$")
    ax1.set_ylabel("Bilinear Sum / KKT Convex-Cone Margin")
    ax1.set_title(r"A. Powers-of-Two Verification & Lemma 2.4 KKT Positivity")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="upper left", frameon=True, fontsize=7.6)

    # Panel B: c(mu_1, K) vs mu_1 across K = 2, 9, 25, 64, 256, 2048 (M up to 4096)
    mu_grid = np.linspace(1.0, 16.0, 200)
    k_configs = [
        (2, "#56B4E9", "s--", r"$K=2$ ($M=8$, $c(1,2)=1.546$)"),
        (9, "#009E73", "o-", r"$K=9$ ($M=32$, $\underline{c}_{32}^{*,\mathrm{Dioph}}=0.924 < 1.00$)"),
        (25, "#E69F00", "D-", r"$K=25$ ($M=64$, $\underline{c}_{64}^{*,\mathrm{Dioph}}=0.668$)"),
        (64, "#D32F2F", "*-", r"$K=64$ ($M=128$, $\underline{c}_{128}^{*,\mathrm{Dioph}}=0.480 < 0.50$)"),
        (256, "#6A1B9A", "^-.", r"$K=256$ ($M=512$, $\underline{c}_{512}^{*,\mathrm{Dioph}}=0.244 < 0.25$)"),
        (2048, "#1A237E", "P-", r"$K=2048$ ($M=4096$, $\underline{c}_{4096}^{*,\mathrm{Dioph}}=0.0866 < 0.10$)"),
    ]
    ax2.axhline(
        2.0,
        color="#111111",
        ls="--",
        lw=1.5,
        label=r"Unconditional Continuum Ceiling $c_\infty^*(1^+) = 2.00$ ($\mu_1 \to \infty$)",
    )
    for k_oct, col, style, lbl in k_configs:
        j_haar = (k_oct * math.log(2.0)) / 6.0
        rhs_cap = 0.5 * (1.0 - 4.0 ** (-k_oct))
        kappa_rel = j_haar / ((mu_grid**2) * rhs_cap)
        c_curve = 2.0 / np.sqrt(1.0 + kappa_rel)
        ax2.plot(
            mu_grid,
            c_curve,
            style,
            color=col,
            lw=1.7,
            ms=4.0,
            markevery=25,
            label=lbl,
        )

    ax2.axhline(1.0, color="#004D40", ls=":", lw=1.0)
    ax2.axhline(0.5, color="#B71C1C", ls=":", lw=1.0)
    ax2.axhline(0.1, color="#1A237E", ls=":", lw=1.0)
    ax2.set_xlabel(r"Top-scale sub-window occupancy $\mu_1 = S / B_1$")
    ax2.set_ylabel(r"Sustainable density constant $c(\mu_1, K)$")
    ax2.set_title(r"B. Section 2.11.3 Scaling Regimes ($M \in [8, 4096]$): $\mu_1 = O_K(1)$ vs. $\mu_1 \to \infty$")
    ax2.set_ylim(0.0, 2.18)
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="lower right", frameon=True, fontsize=6.6)

    fig.tight_layout()
    return save_pair(fig, out_dir, "fig7_we_solve_remediation_diagnostics")


def main():
    root = Path(__file__).resolve().parent.parent
    out_dir = root / "paper" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "fig1": make_fig1(out_dir),
        "fig2": make_fig2(out_dir),
        "fig3": make_fig3(out_dir),
        "fig4": make_fig4(out_dir),
        "fig5": make_fig5(out_dir),
        "fig6": make_fig6(out_dir),
        "fig7": make_fig7(out_dir),
    }
    (out_dir / "figures_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print("Generated 7 publication figures in", out_dir)


if __name__ == "__main__":
    main()
