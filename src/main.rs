//! ============================================================================
//! `main.rs` — CLI entrypoint for `erdos1191`: GPU-accelerated quantitative
//! analysis of Erdős Problem #1191 ($1,000 open prize: lower-envelope fluctuations
//! $\liminf_{x\to\infty} A(x)\sqrt{(\ln x)/x} = 0$ of infinite Sidon sets).
//! ============================================================================

pub mod exact;
pub mod gpu;
pub mod records;
pub mod verify_cpu;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "erdos1191",
    version,
    about = "GPU-accelerated quantitative analysis of Erdős Problem #1191 ($1,000 open prize: liminf A(x) sqrt(log x / x) = 0 for infinite Sidon sets)"
)]
struct Cli {
    /// CUDA device ordinal
    #[arg(long, global = true, default_value_t = 0)]
    device: usize,

    /// NVRTC target architecture
    #[arg(long, global = true, default_value = "compute_86")]
    arch: String,

    /// Force CPU-only execution
    #[arg(long, global = true, default_value_t = false)]
    no_gpu: bool,

    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Display problem status, Theorem 2.1 Harmonic-Difference & Dyadic Budget certificates, and Mian–Chowla vs. Singer-Block comparison
    Info,

    /// Verify Sidon ($B_2$) validity and compute Theorem 2.1(a) Harmonic-Integral Conservation + Theorem 2.1(b) Dyadic Budget
    Verify {
        /// Comma-separated positive integers
        #[arg(long)]
        set: Option<String>,

        /// Path to file containing a candidate Sidon set
        #[arg(long)]
        set_file: Option<PathBuf>,

        /// Horizon N (defaults to max(set))
        #[arg(long)]
        n: Option<usize>,

        /// Lower window threshold x_0 >= 2 for min_{x in [x_0, N]} A(x) sqrt(ln x / x)
        #[arg(long, default_value_t = 4)]
        x0: usize,
    },

    /// Exhaustively determine F_2(N) (maximum Sidon subset of {1, ..., N}, OEIS A143824 / A003022)
    Exact {
        /// Interval size N (1 <= N <= 63)
        #[arg(long, default_value_t = 26)]
        n: usize,

        /// CPU prefix expansion limit before GPU DFS handoff
        #[arg(long, default_value_t = 6)]
        split_x: u32,

        /// Drop OEIS interior F_2(w) table in favor of elementary Erdős–Turán bounds
        #[arg(long, default_value_t = false)]
        unconditional: bool,

        /// Force CPU-only branch-and-bound
        #[arg(long, default_value_t = false)]
        cpu_only: bool,
    },

    /// Solve the Maximin Lower-Envelope Sidon Problem: max_{A subseteq [N] Sidon} min_{x_0 <= x <= N} A(x) sqrt(ln x / x)
    Maximin {
        /// Horizon N (2 <= N <= 63)
        #[arg(long, default_value_t = 32)]
        n: usize,

        /// Lower threshold x_0 (2 <= x_0 <= N)
        #[arg(long, default_value_t = 4)]
        x0: usize,

        /// CPU prefix expansion limit before GPU DFS handoff
        #[arg(long, default_value_t = 6)]
        split_x: u32,

        /// Force CPU-only branch-and-bound
        #[arg(long, default_value_t = false)]
        cpu_only: bool,
    },

    /// Generate k terms of the Mian–Chowla greedy Sidon sequence (OEIS A005282) and evaluate its Erdős #1191 decay
    MianChowla {
        /// Number of terms k to generate
        #[arg(long, default_value_t = 128)]
        terms: usize,
    },

    /// Benchmark GPU vs. CPU O(1) bit-shift Sidon branch-and-bound and Dyadic Budget batch evaluation
    Bench,
}

fn parse_set_str(raw: &str) -> Result<Vec<u64>> {
    let trimmed = raw.trim();
    if trimmed.starts_with('[') {
        let mut v: Vec<u64> = serde_json::from_str(trimmed)?;
        v.sort_unstable();
        v.dedup();
        return Ok(v);
    }
    if trimmed.starts_with('{') {
        #[derive(serde::Deserialize)]
        struct Obj {
            set: Vec<u64>,
        }
        let mut o: Obj = serde_json::from_str(trimmed)?;
        o.set.sort_unstable();
        o.set.dedup();
        return Ok(o.set);
    }
    let mut out = Vec::new();
    for tok in trimmed.split(|c: char| c == ',' || c.is_whitespace()) {
        if !tok.is_empty() {
            out.push(
                tok.parse::<u64>()
                    .with_context(|| format!("invalid integer token '{tok}'"))?,
            );
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

fn cmd_info() {
    println!("=== Erdős Problem #1191: Lower-Envelope Fluctuations of Infinite Sidon Sets ===");
    println!("Prize: $1,000 (Paul Erdős; OPEN)");
    println!("Conjecture: Every infinite Sidon set A ⊆ N satisfies lim inf_{{x->∞}} A(x) * sqrt((ln x) / x) = 0.");
    println!();
    println!("Known Milestones:");
    println!("  * Finite Horizon (Singer 1938; Erdős–Turán 1941): F_2(N) = (1 + o(1)) sqrt(N) (exact for N <= 586 via OEIS A003022)");
    println!("  * Greedy Infinite Sidon (Mian–Chowla 1944, OEIS A005282): A(x) ~ c x^(1/3) => Φ_MC(x) ~ c x^(-1/6) sqrt(ln x) -> 0");
    println!("  * Infinite Upper Density (limsup): Erdős (1955) >= 1/2; Krückeberg (1961) >= 1/sqrt(2); Cilleruelo (2014) = 1");
    println!("  * Simultaneous Power Growth (Ruzsa 1998): exists infinite Sidon A with A(x) = x^(sqrt(2)-1+o(1)) ≈ x^0.4142 for all x");
    println!("  * Infinite Lower Density (Erdős 1955, in Stöhr 1955): lim inf_{{x->∞}} A(x) * sqrt((ln x) / x) <= C < ∞");
    println!();
    println!("Theorem 2.1(b) Exact Dyadic Difference-Budget Verification on Mian–Chowla (OEIS A005282, 128 terms):");
    let mc = verify_cpu::mian_chowla_greedy_cpu(128);
    println!(
        "  {:>3} | {:>8} | {:>6} | {:>12} | {:>12} | {:>12} | {:>10} | {:>10}",
        "M", "2^M", "A(2^M)", "LHS Dyadic", "Diff-Bound", "RHS (M-1)/2", "min Φ(2^k)", "Φ_MC(2^M)"
    );
    println!("  ----+----------+--------+--------------+--------------+--------------+------------+-----------");
    for m in [3usize, 5, 7, 9, 11, 13, 15] {
        let rep = records::evaluate_dyadic_budget(&mc, m);
        let a_2m = mc.iter().filter(|&&x| x <= rep.horizon_n).count();
        let phi_2m = records::erdos1191_phi(a_2m, rep.horizon_n as usize);
        println!(
            "  {:>3} | {:>8} | {:>6} | {:>12.6} | {:>12.6} | {:>12.6} | {:>10.5} | {:>10.5}",
            m,
            rep.horizon_n,
            a_2m,
            rep.lhs_dyadic_sum,
            rep.exact_diff_weighted_sum,
            rep.rhs_universal_ceiling,
            rep.min_phi_dyadic,
            phi_2m
        );
    }
    println!();
    println!("Theorem 3.1 Burst–Desert Dichotomy on Krückeberg–Cilleruelo Dyadic Singer-Block Sidon Set:");
    let (kb, ranges) = records::kruckeberg_singer_block_sidon(&[2, 3, 4, 5, 7, 8], 3);
    println!(
        "  {:>5} | {:>4} | {:>18} | {:>14} | {:>14}",
        "Block", "q_r", "Span [B_r, E_r]", "Desert Φ(B_r-1)", "Burst Φ(E_r)"
    );
    println!("  ------+------+--------------------+----------------+---------------");
    for (idx, (&q, &(start, end))) in [2usize, 3, 4, 5, 7, 8]
        .iter()
        .zip(ranges.iter())
        .enumerate()
    {
        let desert_x = (start.saturating_sub(1)).max(2) as usize;
        let desert_cnt = kb.iter().filter(|&&v| (v as usize) <= desert_x).count();
        let burst_cnt = kb.iter().filter(|&&v| v <= end).count();
        let phi_desert = records::erdos1191_phi(desert_cnt, desert_x);
        let phi_burst = records::erdos1191_phi(burst_cnt, end as usize);
        println!(
            "  {:>5} | {:>4} | [{:>7}, {:>7}] | {:>14.5} | {:>14.5}",
            idx + 1,
            q,
            start,
            end,
            phi_desert,
            phi_burst
        );
    }
    println!();
    println!("Theorem 2.2 Exact Cross-Shell Bilinear Difference Conservation on Mian–Chowla (OEIS A005282):");
    println!(
        "  {:>3} | {:>8} | {:>10} | {:>10} | {:>14} | {:>14} | {:>14}",
        "M", "2^M", "IntraPair", "InterPair", "FarInter 2Shell", "Intra-Only LHS", "Bilinear LHS"
    );
    println!("  ----+----------+------------+------------+----------------+----------------+---------------");
    for m in [4usize, 6, 8, 10, 12, 13] {
        let rep = records::evaluate_cross_shell_bilinear(&mc, m);
        println!(
            "  {:>3} | {:>8} | {:>10} | {:>10} | {:>6}/{:<7} | {:>14.6} | {:>14.6}",
            m,
            rep.horizon_n,
            rep.intra_pairs_total,
            rep.inter_pairs_total,
            rep.far_inter_in_two_shells,
            rep.far_inter_pairs_total,
            rep.intra_only_dyadic_lhs,
            rep.bilinear_dyadic_lhs
        );
    }
    println!();
    println!(
        "Theorem 3.2 Cross-Shell Bilinear QP Envelope vs. Intra-Shell-Only Relaxation (k_0 = 3):"
    );
    println!(
        "  {:>4} | {:>18} | {:>20} | {:>16}",
        "M", "Intra-Only c_M^(1)", "Cross-Shell c_M^(QP)", "Reduction Factor"
    );
    println!("  -----+--------------------+----------------------+-----------------");
    for m in [8usize, 16, 24, 32, 48, 64] {
        let qp = records::solve_cross_shell_qp_envelope(m, 3);
        println!(
            "  {:>4} | {:>18.6} | {:>20.6} | {:>15.4}x",
            m, qp.c_intra_only_upper, qp.c_bilinear_cross_shell_upper, qp.reduction_factor
        );
    }
    println!();
    println!("Theorem 2.3 Exact Sub-Window Short-Difference Projection & Crowding on Mian–Chowla (M = 13):");
    let sub_rep = records::evaluate_subwindow_projection(&mc, 13);
    println!(
        "  {:>3} | {:>7} | {:>8} | {:>12} | {:>12} | {:>12} | {:>12}",
        "r", "2^r-1", "C(A,2)", "PigeonholeLB", "ExactSubWin", "ActualShortD", "Cap 2^r - 1"
    );
    println!(
        "  ----+---------+----------+--------------+--------------+--------------+-------------"
    );
    for row in sub_rep.scale_rows.iter().step_by(2) {
        println!(
            "  {:>3} | {:>7} | {:>8} | {:>12} | {:>12} | {:>12} | {:>12}",
            row.r_scale,
            row.universal_cap_2r_minus_1,
            row.prefix_pairs_binom,
            row.total_pigeonhole_lhs,
            row.total_exact_subwindow_lhs,
            row.actual_short_diffs_rhs,
            row.universal_cap_2r_minus_1
        );
    }
    println!();
    println!("Theorem 2.4 & 3.3 Multi-Resolution & Sliding-Window Fejér QP Envelopes across Phases 1–4 (k_0 = 3):");
    println!(
        "  {:>4} | {:>12} | {:>12} | {:>12} | {:>12} | {:>14} | {:>9} | {:>9}",
        "M",
        "Phase1 Intra",
        "Phase2 Bilin",
        "Phase3 MultR",
        "Phase4 Fejér",
        "Tail c_M*(M/2)",
        "P3 vs P1",
        "P4 vs P1"
    );
    println!("  -----+--------------+--------------+--------------+--------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 24, 32, 48, 64] {
        let mr = records::solve_multires_subwindow_qp_envelope(m, 3);
        println!(
            "  {:>4} | {:>12.6} | {:>12.6} | {:>12.6} | {:>12.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            mr.c_intra_only_upper,
            mr.c_bilinear_cross_shell_upper,
            mr.c_multires_subwindow_upper,
            mr.c_fejer_multires_upper,
            mr.c_monotone_subwindow_star,
            mr.phase3_vs_phase1_factor,
            mr.phase4_vs_phase1_factor
        );
    }
    println!(
        "  Asymptotic Pinched Limit (M -> ∞, B = 2): c_∞*(2) = (1 + √2)√(ln 2) = {:.6}",
        (1.0 + std::f64::consts::SQRT_2) * std::f64::consts::LN_2.sqrt()
    );
    println!();
    println!("Theorem 2.5 Direct Prefix-to-Block Unsplit Fejér Coupling & Base-B Variational Bound (B -> 1+):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "Phase4 Fejér",
        "Phase4 Tail*",
        "P5 Joint(1.01)",
        "P5 Tail*(1.01)",
        "P5 vs P1",
        "P5 vs P4"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 24, 32, 48, 64] {
        let p5 = records::solve_base_b_joint_fejer_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p5.phase4_full_fejer,
            p5.phase4_tail_star,
            p5.c_full_joint_fejer,
            p5.c_tail_octave_joint,
            p5.phase5_vs_phase1_factor,
            p5.phase5_vs_phase4_factor
        );
    }
    println!();
    println!("Theorem 2.5(b) Exact Hyperbolic-Tangent Shell-Efficiency Law Λ(B) = tanh((ln B)/4) / ln B:");
    println!(
        "  {:>10} | {:>14} | {:>14} | {:>16}",
        "Base B", "Λ(B)", "c_∞*(B)", "Gain vs B = 2"
    );
    println!("  -----------+----------------+----------------+-----------------");
    let c_b2 = (1.0 + std::f64::consts::SQRT_2) * std::f64::consts::LN_2.sqrt();
    for &(label, b) in &[
        ("4.0000", 4.0f64),
        ("2.0000", 2.0),
        ("√2≈1.4142", std::f64::consts::SQRT_2),
        ("1.2500", 1.25),
        ("1.1000", 1.10),
        ("1.0100", 1.01),
    ] {
        let env = records::solve_base_b_joint_fejer_envelope(64, 3, b);
        println!(
            "  {:>10} | {:>14.8} | {:>14.6} | {:>+15.6}",
            label,
            env.lambda_b,
            env.c_infty_star_b,
            env.c_infty_star_b - c_b2
        );
    }
    println!(
        "  {:>10} | {:>14.8} | {:>14.6} | {:>+15.6}",
        "B -> 1+",
        0.25f64,
        2.0f64,
        2.0f64 - c_b2
    );
    println!();
    println!("Theorem 2.6 Exact Coupled Lag-Amplified & Discrete Integrality Envelope (Breaking c_∞*(1+) = 2.000000):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P5 Joint(1.01)",
        "P5 Tail*(1.01)",
        "P6 Coup(1.01)",
        "P6 Tail*(1.01)",
        "P6 vs P1",
        "P6 vs P5"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64] {
        let p6 = records::solve_coupled_lag_integrality_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p6.phase5_full_fejer,
            p6.phase5_tail_star,
            p6.c_coupled_full_fejer,
            p6.c_coupled_tail_star,
            p6.phase6_vs_phase1_factor,
            p6.phase6_vs_phase5_factor
        );
    }
    let p6_64 = records::solve_coupled_lag_integrality_envelope(64, 3, 1.01);
    println!(
        "  Asymptotic Coupled Limit (M -> ∞, B -> 1+): c_∞*,coupled(1+) = 2√(6/7) = {:.6} (c_∞*,coupled(1.01) = {:.6}, Δ_tot = {:.8})",
        2.0f64 * (6.0f64 / 7.0f64).sqrt(),
        p6_64.c_infty_coupled_limit,
        p6_64.delta_tot
    );
    println!();
    println!("Theorem 2.7 Higher-Order Cesàro–Fejér Order-2 & Upper-Octave Convex Envelope (Breaking c_∞*(1+) < 1.800000):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P6 Coup(1.01)",
        "P6 Tail*(1.01)",
        "P7 CF2(1.01)",
        "P7 Tail*(1.01)",
        "P7 vs P1",
        "P7 vs P6"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64] {
        let p7 = records::solve_cesaro_fejer_order2_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p7.phase6_full_fejer,
            p7.phase6_tail_star,
            p7.c_cf2_full_fejer,
            p7.c_cf2_tail_star,
            p7.phase7_vs_phase1_factor,
            p7.phase7_vs_phase6_factor
        );
    }
    let p7_64 = records::solve_cesaro_fejer_order2_envelope(64, 3, 1.01);
    println!(
        "  Asymptotic Order-2 Limit (M -> ∞, B -> 1+): c_∞*,CF2(1+) = 4√5/5 = {:.6} < 1.80 (c_∞*,CF2(1.01) = {:.6}, κ_int^(2) = {:.4})",
        2.0f64 * (4.0f64 / 5.0f64).sqrt(),
        p7_64.c_infty_cf2_limit,
        p7_64.kappa_int_order2
    );
    println!();
    println!("Theorem 2.8 Octave-Matched Order p*(s) ~ 2^s ln 2 & Band-Restricted Cesàro–Fejér Envelope (Breaking c_∞*(1+) < 1.770000):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P7 CF2(1.01)",
        "P7 Tail*(1.01)",
        "P8 Band(1.01)",
        "P8 Tail*(1.01)",
        "P8 vs P1",
        "P8 vs P7"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64] {
        let p8 = records::solve_band_restricted_cesaro_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p8.phase7_full_fejer,
            p8.phase7_tail_star,
            p8.c_band_full_fejer,
            p8.c_band_tail_star,
            p8.phase8_vs_phase1_factor,
            p8.phase8_vs_phase7_factor
        );
    }
    let p8_64 = records::solve_band_restricted_cesaro_envelope(64, 3, 1.01);
    println!(
        "  Asymptotic Band-Restricted Limit (M -> ∞, B -> 1+): c_∞*,band(1+) = 2√7/3 = {:.6} < 1.77 (c_∞*,band(1.01) = {:.6}, κ_int^band(1/2) = {:.6})",
        2.0f64 * (7.0f64 / 9.0f64).sqrt(),
        p8_64.c_infty_band_limit,
        p8_64.kappa_int_band_half
    );
    println!();
    println!("Theorem 2.9 Scale-Invariant Block-Window Conservation & Multi-Octave Integrality Cascade (Breaking c_∞*(1+) < 1.640000 & c_16* < 2.00):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P8 Tail*(1.01)",
        "P9 K=2*(1.01)",
        "P9 K=3*(1.01)",
        "P9 Mult*(1.01)",
        "P9 Full(1.01)",
        "P9 vs P1"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64] {
        let p9 = records::solve_multioctave_cesaro_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>9.6} | {:>8.4}x",
            m,
            p9.phase8_tail_star,
            p9.c_k2_tail_star,
            p9.c_k3_tail_star,
            p9.c_multioct_tail_star,
            p9.c_multioct_full_fejer,
            p9.phase9_vs_phase1_factor
        );
    }
    let p9_64 = records::solve_multioctave_cesaro_envelope(64, 3, 1.01);
    println!(
        "  Asymptotic Multi-Octave Limit (M -> ∞, K -> ∞, B -> 1+): c_∞*,multi-oct(1+) = 2√6/3 = {:.6} < 1.64 (c_∞*,multi-oct(1.01) = {:.6}, κ_int^band(0+) = {:.6})",
        2.0f64 * (2.0f64 / 3.0f64).sqrt(),
        p9_64.c_infty_multioct_limit,
        p9_64.kappa_int_multioct_limit
    );
    println!();
    println!("Theorem 2.10 Power-Law & Scale-Invariant Haar-Measure Multi-Octave Cascade (Breaking 1.00 & c_∞*,Haar(1+; K) -> 0):");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P9 Mult*(1.01)",
        "P10 K=3*(1.01)",
        "P10 Haar*(1.01)",
        "P10 Full(1.01)",
        "P10 vs P1",
        "P10 vs P9"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64] {
        let p10 = records::solve_haar_multioctave_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p10.phase9_tail_star,
            p10.c_haar_k3_tail_star,
            p10.c_haar_adapt_tail_star,
            p10.c_haar_adapt_full_fejer,
            p10.phase10_adapt_vs_phase1_factor,
            p10.phase10_adapt_vs_phase9_factor
        );
    }
    let p10_64 = records::solve_haar_multioctave_envelope(64, 3, 1.01);
    println!(
        "  Power-Law Unit Barrier: r = -15/17 => κ_int = 3.000000, c_∞*,(-15/17)(1+) = 1.000000; r = -9/10 => c_∞*,(-9/10)(1+) = 2√2/3 = {:.6} < 1.00",
        2.0f64 * std::f64::consts::SQRT_2 / 3.0
    );
    println!(
        "  Haar Cascade (r -> -1+): c_∞*,Haar(1+; K=3) = {:.6}, c_∞*,Haar(1+; K=13) = {:.6} < 1.00, lim_{{K->∞}} c_∞*,Haar(1+; K) = 0",
        p10_64.c_infty_haar_k3_limit,
        p10_64.c_infty_haar_k13_limit
    );
    println!();
    println!("Theorem 2.11 Unconditional Diophantine Linear-Dilation Sawtooth Identity & Complete Multi-Octave Haar Closure:");
    println!(
        "  {:>4} | {:>14} | {:>14} | {:>14} | {:>14} | {:>9} | {:>9}",
        "M",
        "P10 Haar*(1.01)",
        "P11 K=3*(1.01)",
        "P11 Dioph*(1.01)",
        "P11 Full(1.01)",
        "P11 vs P1",
        "P11 vs P10"
    );
    println!("  -----+----------------+----------------+----------------+----------------+-----------+----------");
    for m in [8usize, 12, 16, 20, 24, 32, 48, 64, 96, 128] {
        let p11 = records::solve_unconditional_diophantine_haar_envelope(m, 3, 1.01);
        println!(
            "  {:>4} | {:>14.6} | {:>14.6} | {:>14.6} | {:>14.6} | {:>8.4}x | {:>8.4}x",
            m,
            p11.phase10_tail_star,
            p11.c_dioph_k3_tail_star,
            p11.c_dioph_adapt_tail_star,
            p11.c_dioph_adapt_full_fejer,
            p11.phase11_adapt_vs_phase1_factor,
            p11.phase11_adapt_vs_phase10_factor
        );
    }
    let p11_128 = records::solve_unconditional_diophantine_haar_envelope(128, 3, 1.01);
    println!(
        "  Glaisher–Kinkelin / Euler–Mascheroni Sawtooth Constants: C_GK = 2 ln A_GK = {:.10}, C_LF = 2 ln A_GK - 1/2 = {:.10}",
        records::C_GK_GLAISHER_SAWTOOTH,
        records::C_LF_EULER_MASCHERONI_SAWTOOTH
    );
    println!(
        "  Unconditional Diophantine Closure: c_∞*,Dioph(1+; K=3) = {:.6}, c_∞*,Dioph(1+; K=9) = {:.6} < 1.00, c_∞*,Dioph(1+; K=64) = {:.6} < 0.50 -> 0",
        p11_128.c_infty_dioph_k3_limit,
        p11_128.c_infty_dioph_k9_limit,
        p11_128.c_infty_dioph_k64_limit
    );
}

fn cmd_verify(
    set_arg: Option<String>,
    set_file: Option<PathBuf>,
    n_opt: Option<usize>,
    x0: usize,
) -> Result<()> {
    let raw = match (set_arg, set_file) {
        (Some(s), None) => s,
        (None, Some(p)) => {
            fs::read_to_string(&p).with_context(|| format!("failed to read {}", p.display()))?
        }
        _ => bail!("specify exactly one of --set or --set-file"),
    };
    let s = parse_set_str(&raw)?;
    if s.is_empty() {
        bail!("candidate set is empty");
    }
    let n = n_opt.unwrap_or_else(|| *s.last().unwrap() as usize).max(2);
    let rep = verify_cpu::harmonic_conservation_report_cpu(&s, n, x0);

    println!(
        "Set (|A ∩ [1,{n}]| = {}, N = {n}): Sidon (B_2) = {} (difference collisions = {})",
        rep.card,
        if rep.is_sidon { "YES ✓" } else { "NO ✗" },
        rep.violations
    );
    println!(
        "Theorem 2.1(a) Harmonic-Integral Conservation Identity on [1, {n}]:\n  Prefix Energy E_A({n}) = {:.8} (Integral = {:.8} + Boundary = {:.8}; sum (j-1)/a_j = {:.8})\n  + Translation Defect Δ_A({n}) = {:.8}\n  = Harmonic Diff Budget B_A({n}) = {:.8} <= H_{n} - 1 = {:.8} (residual = {:.2e})",
        rep.prefix_energy_e,
        rep.integral_term,
        rep.boundary_term,
        rep.sum_j_minus_1_over_aj,
        rep.translation_defect_delta,
        rep.harmonic_diff_budget_b,
        rep.harmonic_ceiling_h_n_minus_1,
        rep.conservation_residual
    );
    println!(
        "Erdős #1191 Envelope on [{}, {n}]: Φ_A({n}) = {:.5}, min Φ_A(x) = {:.5} (at x = {}), max Φ_A(x) = {:.5} (at x = {})",
        x0.clamp(2, n),
        rep.phi_at_n,
        rep.min_phi_from_x0,
        rep.min_phi_x,
        rep.max_phi_from_x0,
        rep.max_phi_x
    );
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Command::Info => {
            cmd_info();
            Ok(())
        }
        Command::Verify {
            set,
            set_file,
            n,
            x0,
        } => cmd_verify(set, set_file, n, x0),
        Command::Exact {
            n,
            split_x,
            unconditional,
            cpu_only,
        } => {
            let gpu_opt = if cli.no_gpu || cpu_only {
                None
            } else {
                gpu::Gpu::init(cli.device, &cli.arch).ok()
            };
            let t0 = Instant::now();
            let res = exact::solve_f2_sidon_exact(n, split_x, unconditional, gpu_opt.as_ref())?;
            let dt = t0.elapsed();
            println!(
                "F_2({n}) = {} (min_{{x>=2}} Φ_A(x) = {:.5} at x = {}) | nodes = {} | elapsed = {:.3}s | engine = {}",
                res.max_card,
                res.min_phi,
                res.min_phi_x,
                res.total_nodes,
                dt.as_secs_f64(),
                if gpu_opt.is_some() { "CPU+GPU" } else { "CPU" }
            );
            println!("Witness: {:?}", res.witness);
            Ok(())
        }
        Command::Maximin {
            n,
            x0,
            split_x,
            cpu_only,
        } => {
            let gpu_opt = if cli.no_gpu || cpu_only {
                None
            } else {
                gpu::Gpu::init(cli.device, &cli.arch).ok()
            };
            let t0 = Instant::now();
            let res = exact::solve_maximin_lower_envelope_exact(n, x0, split_x, gpu_opt.as_ref())?;
            let dt = t0.elapsed();
            println!(
                "Maximin Sidon μ*({x0}, {n}) = {:.6} (at x* = {}, |A*| = {}) | nodes = {} | elapsed = {:.3}s | engine = {}",
                res.min_phi,
                res.min_phi_x,
                res.max_card,
                res.total_nodes,
                dt.as_secs_f64(),
                if gpu_opt.is_some() { "CPU+GPU" } else { "CPU" }
            );
            println!("Maximin Witness A*: {:?}", res.witness);
            Ok(())
        }
        Command::MianChowla { terms } => {
            let t0 = Instant::now();
            let mc = verify_cpu::mian_chowla_greedy_cpu(terms);
            let dt = t0.elapsed();
            let last = *mc.last().unwrap_or(&1) as usize;
            let phi_last = records::erdos1191_phi(mc.len(), last.max(2));
            println!(
                "Generated {} Mian–Chowla (OEIS A005282) terms in {:.4}s: a_{} = {}, Φ_MC(a_{}) = {:.6}",
                mc.len(),
                dt.as_secs_f64(),
                mc.len(),
                last,
                mc.len(),
                phi_last
            );
            println!("First 20 terms: {:?}", &mc[..mc.len().min(20)]);
            Ok(())
        }
        Command::Bench => {
            let gpu = gpu::Gpu::init(cli.device, &cli.arch)?;
            println!(
                "=== erdos1191 GPU vs. CPU Benchmark on {} ===",
                gpu.device_name
            );
            let n = 35usize;
            let t_cpu = Instant::now();
            let res_cpu = exact::solve_f2_sidon_exact(n, 7, false, None)?;
            let dt_cpu = t_cpu.elapsed();

            let t_gpu = Instant::now();
            let res_gpu = exact::solve_f2_sidon_exact(n, 7, false, Some(&gpu))?;
            let dt_gpu = t_gpu.elapsed();
            assert_eq!(res_cpu.max_card, res_gpu.max_card);
            assert_eq!(res_cpu.total_nodes, res_gpu.total_nodes);
            println!(
                "Exact F_2({n}) = {} ({} nodes): CPU = {:.3}s, GPU = {:.3}s",
                res_gpu.max_card,
                res_gpu.total_nodes,
                dt_cpu.as_secs_f64(),
                dt_gpu.as_secs_f64()
            );
            Ok(())
        }
    }
}
