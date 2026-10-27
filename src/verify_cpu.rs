//! ============================================================================
//! `verify_cpu.rs` — Exact Sidon ($B_2$) verifiers, 64-bit bit-shift difference
//! masks, fast bitset Mian–Chowla greedy generator (OEIS A005282), and exact
//! Harmonic-Integral Conservation & Dyadic Difference-Budget evaluators for
//! Erdős Problem #1191 ($1,000 open prize: $\liminf_{x\to\infty} A(x)\sqrt{(\ln x)/x} = 0$).
//! ============================================================================
//!
//! MATHEMATICAL & ALGORITHMIC SUMMARY
//! ----------------------------------
//! 1. **$O(1)$-Per-Element 64-Bit Bit-Shift Sidon Contract:**
//!    For a subset $A \subseteq \{1, \dots, 63\}$ represented by a 64-bit word
//!    `chosen_mask` (bit `x` set iff `x ∈ A`), maintain:
//!      * `rev_mask`: bit `63 - a` set for each `a ∈ A`,
//!      * `diff_mask`: bit `d` ($1 \le d \le 62$) set for each positive difference
//!        $d = b - a$ ($a < b$ in $A$).
//!
//!    When testing a candidate $x > \max(A)$ in $\{1, \dots, 63\}$, the bitmask
//!    of all new positive differences $\{x - a : a \in A\}$ is simply:
//!
//!        new_diffs = rev_mask >> (63 - x).
//!
//!    Since all $a \in A$ are distinct, the integers $x - a$ are pairwise distinct,
//!    so $A \cup \{x\}$ is a Sidon set if and only if `(new_diffs & diff_mask) == 0`!
//!    This contract is shared bit-for-bit with CUDA kernels `k_dfs_maximin_sidon`
//!    and `k_eval_sidon_dyadic_budget` in `src/kernels.cu`.
//!
//! 2. **Theorem 2.1(a) (Exact Harmonic-Integral Conservation Identity):**
//!    For any Sidon set $A = \{a_1 < a_2 < \dots\} \subseteq \mathbb{N}$ and any
//!    horizon $N \ge 2$, let $m = A(N) = |A \cap [1, N]|$ and
//!    $\mathcal{D}_A(N) = \{a_j - a_i : 1 \le i < j \le m\}$. Since $A(t)$ is
//!    constant on $[a_j, a_{j+1})$, integrating $A(t)(A(t)-1)/(2 t^2)$ over
//!    $t \in [1, N]$ yields the exact identity:
//!
//!        \mathcal{E}_A(N) := \int_1^N \frac{A(t)(A(t) - 1)}{2 t^2} dt + \frac{m(m - 1)}{2 N}
//!                          = \sum_{j=2}^m \frac{j - 1}{a_j}.
//!
//!    Adding the non-negative Translation Defect
//!
//!        \Delta_A(N) := \sum_{1 \le i < j \le m} \frac{a_i - 1}{(a_j - a_i + 1) a_j} >= 0
//!
//!    gives the exact Harmonic Difference Budget:
//!
//!        \mathcal{E}_A(N) + \Delta_A(N) = \mathcal{B}_A(N) := \sum_{d \in \mathcal{D}_A(N)} \frac{1}{d + 1} <= H_N - 1.

use crate::records::{erdos1191_phi, harmonic_number};

/// Exact Sidon ($B_2$) test for a strictly increasing slice `set`:
/// returns `true` iff all positive differences $a_j - a_i$ ($i < j$) are distinct.
pub fn is_sidon_cpu(set: &[u64]) -> bool {
    let m = set.len();
    if m <= 2 {
        return true;
    }
    let max_diff = (set[m - 1] - set[0]) as usize;
    if max_diff <= 1_048_576 {
        let mut seen = vec![0u64; (max_diff >> 6) + 1];
        for i in 0..m {
            let ai = set[i];
            for &aj in &set[(i + 1)..m] {
                let d = (aj - ai) as usize;
                let word = d >> 6;
                let bit = 1u64 << (d & 63);
                if (seen[word] & bit) != 0 {
                    return false;
                }
                seen[word] |= bit;
            }
        }
        true
    } else {
        let mut diffs = Vec::with_capacity(m * (m - 1) / 2);
        for i in 0..m {
            for j in (i + 1)..m {
                diffs.push(set[j] - set[i]);
            }
        }
        diffs.sort_unstable();
        diffs.windows(2).all(|w| w[0] != w[1])
    }
}

/// Count the number of duplicate positive difference collisions
/// $|\{(i,j,k,l) : i < j, k < l, (i,j) < (k,l), a_j - a_i = a_l - a_k\}|$ in `set`.
pub fn sidon_violation_count_cpu(set: &[u64]) -> u64 {
    let m = set.len();
    if m <= 2 {
        return 0;
    }
    let mut diffs = Vec::with_capacity(m * (m - 1) / 2);
    for i in 0..m {
        for j in (i + 1)..m {
            diffs.push(set[j] - set[i]);
        }
    }
    diffs.sort_unstable();
    let mut collisions = 0u64;
    let mut idx = 0;
    while idx < diffs.len() {
        let mut end = idx + 1;
        while end < diffs.len() && diffs[end] == diffs[idx] {
            end += 1;
        }
        let c = (end - idx) as u64;
        if c > 1 {
            collisions += c * (c - 1) / 2;
        }
        idx = end;
    }
    collisions
}

/// Fast 64-bit bit-shift Sidon test and positive difference-mask computation on
/// subsets of $\{1, \dots, 63\}$.
///
/// Returns `(is_sidon, diff_mask)`, where bit $d \in \{1, \dots, 62\}$ of `diff_mask`
/// is 1 iff $d \in (A - A) \cap \mathbb{N}$.
/// Contract shared with CUDA `k_dfs_maximin_sidon` and `k_eval_sidon_dyadic_budget`.
#[inline]
pub fn sidon_diff_mask64(mask: u64) -> (bool, u64) {
    let mut m = mask & !1u64; // ignore bit 0
    let mut rev_mask = 0u64;
    let mut diff_mask = 0u64;
    let mut ok = true;
    while m != 0 {
        let x = m.trailing_zeros();
        m &= m - 1;
        let new_diffs = rev_mask >> (63 - x);
        if (new_diffs & diff_mask) != 0 {
            ok = false;
        }
        diff_mask |= new_diffs;
        rev_mask |= 1u64 << (63 - x);
    }
    (ok, diff_mask)
}

/// Fast 64-bit bitmask test for Sidon ($B_2$) validity on subsets of $\{1, \dots, 63\}$.
#[inline]
pub fn is_sidon_mask64(mask: u64) -> bool {
    sidon_diff_mask64(mask).0
}

/// CPU twin of `k_eval_sidon_dyadic_budget` (`src/kernels.cu`).
/// Given a 64-bit bitmask `mask` representing $A \subseteq \{1, \dots, 2^M\}$ ($1 \le M \le 5$,
/// or $M = 6$ restricted to $\{1, \dots, 63\}$), returns:
///   * `is_sidon`: `1` if `mask` is a valid Sidon set, `0` otherwise,
///   * `diff_mask`: 64-bit positive difference mask,
///   * `scaled_dyadic_lhs`: exact integer $\sum_{k=1}^M 2^{M - k} m_k(m_k - 1)$.
#[inline]
pub fn eval_sidon_dyadic_budget_cpu(mask: u64, m_scales: u32) -> (u32, u64, u64) {
    let (ok, diff_mask) = sidon_diff_mask64(mask);
    let mut scaled_lhs = 0u64;
    for k in 1..=m_scales {
        let lo = (1u32 << (k - 1)) + 1;
        let hi = (1u32 << k).min(63);
        let shell_mask = if lo > hi {
            0u64
        } else if hi == 63 {
            (!0u64) << lo
        } else {
            ((1u64 << (hi + 1)) - 1) & !((1u64 << lo) - 1)
        };
        let mk = (mask & shell_mask).count_ones() as u64;
        scaled_lhs += (mk * mk.saturating_sub(1)) << (m_scales - k);
    }
    (u32::from(ok), diff_mask, scaled_lhs)
}

/// CPU twin of `k_profile_sidon_cross_shell` (`src/kernels.cu`) for Phase 2
/// (Theorem 2.2 & Theorem 3.2).
/// Given a 64-bit bitmask `mask` representing $A \subseteq \{1, \dots, 2^M\}$
/// ($1 \le M \le 5$, or $M = 6$ restricted to $\{1, \dots, 63\}$), returns:
///   * `is_sidon`: `1` if `mask` is a valid Sidon set, `0` otherwise,
///   * `intra_pairs`: $\sum_{k=1}^M \binom{m_k}{2}$,
///   * `inter_pairs`: $\sum_{k=1}^M m_k A(2^{k-1})$,
///   * `scaled_bilinear_lhs`: $\sum_{k=1}^M 2^{M - k} [m_k(m_k - 1) + 2 m_k A(2^{k-1})]$.
#[inline]
pub fn profile_sidon_cross_shell_cpu(mask: u64, m_scales: u32) -> (u32, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, _diff_mask) = sidon_diff_mask64(clean);
    let m0 = (clean >> 1) & 1u64;
    let mut prefix_prev = m0;
    let mut intra_pairs = 0u64;
    let mut inter_pairs = 0u64;
    let mut scaled_bilinear = 0u64;

    for k in 1..=m_scales {
        let lo = (1u32 << (k - 1)) + 1;
        let hi = (1u32 << k).min(63);
        let shell_mask = if lo > hi {
            0u64
        } else if hi == 63 {
            (!0u64) << lo
        } else {
            ((1u64 << (hi + 1)) - 1) & !((1u64 << lo) - 1)
        };
        let mk = (clean & shell_mask).count_ones() as u64;
        let binom_mk = mk * mk.saturating_sub(1) / 2;
        let cross_mk = mk * prefix_prev;
        intra_pairs += binom_mk;
        inter_pairs += cross_mk;
        scaled_bilinear += (mk * mk.saturating_sub(1) + 2 * cross_mk) << (m_scales - k);
        prefix_prev += mk;
    }

    (u32::from(ok), intra_pairs, inter_pairs, scaled_bilinear)
}

/// CPU twin of `k_eval_sidon_subwindow_spectrum` (`src/kernels.cu`) for Phase 3
/// (Theorem 2.3 & Theorem 3.3).
/// Given a 64-bit bitmask `mask` representing $A \subseteq \{1, \dots, 2^M\}$
/// and a sub-window scale `r_target` ($1 \le r_{\mathrm{target}} \le M$), returns:
///   * `is_sidon`: `1` if `mask` is a valid Sidon set, `0` otherwise,
///   * `pigeonhole_lb`: $\binom{A(2^r)}{2} + \sum_{k=r+1}^M L(m_k, 2^{k-1-r})$,
///   * `exact_subwindow_lhs`: $\binom{A(2^r)}{2} + \sum_{k=r+1}^M \sum_{b=1}^{2^{k-1-r}} \binom{|A \cap J_{k,r,b}|}{2}$,
///   * `short_diffs_rhs`: $|\mathcal{D}_A(2^M) \cap [1, 2^r - 1]|$.
#[inline]
pub fn eval_sidon_subwindow_spectrum_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
) -> (u32, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, diff_mask) = sidon_diff_mask64(clean);
    let r_hi = (1u32 << r_target).min(63);
    let prefix_mask = if r_hi == 63 {
        !1u64
    } else {
        ((1u64 << (r_hi + 1)) - 1) & !1u64
    };
    let a_2r = (clean & prefix_mask).count_ones() as u64;
    let base_pairs = a_2r * a_2r.saturating_sub(1) / 2;

    let mut pigeon_sum = base_pairs;
    let mut exact_sub_sum = base_pairs;
    let w = 1u32 << r_target;

    for k in (r_target + 1)..=m_scales {
        let num_bins = 1u32 << (k - 1 - r_target);
        let shell_start = (1u32 << (k - 1)) + 1;
        let mut mk = 0u64;
        for b in 0..num_bins {
            let lo = shell_start + b * w;
            let hi = (lo + w - 1).min(63);
            if lo > 63 {
                break;
            }
            let sub_mask = if hi == 63 {
                (!0u64) << lo
            } else {
                ((1u64 << (hi + 1)) - 1) & !((1u64 << lo) - 1)
            };
            let cb = (clean & sub_mask).count_ones() as u64;
            mk += cb;
            exact_sub_sum += cb * cb.saturating_sub(1) / 2;
        }
        let b_u64 = num_bins as u64;
        let q = mk / b_u64;
        let t = mk % b_u64;
        let lb = b_u64 * (q * q.saturating_sub(1) / 2) + t * q;
        pigeon_sum += lb;
    }

    let max_d = ((1u32 << r_target) - 1).min(62);
    let d_mask_win = ((1u64 << (max_d + 1)) - 1) & !1u64;
    let short_diffs = (diff_mask & d_mask_win).count_ones() as u64;

    (u32::from(ok), pigeon_sum, exact_sub_sum, short_diffs)
}

#[inline]
fn eval_refined_block_cpu(
    clean_mask: u64,
    n_hi: u32,
    w: u32,
    u_start: u32,
    u_end: u32,
    levels: u32,
) -> u64 {
    if u_start > u_end {
        return 0;
    }
    let mut segs = vec![(u_start, u_end)];
    let eff_levels = levels.min(4);
    for _ in 0..eff_levels {
        let mut next_segs = Vec::with_capacity(segs.len() * 2);
        for (a, b) in segs {
            if b > a {
                let mid = (a + b) >> 1;
                next_segs.push((a, mid));
                next_segs.push((mid + 1, b));
            } else {
                next_segs.push((a, b));
            }
        }
        segs = next_segs;
    }

    let mut total_lb = 0u64;
    for (a, b) in segs {
        let mut s_blk = 0u64;
        for u in a..=b {
            let lo_x = if u >= w { u - w + 1 } else { 1 };
            let hi_x = u.min(n_hi);
            if lo_x <= hi_x {
                let hi_m = if hi_x >= 63 {
                    !0u64
                } else {
                    (1u64 << (hi_x + 1)) - 1
                };
                let lo_m = (1u64 << lo_x) - 1;
                s_blk += (clean_mask & (hi_m & !lo_m)).count_ones() as u64;
            }
        }
        let bins = (b - a + 1) as u64;
        if bins > 0 && s_blk > bins {
            let q = s_blk / bins;
            let t = s_blk % bins;
            total_lb += bins * (q * q.saturating_sub(1) / 2) + t * q;
        }
    }
    total_lb
}

/// CPU twin of `k_eval_base_b_joint_fejer` (`src/kernels.cu`) for Phase 5
/// (Theorem 2.5: Direct Prefix-to-Block Unsplit Fejér Coupling & Base-$B = 2^{1/L}$
/// Window-Block Refinement Hierarchy).
///
/// Given a 64-bit bitmask `mask` representing $A \subseteq \{1, \dots, 2^M\}$,
/// a window scale `r_target` ($W = 2^{r_{\mathrm{target}}}$), and `refinement_levels`
/// ($\ell \in \{0, \dots, 4\}$, where $L = 2^\ell$ and $B = 2^{1/L}$), returns:
///   * `is_sidon`: `1` if `mask` is a valid Sidon set, `0` otherwise,
///   * `joint_l1_lhs`: $\mathcal{F}_W^{\mathrm{joint}}(A; B=2)$ on unsplit blocks,
///   * `joint_base_b_lhs`: $\mathcal{F}_W^{\mathrm{joint}}(A; B=2^{1/L})$ after $L$-fold subdivision,
///   * `exact_fejer_rhs`: $\sum_{d \in \mathcal{D}_A(2^M) \cap [1, W - 1]} (W - d) \le \binom{W}{2}$.
#[inline]
pub fn eval_base_b_joint_fejer_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, diff_mask) = sidon_diff_mask64(clean);
    let n_hi = (1u32 << m_scales).min(63);
    let w = 1u32 << r_target;
    let t_r = (1u32 << r_target).min(n_hi);

    let mut l1_sum = eval_refined_block_cpu(clean, n_hi, w, 1, t_r + w - 1, 0);
    let mut lb_sum = eval_refined_block_cpu(clean, n_hi, w, 1, t_r + w - 1, refinement_levels);

    for k in (r_target + 1)..=m_scales {
        let t_prev = 1u32 << (k - 1);
        let t_curr = (1u32 << k).min(n_hi);
        if t_prev >= n_hi {
            break;
        }
        let u_start = t_prev + w;
        let u_end = t_curr + w - 1;
        l1_sum += eval_refined_block_cpu(clean, n_hi, w, u_start, u_end, 0);
        lb_sum += eval_refined_block_cpu(clean, n_hi, w, u_start, u_end, refinement_levels);
    }

    let mut fejer_rhs = 0u64;
    for d in 1..w.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            fejer_rhs += (w - d) as u64;
        }
    }

    (u32::from(ok), l1_sum, lb_sum, fejer_rhs)
}

#[inline]
fn eval_refined_block_integrality_cpu(
    clean_mask: u64,
    n_hi: u32,
    w: u32,
    u_start: u32,
    u_end: u32,
    refinement_levels: u32,
) -> (u64, u64, u64) {
    if u_start > u_end {
        return (0, 0, 0);
    }
    let mut num_sub = 1u32 << refinement_levels;
    let total_len = u_end - u_start + 1;
    if num_sub > total_len {
        num_sub = total_len;
    }
    let mut cont_acc = 0u64;
    let mut int_acc = 0u64;
    let mut l_acc = 0u64;
    for s in 0..num_sub {
        let sub_a = u_start + (s * total_len) / num_sub;
        let sub_b = u_start + ((s + 1) * total_len) / num_sub - 1;
        if sub_a > sub_b {
            continue;
        }
        let bins = (sub_b - sub_a + 1) as u64;
        let mut s_val = 0u64;
        for u in sub_a..=sub_b {
            let hi = u.min(n_hi);
            let lo = u.saturating_sub(w);
            let hi_mask = if hi == 63 {
                !1u64
            } else {
                ((1u64 << (hi + 1)) - 1) & !1u64
            };
            let lo_mask = if lo == 0 {
                0u64
            } else {
                ((1u64 << (lo + 1)) - 1) & !1u64
            };
            let win_mask = clean_mask & (hi_mask & !lo_mask);
            s_val += win_mask.count_ones() as u64;
        }
        if s_val > bins {
            let q = s_val / bins;
            let t = s_val % bins;
            cont_acc += s_val * (s_val - bins);
            int_acc += t * (bins - t);
            l_acc += bins * (q * q.saturating_sub(1) / 2) + t * q;
        }
    }
    (cont_acc, int_acc, l_acc)
}

/// CPU twin of `k_eval_coupled_lag_integrality_fejer` (`src/kernels.cu`) for Phase 6
/// (Theorem 2.6: Exact Discrete Legendre–Fenchel Integrality Excess & Coupled Unsplit Fejér).
///
/// Returns `(is_sidon, scaled_cont_parabola, scaled_integrality_excess, exact_pigeonhole_lhs, exact_fejer_rhs)`.
#[inline]
pub fn eval_coupled_lag_integrality_fejer_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, diff_mask) = sidon_diff_mask64(clean);
    let n_hi = (1u32 << m_scales).min(63);
    let w = 1u32 << r_target;
    let t_r = (1u32 << r_target).min(n_hi);

    let (mut cont_sum, mut int_sum, mut l_sum) =
        eval_refined_block_integrality_cpu(clean, n_hi, w, 1, t_r + w - 1, refinement_levels);

    for k in (r_target + 1)..=m_scales {
        let t_prev = 1u32 << (k - 1);
        let t_curr = (1u32 << k).min(n_hi);
        if t_prev >= n_hi {
            break;
        }
        let u_start = t_prev + w;
        let u_end = t_curr + w - 1;
        let (c_k, i_k, l_k) =
            eval_refined_block_integrality_cpu(clean, n_hi, w, u_start, u_end, refinement_levels);
        cont_sum += c_k;
        int_sum += i_k;
        l_sum += l_k;
    }

    let mut fejer_rhs = 0u64;
    for d in 1..w.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            fejer_rhs += (w - d) as u64;
        }
    }

    (u32::from(ok), cont_sum, int_sum, l_sum, fejer_rhs)
}

/// CPU twin of `k_eval_cesaro_fejer_order2` (`src/kernels.cu`) for Phase 7
/// (Theorem 2.7: Higher-Order Cesàro–Fejér Polynomial Kernels $\psi_p(d) = \binom{W - d + p - 1}{p}$
/// via Iterated Sliding-Window Sums).
///
/// Returns `(is_sidon, cf2_full_pigeonhole_lhs, cf2_oct_pigeonhole_lhs, cf2_full_exact_rhs, cf2_oct_exact_rhs)`.
#[inline]
pub fn eval_cesaro_fejer_order2_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, diff_mask) = sidon_diff_mask64(clean);
    let n_hi = (1u32 << m_scales).min(63);
    let w_cap = 1u32 << r_target;
    let w_half = w_cap >> 1;
    let t_r = (1u32 << r_target).min(n_hi);

    let mut full_lhs = 0u64;
    let mut oct_lhs = 0u64;

    for w in 1..=w_cap {
        let (_, _, mut l_w) =
            eval_refined_block_integrality_cpu(clean, n_hi, w, 1, t_r + w - 1, refinement_levels);
        for k in (r_target + 1)..=m_scales {
            let t_prev = 1u32 << (k - 1);
            let t_curr = (1u32 << k).min(n_hi);
            if t_prev >= n_hi {
                break;
            }
            let u_start = t_prev + w;
            let u_end = t_curr + w - 1;
            let (_, _, l_k) = eval_refined_block_integrality_cpu(
                clean,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
            );
            l_w += l_k;
        }
        full_lhs += l_w;
        if w > w_half {
            oct_lhs += l_w;
        }
    }

    let mut full_rhs = 0u64;
    let mut half_rhs = 0u64;
    for d in 1..w_cap.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            let rem_full = (w_cap - d) as u64;
            full_rhs += rem_full * (rem_full + 1) / 2;
            if d < w_half {
                let rem_half = (w_half - d) as u64;
                half_rhs += rem_half * (rem_half + 1) / 2;
            }
        }
    }

    (
        u32::from(ok),
        full_lhs,
        oct_lhs,
        full_rhs,
        full_rhs - half_rhs,
    )
}

/// CPU twin of `k_eval_band_restricted_cesaro_fejer` (`src/kernels.cu`) for Phase 8
/// (Theorem 2.8: Band-Restricted Upper-Octave Cesàro–Fejér Identity $\mathcal{CF}_{W, W_0}(A)$).
///
/// Returns `(is_sidon, band_pigeonhole_lhs, band_scaled_integrality, band_exact_rhs, band_binom_cap)`.
#[inline]
pub fn eval_band_restricted_cesaro_fejer_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    w0_cutoff: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let clean = mask & !1u64;
    let (ok, diff_mask) = sidon_diff_mask64(clean);
    let n_hi = (1u32 << m_scales).min(63);
    let w_cap = 1u32 << r_target;
    let w0 = w0_cutoff.min(w_cap - 1);
    let t_r = (1u32 << r_target).min(n_hi);

    let mut band_lhs = 0u64;
    let mut band_int = 0u64;

    for w in (w0 + 1)..=w_cap {
        let (_, mut i_w, mut l_w) =
            eval_refined_block_integrality_cpu(clean, n_hi, w, 1, t_r + w - 1, refinement_levels);
        for k in (r_target + 1)..=m_scales {
            let t_prev = 1u32 << (k - 1);
            let t_curr = (1u32 << k).min(n_hi);
            if t_prev >= n_hi {
                break;
            }
            let u_start = t_prev + w;
            let u_end = t_curr + w - 1;
            let (_, i_k, l_k) = eval_refined_block_integrality_cpu(
                clean,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
            );
            i_w += i_k;
            l_w += l_k;
        }
        band_lhs += l_w;
        band_int += i_w;
    }

    let mut full_rhs = 0u64;
    let mut sub_rhs = 0u64;
    for d in 1..w_cap.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            let rem_full = (w_cap - d) as u64;
            full_rhs += rem_full * (rem_full + 1) / 2;
            if d < w0 {
                let rem_w0 = (w0 - d) as u64;
                sub_rhs += rem_w0 * (rem_w0 + 1) / 2;
            }
        }
    }

    let w_u64 = w_cap as u64;
    let w0_u64 = w0 as u64;
    let cap_full = (w_u64 - 1) * w_u64 * (w_u64 + 1) / 6;
    let cap_w0 = if w0_u64 >= 2 {
        (w0_u64 - 1) * w0_u64 * (w0_u64 + 1) / 6
    } else {
        0
    };

    (
        u32::from(ok),
        band_lhs,
        band_int,
        full_rhs - sub_rhs,
        cap_full - cap_w0,
    )
}

/// CPU twin of `k_eval_multioctave_cesaro_fejer` (`src/kernels.cu`) for Phase 9
/// (Theorem 2.9: Scale-Invariant Block-Window Conservation & $K$-Octave Integrality Cascade
/// on $w \in [\lfloor 2^{-K} W \rfloor + 1, W]$).
///
/// Returns `(is_sidon, multioct_pigeonhole_lhs, multioct_scaled_integrality, multioct_exact_rhs, multioct_binom_cap)`.
#[inline]
pub fn eval_multioctave_cesaro_fejer_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let w_cap = 1u32 << r_target;
    let shift_k = k_octaves.min(r_target);
    let w0 = if shift_k >= r_target {
        0
    } else {
        (w_cap >> shift_k).min(w_cap - 1)
    };
    eval_band_restricted_cesaro_fejer_cpu(mask, m_scales, r_target, w0, refinement_levels)
}

/// CPU twin of `k_eval_haar_multioctave_cesaro_fejer` (`src/kernels.cu`) for Phase 10
/// (Theorem 2.10: Scale-Invariant Haar-Measure Multi-Octave Cesàro–Fejér Cascade with
/// integer Haar weights $\lfloor W / w \rfloor$ on $w \in [\lfloor 2^{-K} W \rfloor + 1, W]$).
///
/// Returns `(is_sidon, haar_pigeonhole_lhs, haar_scaled_integrality, haar_exact_rhs, haar_binom_cap)`.
#[inline]
pub fn eval_haar_multioctave_cesaro_fejer_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let mask = mask & !1u64;
    let mut m = mask;
    let mut rev_mask = 0u64;
    let mut diff_mask = 0u64;
    let mut ok = true;

    while m != 0 {
        let x = m.trailing_zeros();
        m &= m - 1;
        let new_diffs = rev_mask >> (63 - x);
        if (new_diffs & diff_mask) != 0 {
            ok = false;
        }
        diff_mask |= new_diffs;
        rev_mask |= 1u64 << (63 - x);
    }

    let n_hi = (1u32 << m_scales).min(63);
    let w_cap = 1u32 << r_target;
    let shift_k = k_octaves.min(r_target);
    let w0 = if shift_k >= r_target {
        0
    } else {
        (w_cap >> shift_k).min(w_cap - 1)
    };
    let t_r = (1u32 << r_target).min(n_hi);

    let mut haar_lhs = 0u64;
    let mut haar_int = 0u64;
    let mut haar_cap = 0u64;

    for w in (w0 + 1)..=w_cap {
        let wt = (w_cap / w) as u64;
        let (_, mut i_w, mut l_w) =
            eval_refined_block_integrality_cpu(mask, n_hi, w, 1, t_r + w - 1, refinement_levels);
        for k in (r_target + 1)..=m_scales {
            let t_prev = 1u32 << (k - 1);
            let t_curr = (1u32 << k).min(n_hi);
            if t_prev >= n_hi {
                break;
            }
            let u_start = t_prev + w;
            let u_end = t_curr + w - 1;
            let (_, i_k, l_k) = eval_refined_block_integrality_cpu(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
            );
            i_w += i_k;
            l_w += l_k;
        }
        haar_lhs += wt * l_w;
        haar_int += wt * i_w;
        let w_u64 = w as u64;
        haar_cap += wt * (w_u64 * (w_u64 - 1) / 2);
    }

    let mut haar_rhs = 0u64;
    for d in 1..w_cap.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            let w_start = (w0 + 1).max(d + 1);
            let mut psi_d = 0u64;
            for w in w_start..=w_cap {
                let wt = (w_cap / w) as u64;
                psi_d += wt * ((w - d) as u64);
            }
            haar_rhs += psi_d;
        }
    }

    (u32::from(ok), haar_lhs, haar_int, haar_rhs, haar_cap)
}

/// CPU twin of `k_eval_unconditional_diophantine_haar_closure` (`src/kernels.cu`) for Phase 11
/// (Theorem 2.11: Unconditional Diophantine Scaled-Reciprocal Haar-Measure Multi-Octave
/// Cesàro–Fejér Closure with integer weights $\lfloor W^2 / w \rfloor$ on $w \in [\lfloor 2^{-K} W \rfloor + 1, W]$).
///
/// Returns `(is_sidon, dioph_pigeonhole_lhs, dioph_scaled_integrality, dioph_exact_rhs, dioph_binom_cap)`.
#[inline]
pub fn eval_unconditional_diophantine_haar_closure_cpu(
    mask: u64,
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> (u32, u64, u64, u64, u64) {
    let mask = mask & !1u64;
    let mut m = mask;
    let mut rev_mask = 0u64;
    let mut diff_mask = 0u64;
    let mut ok = true;

    while m != 0 {
        let x = m.trailing_zeros();
        m &= m - 1;
        let new_diffs = rev_mask >> (63 - x);
        if (new_diffs & diff_mask) != 0 {
            ok = false;
        }
        diff_mask |= new_diffs;
        rev_mask |= 1u64 << (63 - x);
    }

    let n_hi = (1u32 << m_scales).min(63);
    let w_cap = 1u32 << r_target;
    let w_sq = (w_cap as u64) * (w_cap as u64);
    let shift_k = k_octaves.min(r_target);
    let w0 = if shift_k >= r_target {
        0
    } else {
        (w_cap >> shift_k).min(w_cap - 1)
    };
    let t_r = (1u32 << r_target).min(n_hi);

    let mut dioph_lhs = 0u64;
    let mut dioph_int = 0u64;
    let mut dioph_cap = 0u64;

    for w in (w0 + 1)..=w_cap {
        let wt = w_sq / (w as u64);
        let (_, mut i_w, mut l_w) =
            eval_refined_block_integrality_cpu(mask, n_hi, w, 1, t_r + w - 1, refinement_levels);
        for k in (r_target + 1)..=m_scales {
            let t_prev = 1u32 << (k - 1);
            let t_curr = (1u32 << k).min(n_hi);
            if t_prev >= n_hi {
                break;
            }
            let u_start = t_prev + w;
            let u_end = t_curr + w - 1;
            let (_, i_k, l_k) = eval_refined_block_integrality_cpu(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
            );
            i_w += i_k;
            l_w += l_k;
        }
        dioph_lhs += wt * l_w;
        dioph_int += wt * i_w;
        let w_u64 = w as u64;
        dioph_cap += wt * (w_u64 * (w_u64 - 1) / 2);
    }

    let mut dioph_rhs = 0u64;
    for d in 1..w_cap.min(63) {
        if (diff_mask & (1u64 << d)) != 0 {
            let w_start = (w0 + 1).max(d + 1);
            let mut psi_d = 0u64;
            for w in w_start..=w_cap {
                let wt = w_sq / (w as u64);
                psi_d += wt * ((w - d) as u64);
            }
            dioph_rhs += psi_d;
        }
    }

    (u32::from(ok), dioph_lhs, dioph_int, dioph_rhs, dioph_cap)
}

/// Generate the first `num_terms` terms of the Mian–Chowla greedy infinite Sidon
/// sequence (OEIS `A005282`) using a dynamically growing difference bitset.
pub fn mian_chowla_greedy_cpu(num_terms: usize) -> Vec<u64> {
    if num_terms == 0 {
        return Vec::new();
    }
    let mut terms = Vec::with_capacity(num_terms);
    terms.push(1u64);
    if num_terms == 1 {
        return terms;
    }

    // Estimate upper bound on a_{num_terms} ~ 0.035 * n^3 + 2048
    let cap_diff = (num_terms * num_terms * num_terms / 16).max(16_384);
    let mut used_diffs = vec![0u64; (cap_diff >> 6) + 1];

    let is_diff_used = |words: &[u64], d: usize| -> bool {
        let w = d >> 6;
        if w >= words.len() {
            false
        } else {
            (words[w] & (1u64 << (d & 63))) != 0
        }
    };

    let set_diff_used = |words: &mut Vec<u64>, d: usize| {
        let w = d >> 6;
        if w >= words.len() {
            words.resize(w + 1024, 0u64);
        }
        words[w] |= 1u64 << (d & 63);
    };

    let mut cand = 2u64;
    while terms.len() < num_terms {
        let mut ok = true;
        // Iterate from largest existing term to smallest (smallest differences first,
        // which are the most densely populated and fail fastest!)
        for &a in terms.iter().rev() {
            let d = (cand - a) as usize;
            if is_diff_used(&used_diffs, d) {
                ok = false;
                break;
            }
        }
        if ok {
            for &a in &terms {
                let d = (cand - a) as usize;
                set_diff_used(&mut used_diffs, d);
            }
            terms.push(cand);
        }
        cand += 1;
    }

    terms
}

/// Complete report for Theorem 2.1(a) (Harmonic-Integral Conservation Identity)
/// and the Erdős Problem #1191 functional $\Phi_A(x) = A(x)\sqrt{(\ln x)/x}$ on $[1, N]$.
#[derive(Debug, Clone)]
pub struct HarmonicConservationReport {
    pub n: usize,
    pub card: usize,
    pub is_sidon: bool,
    pub violations: u64,
    pub integral_term: f64,
    pub boundary_term: f64,
    pub prefix_energy_e: f64,
    pub sum_j_minus_1_over_aj: f64,
    pub translation_defect_delta: f64,
    pub harmonic_diff_budget_b: f64,
    pub harmonic_ceiling_h_n_minus_1: f64,
    pub conservation_residual: f64,
    pub phi_at_n: f64,
    pub min_phi_from_x0: f64,
    pub min_phi_x: usize,
    pub max_phi_from_x0: f64,
    pub max_phi_x: usize,
}

/// Compute the exact Theorem 2.1(a) Harmonic-Integral Conservation decomposition
/// and Erdős #1191 envelope statistics over $x \in [x_0, n]$ for a sorted set `set`.
pub fn harmonic_conservation_report_cpu(
    set: &[u64],
    n: usize,
    x0: usize,
) -> HarmonicConservationReport {
    assert!(n >= 2);
    let in_n: Vec<u64> = set
        .iter()
        .copied()
        .filter(|&x| x >= 1 && (x as usize) <= n)
        .collect();
    let m = in_n.len();
    let is_sidon = is_sidon_cpu(&in_n);
    let violations = sidon_violation_count_cpu(&in_n);

    // 1. Exact piecewise-constant integral \int_1^N \frac{A(t)(A(t) - 1)}{2 t^2} dt:
    //    On [a_k, a_{k+1}) (with a_{m+1} := N), A(t) = k, and \int_{u}^{v} dt/t^2 = 1/u - 1/v.
    let mut integral_term = 0.0f64;
    if m >= 2 {
        for k in 2..=m {
            let u = in_n[k - 1] as f64;
            let v = if k < m { in_n[k] as f64 } else { n as f64 };
            if v > u {
                let binom_k = ((k * (k - 1)) as f64) * 0.5;
                integral_term += binom_k * (1.0 / u - 1.0 / v);
            }
        }
    }
    let boundary_term = ((m * m.saturating_sub(1)) as f64) / (2.0 * (n as f64));
    let prefix_energy_e = integral_term + boundary_term;

    // 2. Direct sum \sum_{j=2}^m (j - 1) / a_j
    let mut sum_j_minus_1_over_aj = 0.0f64;
    for j in 2..=m {
        sum_j_minus_1_over_aj += ((j - 1) as f64) / (in_n[j - 1] as f64);
    }

    // 3. Translation Defect \Delta_A(N) and Harmonic Difference Budget \mathcal{B}_A(N)
    let mut translation_defect_delta = 0.0f64;
    let mut harmonic_diff_budget_b = 0.0f64;
    for i in 0..m {
        let ai = in_n[i] as f64;
        for &aj_u64 in &in_n[(i + 1)..m] {
            let aj = aj_u64 as f64;
            let d = aj - ai;
            translation_defect_delta += (ai - 1.0) / ((d + 1.0) * aj);
            harmonic_diff_budget_b += 1.0 / (d + 1.0);
        }
    }

    let harmonic_ceiling_h_n_minus_1 = harmonic_number(n) - 1.0;
    let conservation_residual =
        (prefix_energy_e + translation_defect_delta - harmonic_diff_budget_b).abs();

    // 4. Evaluate \Phi_A(x) = A(x) \sqrt{(\ln x)/x} across [x0, n]
    let x_start = x0.clamp(2, n);
    let mut min_phi_from_x0 = f64::INFINITY;
    let mut min_phi_x = x_start;
    let mut max_phi_from_x0 = 0.0f64;
    let mut max_phi_x = x_start;

    let mut ptr = 0usize;
    while ptr < m && (in_n[ptr] as usize) < x_start {
        ptr += 1;
    }
    for x in x_start..=n {
        while ptr < m && (in_n[ptr] as usize) <= x {
            ptr += 1;
        }
        let phi = erdos1191_phi(ptr, x);
        if phi < min_phi_from_x0 {
            min_phi_from_x0 = phi;
            min_phi_x = x;
        }
        if phi > max_phi_from_x0 {
            max_phi_from_x0 = phi;
            max_phi_x = x;
        }
    }
    let phi_at_n = erdos1191_phi(m, n);

    HarmonicConservationReport {
        n,
        card: m,
        is_sidon,
        violations,
        integral_term,
        boundary_term,
        prefix_energy_e,
        sum_j_minus_1_over_aj,
        translation_defect_delta,
        harmonic_diff_budget_b,
        harmonic_ceiling_h_n_minus_1,
        conservation_residual,
        phi_at_n,
        min_phi_from_x0,
        min_phi_x,
        max_phi_from_x0,
        max_phi_x,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::A005282_MIAN_CHOWLA_64;

    #[test]
    fn bit_shift_sidon_mask64_matches_reference_verifier() {
        // Test on Mian-Chowla prefix inside {1..63}
        let mc_small: Vec<u64> = A005282_MIAN_CHOWLA_64
            .iter()
            .copied()
            .filter(|&x| x <= 63)
            .collect();
        let mut mask = 0u64;
        for &x in &mc_small {
            mask |= 1u64 << x;
        }
        assert!(is_sidon_cpu(&mc_small));
        assert!(is_sidon_mask64(mask));

        // Adding 3 to {1, 2, 4} creates duplicate difference 3-2 = 2-1 = 1
        let bad_mask = (1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4);
        assert!(!is_sidon_mask64(bad_mask));
        assert!(!is_sidon_cpu(&[1, 2, 3, 4]));
    }

    #[test]
    fn theorem_2_1a_harmonic_conservation_identity_is_exact() {
        let mc = mian_chowla_greedy_cpu(32);
        for &n in &[20usize, 50, 100, 250, 500, 1000] {
            let rep = harmonic_conservation_report_cpu(&mc, n, 4);
            assert!(rep.is_sidon);
            assert_eq!(rep.violations, 0);
            assert!(
                (rep.prefix_energy_e - rep.sum_j_minus_1_over_aj).abs() < 1e-12,
                "integral identity mismatch at n={n}"
            );
            assert!(
                rep.conservation_residual < 1e-12,
                "conservation residual {} at n={n}",
                rep.conservation_residual
            );
            assert!(
                rep.harmonic_diff_budget_b <= rep.harmonic_ceiling_h_n_minus_1 + 1e-12,
                "harmonic ceiling exceeded at n={n}"
            );
        }
    }
}
