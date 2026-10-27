//! ============================================================================
//! `exact.rs` — Hybrid CPU + GPU exact branch-and-bound solver for finite Sidon
//! capacity $F_2(N)$ (OEIS A143824 / A003022) and the Maximin Lower-Envelope
//! Sidon Problem $\mu^*(x_0, N) = \max_{A \subseteq [N] \text{ Sidon}} \min_{x_0 \le x \le N} A(x)\sqrt{(\ln x)/x}$
//! for Erdős Problem #1191 ($1,000 open prize).
//! ============================================================================

use crate::gpu::Gpu;
use crate::records::{erdos1191_phi, f2_sidon_exact, sidon_elementary_upper_bound};
use crate::verify_cpu::{is_sidon_cpu, is_sidon_mask64};
use anyhow::{bail, Result};
use cudarc::driver::PushKernelArg;
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct SidonPrefixState {
    pub chosen_mask: u64,
    pub rev_mask: u64,
    pub diff_mask: u64,
    pub last: u32,
    pub count: u32,
}

#[derive(Debug, Clone)]
pub struct ExactSidonResult {
    pub n: usize,
    pub max_card: u32,
    pub min_phi: f64,
    pub min_phi_x: usize,
    pub witness: Vec<u64>,
    pub total_nodes: u64,
}

/// Build the non-circular subadditive tail-capacity table `sidon_cap[0..=n]` for
/// interior intervals of width `w <= n`. For `w < n`, uses `f2_sidon_exact(w)`
/// (or falls back to the elementary Erdős–Turán difference-counting bound
/// $\lfloor (1 + \sqrt{1 + 8(w - 1)})/2 \rfloor$ when `unconditional` is set and
/// `w > 12`); never assumes $F_2(n)$ at $w = n$.
pub fn build_sidon_cap(n: usize, unconditional: bool) -> Vec<u32> {
    let mut cap = vec![0u32; n + 1];
    for (w, slot) in cap.iter_mut().enumerate().take(n + 1).skip(1) {
        if w < n && (!unconditional || w <= 12) {
            if let Some(f) = f2_sidon_exact(w) {
                *slot = f as u32;
                continue;
            }
        }
        *slot = sidon_elementary_upper_bound(w) as u32;
    }
    if n >= 1 {
        cap[n] = (cap[n - 1] + 1).min(sidon_elementary_upper_bound(n) as u32);
    }
    cap
}

/// Decode a 64-bit bitmask into a sorted `Vec<u64>` of elements in $\{1, \dots, 63\}$.
pub fn mask_to_vec(mask: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut m = mask & !1u64;
    while m != 0 {
        let x = m.trailing_zeros() as u64;
        m &= m - 1;
        out.push(x);
    }
    out
}

/// Encode a slice of elements in $\{1, \dots, 63\}$ into a 64-bit bitmask.
pub fn vec_to_mask(set: &[u64]) -> u64 {
    let mut mask = 0u64;
    for &x in set {
        if (1..=63).contains(&x) {
            mask |= 1u64 << x;
        }
    }
    mask
}

/// Expand valid Sidon prefixes up to element `split_x` (inclusive) that satisfy
/// the lower-envelope floor `req_prefix[1..=split_x]` and can still reach
/// `eff_target = target_card.max(req_prefix[n])` within `{1, ..., n}`.
pub fn expand_sidon_prefixes(
    n: u32,
    split_x: u32,
    target_card: u32,
    sidon_cap: &[u32],
    req_prefix: &[u32],
) -> Vec<SidonPrefixState> {
    let eff_target = target_card.max(req_prefix[n as usize]);
    let mut current = vec![SidonPrefixState {
        chosen_mask: 0,
        rev_mask: 0,
        diff_mask: 0,
        last: 0,
        count: 0,
    }];
    let limit = split_x.min(n);
    for x in 1..=limit {
        let mut next_level = Vec::with_capacity(current.len() * 2);
        let rem_cap = sidon_cap[(n - x) as usize];
        let req_at_x = req_prefix[x as usize];
        for st in current {
            // Branch 1: include x (if O(1) bit-shift Sidon test passes and can reach eff_target)
            let new_diffs = st.rev_mask >> (63 - x);
            if (new_diffs & st.diff_mask) == 0
                && st.count + 1 >= req_at_x
                && st.count + 1 + rem_cap >= eff_target
            {
                next_level.push(SidonPrefixState {
                    chosen_mask: st.chosen_mask | (1u64 << x),
                    rev_mask: st.rev_mask | (1u64 << (63 - x)),
                    diff_mask: st.diff_mask | new_diffs,
                    last: x,
                    count: st.count + 1,
                });
            }
            // Branch 2: skip x (if st.count >= req_prefix[x] and remaining tail can reach eff_target)
            if st.count >= req_at_x && st.count + rem_cap >= eff_target {
                next_level.push(SidonPrefixState {
                    chosen_mask: st.chosen_mask,
                    rev_mask: st.rev_mask,
                    diff_mask: st.diff_mask,
                    last: x,
                    count: st.count,
                });
            }
        }
        current = next_level;
    }
    current
}

/// CPU twin of `k_dfs_maximin_sidon` (`src/kernels.cu`).
/// Must maintain identical traversal order, pruning, and node-count semantics.
pub fn extend_sidon_dfs_cpu(
    prefix: SidonPrefixState,
    n: u32,
    target_card: u32,
    sidon_cap: &[u32],
    req_prefix: &[u32],
) -> (Option<u64>, u64) {
    let eff_target = target_card.max(req_prefix[n as usize]);
    if prefix.count >= eff_target {
        return (Some(prefix.chosen_mask), 0);
    }

    let mut stk_chosen = [0u64; 64];
    let mut stk_rev = [0u64; 64];
    let mut stk_diff = [0u64; 64];
    let mut stk_next_x = [0u32; 64];
    let mut stk_count = [0u32; 64];

    stk_chosen[0] = prefix.chosen_mask;
    stk_rev[0] = prefix.rev_mask;
    stk_diff[0] = prefix.diff_mask;
    stk_next_x[0] = prefix.last + 1;
    stk_count[0] = prefix.count;

    let mut sp: isize = 0;
    let mut nodes = 0u64;
    let mut first_witness: Option<u64> = None;

    while sp >= 0 {
        let idx = sp as usize;
        let x = stk_next_x[idx];
        let cnt = stk_count[idx];

        if x > 1 && cnt < req_prefix[(x - 1) as usize] {
            sp -= 1;
            continue;
        }
        if x > n {
            sp -= 1;
            continue;
        }
        stk_next_x[idx] = x + 1;

        if cnt + sidon_cap[(n - x + 1) as usize] < eff_target {
            sp -= 1;
            continue;
        }
        if cnt + 1 + sidon_cap[(n - x) as usize] < eff_target {
            continue;
        }

        let rev = stk_rev[idx];
        let diff = stk_diff[idx];
        let new_diffs = rev >> (63 - x);
        if (new_diffs & diff) != 0 {
            continue;
        }

        nodes += 1;
        let next_chosen = stk_chosen[idx] | (1u64 << x);
        let next_cnt = cnt + 1;

        if next_cnt >= eff_target {
            if first_witness.is_none() {
                first_witness = Some(next_chosen);
            }
            continue;
        }

        sp += 1;
        let next_idx = sp as usize;
        stk_chosen[next_idx] = next_chosen;
        stk_rev[next_idx] = rev | (1u64 << (63 - x));
        stk_diff[next_idx] = diff | new_diffs;
        stk_next_x[next_idx] = x + 1;
        stk_count[next_idx] = next_cnt;
    }

    (first_witness, nodes)
}

/// Launch `k_dfs_maximin_sidon` on a batch of prefixes on the GPU.
pub fn gpu_extend_sidon_batch(
    gpu: &Gpu,
    prefixes: &[SidonPrefixState],
    n: u32,
    target_card: u32,
    sidon_cap: &[u32],
    req_prefix: &[u32],
) -> Result<(Option<u64>, u64)> {
    if prefixes.is_empty() {
        return Ok((None, 0));
    }
    let num = prefixes.len();
    let h_chosen: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
    let h_rev: Vec<u64> = prefixes.iter().map(|p| p.rev_mask).collect();
    let h_diff: Vec<u64> = prefixes.iter().map(|p| p.diff_mask).collect();
    let h_last: Vec<u32> = prefixes.iter().map(|p| p.last).collect();
    let h_count: Vec<u32> = prefixes.iter().map(|p| p.count).collect();

    let d_chosen = gpu.stream.clone_htod(&h_chosen)?;
    let d_rev = gpu.stream.clone_htod(&h_rev)?;
    let d_diff = gpu.stream.clone_htod(&h_diff)?;
    let d_last = gpu.stream.clone_htod(&h_last)?;
    let d_count = gpu.stream.clone_htod(&h_count)?;
    let d_cap = gpu.stream.clone_htod(sidon_cap)?;
    let d_req = gpu.stream.clone_htod(req_prefix)?;
    let mut d_found = gpu.stream.alloc_zeros::<u32>(1)?;
    let mut d_witness = gpu.stream.alloc_zeros::<u64>(1)?;
    let mut d_nodes = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu.stream.launch_builder(&gpu.fn_dfs_maximin_sidon);
    builder.arg(&d_chosen);
    builder.arg(&d_rev);
    builder.arg(&d_diff);
    builder.arg(&d_last);
    builder.arg(&d_count);
    builder.arg(&num_u32);
    builder.arg(&n);
    builder.arg(&target_card);
    builder.arg(&d_cap);
    builder.arg(&d_req);
    builder.arg(&mut d_found);
    builder.arg(&mut d_witness);
    builder.arg(&mut d_nodes);
    unsafe { builder.launch(cfg) }?;

    let h_found: Vec<u32> = gpu.stream.clone_dtoh(&d_found)?;
    let h_witness: Vec<u64> = gpu.stream.clone_dtoh(&d_witness)?;
    let h_nodes: Vec<u64> = gpu.stream.clone_dtoh(&d_nodes)?;

    let total_nodes: u64 = h_nodes.into_iter().sum();
    let witness = if h_found[0] != 0 {
        Some(h_witness[0])
    } else {
        None
    };
    Ok((witness, total_nodes))
}

/// Launch `k_eval_sidon_dyadic_budget` on a batch of 64-bit masks on the GPU.
pub fn gpu_eval_dyadic_budget_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
) -> Result<Vec<(u32, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_diff = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_lhs = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu.stream.launch_builder(&gpu.fn_eval_sidon_dyadic_budget);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_diff);
    builder.arg(&mut d_lhs);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_diff: Vec<u64> = gpu.stream.clone_dtoh(&d_diff)?;
    let h_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_lhs)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_diff[i], h_lhs[i]));
    }
    Ok(out)
}

/// Launch `k_profile_sidon_cross_shell` on a batch of 64-bit masks on the GPU
/// (Phase 2: Theorem 2.2 & Theorem 3.2).
pub fn gpu_profile_cross_shell_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
) -> Result<Vec<(u32, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_intra = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_inter = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_bilinear = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu.stream.launch_builder(&gpu.fn_profile_sidon_cross_shell);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_intra);
    builder.arg(&mut d_inter);
    builder.arg(&mut d_bilinear);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_intra: Vec<u64> = gpu.stream.clone_dtoh(&d_intra)?;
    let h_inter: Vec<u64> = gpu.stream.clone_dtoh(&d_inter)?;
    let h_bilinear: Vec<u64> = gpu.stream.clone_dtoh(&d_bilinear)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_intra[i], h_inter[i], h_bilinear[i]));
    }
    Ok(out)
}

/// Launch `k_eval_sidon_subwindow_spectrum` on a batch of 64-bit masks on the GPU
/// (Phase 3: Theorem 2.3 & Theorem 3.3).
pub fn gpu_eval_subwindow_spectrum_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
) -> Result<Vec<(u32, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_pigeon = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_exact = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_short = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_sidon_subwindow_spectrum);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_pigeon);
    builder.arg(&mut d_exact);
    builder.arg(&mut d_short);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_pigeon: Vec<u64> = gpu.stream.clone_dtoh(&d_pigeon)?;
    let h_exact: Vec<u64> = gpu.stream.clone_dtoh(&d_exact)?;
    let h_short: Vec<u64> = gpu.stream.clone_dtoh(&d_short)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_pigeon[i], h_exact[i], h_short[i]));
    }
    Ok(out)
}

/// Launch `k_eval_base_b_joint_fejer` on a batch of 64-bit masks on the GPU
/// (Phase 5: Theorem 2.5 — Direct Prefix-to-Block Unsplit Fejér Coupling &
/// Base-$B = 2^{1/L}$ Window-Block Refinement Hierarchy).
pub fn gpu_eval_base_b_joint_fejer_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> Result<Vec<(u32, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_l1 = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_lb = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_rhs = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu.stream.launch_builder(&gpu.fn_eval_base_b_joint_fejer);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_l1);
    builder.arg(&mut d_lb);
    builder.arg(&mut d_rhs);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_l1: Vec<u64> = gpu.stream.clone_dtoh(&d_l1)?;
    let h_lb: Vec<u64> = gpu.stream.clone_dtoh(&d_lb)?;
    let h_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_rhs)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_l1[i], h_lb[i], h_rhs[i]));
    }
    Ok(out)
}

/// Output tuple `(is_sidon, scaled_continuous_parabola, scaled_integrality_excess, exact_refined_pigeonhole_lb, exact_fejer_rhs)`
/// of `k_eval_coupled_lag_integrality_fejer`.
pub type CoupledLagFejerRow = (u32, u64, u64, u64, u64);

/// Launch `k_eval_coupled_lag_integrality_fejer` on a batch of 64-bit masks on the GPU
/// (Phase 6: Theorem 2.6 — Exact Discrete Legendre–Fenchel Integrality Excess &
/// Coupled Unsplit Fejér Evaluator).
pub fn gpu_eval_coupled_lag_integrality_fejer_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> Result<Vec<CoupledLagFejerRow>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_cont = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_int = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_l = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_rhs = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_coupled_lag_integrality_fejer);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_cont);
    builder.arg(&mut d_int);
    builder.arg(&mut d_l);
    builder.arg(&mut d_rhs);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_cont: Vec<u64> = gpu.stream.clone_dtoh(&d_cont)?;
    let h_int: Vec<u64> = gpu.stream.clone_dtoh(&d_int)?;
    let h_l: Vec<u64> = gpu.stream.clone_dtoh(&d_l)?;
    let h_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_rhs)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_cont[i], h_int[i], h_l[i], h_rhs[i]));
    }
    Ok(out)
}

/// Output tuple `(is_sidon, cf2_full_pigeonhole_lhs, cf2_oct_pigeonhole_lhs, cf2_full_exact_rhs, cf2_oct_exact_rhs)`
/// of `k_eval_cesaro_fejer_order2`.
pub type CesaroFejerOrder2Row = (u32, u64, u64, u64, u64);

/// Launch `k_eval_cesaro_fejer_order2` on a batch of 64-bit masks on the GPU
/// (Phase 7: Theorem 2.7 — Higher-Order Cesàro–Fejér Order-2 & Upper-Octave
/// Iterated Sliding-Window Evaluator).
pub fn gpu_eval_cesaro_fejer_order2_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    refinement_levels: u32,
) -> Result<Vec<CesaroFejerOrder2Row>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_full_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_oct_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_full_rhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_oct_rhs = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu.stream.launch_builder(&gpu.fn_eval_cesaro_fejer_order2);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_full_lhs);
    builder.arg(&mut d_oct_lhs);
    builder.arg(&mut d_full_rhs);
    builder.arg(&mut d_oct_rhs);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_full_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_full_lhs)?;
    let h_oct_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_oct_lhs)?;
    let h_full_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_full_rhs)?;
    let h_oct_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_oct_rhs)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((
            h_ok[i],
            h_full_lhs[i],
            h_oct_lhs[i],
            h_full_rhs[i],
            h_oct_rhs[i],
        ));
    }
    Ok(out)
}

/// Output tuple `(is_sidon, band_pigeonhole_lhs, band_scaled_integrality, band_exact_rhs, band_binom_cap)`
/// of `k_eval_band_restricted_cesaro_fejer`.
pub type BandRestrictedCesaroRow = (u32, u64, u64, u64, u64);

/// Launch `k_eval_band_restricted_cesaro_fejer` on a batch of 64-bit masks on the GPU
/// (Phase 8: Theorem 2.8 — Band-Restricted Upper-Octave Cesàro–Fejér Evaluator).
pub fn gpu_eval_band_restricted_cesaro_fejer_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    w0_cutoff: u32,
    refinement_levels: u32,
) -> Result<Vec<BandRestrictedCesaroRow>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_band_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_band_int = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_band_rhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_band_cap = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_band_restricted_cesaro_fejer);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&w0_cutoff);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_band_lhs);
    builder.arg(&mut d_band_int);
    builder.arg(&mut d_band_rhs);
    builder.arg(&mut d_band_cap);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_band_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_band_lhs)?;
    let h_band_int: Vec<u64> = gpu.stream.clone_dtoh(&d_band_int)?;
    let h_band_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_band_rhs)?;
    let h_band_cap: Vec<u64> = gpu.stream.clone_dtoh(&d_band_cap)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((
            h_ok[i],
            h_band_lhs[i],
            h_band_int[i],
            h_band_rhs[i],
            h_band_cap[i],
        ));
    }
    Ok(out)
}

/// Launch `k_eval_multioctave_cesaro_fejer` on the GPU across a batch of 64-bit masks.
#[allow(clippy::type_complexity)]
pub fn gpu_eval_multioctave_cesaro_fejer_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> Result<Vec<(u32, u64, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_int = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_rhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_cap = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_multioctave_cesaro_fejer);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&k_octaves);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_lhs);
    builder.arg(&mut d_int);
    builder.arg(&mut d_rhs);
    builder.arg(&mut d_cap);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_lhs)?;
    let h_int: Vec<u64> = gpu.stream.clone_dtoh(&d_int)?;
    let h_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_rhs)?;
    let h_cap: Vec<u64> = gpu.stream.clone_dtoh(&d_cap)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_lhs[i], h_int[i], h_rhs[i], h_cap[i]));
    }
    Ok(out)
}

/// Launch `k_eval_haar_multioctave_cesaro_fejer` on the GPU across a batch of 64-bit masks.
#[allow(clippy::type_complexity)]
pub fn gpu_eval_haar_multioctave_cesaro_fejer_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> Result<Vec<(u32, u64, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_int = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_rhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_cap = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_haar_multioctave_cesaro_fejer);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&k_octaves);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_lhs);
    builder.arg(&mut d_int);
    builder.arg(&mut d_rhs);
    builder.arg(&mut d_cap);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_lhs)?;
    let h_int: Vec<u64> = gpu.stream.clone_dtoh(&d_int)?;
    let h_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_rhs)?;
    let h_cap: Vec<u64> = gpu.stream.clone_dtoh(&d_cap)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_lhs[i], h_int[i], h_rhs[i], h_cap[i]));
    }
    Ok(out)
}

/// Launch `k_eval_unconditional_diophantine_haar_closure` on the GPU across a batch of 64-bit masks.
#[allow(clippy::type_complexity)]
pub fn gpu_eval_unconditional_diophantine_haar_closure_batch(
    gpu: &Gpu,
    masks: &[u64],
    m_scales: u32,
    r_target: u32,
    k_octaves: u32,
    refinement_levels: u32,
) -> Result<Vec<(u32, u64, u64, u64, u64)>> {
    if masks.is_empty() {
        return Ok(Vec::new());
    }
    let num = masks.len();
    let d_masks = gpu.stream.clone_htod(masks)?;
    let mut d_ok = gpu.stream.alloc_zeros::<u32>(num)?;
    let mut d_lhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_int = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_rhs = gpu.stream.alloc_zeros::<u64>(num)?;
    let mut d_cap = gpu.stream.alloc_zeros::<u64>(num)?;

    let num_u32 = num as u32;
    let cfg = gpu.launch_cfg(num_u32);
    let mut builder = gpu
        .stream
        .launch_builder(&gpu.fn_eval_unconditional_diophantine_haar_closure);
    builder.arg(&d_masks);
    builder.arg(&num_u32);
    builder.arg(&m_scales);
    builder.arg(&r_target);
    builder.arg(&k_octaves);
    builder.arg(&refinement_levels);
    builder.arg(&mut d_ok);
    builder.arg(&mut d_lhs);
    builder.arg(&mut d_int);
    builder.arg(&mut d_rhs);
    builder.arg(&mut d_cap);
    unsafe { builder.launch(cfg) }?;

    let h_ok: Vec<u32> = gpu.stream.clone_dtoh(&d_ok)?;
    let h_lhs: Vec<u64> = gpu.stream.clone_dtoh(&d_lhs)?;
    let h_int: Vec<u64> = gpu.stream.clone_dtoh(&d_int)?;
    let h_rhs: Vec<u64> = gpu.stream.clone_dtoh(&d_rhs)?;
    let h_cap: Vec<u64> = gpu.stream.clone_dtoh(&d_cap)?;

    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push((h_ok[i], h_lhs[i], h_int[i], h_rhs[i], h_cap[i]));
    }
    Ok(out)
}

/// Compute the minimum Erdős #1191 ratio $\min_{x_0 \le x \le n} A(x)\sqrt{(\ln x)/x}$
/// and minimizing coordinate $x^*$ for a sorted Sidon set `set` on $[1, n]$.
pub fn min_phi_on_window(set: &[u64], x0: usize, n: usize) -> (f64, usize) {
    let x_start = x0.clamp(2, n);
    let mut min_phi = f64::INFINITY;
    let mut min_x = x_start;
    let mut ptr = 0usize;
    while ptr < set.len() && (set[ptr] as usize) < x_start {
        ptr += 1;
    }
    for x in x_start..=n {
        while ptr < set.len() && (set[ptr] as usize) <= x {
            ptr += 1;
        }
        let phi = erdos1191_phi(ptr, x);
        if phi < min_phi {
            min_phi = phi;
            min_x = x;
        }
    }
    (min_phi, min_x)
}

/// Exhaustively determine the finite Sidon extremal value $F_2(n)$ (maximum size
/// of a Sidon subset of $\{1, \dots, n\}$, $1 \le n \le 63$) using $O(1)$ bit-shift
/// branch-and-bound.
pub fn solve_f2_sidon_exact(
    n: usize,
    split_x: u32,
    unconditional: bool,
    gpu: Option<&Gpu>,
) -> Result<ExactSidonResult> {
    if n == 0 || n > 63 {
        bail!("solve_f2_sidon_exact requires 1 <= n <= 63 (got {n})");
    }
    let sidon_cap = build_sidon_cap(n, unconditional);
    let req_prefix = vec![0u32; n + 1];
    let n_u32 = n as u32;

    // Seed with the Mian–Chowla greedy prefix inside {1..=n}
    let mc: Vec<u64> = crate::records::A005282_MIAN_CHOWLA_64
        .iter()
        .copied()
        .filter(|&x| (x as usize) <= n)
        .collect();
    assert!(is_sidon_cpu(&mc));
    let mut best_card = mc.len() as u32;
    let mut best_witness = mc;
    let mut grand_nodes = 0u64;

    let mut target = best_card + 1;
    while target <= sidon_cap[n] {
        let prefixes =
            expand_sidon_prefixes(n_u32, split_x.min(n_u32), target, &sidon_cap, &req_prefix);
        let (witness_opt, nodes) = if let Some(g) = gpu {
            gpu_extend_sidon_batch(g, &prefixes, n_u32, target, &sidon_cap, &req_prefix)?
        } else {
            let results: Vec<(Option<u64>, u64)> = prefixes
                .par_iter()
                .map(|&p| extend_sidon_dfs_cpu(p, n_u32, target, &sidon_cap, &req_prefix))
                .collect();
            let tot: u64 = results.iter().map(|r| r.1).sum();
            let wit = results.into_iter().find_map(|r| r.0);
            (wit, tot)
        };
        grand_nodes += nodes;

        if let Some(mask) = witness_opt {
            assert!(
                is_sidon_mask64(mask),
                "witness mask {mask:#x} failed 64-bit Sidon check!"
            );
            let vec_w = mask_to_vec(mask);
            assert!(
                vec_w.len() as u32 >= target && is_sidon_cpu(&vec_w),
                "witness {vec_w:?} failed independent CPU Sidon check!"
            );
            best_card = vec_w.len() as u32;
            best_witness = vec_w;
            target = best_card + 1;
        } else {
            break;
        }
    }

    let (min_phi, min_phi_x) = min_phi_on_window(&best_witness, 2, n);
    Ok(ExactSidonResult {
        n,
        max_card: best_card,
        min_phi,
        min_phi_x,
        witness: best_witness,
        total_nodes: grand_nodes,
    })
}

/// Exhaustively solve the **Maximin Lower-Envelope Sidon Problem** on $\{1, \dots, n\}$:
///
///     \mu^*(x_0, n) := \max_{A \subseteq \{1, \dots, n\} \text{ Sidon}} \min_{x \in \{x_0, \dots, n\}} A(x)\sqrt{(\ln x)/x}.
///
/// Because $A(x)$ is integer-valued, $\mu^*(x_0, n)$ belongs to the finite set of
/// critical thresholds $\tau = k \sqrt{(\ln x)/x}$ ($x \in [x_0, n]$, $1 \le k \le F_2(x)$),
/// and testing $\min_{x \in [x_0, n]} \Phi_A(x) \ge \tau$ reduces to the integer-exact
/// prefix-floor constraint $A(x) \ge \lceil \tau \sqrt{x / \ln x} \rceil$.
pub fn solve_maximin_lower_envelope_exact(
    n: usize,
    x0: usize,
    split_x: u32,
    gpu: Option<&Gpu>,
) -> Result<ExactSidonResult> {
    if !(2..=63).contains(&n) {
        bail!("solve_maximin_lower_envelope_exact requires 2 <= n <= 63 (got {n})");
    }
    let x_start = x0.clamp(2, n);
    let sidon_cap = build_sidon_cap(n, false);
    let n_u32 = n as u32;

    // Seed with Mian–Chowla prefix inside {1..=n}
    let mc: Vec<u64> = crate::records::A005282_MIAN_CHOWLA_64
        .iter()
        .copied()
        .filter(|&x| (x as usize) <= n)
        .collect();
    let (mut best_phi, mut best_phi_x) = min_phi_on_window(&mc, x_start, n);
    let mut best_witness = mc;
    let mut grand_nodes = 0u64;

    // Collect all candidate critical thresholds \tau > best_phi
    let mut candidates = Vec::new();
    for (x, &cap_x) in sidon_cap.iter().enumerate().take(n + 1).skip(x_start) {
        for k in 1..=cap_x {
            let tau = erdos1191_phi(k as usize, x);
            if tau > best_phi + 1e-12 {
                candidates.push(tau);
            }
        }
    }
    candidates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    candidates.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

    // Binary search over discrete critical thresholds
    let mut lo = 0usize;
    let mut hi = candidates.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        let tau = candidates[mid];

        // Build non-decreasing integer prefix-floor array req_prefix[0..=n]
        let mut req_prefix = vec![0u32; n + 1];
        let mut possible = true;
        for x in x_start..=n {
            let raw = (tau * ((x as f64) / (x as f64).ln()).sqrt() - 1e-12)
                .ceil()
                .max(0.0) as u32;
            let req = raw.max(req_prefix[x - 1]);
            if req > sidon_cap[x] {
                possible = false;
                break;
            }
            req_prefix[x] = req;
        }

        if !possible {
            hi = mid;
            continue;
        }

        let target = req_prefix[n];
        let prefixes =
            expand_sidon_prefixes(n_u32, split_x.min(n_u32), target, &sidon_cap, &req_prefix);
        let (witness_opt, nodes) = if let Some(g) = gpu {
            gpu_extend_sidon_batch(g, &prefixes, n_u32, target, &sidon_cap, &req_prefix)?
        } else {
            let results: Vec<(Option<u64>, u64)> = prefixes
                .par_iter()
                .map(|&p| extend_sidon_dfs_cpu(p, n_u32, target, &sidon_cap, &req_prefix))
                .collect();
            let tot: u64 = results.iter().map(|r| r.1).sum();
            let wit = results.into_iter().find_map(|r| r.0);
            (wit, tot)
        };
        grand_nodes += nodes;

        if let Some(mask) = witness_opt {
            assert!(is_sidon_mask64(mask));
            let vec_w = mask_to_vec(mask);
            assert!(is_sidon_cpu(&vec_w));
            let (phi_w, phi_x) = min_phi_on_window(&vec_w, x_start, n);
            assert!(phi_w + 1e-10 >= tau);
            best_phi = phi_w;
            best_phi_x = phi_x;
            best_witness = vec_w;
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    // Maximize total cardinality among Sidon sets achieving the optimal lower envelope best_phi
    let mut final_req = vec![0u32; n + 1];
    for x in x_start..=n {
        let raw = (best_phi * ((x as f64) / (x as f64).ln()).sqrt() - 1e-12)
            .ceil()
            .max(0.0) as u32;
        final_req[x] = raw.max(final_req[x - 1]);
    }
    let mut next_card = best_witness.len() as u32 + 1;
    while next_card <= sidon_cap[n] {
        let prefixes =
            expand_sidon_prefixes(n_u32, split_x.min(n_u32), next_card, &sidon_cap, &final_req);
        let (witness_opt, nodes) = if let Some(g) = gpu {
            gpu_extend_sidon_batch(g, &prefixes, n_u32, next_card, &sidon_cap, &final_req)?
        } else {
            let results: Vec<(Option<u64>, u64)> = prefixes
                .par_iter()
                .map(|&p| extend_sidon_dfs_cpu(p, n_u32, next_card, &sidon_cap, &final_req))
                .collect();
            let tot: u64 = results.iter().map(|r| r.1).sum();
            let wit = results.into_iter().find_map(|r| r.0);
            (wit, tot)
        };
        grand_nodes += nodes;
        if let Some(mask) = witness_opt {
            let vec_w = mask_to_vec(mask);
            assert!(is_sidon_cpu(&vec_w));
            let (phi_w, phi_x) = min_phi_on_window(&vec_w, x_start, n);
            best_phi = phi_w;
            best_phi_x = phi_x;
            best_witness = vec_w;
            next_card += 1;
        } else {
            break;
        }
    }

    Ok(ExactSidonResult {
        n,
        max_card: best_witness.len() as u32,
        min_phi: best_phi,
        min_phi_x: best_phi_x,
        witness: best_witness,
        total_nodes: grand_nodes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify_cpu::{
        eval_band_restricted_cesaro_fejer_cpu, eval_base_b_joint_fejer_cpu,
        eval_cesaro_fejer_order2_cpu, eval_coupled_lag_integrality_fejer_cpu,
        eval_haar_multioctave_cesaro_fejer_cpu, eval_multioctave_cesaro_fejer_cpu,
        eval_sidon_dyadic_budget_cpu, eval_sidon_subwindow_spectrum_cpu,
        eval_unconditional_diophantine_haar_closure_cpu, profile_sidon_cross_shell_cpu,
    };

    #[test]
    fn exact_sidon_search_reproduces_oeis_a003022_through_26() {
        for &(n, expected_f2) in &[
            (4usize, 3u32),
            (6, 3),
            (7, 4),
            (11, 4),
            (12, 5),
            (18, 6),
            (26, 7),
        ] {
            let res = solve_f2_sidon_exact(n, 4, false, None).unwrap();
            assert_eq!(res.max_card, expected_f2, "mismatch at n = {n}");
            assert!(is_sidon_cpu(&res.witness));
        }
    }

    #[test]
    fn unconditional_sidon_search_matches_conditional_n18() {
        let r_cond = solve_f2_sidon_exact(18, 4, false, None).unwrap();
        let r_uncond = solve_f2_sidon_exact(18, 4, true, None).unwrap();
        assert_eq!(r_cond.max_card, r_uncond.max_card);
        assert_eq!(r_cond.max_card, 6);
    }

    #[test]
    fn maximin_lower_envelope_strictly_improves_over_mian_chowla_n25() {
        let n = 25usize;
        let x0 = 7usize;
        let mc_n: Vec<u64> = crate::records::A005282_MIAN_CHOWLA_64
            .iter()
            .copied()
            .filter(|&x| (x as usize) <= n)
            .collect();
        let (mc_phi, _) = min_phi_on_window(&mc_n, x0, n);
        let res = solve_maximin_lower_envelope_exact(n, x0, 4, None).unwrap();
        assert!(is_sidon_cpu(&res.witness));
        assert!(
            res.min_phi >= mc_phi - 1e-12,
            "maximin {} < Mian-Chowla {}",
            res.min_phi,
            mc_phi
        );
        assert!(res.min_phi > mc_phi + 0.05);
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_sidon_dfs_matches_cpu_replica_n26() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let n = 26usize;
        let sidon_cap = build_sidon_cap(n, false);
        let req_prefix = vec![0u32; n + 1];
        for target in 5..=8 {
            let prefixes = expand_sidon_prefixes(n as u32, 5, target, &sidon_cap, &req_prefix);
            let (_, gpu_nodes) =
                gpu_extend_sidon_batch(&gpu, &prefixes, n as u32, target, &sidon_cap, &req_prefix)
                    .unwrap();
            let cpu_nodes: u64 = prefixes
                .iter()
                .map(|&p| extend_sidon_dfs_cpu(p, n as u32, target, &sidon_cap, &req_prefix).1)
                .sum();
            assert_eq!(
                gpu_nodes, cpu_nodes,
                "CPU/GPU node count divergence at n={n}, target={target}"
            );
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_dyadic_budget_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        // Also add non-Sidon masks to test both ok=1 and ok=0 paths
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        let gpu_out = gpu_eval_dyadic_budget_batch(&gpu, &masks, 5).unwrap();
        for (i, &m) in masks.iter().enumerate() {
            let cpu_tuple = eval_sidon_dyadic_budget_cpu(m, 5);
            assert_eq!(gpu_out[i], cpu_tuple, "mismatch at index {i}, mask {m:#x}");
        }
    }

    #[test]
    fn cpu_cross_shell_profiler_matches_bilinear_conservation() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 8, 4, &sidon_cap, &req);
        for p in prefixes {
            let (ok, intra, inter, scaled_bil) = profile_sidon_cross_shell_cpu(p.chosen_mask, 5);
            assert_eq!(ok, 1);
            let k_tot = p.chosen_mask.count_ones() as u64;
            assert_eq!(intra + inter, k_tot * k_tot.saturating_sub(1) / 2);
            let (_, _, scaled_intra) = eval_sidon_dyadic_budget_cpu(p.chosen_mask, 5);
            assert!(scaled_bil >= scaled_intra);
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_cross_shell_profiler_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        let gpu_out = gpu_profile_cross_shell_batch(&gpu, &masks, 5).unwrap();
        for (i, &m) in masks.iter().enumerate() {
            let cpu_tuple = profile_sidon_cross_shell_cpu(m, 5);
            assert_eq!(
                gpu_out[i], cpu_tuple,
                "cross-shell mismatch at index {i}, mask {m:#x}"
            );
        }
    }

    #[test]
    fn cpu_subwindow_spectrum_satisfies_theorem_2_3_chain() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 1..=5 {
                let (ok, pigeon_lb, exact_sub, short_diffs) =
                    eval_sidon_subwindow_spectrum_cpu(p.chosen_mask, 5, r);
                assert_eq!(ok, 1);
                assert!(
                    pigeon_lb <= exact_sub,
                    "pigeon_lb {pigeon_lb} > exact_sub {exact_sub} at r={r}"
                );
                assert!(
                    exact_sub <= short_diffs,
                    "exact_sub {exact_sub} > short_diffs {short_diffs} at r={r}"
                );
                assert!(
                    short_diffs < (1u64 << r),
                    "short_diffs {short_diffs} > 2^r - 1 at r={r}"
                );
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_subwindow_spectrum_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 1..=5 {
            let gpu_out = gpu_eval_subwindow_spectrum_batch(&gpu, &masks, 5, r).unwrap();
            for (i, &m) in masks.iter().enumerate() {
                let cpu_tuple = eval_sidon_subwindow_spectrum_cpu(m, 5, r);
                assert_eq!(
                    gpu_out[i], cpu_tuple,
                    "subwindow spectrum mismatch at r={r}, index {i}, mask {m:#x}"
                );
            }
        }
    }

    #[test]
    fn cpu_base_b_joint_fejer_satisfies_theorem_2_5_refinement_chain() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 1..=5 {
                let w = 1u64 << r;
                let cap_binom = w * (w - 1) / 2;
                let (ok1, l1, lb1, rhs1) = eval_base_b_joint_fejer_cpu(p.chosen_mask, 5, r, 1);
                let (ok2, _, lb2, _) = eval_base_b_joint_fejer_cpu(p.chosen_mask, 5, r, 2);
                let (ok3, _, lb3, _) = eval_base_b_joint_fejer_cpu(p.chosen_mask, 5, r, 3);
                assert_eq!(ok1, 1);
                assert_eq!(ok2, 1);
                assert_eq!(ok3, 1);
                assert!(l1 <= lb1 && lb1 <= lb2 && lb2 <= lb3 && lb3 <= rhs1 && rhs1 <= cap_binom);
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_base_b_joint_fejer_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 1..=5 {
            for lv in [1u32, 2, 3] {
                let gpu_out = gpu_eval_base_b_joint_fejer_batch(&gpu, &masks, 5, r, lv).unwrap();
                for (i, &m) in masks.iter().enumerate() {
                    let cpu_tuple = eval_base_b_joint_fejer_cpu(m, 5, r, lv);
                    assert_eq!(
                        gpu_out[i], cpu_tuple,
                        "base-B joint Fejér mismatch at r={r}, lv={lv}, index {i}, mask {m:#x}"
                    );
                }
            }
        }
    }

    #[test]
    fn cpu_coupled_lag_integrality_fejer_satisfies_theorem_2_6b() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 1..=5 {
                let w = 1u64 << r;
                let cap_binom = w * (w - 1) / 2;
                for lv in [0u32, 1, 2] {
                    let (ok, cont, int_ex, l_sum, rhs) =
                        eval_coupled_lag_integrality_fejer_cpu(p.chosen_mask, 5, r, lv);
                    assert_eq!(ok, 1);
                    assert!(l_sum <= rhs && rhs <= cap_binom);
                    // Both scaled continuous parabola and scaled integrality excess are even integers
                    // whose sum equals 2 * sum_b B_b * L(S_b, B_b) >= 2 * l_sum
                    assert_eq!((cont + int_ex) % 2, 0);
                    assert!(cont + int_ex >= 2 * l_sum);
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_coupled_lag_integrality_fejer_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 1..=5 {
            for lv in [0u32, 1, 2] {
                let gpu_out =
                    gpu_eval_coupled_lag_integrality_fejer_batch(&gpu, &masks, 5, r, lv).unwrap();
                for (i, &m) in masks.iter().enumerate() {
                    let cpu_tuple = eval_coupled_lag_integrality_fejer_cpu(m, 5, r, lv);
                    assert_eq!(
                        gpu_out[i], cpu_tuple,
                        "coupled lag integrality mismatch at r={r}, lv={lv}, index {i}, mask {m:#x}"
                    );
                }
            }
        }
    }

    #[test]
    fn cpu_cesaro_fejer_order2_satisfies_theorem_2_7_chains() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 1..=5 {
                let w = 1u64 << r;
                let w_half = w >> 1;
                let full_binom_cap = (w - 1) * w * (w + 1) / 6;
                let half_binom_cap = (w_half.saturating_sub(1)) * w_half * (w_half + 1) / 6;
                let oct_binom_cap = full_binom_cap - half_binom_cap;
                for lv in [0u32, 1, 2] {
                    let (ok, full_lhs, oct_lhs, full_rhs, oct_rhs) =
                        eval_cesaro_fejer_order2_cpu(p.chosen_mask, 5, r, lv);
                    assert_eq!(ok, 1);
                    assert!(full_lhs <= full_rhs && full_rhs <= full_binom_cap);
                    assert!(oct_lhs <= oct_rhs && oct_rhs <= oct_binom_cap);
                    assert!(oct_lhs <= full_lhs);
                    assert!(oct_rhs <= full_rhs);
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_cesaro_fejer_order2_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 1..=5 {
            for lv in [0u32, 1, 2] {
                let gpu_out = gpu_eval_cesaro_fejer_order2_batch(&gpu, &masks, 5, r, lv).unwrap();
                for (i, &m) in masks.iter().enumerate() {
                    let cpu_tuple = eval_cesaro_fejer_order2_cpu(m, 5, r, lv);
                    assert_eq!(
                        gpu_out[i], cpu_tuple,
                        "Cesàro–Fejér order-2 mismatch at r={r}, lv={lv}, index {i}, mask {m:#x}"
                    );
                }
            }
        }
    }

    #[test]
    fn cpu_band_restricted_cesaro_satisfies_theorem_2_8_chains() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 2..=5 {
                let w = 1u32 << r;
                let w_half = w >> 1;
                let w_tq = (3 * w) >> 2;
                for lv in [0u32, 1, 2] {
                    let (ok_h, lhs_h, int_h, rhs_h, cap_h) =
                        eval_band_restricted_cesaro_fejer_cpu(p.chosen_mask, 5, r, w_half, lv);
                    let (ok_q, lhs_q, int_q, rhs_q, cap_q) =
                        eval_band_restricted_cesaro_fejer_cpu(p.chosen_mask, 5, r, w_tq, lv);
                    assert_eq!(ok_h, 1);
                    assert_eq!(ok_q, 1);
                    assert!(lhs_h <= rhs_h && rhs_h <= cap_h);
                    assert!(lhs_q <= rhs_q && rhs_q <= cap_q);
                    assert!(lhs_q <= lhs_h);
                    assert!(int_q <= int_h);
                    assert!(rhs_q <= rhs_h);
                    assert!(cap_q < cap_h);
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_band_restricted_cesaro_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 2..=5 {
            let w = 1u32 << r;
            for w0 in [w >> 1, (3 * w) >> 2] {
                for lv in [0u32, 1, 2] {
                    let gpu_out =
                        gpu_eval_band_restricted_cesaro_fejer_batch(&gpu, &masks, 5, r, w0, lv)
                            .unwrap();
                    for (i, &m) in masks.iter().enumerate() {
                        let cpu_tuple = eval_band_restricted_cesaro_fejer_cpu(m, 5, r, w0, lv);
                        assert_eq!(
                            gpu_out[i], cpu_tuple,
                            "band-restricted Cesàro–Fejér mismatch at r={r}, w0={w0}, lv={lv}, index {i}, mask {m:#x}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cpu_multioctave_cesaro_satisfies_theorem_2_9_cascade() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 3..=5 {
                let w = 1u32 << r;
                for lv in [0u32, 1, 2] {
                    let (ok1, lhs1, int1, rhs1, cap1) =
                        eval_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 1, lv);
                    let (ok2, lhs2, int2, rhs2, cap2) =
                        eval_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 2, lv);
                    let (ok3, lhs3, int3, rhs3, cap3) =
                        eval_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 3, lv);
                    let band_k1 =
                        eval_band_restricted_cesaro_fejer_cpu(p.chosen_mask, 5, r, w >> 1, lv);
                    assert_eq!((ok1, lhs1, int1, rhs1, cap1), band_k1);
                    assert_eq!(ok1, 1);
                    assert_eq!(ok2, 1);
                    assert_eq!(ok3, 1);
                    assert!(lhs1 <= rhs1 && rhs1 <= cap1);
                    assert!(lhs2 <= rhs2 && rhs2 <= cap2);
                    assert!(lhs3 <= rhs3 && rhs3 <= cap3);
                    assert!(lhs1 <= lhs2 && lhs2 <= lhs3);
                    assert!(int1 <= int2 && int2 <= int3);
                    assert!(rhs1 <= rhs2 && rhs2 <= rhs3);
                    assert!(cap1 < cap2 && cap2 < cap3);
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_multioctave_cesaro_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 3..=5 {
            for k_oct in [1u32, 2, 3] {
                for lv in [0u32, 1, 2] {
                    let gpu_out =
                        gpu_eval_multioctave_cesaro_fejer_batch(&gpu, &masks, 5, r, k_oct, lv)
                            .unwrap();
                    for (i, &m) in masks.iter().enumerate() {
                        let cpu_tuple = eval_multioctave_cesaro_fejer_cpu(m, 5, r, k_oct, lv);
                        assert_eq!(
                            gpu_out[i], cpu_tuple,
                            "multi-octave Cesàro–Fejér mismatch at r={r}, K={k_oct}, lv={lv}, index {i}, mask {m:#x}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cpu_haar_multioctave_cesaro_satisfies_theorem_2_10_cascade() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 3..=5 {
                for lv in [0u32, 1, 2] {
                    let (ok1, lhs1, int1, rhs1, cap1) =
                        eval_haar_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 1, lv);
                    let (ok2, lhs2, int2, rhs2, cap2) =
                        eval_haar_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 2, lv);
                    let (ok3, lhs3, int3, rhs3, cap3) =
                        eval_haar_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 3, lv);
                    let unwt_k1 = eval_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 1, lv);
                    let (_, unwt_lhs3, unwt_int3, unwt_rhs3, unwt_cap3) =
                        eval_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 3, lv);
                    assert_eq!((ok1, lhs1, int1, rhs1, cap1), unwt_k1);
                    assert_eq!(ok1, 1);
                    assert_eq!(ok2, 1);
                    assert_eq!(ok3, 1);
                    assert!(lhs1 <= rhs1 && rhs1 <= cap1);
                    assert!(lhs2 <= rhs2 && rhs2 <= cap2);
                    assert!(lhs3 <= rhs3 && rhs3 <= cap3);
                    assert!(lhs1 <= lhs2 && lhs2 <= lhs3);
                    assert!(int1 <= int2 && int2 <= int3);
                    assert!(rhs1 <= rhs2 && rhs2 <= rhs3);
                    assert!(cap1 < cap2 && cap2 < cap3);
                    assert!(
                        lhs3 >= unwt_lhs3
                            && int3 >= unwt_int3
                            && rhs3 >= unwt_rhs3
                            && cap3 > unwt_cap3
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_haar_multioctave_cesaro_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 3..=5 {
            for k_oct in [1u32, 2, 3] {
                for lv in [0u32, 1, 2] {
                    let gpu_out =
                        gpu_eval_haar_multioctave_cesaro_fejer_batch(&gpu, &masks, 5, r, k_oct, lv)
                            .unwrap();
                    for (i, &m) in masks.iter().enumerate() {
                        let cpu_tuple = eval_haar_multioctave_cesaro_fejer_cpu(m, 5, r, k_oct, lv);
                        assert_eq!(
                            gpu_out[i], cpu_tuple,
                            "Haar multi-octave Cesàro–Fejér mismatch at r={r}, K={k_oct}, lv={lv}, index {i}, mask {m:#x}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cpu_unconditional_diophantine_haar_satisfies_theorem_2_11_closure() {
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 9, 4, &sidon_cap, &req);
        for p in prefixes {
            for r in 3..=5 {
                let w_cap = (1u32 << r) as u64;
                for lv in [0u32, 1, 2] {
                    let (ok1, lhs1, int1, rhs1, cap1) =
                        eval_unconditional_diophantine_haar_closure_cpu(p.chosen_mask, 5, r, 1, lv);
                    let (ok2, lhs2, int2, rhs2, cap2) =
                        eval_unconditional_diophantine_haar_closure_cpu(p.chosen_mask, 5, r, 2, lv);
                    let (ok3, lhs3, int3, rhs3, cap3) =
                        eval_unconditional_diophantine_haar_closure_cpu(p.chosen_mask, 5, r, 3, lv);
                    let (_, h_lhs3, h_int3, h_rhs3, h_cap3) =
                        eval_haar_multioctave_cesaro_fejer_cpu(p.chosen_mask, 5, r, 3, lv);
                    assert_eq!(ok1, 1);
                    assert_eq!(ok2, 1);
                    assert_eq!(ok3, 1);
                    assert!(lhs1 <= rhs1 && rhs1 <= cap1);
                    assert!(lhs2 <= rhs2 && rhs2 <= cap2);
                    assert!(lhs3 <= rhs3 && rhs3 <= cap3);
                    assert!(lhs1 <= lhs2 && lhs2 <= lhs3);
                    assert!(int1 <= int2 && int2 <= int3);
                    assert!(rhs1 <= rhs2 && rhs2 <= rhs3);
                    assert!(cap1 < cap2 && cap2 < cap3);
                    // Because floor(W^2 / w) >= W * floor(W / w), the scaled-reciprocal
                    // Diophantine Haar operator strictly dominates W times Phase 10's floor(W / w):
                    assert!(lhs3 >= w_cap * h_lhs3);
                    assert!(int3 >= w_cap * h_int3);
                    assert!(rhs3 >= w_cap * h_rhs3);
                    assert!(cap3 > w_cap * h_cap3);
                }
            }
        }
    }

    #[test]
    #[ignore = "requires CUDA GPU (compute_86)"]
    fn gpu_unconditional_diophantine_haar_matches_cpu_replica() {
        let gpu = Gpu::init(0, "compute_86").expect("GPU init failed");
        let sidon_cap = build_sidon_cap(32, false);
        let req = vec![0u32; 33];
        let prefixes = expand_sidon_prefixes(32, 10, 4, &sidon_cap, &req);
        let mut masks: Vec<u64> = prefixes.iter().map(|p| p.chosen_mask).collect();
        masks.push((1u64 << 1) | (1u64 << 2) | (1u64 << 3) | (1u64 << 4));
        masks.push(
            (1u64 << 1)
                | (1u64 << 2)
                | (1u64 << 5)
                | (1u64 << 11)
                | (1u64 << 19)
                | (1u64 << 24)
                | (1u64 << 26),
        );

        for r in 3..=5 {
            for k_oct in [1u32, 2, 3] {
                for lv in [0u32, 1, 2] {
                    let gpu_out = gpu_eval_unconditional_diophantine_haar_closure_batch(
                        &gpu, &masks, 5, r, k_oct, lv,
                    )
                    .unwrap();
                    for (i, &m) in masks.iter().enumerate() {
                        let cpu_tuple =
                            eval_unconditional_diophantine_haar_closure_cpu(m, 5, r, k_oct, lv);
                        assert_eq!(
                            gpu_out[i], cpu_tuple,
                            "Diophantine Haar multi-octave closure mismatch at r={r}, K={k_oct}, lv={lv}, index {i}, mask {m:#x}"
                        );
                    }
                }
            }
        }
    }
}
