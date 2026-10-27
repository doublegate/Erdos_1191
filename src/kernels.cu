// ============================================================================
// `kernels.cu` — CUDA C kernels for Erdős Problem #1191 ($1,000 open prize:
// lower-envelope fluctuations liminf A(x) sqrt(log x / x) of infinite Sidon sets):
//   1. `k_dfs_maximin_sidon`: Exhaustive O(1) bit-shift branch-and-bound solver
//      for both the finite Sidon extremal function F_2(N) (OEIS A143824 / A003022)
//      and the Maximin Lower-Envelope Sidon Problem
//          \mu^*(x_0, N) = \max_{A \subseteq [N] \text{ Sidon}} \min_{x_0 \le x \le N} A(x)\sqrt{(\ln x)/x}
//      via integer prefix-floor constraints `req_prefix[1..=N]` and non-circular
//      subadditive tail capacities `sidon_cap[0..=N]`.
//      Bit-exact twin of `extend_sidon_dfs_cpu` (`src/exact.rs`).
//   2. `k_eval_sidon_dyadic_budget`: Parallel batch evaluation of Sidon ($B_2$)
//      validity, exact 64-bit positive difference bitmasks, and Theorem 2.1(b)
//      scaled Dyadic Difference-Budget sums \sum_{k=1}^M 2^{M-k} m_k(m_k - 1).
//      Bit-exact twin of `eval_sidon_dyadic_budget_cpu` (`src/verify_cpu.rs`).
// ============================================================================
//
// ENCODING CONTRACTS (shared with `src/verify_cpu.rs` and `src/exact.rs` —
// DO NOT CHANGE UNILATERALLY):
//   * Set bitmask layout (`u64`): element `x ∈ {1, ..., 63}` is represented by
//     bit `(1ULL << x)`. Bit 0 is always 0.
//   * Reversed bitmask (`rev_mask`): bit `(63 - a)` is 1 for every `a ∈ A`.
//   * O(1) bit-shift Sidon test: for candidate `x > max(A)`, the mask of all
//     new positive differences `{x - a : a ∈ A}` is `new_diffs = rev_mask >> (63 - x)`.
//     `A ∪ {x}` is Sidon iff `(new_diffs & diff_mask) == 0`.
//   * Node-count definition: a DFS node is counted whenever an element `x` is
//     legally added (passing the prefix-floor check, tail-capacity check, and
//     `(new_diffs & diff_mask) == 0` Sidon check).

typedef unsigned int       u32;
typedef unsigned long long u64;

extern "C" __global__ void k_dfs_maximin_sidon(
    const u64* __restrict__ prefix_chosen,
    const u64* __restrict__ prefix_rev,
    const u64* __restrict__ prefix_diff,
    const u32* __restrict__ prefix_last,
    const u32* __restrict__ prefix_count,
    u32 num_prefixes,
    u32 n,
    u32 target_card,
    const u32* __restrict__ sidon_cap,
    const u32* __restrict__ req_prefix,
    u32* __restrict__ found_flag,
    u64* __restrict__ witness_mask_out,
    u64* __restrict__ node_counts_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_prefixes) {
        return;
    }

    u64 init_chosen = prefix_chosen[tid];
    u64 init_rev    = prefix_rev[tid];
    u64 init_diff   = prefix_diff[tid];
    u32 init_last   = prefix_last[tid];
    u32 init_count  = prefix_count[tid];

    u32 eff_target = target_card >= req_prefix[n] ? target_card : req_prefix[n];

    if (init_count >= eff_target) {
        if (atomicCAS(found_flag, 0U, 1U) == 0U) {
            *witness_mask_out = init_chosen;
        }
        node_counts_out[tid] = 0ULL;
        return;
    }

    // Explicit DFS stack up to depth 64
    u64 stk_chosen[64];
    u64 stk_rev[64];
    u64 stk_diff[64];
    u32 stk_next_x[64];
    u32 stk_count[64];

    stk_chosen[0] = init_chosen;
    stk_rev[0]    = init_rev;
    stk_diff[0]   = init_diff;
    stk_next_x[0] = init_last + 1U;
    stk_count[0]  = init_count;

    int sp = 0;
    u64 nodes = 0ULL;

    while (sp >= 0) {
        u32 x = stk_next_x[sp];
        u32 cnt = stk_count[sp];

        // Prefix floor check at (x - 1): since all elements <= x - 1 have been
        // decided, A(x - 1) == cnt. If cnt < req_prefix[x - 1], this branch
        // violates the lower-envelope constraint and must be pruned immediately.
        if (x > 1U && cnt < req_prefix[x - 1U]) {
            sp--;
            continue;
        }

        if (x > n) {
            sp--;
            continue;
        }

        // Advance cursor at current stack frame
        stk_next_x[sp] = x + 1U;

        // Subadditive tail-capacity prune: even [x, n] of width (n - x + 1)
        // can contribute at most sidon_cap[n - x + 1] additional elements.
        if (cnt + sidon_cap[n - x + 1U] < eff_target) {
            sp--;
            continue;
        }
        if (cnt + 1U + sidon_cap[n - x] < eff_target) {
            continue;
        }

        // O(1) bit-shift Sidon test
        u64 rev = stk_rev[sp];
        u64 diff = stk_diff[sp];
        u64 new_diffs = rev >> (63U - x);
        if ((new_diffs & diff) != 0ULL) {
            continue;
        }

        // Legal Sidon addition of x
        nodes += 1ULL;
        u64 next_chosen = stk_chosen[sp] | (1ULL << x);
        u32 next_cnt = cnt + 1U;

        if (next_cnt >= eff_target) {
            if (atomicCAS(found_flag, 0U, 1U) == 0U) {
                *witness_mask_out = next_chosen;
            }
            continue;
        }

        sp++;
        stk_chosen[sp] = next_chosen;
        stk_rev[sp]    = rev | (1ULL << (63U - x));
        stk_diff[sp]   = diff | new_diffs;
        stk_next_x[sp] = x + 1U;
        stk_count[sp]  = next_cnt;
    }

    node_counts_out[tid] = nodes;
}

extern "C" __global__ void k_eval_sidon_dyadic_budget(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ diff_mask_out,
    u64* __restrict__ scaled_lhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u64 scaled_lhs = 0ULL;
    for (u32 k = 1U; k <= m_scales; k++) {
        u32 lo = (1U << (k - 1U)) + 1U;
        u32 hi = (1U << k);
        if (hi > 63U) {
            hi = 63U;
        }
        u64 shell_mask = 0ULL;
        if (lo <= hi) {
            u64 hi_mask = (hi == 63U) ? (~0ULL) : ((1ULL << (hi + 1U)) - 1ULL);
            u64 lo_mask = (1ULL << lo) - 1ULL;
            shell_mask = hi_mask & ~lo_mask;
        }
        u64 mk = (u64)__popcll(mask & shell_mask);
        if (mk >= 2ULL) {
            scaled_lhs += (mk * (mk - 1ULL)) << (m_scales - k);
        }
    }

    is_sidon_out[tid]   = ok;
    diff_mask_out[tid]  = diff_mask;
    scaled_lhs_out[tid] = scaled_lhs;
}

// ============================================================================
// 3. `k_profile_sidon_cross_shell`: Phase 2 (Theorem 2.2 & Theorem 3.2)
//    Parallel batch evaluation of exact intra-shell vs. inter-shell difference
//    saturation and the Exact Bilinear Cross-Shell Quadratic Form LHS:
//        \sum_{k=1}^M 2^{M-k} [ m_k(m_k - 1) + 2 m_k A(2^{k-1}) ]
//    alongside the unweighted total pair count:
//        \sum_{k=1}^M \binom{m_k}{2} + \sum_{1 \le j < k \le M} m_j m_k
//        = \binom{A(2^M)}{2} - m_0 (A(2^M) - m_0).
//    Bit-exact twin of `profile_sidon_cross_shell_cpu` (`src/verify_cpu.rs`).
// ============================================================================
extern "C" __global__ void k_profile_sidon_cross_shell(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ intra_pairs_out,
    u64* __restrict__ inter_pairs_out,
    u64* __restrict__ scaled_bilinear_lhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u64 m0 = (mask >> 1U) & 1ULL;
    u64 prefix_prev = m0;
    u64 intra_pairs = 0ULL;
    u64 inter_pairs = 0ULL;
    u64 scaled_bilinear = 0ULL;

    for (u32 k = 1U; k <= m_scales; k++) {
        u32 lo = (1U << (k - 1U)) + 1U;
        u32 hi = (1U << k);
        if (hi > 63U) {
            hi = 63U;
        }
        u64 shell_mask = 0ULL;
        if (lo <= hi) {
            u64 hi_mask = (hi == 63U) ? (~0ULL) : ((1ULL << (hi + 1U)) - 1ULL);
            u64 lo_mask = (1ULL << lo) - 1ULL;
            shell_mask = hi_mask & ~lo_mask;
        }
        u64 mk = (u64)__popcll(mask & shell_mask);
        u64 binom_mk = (mk >= 2ULL) ? (mk * (mk - 1ULL) / 2ULL) : 0ULL;
        u64 cross_mk = mk * prefix_prev;
        intra_pairs += binom_mk;
        inter_pairs += cross_mk;
        scaled_bilinear += (mk * (mk >= 1ULL ? (mk - 1ULL) : 0ULL) + 2ULL * cross_mk)
                           << (m_scales - k);
        prefix_prev += mk;
    }

    is_sidon_out[tid]            = ok;
    intra_pairs_out[tid]         = intra_pairs;
    inter_pairs_out[tid]         = inter_pairs;
    scaled_bilinear_lhs_out[tid] = scaled_bilinear;
}

// ============================================================================
// 4. `k_eval_sidon_subwindow_spectrum`: Phase 3 (Theorem 2.3 & Theorem 3.3)
//    Parallel batch evaluation of the Exact Sub-Window Short-Difference
//    Projection & Multi-Resolution Packing inequalities across sub-window
//    scales r = 1, ..., m_scales:
//      * `subwindow_pairs_out[tid]`: exact count of intra-sub-interval pairs
//        \sum_{k=r+1}^{m_scales} \sum_{b=1}^{2^{k-1-r}} \binom{|A \cap J_{k,r,b}|}{2}
//        at scale r = r_target, each of which generates a distinct difference
//        d = y - x <= 2^{r_target} - 1 disjoint from \mathcal{D}_A(2^{r_target}).
//      * `short_diffs_out[tid]`: exact count of positive differences in
//        \mathcal{D}_A(2^M) \cap [1, 2^{r_target} - 1], computed in O(1) via
//        `__popcll(diff_mask & ((1ULL << (1U << r_target)) - 2ULL))`.
//      * `pigeonhole_lb_out[tid]`: exact integer pigeonhole lower bound
//        \binom{A(2^{r_target})}{2} + \sum_{k=r_target+1}^{m_scales} L(m_k, 2^{k-1-r_target}),
//        where L(m, B) = B \binom{q}{2} + t q (m = q B + t, 0 <= t < B).
//    Bit-exact twin of `eval_sidon_subwindow_spectrum_cpu` (`src/verify_cpu.rs`).
// ============================================================================
extern "C" __global__ void k_eval_sidon_subwindow_spectrum(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ pigeonhole_lb_out,
    u64* __restrict__ exact_subwindow_lhs_out,
    u64* __restrict__ short_diffs_rhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    // Prefix count A(2^r_target) in [1, 2^r_target]
    u32 r_hi = (1U << r_target);
    if (r_hi > 63U) {
        r_hi = 63U;
    }
    u64 prefix_mask = (r_hi == 63U) ? (~1ULL) : (((1ULL << (r_hi + 1U)) - 1ULL) & ~1ULL);
    u64 a_2r = (u64)__popcll(mask & prefix_mask);
    u64 base_pairs = (a_2r >= 2ULL) ? (a_2r * (a_2r - 1ULL) / 2ULL) : 0ULL;

    u64 pigeon_sum = base_pairs;
    u64 exact_sub_sum = base_pairs;
    u32 w = (1U << r_target);

    for (u32 k = r_target + 1U; k <= m_scales; k++) {
        u32 num_bins = 1U << (k - 1U - r_target);
        u32 shell_start = (1U << (k - 1U)) + 1U;
        u64 mk = 0ULL;
        for (u32 b = 0U; b < num_bins; b++) {
            u32 lo = shell_start + b * w;
            u32 hi = lo + w - 1U;
            if (lo > 63U) {
                break;
            }
            if (hi > 63U) {
                hi = 63U;
            }
            u64 hi_m = (hi == 63U) ? (~0ULL) : ((1ULL << (hi + 1U)) - 1ULL);
            u64 lo_m = (1ULL << lo) - 1ULL;
            u64 cb = (u64)__popcll(mask & (hi_m & ~lo_m));
            mk += cb;
            if (cb >= 2ULL) {
                exact_sub_sum += cb * (cb - 1ULL) / 2ULL;
            }
        }
        u64 B = (u64)num_bins;
        u64 q = mk / B;
        u64 t = mk % B;
        u64 lb = B * (q >= 2ULL ? (q * (q - 1ULL) / 2ULL) : 0ULL) + t * q;
        pigeon_sum += lb;
    }

    // Count exact positive differences in [1, 2^r_target - 1]
    u32 max_d = (1U << r_target) - 1U;
    if (max_d > 62U) {
        max_d = 62U;
    }
    u64 d_mask_win = ((1ULL << (max_d + 1U)) - 1ULL) & ~1ULL;
    u64 short_diffs = (u64)__popcll(diff_mask & d_mask_win);

    is_sidon_out[tid]            = ok;
    pigeonhole_lb_out[tid]       = pigeon_sum;
    exact_subwindow_lhs_out[tid] = exact_sub_sum;
    short_diffs_rhs_out[tid]     = short_diffs;
}

// ============================================================================
// 5. `k_eval_base_b_joint_fejer`: Phase 5 (Theorem 2.5)
//    Parallel batch evaluation of the Direct Prefix-to-Block Unsplit Fejér
//    Coupling & Base-B (B = 2^{1/L}, L = 2^refinement_levels) Window-Block
//    Refinement Hierarchy across 64-bit Sidon bitmasks:
//      * `joint_l1_lhs_out[tid]`: F_W^{joint}(A; B=2) on the unsplit window
//        partition U_{<=r}^{(W)} = [1, t_r + W - 1] and U_k^{(W)} = [t_{k-1} + W, t_k + W - 1].
//      * `joint_base_b_lhs_out[tid]`: F_W^{joint}(A; B=2^{1/L}) after subdividing
//        each window block via `refinement_levels` binary halvings (L = 2^refinement_levels).
//      * `exact_fejer_rhs_out[tid]`: exact Fejér short-difference sum
//        \sum_{d \in \mathcal{D}_A(N) \cap [1, W - 1]} (W - d) <= \binom{W}{2}.
//    Bit-exact twin of `eval_base_b_joint_fejer_cpu` (`src/verify_cpu.rs`).
// ============================================================================
__device__ __forceinline__ u64 d_exact_pigeonhole_lb(u64 balls, u64 bins) {
    if (bins == 0ULL || balls <= bins) {
        return 0ULL;
    }
    u64 q = balls / bins;
    u64 t = balls % bins;
    return bins * (q >= 2ULL ? (q * (q - 1ULL) / 2ULL) : 0ULL) + t * q;
}

__device__ __forceinline__ u64 d_eval_refined_block(
    u64 mask,
    u32 n_hi,
    u32 w,
    u32 u_start,
    u32 u_end,
    u32 levels
) {
    if (u_start > u_end) {
        return 0ULL;
    }
    // Iterative binary subdivision of [u_start, u_end] up to depth `levels` (<= 4, i.e. 16 sub-blocks)
    u32 seg_lo[16];
    u32 seg_hi[16];
    seg_lo[0] = u_start;
    seg_hi[0] = u_end;
    u32 num_segs = 1U;
    u32 eff_levels = levels > 4U ? 4U : levels;

    for (u32 lv = 0U; lv < eff_levels; lv++) {
        u32 next_lo[16];
        u32 next_hi[16];
        u32 next_cnt = 0U;
        for (u32 i = 0U; i < num_segs; i++) {
            u32 a = seg_lo[i];
            u32 b = seg_hi[i];
            u32 mid = (a + b) >> 1U;
            if (b > a) {
                next_lo[next_cnt] = a;
                next_hi[next_cnt] = mid;
                next_cnt++;
                next_lo[next_cnt] = mid + 1U;
                next_hi[next_cnt] = b;
                next_cnt++;
            } else {
                next_lo[next_cnt] = a;
                next_hi[next_cnt] = b;
                next_cnt++;
            }
        }
        for (u32 i = 0U; i < next_cnt; i++) {
            seg_lo[i] = next_lo[i];
            seg_hi[i] = next_hi[i];
        }
        num_segs = next_cnt;
    }

    u64 total_lb = 0ULL;
    for (u32 i = 0U; i < num_segs; i++) {
        u32 a = seg_lo[i];
        u32 b = seg_hi[i];
        u64 s_blk = 0ULL;
        for (u32 u = a; u <= b; u++) {
            u32 lo_x = (u >= w) ? (u - w + 1U) : 1U;
            u32 hi_x = (u <= n_hi) ? u : n_hi;
            if (lo_x <= hi_x) {
                u64 hi_m = (hi_x >= 63U) ? (~0ULL) : ((1ULL << (hi_x + 1U)) - 1ULL);
                u64 lo_m = (1ULL << lo_x) - 1ULL;
                s_blk += (u64)__popcll(mask & (hi_m & ~lo_m));
            }
        }
        total_lb += d_exact_pigeonhole_lb(s_blk, (u64)(b - a + 1U));
    }
    return total_lb;
}

extern "C" __global__ void k_eval_base_b_joint_fejer(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ joint_l1_lhs_out,
    u64* __restrict__ joint_base_b_lhs_out,
    u64* __restrict__ exact_fejer_rhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w = (1U << r_target);
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    // Prefix block U_{<=r}^{(W)} = [1, t_r + W - 1]
    u64 l1_sum = d_eval_refined_block(mask, n_hi, w, 1U, t_r + w - 1U, 0U);
    u64 lb_sum = d_eval_refined_block(mask, n_hi, w, 1U, t_r + w - 1U, refinement_levels);

    // Unsplit shell blocks U_k^{(W)} = [t_{k-1} + W, t_k + W - 1] for k = r_target + 1 .. m_scales
    for (u32 k = r_target + 1U; k <= m_scales; k++) {
        u32 t_prev = (1U << (k - 1U));
        u32 t_curr = (1U << k);
        if (t_prev >= n_hi) {
            break;
        }
        if (t_curr > n_hi) {
            t_curr = n_hi;
        }
        u32 u_start = t_prev + w;
        u32 u_end = t_curr + w - 1U;
        l1_sum += d_eval_refined_block(mask, n_hi, w, u_start, u_end, 0U);
        lb_sum += d_eval_refined_block(mask, n_hi, w, u_start, u_end, refinement_levels);
    }

    // Exact Fejér short-difference RHS: \sum_{d \in D_A(N) \cap [1, W - 1]} (W - d)
    u64 fejer_rhs = 0ULL;
    for (u32 d = 1U; d < w && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            fejer_rhs += (u64)(w - d);
        }
    }

    is_sidon_out[tid]         = ok;
    joint_l1_lhs_out[tid]     = l1_sum;
    joint_base_b_lhs_out[tid] = lb_sum;
    exact_fejer_rhs_out[tid]  = fejer_rhs;
}

// Helper for Kernel 6: Evaluate scaled continuous parabola S*(S - B), scaled discrete
// Legendre–Fenchel integrality excess t*(B - t), and exact pigeonhole bound L(S, B)
// on [u_start, u_end] subdivided into L = 2^refinement_levels sub-blocks.
__device__ __forceinline__ void d_eval_refined_block_integrality(
    u64 mask,
    u32 n_hi,
    u32 w,
    u32 u_start,
    u32 u_end,
    u32 refinement_levels,
    u64* cont_acc,
    u64* int_acc,
    u64* l_acc
) {
    if (u_start > u_end) {
        return;
    }
    u32 num_sub = 1U << refinement_levels;
    u32 total_len = u_end - u_start + 1U;
    if (num_sub > total_len) {
        num_sub = total_len;
    }
    for (u32 s = 0U; s < num_sub; s++) {
        u32 sub_a = u_start + (s * total_len) / num_sub;
        u32 sub_b = u_start + ((s + 1U) * total_len) / num_sub - 1U;
        if (sub_a > sub_b) {
            continue;
        }
        u64 bins = (u64)(sub_b - sub_a + 1U);
        u64 s_val = 0ULL;
        for (u32 u = sub_a; u <= sub_b; u++) {
            u32 hi = (u <= n_hi) ? u : n_hi;
            u32 lo = (u > w) ? (u - w) : 0U;
            u64 hi_mask = (hi == 63U) ? ~1ULL : (((1ULL << (hi + 1U)) - 1ULL) & ~1ULL);
            u64 lo_mask = (lo == 0U) ? 0ULL : (((1ULL << (lo + 1U)) - 1ULL) & ~1ULL);
            u64 win_mask = mask & (hi_mask & ~lo_mask);
            s_val += (u64)__popcll((long long)win_mask);
        }
        if (s_val > bins) {
            u64 t_val = s_val % bins;
            *cont_acc += s_val * (s_val - bins);
            *int_acc += t_val * (bins - t_val);
            *l_acc += d_exact_pigeonhole_lb(s_val, bins);
        }
    }
}

// ----------------------------------------------------------------------------
// Kernel 6: Phase 6 (Theorem 2.6) Exact Discrete Legendre–Fenchel Integrality
// Excess & Coupled Unsplit Fejér Evaluator
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_coupled_lag_integrality_fejer(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ scaled_cont_parabola_out,
    u64* __restrict__ scaled_integrality_excess_out,
    u64* __restrict__ exact_pigeonhole_lhs_out,
    u64* __restrict__ exact_fejer_rhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w = (1U << r_target);
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 cont_sum = 0ULL;
    u64 int_sum = 0ULL;
    u64 l_sum = 0ULL;

    // Prefix block U_{<=r}^{(W)} = [1, t_r + W - 1]
    d_eval_refined_block_integrality(
        mask,
        n_hi,
        w,
        1U,
        t_r + w - 1U,
        refinement_levels,
        &cont_sum,
        &int_sum,
        &l_sum
    );

    // Unsplit shell blocks U_k^{(W)} = [t_{k-1} + W, t_k + W - 1] for k = r_target + 1 .. m_scales
    for (u32 k = r_target + 1U; k <= m_scales; k++) {
        u32 t_prev = (1U << (k - 1U));
        u32 t_curr = (1U << k);
        if (t_prev >= n_hi) {
            break;
        }
        if (t_curr > n_hi) {
            t_curr = n_hi;
        }
        u32 u_start = t_prev + w;
        u32 u_end = t_curr + w - 1U;
        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            u_start,
            u_end,
            refinement_levels,
            &cont_sum,
            &int_sum,
            &l_sum
        );
    }

    // Exact Fejér short-difference RHS: \sum_{d \in D_A(N) \cap [1, W - 1]} (W - d)
    u64 fejer_rhs = 0ULL;
    for (u32 d = 1U; d < w && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            fejer_rhs += (u64)(w - d);
        }
    }

    is_sidon_out[tid]                  = ok;
    scaled_cont_parabola_out[tid]      = cont_sum;
    scaled_integrality_excess_out[tid] = int_sum;
    exact_pigeonhole_lhs_out[tid]      = l_sum;
    exact_fejer_rhs_out[tid]           = fejer_rhs;
}

// ----------------------------------------------------------------------------
// Kernel 7: Phase 7 (Theorem 2.7) Higher-Order Cesàro–Fejér Order-2 &
// Upper-Octave Iterated Sliding-Window Evaluator
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_cesaro_fejer_order2(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ cf2_full_pigeonhole_lhs_out,
    u64* __restrict__ cf2_oct_pigeonhole_lhs_out,
    u64* __restrict__ cf2_full_exact_rhs_out,
    u64* __restrict__ cf2_oct_exact_rhs_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w_cap = (1U << r_target);
    u32 w_half = (w_cap >> 1U);
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 full_lhs = 0ULL;
    u64 oct_lhs = 0ULL;

    for (u32 w = 1U; w <= w_cap; w++) {
        u64 c_w = 0ULL;
        u64 i_w = 0ULL;
        u64 l_w = 0ULL;

        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            1U,
            t_r + w - 1U,
            refinement_levels,
            &c_w,
            &i_w,
            &l_w
        );

        for (u32 k = r_target + 1U; k <= m_scales; k++) {
            u32 t_prev = (1U << (k - 1U));
            u32 t_curr = (1U << k);
            if (t_prev >= n_hi) {
                break;
            }
            if (t_curr > n_hi) {
                t_curr = n_hi;
            }
            u32 u_start = t_prev + w;
            u32 u_end = t_curr + w - 1U;
            d_eval_refined_block_integrality(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
                &c_w,
                &i_w,
                &l_w
            );
        }

        full_lhs += l_w;
        if (w > w_half) {
            oct_lhs += l_w;
        }
    }

    // Exact Order-2 Cesàro–Fejér difference RHS: \sum_{d \in D_A(N) \cap [1, W - 1]} \binom{W - d + 1}{2}
    // and Upper-Octave difference RHS: \binom{W - d + 1}{2} - \binom{\max(W/2 - d + 1, 0)}{2}
    u64 full_rhs = 0ULL;
    u64 half_rhs = 0ULL;
    for (u32 d = 1U; d < w_cap && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            u64 rem_full = (u64)(w_cap - d);
            full_rhs += (rem_full * (rem_full + 1ULL)) >> 1U;
            if (d < w_half) {
                u64 rem_half = (u64)(w_half - d);
                half_rhs += (rem_half * (rem_half + 1ULL)) >> 1U;
            }
        }
    }

    is_sidon_out[tid]                = ok;
    cf2_full_pigeonhole_lhs_out[tid] = full_lhs;
    cf2_oct_pigeonhole_lhs_out[tid]  = oct_lhs;
    cf2_full_exact_rhs_out[tid]      = full_rhs;
    cf2_oct_exact_rhs_out[tid]       = full_rhs - half_rhs;
}

// ----------------------------------------------------------------------------
// Kernel 8: Phase 8 (Theorem 2.8) Band-Restricted Upper-Octave Cesàro–Fejér
// Evaluator on w in [W_0 + 1, W]
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_band_restricted_cesaro_fejer(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 w0_cutoff,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ band_pigeonhole_lhs_out,
    u64* __restrict__ band_scaled_integrality_out,
    u64* __restrict__ band_exact_rhs_out,
    u64* __restrict__ band_binom_cap_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w_cap = (1U << r_target);
    u32 w0 = (w0_cutoff < w_cap) ? w0_cutoff : (w_cap - 1U);
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 band_lhs = 0ULL;
    u64 band_int = 0ULL;

    for (u32 w = w0 + 1U; w <= w_cap; w++) {
        u64 c_w = 0ULL;
        u64 i_w = 0ULL;
        u64 l_w = 0ULL;

        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            1U,
            t_r + w - 1U,
            refinement_levels,
            &c_w,
            &i_w,
            &l_w
        );

        for (u32 k = r_target + 1U; k <= m_scales; k++) {
            u32 t_prev = (1U << (k - 1U));
            u32 t_curr = (1U << k);
            if (t_prev >= n_hi) {
                break;
            }
            if (t_curr > n_hi) {
                t_curr = n_hi;
            }
            u32 u_start = t_prev + w;
            u32 u_end = t_curr + w - 1U;
            d_eval_refined_block_integrality(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
                &c_w,
                &i_w,
                &l_w
            );
        }

        band_lhs += l_w;
        band_int += i_w;
    }

    // Exact Band-Restricted Cesàro–Fejér difference RHS:
    //   \sum_{d \in D_A(N) \cap [1, W - 1]} [\binom{W - d + 1}{2} - \binom{\max(W_0 - d + 1, 0)}{2}]
    u64 full_rhs = 0ULL;
    u64 sub_rhs = 0ULL;
    for (u32 d = 1U; d < w_cap && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            u64 rem_full = (u64)(w_cap - d);
            full_rhs += (rem_full * (rem_full + 1ULL)) >> 1U;
            if (d < w0) {
                u64 rem_w0 = (u64)(w0 - d);
                sub_rhs += (rem_w0 * (rem_w0 + 1ULL)) >> 1U;
            }
        }
    }

    u64 w_u64 = (u64)w_cap;
    u64 w0_u64 = (u64)w0;
    u64 cap_full = ((w_u64 - 1ULL) * w_u64 * (w_u64 + 1ULL)) / 6ULL;
    u64 cap_w0 = (w0_u64 >= 2ULL) ? (((w0_u64 - 1ULL) * w0_u64 * (w0_u64 + 1ULL)) / 6ULL) : 0ULL;

    is_sidon_out[tid]                = ok;
    band_pigeonhole_lhs_out[tid]     = band_lhs;
    band_scaled_integrality_out[tid] = band_int;
    band_exact_rhs_out[tid]          = full_rhs - sub_rhs;
    band_binom_cap_out[tid]          = cap_full - cap_w0;
}

// ----------------------------------------------------------------------------
// Kernel 9: Phase 9 (Theorem 2.9) Multi-Octave K-Octave Integrality Cascade &
// Cesàro–Fejér Evaluator on w in [floor(2^{-K} W) + 1, W]
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_multioctave_cesaro_fejer(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 k_octaves,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ multioct_pigeonhole_lhs_out,
    u64* __restrict__ multioct_scaled_integrality_out,
    u64* __restrict__ multioct_exact_rhs_out,
    u64* __restrict__ multioct_binom_cap_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w_cap = (1U << r_target);
    u32 shift_k = (k_octaves <= r_target) ? k_octaves : r_target;
    u32 w0 = (shift_k >= r_target) ? 0U : (w_cap >> shift_k);
    if (w0 >= w_cap) {
        w0 = w_cap - 1U;
    }
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 multioct_lhs = 0ULL;
    u64 multioct_int = 0ULL;

    for (u32 w = w0 + 1U; w <= w_cap; w++) {
        u64 c_w = 0ULL;
        u64 i_w = 0ULL;
        u64 l_w = 0ULL;

        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            1U,
            t_r + w - 1U,
            refinement_levels,
            &c_w,
            &i_w,
            &l_w
        );

        for (u32 k = r_target + 1U; k <= m_scales; k++) {
            u32 t_prev = (1U << (k - 1U));
            u32 t_curr = (1U << k);
            if (t_prev >= n_hi) {
                break;
            }
            if (t_curr > n_hi) {
                t_curr = n_hi;
            }
            u32 u_start = t_prev + w;
            u32 u_end = t_curr + w - 1U;
            d_eval_refined_block_integrality(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
                &c_w,
                &i_w,
                &l_w
            );
        }

        multioct_lhs += l_w;
        multioct_int += i_w;
    }

    u64 full_rhs = 0ULL;
    u64 sub_rhs = 0ULL;
    for (u32 d = 1U; d < w_cap && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            u64 rem_full = (u64)(w_cap - d);
            full_rhs += (rem_full * (rem_full + 1ULL)) >> 1U;
            if (d < w0) {
                u64 rem_w0 = (u64)(w0 - d);
                sub_rhs += (rem_w0 * (rem_w0 + 1ULL)) >> 1U;
            }
        }
    }

    u64 w_u64 = (u64)w_cap;
    u64 w0_u64 = (u64)w0;
    u64 cap_full = ((w_u64 - 1ULL) * w_u64 * (w_u64 + 1ULL)) / 6ULL;
    u64 cap_w0 = (w0_u64 >= 2ULL) ? (((w0_u64 - 1ULL) * w0_u64 * (w0_u64 + 1ULL)) / 6ULL) : 0ULL;

    is_sidon_out[tid]                    = ok;
    multioct_pigeonhole_lhs_out[tid]     = multioct_lhs;
    multioct_scaled_integrality_out[tid] = multioct_int;
    multioct_exact_rhs_out[tid]          = full_rhs - sub_rhs;
    multioct_binom_cap_out[tid]          = cap_full - cap_w0;
}

// ----------------------------------------------------------------------------
// Kernel 10: Phase 10 (Theorem 2.10) Scale-Invariant Haar-Measure Multi-Octave
// Cesàro–Fejér Evaluator with Integer Haar Weights floor(W / w) on w in [W_0 + 1, W]
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_haar_multioctave_cesaro_fejer(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 k_octaves,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ haar_pigeonhole_lhs_out,
    u64* __restrict__ haar_scaled_integrality_out,
    u64* __restrict__ haar_exact_rhs_out,
    u64* __restrict__ haar_binom_cap_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w_cap = (1U << r_target);
    u32 shift_k = (k_octaves <= r_target) ? k_octaves : r_target;
    u32 w0 = (shift_k >= r_target) ? 0U : (w_cap >> shift_k);
    if (w0 >= w_cap) {
        w0 = w_cap - 1U;
    }
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 haar_lhs = 0ULL;
    u64 haar_int = 0ULL;
    u64 haar_cap = 0ULL;

    for (u32 w = w0 + 1U; w <= w_cap; w++) {
        u64 wt = (u64)(w_cap / w);
        u64 c_w = 0ULL;
        u64 i_w = 0ULL;
        u64 l_w = 0ULL;

        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            1U,
            t_r + w - 1U,
            refinement_levels,
            &c_w,
            &i_w,
            &l_w
        );

        for (u32 k = r_target + 1U; k <= m_scales; k++) {
            u32 t_prev = (1U << (k - 1U));
            u32 t_curr = (1U << k);
            if (t_prev >= n_hi) {
                break;
            }
            if (t_curr > n_hi) {
                t_curr = n_hi;
            }
            u32 u_start = t_prev + w;
            u32 u_end = t_curr + w - 1U;
            d_eval_refined_block_integrality(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
                &c_w,
                &i_w,
                &l_w
            );
        }

        haar_lhs += wt * l_w;
        haar_int += wt * i_w;
        u64 w_u64 = (u64)w;
        haar_cap += wt * ((w_u64 * (w_u64 - 1ULL)) >> 1U);
    }

    u64 haar_rhs = 0ULL;
    for (u32 d = 1U; d < w_cap && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            u32 w_start = (w0 + 1U > d + 1U) ? (w0 + 1U) : (d + 1U);
            u64 psi_d = 0ULL;
            for (u32 w = w_start; w <= w_cap; w++) {
                u64 wt = (u64)(w_cap / w);
                psi_d += wt * (u64)(w - d);
            }
            haar_rhs += psi_d;
        }
    }

    is_sidon_out[tid]                = ok;
    haar_pigeonhole_lhs_out[tid]     = haar_lhs;
    haar_scaled_integrality_out[tid] = haar_int;
    haar_exact_rhs_out[tid]          = haar_rhs;
    haar_binom_cap_out[tid]          = haar_cap;
}

// ----------------------------------------------------------------------------
// Kernel 11: Phase 11 (Theorem 2.11) Unconditional Diophantine Scaled-Reciprocal
// Haar-Measure Multi-Octave Cesàro–Fejér Closure Evaluator with Integer Weights
// floor(W^2 / w) on w in [W_0 + 1, W]
// ----------------------------------------------------------------------------
extern "C" __global__ void k_eval_unconditional_diophantine_haar_closure(
    const u64* __restrict__ masks,
    u32 num_masks,
    u32 m_scales,
    u32 r_target,
    u32 k_octaves,
    u32 refinement_levels,
    u32* __restrict__ is_sidon_out,
    u64* __restrict__ dioph_pigeonhole_lhs_out,
    u64* __restrict__ dioph_scaled_integrality_out,
    u64* __restrict__ dioph_exact_rhs_out,
    u64* __restrict__ dioph_binom_cap_out
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_masks) {
        return;
    }

    u64 mask = masks[tid] & ~1ULL;
    u64 m = mask;
    u64 rev_mask = 0ULL;
    u64 diff_mask = 0ULL;
    u32 ok = 1U;

    while (m != 0ULL) {
        u32 x = (u32)__ffsll((long long)m) - 1U;
        m &= (m - 1ULL);
        u64 new_diffs = rev_mask >> (63U - x);
        if ((new_diffs & diff_mask) != 0ULL) {
            ok = 0U;
        }
        diff_mask |= new_diffs;
        rev_mask |= (1ULL << (63U - x));
    }

    u32 n_hi = (1U << m_scales);
    if (n_hi > 63U) {
        n_hi = 63U;
    }
    u32 w_cap = (1U << r_target);
    u64 w_sq = (u64)w_cap * (u64)w_cap;
    u32 shift_k = (k_octaves <= r_target) ? k_octaves : r_target;
    u32 w0 = (shift_k >= r_target) ? 0U : (w_cap >> shift_k);
    if (w0 >= w_cap) {
        w0 = w_cap - 1U;
    }
    u32 t_r = (1U << r_target);
    if (t_r > n_hi) {
        t_r = n_hi;
    }

    u64 dioph_lhs = 0ULL;
    u64 dioph_int = 0ULL;
    u64 dioph_cap = 0ULL;

    for (u32 w = w0 + 1U; w <= w_cap; w++) {
        u64 wt = w_sq / (u64)w;
        u64 c_w = 0ULL;
        u64 i_w = 0ULL;
        u64 l_w = 0ULL;

        d_eval_refined_block_integrality(
            mask,
            n_hi,
            w,
            1U,
            t_r + w - 1U,
            refinement_levels,
            &c_w,
            &i_w,
            &l_w
        );

        for (u32 k = r_target + 1U; k <= m_scales; k++) {
            u32 t_prev = (1U << (k - 1U));
            u32 t_curr = (1U << k);
            if (t_prev >= n_hi) {
                break;
            }
            if (t_curr > n_hi) {
                t_curr = n_hi;
            }
            u32 u_start = t_prev + w;
            u32 u_end = t_curr + w - 1U;
            d_eval_refined_block_integrality(
                mask,
                n_hi,
                w,
                u_start,
                u_end,
                refinement_levels,
                &c_w,
                &i_w,
                &l_w
            );
        }

        dioph_lhs += wt * l_w;
        dioph_int += wt * i_w;
        u64 w_u64 = (u64)w;
        dioph_cap += wt * ((w_u64 * (w_u64 - 1ULL)) >> 1U);
    }

    u64 dioph_rhs = 0ULL;
    for (u32 d = 1U; d < w_cap && d <= 62U; d++) {
        if ((diff_mask & (1ULL << d)) != 0ULL) {
            u32 w_start = (w0 + 1U > d + 1U) ? (w0 + 1U) : (d + 1U);
            u64 psi_d = 0ULL;
            for (u32 w = w_start; w <= w_cap; w++) {
                u64 wt = w_sq / (u64)w;
                psi_d += wt * (u64)(w - d);
            }
            dioph_rhs += psi_d;
        }
    }

    is_sidon_out[tid]                 = ok;
    dioph_pigeonhole_lhs_out[tid]     = dioph_lhs;
    dioph_scaled_integrality_out[tid] = dioph_int;
    dioph_exact_rhs_out[tid]          = dioph_rhs;
    dioph_binom_cap_out[tid]          = dioph_cap;
}
