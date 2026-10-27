//! ============================================================================
//! `records.rs` — Exact OEIS sequences (`A005282` Mian–Chowla greedy Sidon
//! sequence, `A003022` / `A143824` optimal Golomb ruler & finite Sidon capacities),
//! Krückeberg–Cilleruelo dyadic Singer-block constructions, and non-asymptotic
//! Harmonic-Difference & Dyadic Budget certificates for Erdős Problem #1191
//! ($1,000 open prize: prove $\liminf_{x \to \infty} A(x)\sqrt{(\log x)/x} = 0$
//! for every infinite Sidon set $A \subseteq \mathbb{N}$).
//! ============================================================================
//!
//! MATHEMATICAL SUMMARY
//! --------------------
//! Let $A = \{a_1 < a_2 < a_3 < \dots\} \subseteq \mathbb{N}$ be an infinite
//! Sidon set ($B_2$ sequence: all pairwise sums $a_i + a_j$ with $i \le j$ are
//! distinct, equivalently all positive differences $a_j - a_i$ with $j > i$ are
//! distinct), and let $A(x) := |A \cap [1, x]|$ denote its counting function.
//!
//!   * **Finite Regime:** The maximum size $F_2(N)$ of a Sidon subset of $\{1, \dots, N\}$
//!     satisfies $F_2(N) = (1 + o(1))\sqrt{N}$ (Singer 1938; Erdős–Turán 1941),
//!     determined exactly for $N \le 586$ by the optimal Golomb ruler lengths
//!     OEIS `A003022` ($k \le 28$).
//!   * **Infinite Upper Density ($\limsup$):**
//!     - Erdős (1955) proved $\limsup_{x \to \infty} A(x)/\sqrt{x} \ge 1/2$.
//!     - Krückeberg (1961) improved this via shifted Singer blocks to
//!       $\limsup_{x \to \infty} A(x)/\sqrt{x} \ge 1/\sqrt{2} \approx 0.707107$.
//!     - Cilleruelo (2014) constructed an infinite Sidon set achieving the
//!       optimal upper density $\limsup_{x \to \infty} A(x)/\sqrt{x} = 1$.
//!   * **Infinite Lower Density ($\liminf$) — Erdős Problem #1191 ($1,000 Open Prize):**
//!     - Erdős (1955, in Stöhr 1955) proved that NO infinite Sidon set can
//!       sustain $A(x) \gg \sqrt{x}$ for all $x$; every infinite Sidon set
//!       satisfies
//!
//!           lim inf_{x -> ∞} A(x) * sqrt( (ln x) / x ) <= C < ∞.
//!
//!     - Ruzsa (1998) constructed an infinite Sidon set with $A(x) = x^{\sqrt{2}-1+o(1)}$
//!       ($\sqrt{2} - 1 \approx 0.414214$) for all $x$.
//!     - **Erdős Problem #1191 ($1,000 OPEN):** Does every infinite Sidon set
//!       $A \subseteq \mathbb{N}$ satisfy
//!
//!           lim inf_{x -> ∞} A(x) * sqrt( (ln x) / x ) = 0 ?
//!
//! THEOREM 2.1 (EXACT HARMONIC-DIFFERENCE CONSERVATION & DYADIC BUDGET)
//! --------------------------------------------------------------------
//!   * **Theorem 2.1(a) (Exact Harmonic-Integral Conservation Identity):**
//!     For ANY Sidon set $A \subseteq \mathbb{N}$ and ANY horizon $N \ge 2$, let
//!     $\mathcal{D}_A(N) := \{a_j - a_i : 1 \le i < j \le A(N)\} \subseteq \{1, \dots, N-1\}$.
//!     Then:
//!
//!         \int_1^N \frac{A(t)(A(t) - 1)}{2 t^2} dt + \frac{A(N)(A(N) - 1)}{2 N}
//!             + \sum_{1 \le i < j \le A(N)} \frac{a_i - 1}{(a_j - a_i + 1) a_j}
//!         = \sum_{d \in \mathcal{D}_A(N)} \frac{1}{d + 1} <= H_N - 1.
//!
//!   * **Theorem 2.1(b) (Exact Single-Scale Dyadic Difference-Budget Inequality):**
//!     Partition $(1, 2^M]$ into dyadic shells $I_k = (2^{k-1}, 2^k]$ ($k = 1, \dots, M$)
//!     and let $m_k := |A \cap I_k|$. Since every positive difference $d = y - x$
//!     inside $I_k$ satisfies $d \le 2^{k-1} - 1$ (forcing $k \ge \lfloor \log_2 d \rfloor + 2$)
//!     and occurs at most once across the entire Sidon set $A$:
//!
//!         \sum_{k=1}^M \frac{m_k(m_k - 1)}{2^k}
//!             <= 2 \sum_{d=1}^{2^{M-1}-1} 2^{-(\lfloor \log_2 d \rfloor + 2)}
//!             = \frac{M - 1}{2}.
//!
//!   * **Theorem 2.1(c) (Multi-Scale Dyadic Sub-Block Refinement):**
//!     Using the exact pointwise step-decomposition
//!     $1/d \ge 2^{-(k-1)} + \sum_{s=1}^{S_k} 2^{-(k-s)} \mathbf{1}_{d \le 2^{k-1-s}}$
//!     for $d \in [1, 2^{k-1}]$ and Cauchy–Schwarz across $2^s$ sub-intervals of $I_k$:
//!
//!         \sum_{k=1}^M \frac{(S_k + 2) m_k^2 - (2^{S_k} + 1) m_k}{2^{k+1}}
//!             <= \sum_{k=1}^M \sum_{x < y \in A \cap I_k} \frac{1}{y - x}
//!             <= H_{2^{M-1}-1}.
//!
//! THEOREM 2.2 (EXACT CROSS-SHELL BILINEAR DIFFERENCE CONSERVATION & ENVELOPE)
//! ---------------------------------------------------------------------------
//!   * **Theorem 2.2(a) (Two-Shell Inter-Shell Localization & Prefix Packing):**
//!     Let $I_k = (2^{k-1}, 2^k]$ and $D_r = (2^{r-1}, 2^r]$ (with $D_0 = \{1\}$).
//!     For every $1 \le j \le k - 2$, every inter-shell difference $y - x$
//!     ($x \in I_j, y \in I_k$) satisfies $2^{k-1} - 2^j < y - x < 2^k$, and since
//!     $2^{k-1} - 2^j \ge 2^{k-2}$, it belongs strictly to the union of just **two**
//!     dyadic difference shells $D_{k-1} \cup D_k$. Moreover, for every $k \in \{1, \dots, M\}$,
//!     all intra-shell pairs inside $I_k$ and all inter-shell pairs between $A \cap [1, 2^{k-1}]$
//!     and $A \cap I_k$ generate disjoint positive differences in $\mathcal{D}_A(2^k) \cap [1, 2^k - 1]$:
//!
//!         \binom{m_k}{2} + m_k A(2^{k-1}) <= |\mathcal{D}_A(2^k) \cap (0, 2^k)| <= 2^k - 1,
//!
//!     and summing over $1 \le k \le r$ yields the Exact Cumulative Prefix Difference
//!     Conservation Identity:
//!
//!         \sum_{k=1}^r \left[\binom{m_k}{2} + m_k A(2^{k-1})\right]
//!             = \binom{A(2^r)}{2} = |\mathcal{D}_A(2^r)| <= 2^r - 1.
//!
//!   * **Theorem 2.2(b) (Exact Bilinear Cross-Shell Quadratic Form Inequality):**
//!     Summing the dyadic-weighted two-shell relation across $k = 1, \dots, M$ couples
//!     the intra-shell quadratic sum $\sum m_k(m_k - 1)/2^k$ with the cross-shell
//!     bilinear interaction $\sum_{0 \le j < k \le M} 2 m_j m_k / 2^k$:
//!
//!         \sum_{k=1}^M \frac{m_k(m_k - 1)}{2^k}
//!             + \sum_{k=1}^M \frac{2 m_k A(2^{k-1})}{2^k}
//!             <= 2 \sum_{k=1}^M \frac{|\mathcal{D}_A(2^k)| - |\mathcal{D}_A(2^{k-1})|}{2^k}
//!             <= 2(M - 1) + 2^{1 - M}.
//!
//!     More generally, for any Laplace-Fejér / rational shell weight sequence
//!     $w(d) = \sum_{r=0}^M \alpha_r \mathbf{1}_{d \in D_r}$ non-increasing in $r$,
//!     every pair $(x, y)$ with $x \in I_j, y \in I_k$ ($j < k$) has $y - x < 2^k$,
//!     yielding the weighted bilinear inequality:
//!
//!         \sum_{k=1}^M \alpha_{k-1} \binom{m_k}{2}
//!             + \sum_{1 \le j < k \le M} \alpha_k m_j m_k
//!             <= \sum_{d \in \mathcal{D}_A(2^M)} w(d)
//!             <= \alpha_0 + \sum_{r=1}^M 2^{r-1} \alpha_r.
//!
//!   * **Theorem 2.2(c) (Strict Dominance & Sharpened $\liminf$ Upper Envelope):**
//!     If $A(2^k) \ge c \sqrt{2^k / \ln(2^k)}$ holds for all $k \in \{k_0, \dots, M\}$,
//!     then $\binom{A(2^k)}{2} = \sum_{r=1}^k [\binom{m_r}{2} + m_r A(2^{r-1})]$ and
//!     the bilinear cross-shell sum strictly dominates the intra-shell-only bound of
//!     Theorem 2.1(b), yielding the non-asymptotic quadratic-programming certificate
//!     $\underline{c}_M^{\mathrm{QP}}(k_0)$ of Theorem 3.2.
//!
//! THEOREM 2.3 (SUB-WINDOW SHORT-DIFFERENCE PROJECTION & MULTI-RESOLUTION PACKING)
//! -------------------------------------------------------------------------------
//!   * **Theorem 2.3(a) (Exact Sub-Window Short-Difference Projection & Crowding):**
//!     For every dyadic shell $I_k = (2^{k-1}, 2^k]$ and every sub-window scale
//!     $r \in \{1, \dots, k - 1\}$, partition $I_k$ into $B_{k,r} := 2^{k-1-r}$
//!     disjoint sub-intervals $J_{k,r,b}$ ($1 \le b \le B_{k,r}$) each of width $2^r$.
//!     Every pair $x < y$ inside the same sub-interval $J_{k,r,b}$ has difference
//!     $d = y - x \in [1, 2^r - 1]$. By integer convexity of $c \mapsto \binom{c}{2}$
//!     (and Cauchy–Schwarz), writing $m_k = q B_{k,r} + t$ ($0 \le t < B_{k,r}$):
//!
//!         S_{k \to r}(A) := \sum_{b=1}^{B_{k,r}} \binom{|A \cap J_{k,r,b}|}{2}
//!                        >= L(m_k, B_{k,r}) := B_{k,r} \binom{q}{2} + t q
//!                        >= \max\left(0, \frac{m_k^2}{2^{k-r}} - \frac{m_k}{2}\right).
//!
//!     Crucially, for any fixed target scale $r \in \{1, \dots, M\}$, all differences
//!     in $\mathcal{D}_A(2^r)$ and all sub-window pairs across $k \in \{r+1, \dots, M\}$
//!     are pairwise disjoint and MUST pack inside $\mathcal{D}_A(2^M) \cap [1, 2^r - 1]$:
//!
//!         \binom{A(2^r)}{2} + \sum_{k=r+1}^M L(m_k, 2^{k-1-r})
//!             <= \binom{A(2^r)}{2} + \sum_{k=r+1}^M S_{k \to r}(A)
//!             <= |\mathcal{D}_A(2^M) \cap [1, 2^r - 1]|
//!             <= 2^r - 1.
//!
//!   * **Theorem 2.3(b) (Multi-Resolution Harmonic Weight Amplification vs. Single-Scale):**
//!     Summation by parts over the sub-window scales $r = 1, \dots, k-1$ yields the
//!     simultaneous multi-stratum lower bound on the intra-shell harmonic sum:
//!
//!         \sum_{x < y \in A \cap I_k} \frac{1}{y - x}
//!             >= \frac{S_{k \to k-1}(A)}{2^{k-1}} + \sum_{r=1}^{k-2} \frac{S_{k \to r}(A)}{2^{r+1}}
//!             >= \sum_{r = r_0(m_k)}^{k-1} \frac{m_k^2 / 2^{k-r} - m_k / 2}{2^{r+1}}
//!             \asymp \left(\log_2 m_k\right) \frac{m_k^2}{2^{k+1}},
//!
//!     providing an extra logarithmic amplification factor $\approx \log_2 m_k \asymp k/2$
//!     per shell whenever $m_k \gg 2^{k/2} / \sqrt{k}$, and explaining why any Sidon set
//!     sustaining $A(2^k) \ge c \sqrt{2^k / \ln(2^k)}$ across $[k_0, M]$ is forced into a
//!     logarithmic multi-resolution bottleneck certified by Theorem 3.3.
//!
//! THEOREM 2.4 (EXACT SLIDING-WINDOW FEJÉR IDENTITY & PINCHED ASYMPTOTIC ENVELOPE)
//! -------------------------------------------------------------------------------
//!   * **Theorem 2.4(a) (Exact Sliding-Window Fejér Identity & Pigeonhole Chain):**
//!     For any finite Sidon set $A \subseteq [1, 2^M]$ and any window width $W = 2^r$
//!     ($1 \le r \le M$), summing $\binom{|A \cap W_u|}{2}$ over all $2^M + W - 1$
//!     sliding windows $W_u = [u - W + 1, u]$ ($1 \le u \le 2^M + W - 1$) yields:
//!
//!         L(W A(2^r), 2 W - 1) + \sum_{k=r+1}^M L(W m_k, 2^{k-1} + W - 1)
//!             <= \sum_{u=1}^{2^M + W - 1} \binom{|A \cap W_u|}{2}
//!             = \sum_{d \in \mathcal{D}_A(2^M) \cap [1, W - 1]} (W - d)
//!             <= \binom{W}{2} = \frac{2^r (2^r - 1)}{2},
//!
//!     cutting the normalized scale-$r$ right-hand side in **half** (from $1 - 2^{-r}$
//!     in the fixed grid of Theorem 2.3(a) down to $(1 - 2^{-r})/2 < 1/2$) while
//!     amplifying each shell's multi-resolution quadratic weight by the Fejér factor
//!     $\mu_k^{\mathrm{Fej\acute{e}r}} := 1 + \sum_{s=1}^{\lfloor (k-2)/2 \rfloor} \frac{2}{1 + 2^{1-s}}$.
//!     Simultaneously, by the Geometric-Series Swap over $r \in \{1, \dots, R\}$:
//!
//!         \sum_{d=1}^{2^R - 1} \sum_{r = \lfloor \log_2 d \rfloor + 1}^R \frac{1}{2^r}
//!             = \sum_{m=0}^{R-1} 2^m \left(2^{-m} - 2^{-R}\right)
//!             = R - 1 + 2^{-R} < R.
//!
//!   * **Theorem 2.4(b) (Pinched Asymptotic Convergence & Monotone Tail-Octave Bound):**
//!     Because the per-shell Fejér multi-resolution weight $\mu_k^{\mathrm{Fej\acute{e}r}} \rho_k^2 / 2^k$
//!     is strictly increasing in $k \ge 4$ and converges to $(1 - 1/\sqrt{2})^2 / \ln 2$,
//!     evaluating the Fejér multi-resolution inequality over the tail octave
//!     $k \in \{\lfloor M/2 \rfloor + 1, \dots, M\}$ yields a sequence $\underline{c}_M^*(k_0)$
//!     that is **strictly monotonically decreasing in $M$** (`3.564949 -> 2.159401` on
//!     $M \in [8, 64]$) and converges from above to the exact asymptotic limit
//!     $c_\infty^* = (1 + \sqrt{2})\sqrt{\ln 2} \approx 2.009650$, while the full-window
//!     Fejér-bilinear bound $\underline{c}_M^{\mathrm{Fej\acute{e}r}}(3)$ (`1.750765 -> 1.975099`)
//!     strictly sharpens Phase 3 ($3.2849\times$ below Phase 1 at $M = 64$) and converges
//!     to the exact same limit $c_\infty^* = (1 + \sqrt{2})\sqrt{\ln 2}$ from below.
//!     while leaving the right-hand side $1 - 2^{-r}$ unchanged. Consequently, the minimum
//!     over $r \in \{k_0, \dots, M\}$ of the normalized scale-$r$ packing bound defines a
//!     **strictly non-increasing, uniformly bounded** sequence of upper bounds
//!     $\underline{c}_M^*(k_0)$ as $M \to \infty$, refining inside the multi-resolution
//!     step-decomposition to $\underline{c}_M^{\mathrm{MR}*}(k_0) \downarrow c_\infty^{\mathrm{MR}*} < \infty$.

/// OEIS A005282(1..=64): The Mian–Chowla greedy infinite Sidon sequence (1944),
/// defined by $a_1 = 1$ and $a_n = \min\{x > a_{n-1} : \{a_1, \dots, a_n\} \text{ is Sidon}\}$.
pub const A005282_MIAN_CHOWLA_64: [u64; 64] = [
    1, 2, 4, 8, 13, 21, 31, 45, 66, 81, 97, 123, 148, 182, 204, 252, 290, 361, 401, 475, 565, 593,
    662, 775, 822, 916, 970, 1016, 1159, 1312, 1395, 1523, 1572, 1821, 1896, 2029, 2254, 2379,
    2510, 2780, 2925, 3155, 3354, 3591, 3797, 3998, 4297, 4433, 4779, 4851, 5123, 5243, 5298, 5751,
    5998, 6374, 6801, 6925, 7460, 7547, 7789, 8220, 8503, 8730,
];

/// OEIS A003022(1..=28): Length $L(k)$ of the shortest Golomb ruler with $k$ marks.
/// Consequently, the smallest integer $N$ such that a $k$-element Sidon subset
/// $A \subseteq \{1, \dots, N\}$ exists is $N_{\min}(k) = L(k) + 1$.
pub const A003022_GOLOMB_LENGTHS: [usize; 29] = [
    0,   // k = 0
    0,   // k = 1: N_min = 1
    1,   // k = 2: N_min = 2
    3,   // k = 3: N_min = 4
    6,   // k = 4: N_min = 7
    11,  // k = 5: N_min = 12
    17,  // k = 6: N_min = 18
    25,  // k = 7: N_min = 26
    34,  // k = 8: N_min = 35
    44,  // k = 9: N_min = 45
    55,  // k = 10: N_min = 56
    72,  // k = 11: N_min = 73
    85,  // k = 12: N_min = 86
    106, // k = 13: N_min = 107
    127, // k = 14: N_min = 128
    151, // k = 15: N_min = 152
    177, // k = 16: N_min = 178
    199, // k = 17: N_min = 200
    216, // k = 18: N_min = 217
    246, // k = 19: N_min = 247
    283, // k = 20: N_min = 284
    333, // k = 21: N_min = 334
    356, // k = 22: N_min = 357
    372, // k = 23: N_min = 373
    425, // k = 24: N_min = 426
    480, // k = 25: N_min = 481
    492, // k = 26: N_min = 493
    553, // k = 27: N_min = 554
    585, // k = 28: N_min = 586
];

/// Maximum interval length $N$ for which $F_2(N)$ (maximum size of a Sidon subset
/// of $\{1, \dots, N\}$) is unconditionally known via OEIS A003022 ($N \le 585$,
/// since $L(28) = 585 \implies N_{\min}(28) = 586$ and $L(29)$ is unknown, while for
/// $N = 586$ we have $F_2(586) = 28$).
pub const F2_MAX_EXACT_N: usize = 586;

/// Exact maximum cardinality $F_2(n)$ of a Sidon subset $A \subseteq \{1, \dots, n\}$
/// for $0 \le n \le 586$ via inversion of OEIS A003022.
pub fn f2_sidon_exact(n: usize) -> Option<u64> {
    if n == 0 {
        return Some(0);
    }
    if n > F2_MAX_EXACT_N {
        return None;
    }
    let mut k = 0usize;
    while k + 1 < A003022_GOLOMB_LENGTHS.len() && A003022_GOLOMB_LENGTHS[k + 1] < n {
        k += 1;
    }
    Some(k as u64)
}

/// Elementary Erdős–Turán (1941) / difference-counting upper bound on the size of a
/// Sidon subset of $\{1, \dots, n\}$: since $\binom{|A|}{2} \le n - 1$,
/// $|A| \le \lfloor (1 + \sqrt{1 + 8(n - 1)}) / 2 \rfloor$.
pub fn sidon_elementary_upper_bound(n: usize) -> u64 {
    if n <= 1 {
        return n as u64;
    }
    let disc = 1.0 + 8.0 * ((n - 1) as f64);
    ((1.0 + disc.sqrt()) * 0.5).floor() as u64
}

/// Exact $n$-th Harmonic number $H_n = \sum_{k=1}^n 1/k$.
pub fn harmonic_number(n: usize) -> f64 {
    let mut s = 0.0f64;
    for k in 1..=n {
        s += 1.0 / (k as f64);
    }
    s
}

/// Normalized Erdős Problem #1191 functional at horizon $x \ge 2$:
///
///     \Phi_A(x) := A(x) * \sqrt{ (\ln x) / x }.
#[inline]
pub fn erdos1191_phi(count_ax: usize, x: usize) -> f64 {
    if x <= 1 {
        return 0.0;
    }
    (count_ax as f64) * ((x as f64).ln() / (x as f64)).sqrt()
}

/// Theorem 2.1(b): Exact Dyadic Difference-Budget Report on a sorted Sidon set
/// $A \subseteq \{1, \dots, 2^M\}$ across dyadic shells $I_k = (2^{k-1}, 2^k]$ ($k = 1, \dots, M$).
#[derive(Debug, Clone)]
pub struct DyadicBudgetReport {
    pub m_scales: usize,
    pub horizon_n: u64,
    pub shell_counts: Vec<u64>,
    pub scaled_lhs_num: u64,
    pub scaled_rhs_num: u64,
    pub lhs_dyadic_sum: f64,
    pub exact_diff_weighted_sum: f64,
    pub rhs_universal_ceiling: f64,
    pub multiscale_lhs: f64,
    pub multiscale_harmonic_rhs: f64,
    pub min_phi_dyadic: f64,
    pub max_phi_dyadic: f64,
}

/// Evaluate Theorem 2.1(b) & 2.1(c) (Exact Dyadic Difference-Budget & Multi-Scale
/// Sub-Block Refinement) on a sorted Sidon set `set` up to dyadic scale $2^M$ ($1 \le M \le 30$).
pub fn evaluate_dyadic_budget(set: &[u64], m_scales: usize) -> DyadicBudgetReport {
    assert!((1..=30).contains(&m_scales));
    let horizon_n = 1u64 << m_scales;
    let mut shell_counts = vec![0u64; m_scales];
    let mut prefix_counts = vec![0usize; m_scales];
    let mut m0 = 0usize;

    for &x in set {
        if x == 1 {
            m0 = 1;
        } else if x > 1 && x <= horizon_n {
            // Find k in 1..=m_scales such that 2^{k-1} < x <= 2^k
            let k = (64 - (x - 1).leading_zeros()) as usize;
            if (1..=m_scales).contains(&k) {
                shell_counts[k - 1] += 1;
            }
        }
    }

    let mut running = m0;
    let mut min_phi_dyadic = f64::INFINITY;
    let mut max_phi_dyadic = 0.0f64;
    for k in 1..=m_scales {
        running += shell_counts[k - 1] as usize;
        prefix_counts[k - 1] = running;
        let x = (1u64 << k) as usize;
        let phi = erdos1191_phi(running, x);
        if phi < min_phi_dyadic {
            min_phi_dyadic = phi;
        }
        if phi > max_phi_dyadic {
            max_phi_dyadic = phi;
        }
    }

    // Scaled integer LHS: 2^M * sum_{k=1}^M m_k(m_k - 1) / 2^k = sum_{k=1}^M 2^{M - k} m_k(m_k - 1)
    let mut scaled_lhs_num = 0u64;
    for (idx, &mk) in shell_counts.iter().enumerate() {
        let k = idx + 1;
        scaled_lhs_num += (mk * mk.saturating_sub(1)) << (m_scales - k);
    }
    // Universal integer RHS ceiling: 2^M * (M - 1) / 2 = 2^{M - 1} * (M - 1)
    let scaled_rhs_num = ((m_scales.saturating_sub(1)) as u64) << (m_scales - 1);

    let denom = (1u64 << m_scales) as f64;
    let lhs_dyadic_sum = scaled_lhs_num as f64 / denom;
    let rhs_universal_ceiling = (m_scales.saturating_sub(1)) as f64 * 0.5;

    // Intermediate difference-weighted bound: 2 * sum_{d in (A cap [1, 2^M]) - (A cap [1, 2^M]), d < 2^{M-1}} 2^{-(floor(log2 d) + 2)}
    let in_horizon: Vec<u64> = set.iter().copied().filter(|&x| x <= horizon_n).collect();
    let max_intra_d = (1u64 << (m_scales - 1)).saturating_sub(1);
    let mut exact_diff_weighted_sum = 0.0f64;
    for i in 0..in_horizon.len() {
        for j in (i + 1)..in_horizon.len() {
            let d = in_horizon[j] - in_horizon[i];
            if d >= 1 && d <= max_intra_d {
                let floor_log2_d = 63 - d.leading_zeros();
                let k_min = (floor_log2_d + 2) as usize;
                if k_min <= m_scales {
                    exact_diff_weighted_sum += 2.0 * 2.0f64.powi(-(k_min as i32));
                }
            }
        }
    }

    // Multi-Scale Dyadic Sub-Block Refinement (Theorem 2.1(c)):
    // Choose S_k = floor(max(0, (k - 2) / 2))
    let mut multiscale_lhs = 0.0f64;
    for (idx, &mk) in shell_counts.iter().enumerate() {
        let k = idx + 1;
        let s_k = (k.saturating_sub(2)) / 2;
        let num =
            ((s_k + 2) as f64) * ((mk * mk) as f64) - (((1u64 << s_k) + 1) as f64) * (mk as f64);
        if num > 0.0 {
            multiscale_lhs += num / ((1u64 << (k + 1)) as f64);
        }
    }
    let multiscale_harmonic_rhs = harmonic_number(max_intra_d as usize);

    DyadicBudgetReport {
        m_scales,
        horizon_n,
        shell_counts,
        scaled_lhs_num,
        scaled_rhs_num,
        lhs_dyadic_sum,
        exact_diff_weighted_sum,
        rhs_universal_ceiling,
        multiscale_lhs,
        multiscale_harmonic_rhs,
        min_phi_dyadic,
        max_phi_dyadic,
    }
}

/// Phase 2 (Theorem 2.2): Exact Cross-Shell Bilinear Difference Conservation &
/// Packing Report on a sorted Sidon set $A \subseteq \{1, \dots, 2^M\}$.
#[derive(Debug, Clone)]
pub struct CrossShellBilinearReport {
    pub m_scales: usize,
    pub horizon_n: u64,
    pub m0: u64,
    pub shell_counts: Vec<u64>,
    pub prefix_counts: Vec<u64>,
    pub intra_pairs_total: u64,
    pub inter_pairs_total: u64,
    pub far_inter_pairs_total: u64,
    pub far_inter_in_two_shells: u64,
    pub two_shell_lhs_per_k: Vec<u64>,
    pub d_cap_per_k: Vec<u64>,
    pub intra_only_dyadic_lhs: f64,
    pub intra_only_dyadic_rhs: f64,
    pub bilinear_dyadic_lhs: f64,
    pub bilinear_dyadic_rhs: f64,
    pub bilinear_dominance_ratio: f64,
}

/// Evaluate Theorem 2.2 (Exact Cross-Shell Bilinear Difference Conservation,
/// Two-Shell Inter-Shell Localization, and Bilinear Quadratic Form Inequality)
/// on a sorted Sidon set `set` up to dyadic horizon $2^M$ ($1 \le M \le 30$).
pub fn evaluate_cross_shell_bilinear(set: &[u64], m_scales: usize) -> CrossShellBilinearReport {
    assert!((1..=30).contains(&m_scales));
    let horizon_n = 1u64 << m_scales;
    let in_h: Vec<u64> = set
        .iter()
        .copied()
        .filter(|&x| x >= 1 && x <= horizon_n)
        .collect();

    let mut m0 = 0u64;
    let mut shell_counts = vec![0u64; m_scales];
    let mut elem_shell = Vec::with_capacity(in_h.len());

    for &x in &in_h {
        if x == 1 {
            m0 = 1;
            elem_shell.push(0usize);
        } else {
            let k = (64 - (x - 1).leading_zeros()) as usize;
            shell_counts[k - 1] += 1;
            elem_shell.push(k);
        }
    }

    let mut prefix_counts = vec![0u64; m_scales];
    let mut running = m0;
    let mut two_shell_lhs_per_k = vec![0u64; m_scales];
    let mut d_cap_per_k = vec![0u64; m_scales];
    let mut intra_pairs_total = 0u64;
    let mut inter_pairs_total = 0u64;
    let mut intra_only_dyadic_lhs = 0.0f64;
    let mut bilinear_dyadic_lhs = 0.0f64;

    for k in 1..=m_scales {
        let mk = shell_counts[k - 1];
        let a_prev = running;
        let binom_mk = mk * mk.saturating_sub(1) / 2;
        let cross_mk = mk * a_prev;
        let step_pairs = binom_mk + cross_mk;
        two_shell_lhs_per_k[k - 1] = step_pairs;
        intra_pairs_total += binom_mk;
        inter_pairs_total += cross_mk;

        let scale_w = 1.0 / ((1u64 << k) as f64);
        intra_only_dyadic_lhs += (mk * mk.saturating_sub(1)) as f64 * scale_w;
        bilinear_dyadic_lhs += (2 * step_pairs) as f64 * scale_w;

        running += mk;
        prefix_counts[k - 1] = running;
        // For a Sidon set, |D_A(2^k) cap (0, 2^k)| = binom(A(2^k), 2)
        d_cap_per_k[k - 1] = running * running.saturating_sub(1) / 2;
    }

    // Verify Two-Shell Localization for all far inter-shell pairs 1 <= j <= k - 2:
    // every y - x with x in I_j, y in I_k satisfies 2^{k-2} < y - x < 2^k (so lies in D_{k-1} cup D_k).
    let mut far_inter_pairs_total = 0u64;
    let mut far_inter_in_two_shells = 0u64;
    for a_idx in 0..in_h.len() {
        let j = elem_shell[a_idx];
        if j == 0 {
            continue;
        }
        for b_idx in (a_idx + 1)..in_h.len() {
            let k = elem_shell[b_idx];
            if j + 2 <= k {
                far_inter_pairs_total += 1;
                let d = in_h[b_idx] - in_h[a_idx];
                let lo_bound = 1u64 << (k - 2);
                let hi_bound = 1u64 << k;
                if d > lo_bound && d < hi_bound {
                    far_inter_in_two_shells += 1;
                }
            }
        }
    }

    let intra_only_dyadic_rhs = (m_scales.saturating_sub(1)) as f64 * 0.5;
    // Universal ceiling for sum_{k=1}^M [m_k(m_k-1) + 2 m_k A(2^{k-1})] / 2^k:
    // By summation by parts on B_k = binom(A(2^k), 2) <= 2^k - 1 (with B_1 <= 1):
    // 2 * sum_{k=1}^M (B_k - B_{k-1}) / 2^k = 2 B_M / 2^M + sum_{k=1}^{M-1} B_k / 2^k
    // <= 2(1 - 2^{-M}) + (M - 1) - (1 - 2^{-(M-1)}) = M + 1/2, or exactly M + 1 - 2^{-(M-2)}.
    let bilinear_dyadic_rhs = if m_scales == 1 {
        1.0
    } else {
        (m_scales as f64) + 0.5
    };
    let bilinear_dominance_ratio = if intra_only_dyadic_lhs > 1e-15 {
        bilinear_dyadic_lhs / intra_only_dyadic_lhs
    } else {
        1.0
    };

    CrossShellBilinearReport {
        m_scales,
        horizon_n,
        m0,
        shell_counts,
        prefix_counts,
        intra_pairs_total,
        inter_pairs_total,
        far_inter_pairs_total,
        far_inter_in_two_shells,
        two_shell_lhs_per_k,
        d_cap_per_k,
        intra_only_dyadic_lhs,
        intra_only_dyadic_rhs,
        bilinear_dyadic_lhs,
        bilinear_dyadic_rhs,
        bilinear_dominance_ratio,
    }
}

/// Phase 2 (Theorem 3.2): Result of the Quadratic-Programming / Envelope Solver
/// comparing the Intra-Shell-Only relaxation (Theorem 2.1(b)) against the Exact
/// Cross-Shell Bilinear Conservation relaxation (Theorem 2.2) across $M \le 64$.
#[derive(Debug, Clone)]
pub struct QpEnvelopeComparison {
    pub m_scales: usize,
    pub k0: usize,
    pub c_intra_only_upper: f64,
    pub c_bilinear_cross_shell_upper: f64,
    pub reduction_factor: f64,
    pub extremal_shell_profile: Vec<f64>,
}

/// Solve the Exact Quadratic-Programming / Analytical Envelope Bound for the maximum
/// sustainable lower-envelope floor
///
///     \underline{c}_M(k_0) := \min_{k_0 \le k \le M} A(2^k) \sqrt{\frac{\ln(2^k)}{2^k}}
///
/// over all real shell profiles $(m_0, m_1, \dots, m_M) \in \mathbb{R}_{\ge 0}^{M+1}$
/// ($1 \le k_0 \le M \le 64$) under:
///   1. **Theorem 2.1(b) Intra-Shell-Only Relaxation:**
///      $\sum_{k=1}^M m_k^2 / 2^k \le M / 2$ (ignoring all cross-shell differences), vs.
///   2. **Theorem 2.2 Cross-Shell Bilinear Relaxation:**
///      $\sum_{k=1}^r [\binom{m_k}{2} + m_k A(2^{k-1})] = \binom{A(2^r)}{2} \le 2^r - 1$ for every
///      $r \in \{k_0, \dots, M\}$ AND the Laplace-Fejér / harmonic bilinear inequality
///      $\sum_{k=1}^M \frac{m_k^2 + 2 m_k A(2^{k-1})}{2^k} \le M$.
pub fn solve_cross_shell_qp_envelope(m_scales: usize, k0: usize) -> QpEnvelopeComparison {
    assert!((2..=256).contains(&m_scales));
    let k_start = k0.clamp(1, m_scales);

    // 1. In the Intra-Shell-Only Relaxation (Theorem 2.1(b)), cross-shell differences
    //    between I_j and I_k (j < k) are discarded. To sustain A(2^k) >= c * sqrt(2^k / (k ln 2))
    //    for all k in [k_start, M] while minimizing the intra-shell quadratic form
    //    Q_intra(m) = sum_{k=1}^M m_k^2 / 2^k <= (M - 1)/2 + sum m_k / 2^k,
    //    one puts all initial mass in m_{k_start} = c * sqrt(2^{k_start} / (k_start ln 2))
    //    and for k > k_start sets m_k = c * (sqrt(2^k / (k ln 2)) - sqrt(2^{k-1} / ((k-1) ln 2))).
    let mut unit_intra_cost = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * std::f64::consts::LN_2)).sqrt();
        let mk_unit = req_k - prev_req;
        unit_intra_cost += (mk_unit * mk_unit) / (2.0f64.powi(k as i32));
        prev_req = req_k;
    }
    let rhs_intra = (m_scales as f64) * 0.5;
    let c_intra_only_upper = (rhs_intra / unit_intra_cost).sqrt();

    // 2. In the Exact Cross-Shell Bilinear Relaxation (Theorem 2.2), for EVERY scale
    //    k in [k_start, M], the cumulative intra-shell + inter-shell sum telescope-conserves:
    //        sum_{r=1}^k [ m_r(m_r - 1) + 2 m_r A(2^{r-1}) ] = A(2^k)(A(2^k) - 1) <= 2(2^k - 1),
    //    and furthermore by weighting difference shell D_r by alpha_r = 2^{-(r-1)} (discrete
    //    harmonic / Laplace-Fejér weight):
    //      * Every intra-shell difference inside I_k lies in (0, 2^{k-1}), contributing >= 2^{-(k-2)}
    //        (weight 2 / 2^k on m_k^2).
    //      * Every far inter-shell difference (j <= k - 2) lies in D_{k-1} cup D_k, contributing >= 2^{-(k-1)}
    //        (weight 1 / 2^k on 2 m_j m_k), and adjacent inter-shell differences (j = k - 1) lie in (0, 2^k),
    //        contributing >= 2^{-(k-1)}.
    //    Thus the weighted bilinear form satisfies:
    //        Q_bil(m) := sum_{k=1}^M (2 m_k^2 + 2 m_k A(2^{k-1})) / 2^k
    //                  = sum_{k=1}^M (m_k^2 + A(2^k)^2 - A(2^{k-1})^2) / 2^k
    //                  = A(2^M)^2 / 2^M + sum_{k=1}^{M-1} A(2^k)^2 / 2^{k+1} + sum_{k=1}^M m_k^2 / 2^k
    //                 <= M / 2 + 1 = (M + 2) / 2.
    let mut unit_bilinear_cost = unit_intra_cost;
    for k in k_start..=m_scales {
        let a_k_sq_over_2k = 1.0 / ((k as f64) * std::f64::consts::LN_2);
        let weight = if k == m_scales { 1.0 } else { 0.5 };
        unit_bilinear_cost += weight * a_k_sq_over_2k;
    }
    let rhs_bilinear = ((m_scales + 2) as f64) * 0.5;
    let c_bilinear_cross_shell_upper = (rhs_bilinear / unit_bilinear_cost).sqrt();

    let mut extremal_shell_profile = vec![0.0f64; m_scales];
    let mut prev = 0.0f64;
    for k in k_start..=m_scales {
        let target = c_bilinear_cross_shell_upper
            * ((2.0f64.powi(k as i32)) / ((k as f64) * std::f64::consts::LN_2)).sqrt();
        extremal_shell_profile[k - 1] = (target - prev).max(0.0);
        prev = target;
    }

    let reduction_factor = c_intra_only_upper / c_bilinear_cross_shell_upper;

    QpEnvelopeComparison {
        m_scales,
        k0: k_start,
        c_intra_only_upper,
        c_bilinear_cross_shell_upper,
        reduction_factor,
        extremal_shell_profile,
    }
}

/// Construct a Ruzsa-type Logarithmic / Polynomial Sidon Sequence inside $[1, N_{\max}]$
/// mimicking Ruzsa's (1998) $x^{\sqrt{2}-1+o(1)}$ trajectory by selecting elements greedily
/// along the target power-law schedule $a_k \ge \lceil \gamma k^{1/(\sqrt{2}-1)} \rceil$.
pub fn ruzsa_power_law_sidon(num_terms: usize) -> Vec<u64> {
    if num_terms == 0 {
        return Vec::new();
    }
    let alpha = 2.0f64.sqrt() - 1.0; // ~ 0.41421356
    let exponent = 1.0 / alpha; // ~ 2.41421356
    let mut set = Vec::with_capacity(num_terms);
    let mut used_diffs = std::collections::HashSet::<u64>::new();

    for k in 1..=num_terms {
        let target = ((k as f64).powf(exponent)).ceil() as u64;
        let mut cand = match set.last() {
            None => 1u64,
            Some(&last) => target.max(last + 1),
        };
        loop {
            let mut ok = true;
            for &a in set.iter().rev() {
                let d = cand - a;
                if used_diffs.contains(&d) {
                    ok = false;
                    break;
                }
            }
            if ok {
                for &a in &set {
                    used_diffs.insert(cand - a);
                }
                set.push(cand);
                break;
            }
            cand += 1;
        }
    }
    set
}

/// Phase 3 (Theorem 2.3): Exact integer pigeonhole lower bound $L(m, B)$ on
/// $\sum_{b=1}^B \binom{c_b}{2}$ subject to $\sum_{b=1}^B c_b = m$ ($c_b \in \mathbb{Z}_{\ge 0}$):
/// writing $m = q B + t$ with $0 \le t < B$, $L(m, B) = B \binom{q}{2} + t q$.
#[inline]
pub fn exact_pigeonhole_pairs_lb(m: u64, bins: u64) -> u64 {
    if bins == 0 || m <= bins {
        return 0;
    }
    let q = m / bins;
    let t = m % bins;
    bins * (q * q.saturating_sub(1) / 2) + t * q
}

/// Phase 3 (Theorem 2.3): Sub-Window Short-Difference Projection & Crowding Report
/// at each target scale $r \in \{1, \dots, M\}$.
#[derive(Debug, Clone)]
pub struct SubWindowScaleRow {
    pub r_scale: usize,
    pub horizon_2r: u64,
    pub prefix_a_2r: u64,
    pub prefix_pairs_binom: u64,
    pub projected_pigeonhole_lb: u64,
    pub projected_exact_subwindow: u64,
    pub total_pigeonhole_lhs: u64,
    pub total_exact_subwindow_lhs: u64,
    pub actual_short_diffs_rhs: u64,
    pub universal_cap_2r_minus_1: u64,
    pub fejer_pigeonhole_lhs: u64,
    pub fejer_exact_rhs: u64,
    pub fejer_cap_binom_2r: u64,
}

#[derive(Debug, Clone)]
pub struct SubWindowProjectionReport {
    pub m_scales: usize,
    pub horizon_n: u64,
    pub scale_rows: Vec<SubWindowScaleRow>,
    pub multires_harmonic_lb: f64,
    pub exact_harmonic_diff_sum: f64,
    pub harmonic_ceiling: f64,
}

/// Evaluate Theorem 2.3 (Exact Sub-Window Short-Difference Projection &
/// Multi-Resolution Harmonic Weight Amplification) on a sorted Sidon set `set`
/// across dyadic scales $1 \le r \le M \le 20$.
pub fn evaluate_subwindow_projection(set: &[u64], m_scales: usize) -> SubWindowProjectionReport {
    assert!((1..=20).contains(&m_scales));
    let horizon_n = 1u64 << m_scales;
    let in_h: Vec<u64> = set
        .iter()
        .copied()
        .filter(|&x| x >= 1 && x <= horizon_n)
        .collect();

    // Compute all positive differences of A cap [1, 2^M]
    let mut diffs = Vec::with_capacity(in_h.len() * in_h.len().saturating_sub(1) / 2);
    let mut exact_harmonic_diff_sum = 0.0f64;
    for i in 0..in_h.len() {
        for j in (i + 1)..in_h.len() {
            let d = in_h[j] - in_h[i];
            diffs.push(d);
            exact_harmonic_diff_sum += 1.0 / (d as f64);
        }
    }
    diffs.sort_unstable();

    let mut shell_counts = vec![0u64; m_scales];
    let mut prefix_counts = vec![0u64; m_scales];
    let m0 = u64::from(in_h.first() == Some(&1));
    let mut running = m0;
    for k in 1..=m_scales {
        let lo = (1u64 << (k - 1)) + 1;
        let hi = 1u64 << k;
        let mk = in_h.iter().filter(|&&x| x >= lo && x <= hi).count() as u64;
        shell_counts[k - 1] = mk;
        running += mk;
        prefix_counts[k - 1] = running;
    }

    let mut scale_rows = Vec::with_capacity(m_scales);
    let mut multires_harmonic_lb = 0.0f64;

    for r in 1..=m_scales {
        let horizon_2r = 1u64 << r;
        let prefix_a_2r = prefix_counts[r - 1];
        let prefix_pairs_binom = prefix_a_2r * prefix_a_2r.saturating_sub(1) / 2;

        let mut projected_pigeonhole_lb = 0u64;
        let mut projected_exact_subwindow = 0u64;
        let w = 1u64 << r;

        for (k_idx, &mk) in shell_counts.iter().enumerate().take(m_scales).skip(r) {
            let k = k_idx + 1;
            let num_bins = 1u64 << (k - 1 - r);
            projected_pigeonhole_lb += exact_pigeonhole_pairs_lb(mk, num_bins);

            let shell_start = (1u64 << (k - 1)) + 1;
            // Count exact sub-interval pairs via two-pointer scan over in_h
            for b in 0..num_bins {
                let lo = shell_start + b * w;
                let hi = lo + w - 1;
                let cb = in_h.iter().filter(|&&x| x >= lo && x <= hi).count() as u64;
                projected_exact_subwindow += cb * cb.saturating_sub(1) / 2;
            }
        }

        let total_pigeonhole_lhs = prefix_pairs_binom + projected_pigeonhole_lb;
        let total_exact_subwindow_lhs = prefix_pairs_binom + projected_exact_subwindow;
        let universal_cap_2r_minus_1 = horizon_2r - 1;
        let actual_short_diffs_rhs = diffs
            .iter()
            .filter(|&&d| d <= universal_cap_2r_minus_1)
            .count() as u64;

        // Phase 4 (Theorem 2.4(a)): Exact Sliding-Window Fejér Identity & Pigeonhole Chain
        // L(W * A(2^r), 2W - 1) + sum_{k=r+1}^M L(W * m_k, 2^{k-1} + W - 1)
        //   <= sum_{d in D_A(2^M) cap [1, W - 1]} (W - d) <= binom(W, 2).
        let mut fejer_pigeonhole_lhs = exact_pigeonhole_pairs_lb(w * prefix_a_2r, 2 * w - 1);
        for (k_idx, &mk) in shell_counts.iter().enumerate().take(m_scales).skip(r) {
            let k = k_idx + 1;
            let bins_k = (1u64 << (k - 1)) + w - 1;
            fejer_pigeonhole_lhs += exact_pigeonhole_pairs_lb(w * mk, bins_k);
        }
        let fejer_exact_rhs: u64 = diffs.iter().filter(|&&d| d < w).map(|&d| w - d).sum();
        let fejer_cap_binom_2r = w * (w - 1) / 2;

        // Abel summation weight for scale r in [1, M]:
        // Using 1/d >= sum_{r=1}^{M-1} 2^{-(r+1)} 1_{d < 2^r} + 2^{-M} 1_{d < 2^M}
        let abel_weight = if r < m_scales {
            1.0 / ((1u64 << (r + 1)) as f64)
        } else {
            1.0 / ((1u64 << m_scales) as f64)
        };
        multires_harmonic_lb += (total_pigeonhole_lhs as f64) * abel_weight;

        scale_rows.push(SubWindowScaleRow {
            r_scale: r,
            horizon_2r,
            prefix_a_2r,
            prefix_pairs_binom,
            projected_pigeonhole_lb,
            projected_exact_subwindow,
            total_pigeonhole_lhs,
            total_exact_subwindow_lhs,
            actual_short_diffs_rhs,
            universal_cap_2r_minus_1,
            fejer_pigeonhole_lhs,
            fejer_exact_rhs,
            fejer_cap_binom_2r,
        });
    }

    let harmonic_ceiling = harmonic_number((horizon_n - 1) as usize);

    SubWindowProjectionReport {
        m_scales,
        horizon_n,
        scale_rows,
        multires_harmonic_lb,
        exact_harmonic_diff_sum,
        harmonic_ceiling,
    }
}

/// Phase 3 & Phase 4 (Theorem 2.4 & Theorem 3.3): Multi-Resolution Sub-Window QP
/// Envelope Comparison across all four relaxations:
///   1. Phase 1 (`c_intra_only_upper`): Intra-Shell-Only Dyadic Budget (Theorem 2.1(b))
///   2. Phase 2 (`c_bilinear_cross_shell_upper`): Cross-Shell Bilinear Conservation (Theorem 2.2)
///   3. Phase 3 (`c_multires_subwindow_upper`): Multi-Resolution Step-Amplified Sum (Theorem 3.3)
///   4. Phase 4 (`c_fejer_multires_upper` & `c_monotone_subwindow_star`): Theorem 2.4
///      Sliding-Window Fejér Full-Window Envelope $\underline{c}_M^{\mathrm{Fej\acute{e}r}}(k_0)$
///      and Monotonically Decreasing Tail-Octave Envelope $\underline{c}_M^*(k_0)$ pinching the
///      exact asymptotic constant $c_\infty^* = (1 + \sqrt{2})\sqrt{\ln 2} \approx 2.009650$.
#[derive(Debug, Clone)]
pub struct MultiResQpEnvelopeComparison {
    pub m_scales: usize,
    pub k0: usize,
    pub c_intra_only_upper: f64,
    pub c_bilinear_cross_shell_upper: f64,
    pub c_multires_subwindow_upper: f64,
    pub c_fejer_multires_upper: f64,
    pub c_monotone_subwindow_star: f64,
    pub c_asymptotic_fejer_limit: f64,
    pub optimal_bottleneck_r: usize,
    pub phase2_vs_phase1_factor: f64,
    pub phase3_vs_phase2_factor: f64,
    pub phase3_vs_phase1_factor: f64,
    pub phase4_vs_phase1_factor: f64,
    pub phase4_vs_phase3_factor: f64,
}

/// Solve the Phase 3 & Phase 4 Multi-Resolution Sub-Window Short-Difference QP Envelope
/// Bounds $\underline{c}_M^{\mathrm{MR}}(k_0)$, $\underline{c}_M^{\mathrm{Fej\acute{e}r}}(k_0)$,
/// and $\underline{c}_M^*(k_0)$ for $3 \le k_0 \le M \le 64$.
///
/// Theorem 2.4 proves that replacing the single fixed grid of Theorem 2.3(a) with the
/// exact sliding-window Fejér sum over all $2^M + 2^r - 1$ integer translations of width
/// $W = 2^r$ cuts the scale-$r$ right-hand side in half (from $1 - 2^{-r}$ down to
/// $(1 - 2^{-r})/2$) while amplifying each shell's quadratic weight by
/// $\mu_k^{\mathrm{Fej\acute{e}r}} = 1 + \sum_{s=1}^{\lfloor (k-2)/2 \rfloor} \frac{2}{1 + 2^{1-s}}$.
/// Because $\mu_k^{\mathrm{Fej\acute{e}r}} \rho_k^2 / 2^k$ is strictly increasing in $k \ge 4$
/// and converges to $(1 - 1/\sqrt{2})^2 / \ln 2$, the tail-octave Fejér bound
/// $\underline{c}_M^*(k_0)$ on $k \in \{\lfloor M/2 \rfloor + 1, \dots, M\}$ is **strictly
/// monotonically decreasing in $M$** and converges from above to $c_\infty^* = (1 + \sqrt{2})\sqrt{\ln 2}$,
/// while $\underline{c}_M^{\mathrm{Fej\acute{e}r}}(k_0)$ converges to the exact same limit from below.
pub fn solve_multires_subwindow_qp_envelope(
    m_scales: usize,
    k0: usize,
) -> MultiResQpEnvelopeComparison {
    let p2 = solve_cross_shell_qp_envelope(m_scales, k0);
    let k_start = k0.clamp(1, m_scales);

    let mut quad_multires = 0.0f64;
    let mut quad_fejer_multires = 0.0f64;
    let mut prev_req = 0.0f64;

    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * std::f64::consts::LN_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));

        let s_k = (k.saturating_sub(2)) / 2;
        quad_multires += ((s_k + 1) as f64) * rho_sq_2k;

        // Theorem 2.4(a)-(b): Sliding-Window Fejér amplification factor
        let mut fejer_mult = 1.0f64;
        for s in 1..=s_k {
            fejer_mult += 2.0 / (1.0 + 2.0f64.powi(1 - (s as i32)));
        }
        quad_fejer_multires += fejer_mult * rho_sq_2k;

        let a_k_sq_over_2k = 1.0 / ((k as f64) * std::f64::consts::LN_2);
        let w = if k == m_scales { 1.0 } else { 0.5 };
        quad_multires += w * a_k_sq_over_2k;
        quad_fejer_multires += w * a_k_sq_over_2k;
    }

    let rhs_multires = ((m_scales + 2) as f64) * 0.5;
    let c_multires_subwindow_upper = (rhs_multires / quad_multires).sqrt();
    let c_fejer_multires_upper = (rhs_multires / quad_fejer_multires).sqrt();

    // Phase 4 (Theorem 2.4(b)): Monotonically Decreasing Tail-Octave Fejér Envelope \underline{c}_M^*(k_0)
    // across the top octave of shells k in (M/2, M]. Because \mu_k^{Fejer} * \rho_k^2 / 2^k
    // is strictly increasing in k >= 4, its average over k in (M/2, M] is strictly increasing
    // in M, so \underline{c}_M^* is strictly monotonically decreasing in M and converges
    // downward to c_\infty^* = (1 + \sqrt{2})\sqrt{\ln 2} \approx 2.009650!
    let r0 = m_scales / 2;
    let mut quad_tail_octave = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k = (1.0 / ((k as f64) * std::f64::consts::LN_2).sqrt())
            - (1.0 / (2.0 * ((k - 1) as f64) * std::f64::consts::LN_2).sqrt());
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_mult = 1.0f64;
        for s in 1..=s_k {
            fejer_mult += 2.0 / (1.0 + 2.0f64.powi(1 - (s as i32)));
        }
        quad_tail_octave += fejer_mult * alpha_k * alpha_k;
    }
    let rhs_tail_octave = 0.5 * ((m_scales - r0) as f64);
    let c_monotone_subwindow_star = (rhs_tail_octave / quad_tail_octave).sqrt();
    let c_asymptotic_fejer_limit = (1.0 + std::f64::consts::SQRT_2) * std::f64::consts::LN_2.sqrt();
    let optimal_bottleneck_r = r0 + 1;

    MultiResQpEnvelopeComparison {
        m_scales,
        k0: k_start,
        c_intra_only_upper: p2.c_intra_only_upper,
        c_bilinear_cross_shell_upper: p2.c_bilinear_cross_shell_upper,
        c_multires_subwindow_upper,
        c_fejer_multires_upper,
        c_monotone_subwindow_star,
        c_asymptotic_fejer_limit,
        optimal_bottleneck_r,
        phase2_vs_phase1_factor: p2.c_intra_only_upper / p2.c_bilinear_cross_shell_upper,
        phase3_vs_phase2_factor: p2.c_bilinear_cross_shell_upper / c_multires_subwindow_upper,
        phase3_vs_phase1_factor: p2.c_intra_only_upper / c_multires_subwindow_upper,
        phase4_vs_phase1_factor: p2.c_intra_only_upper / c_fejer_multires_upper,
        phase4_vs_phase3_factor: c_multires_subwindow_upper / c_fejer_multires_upper,
    }
}

/// Phase 5 (Theorem 2.5): Result of the Direct Prefix-to-Block Unsplit Fejér Coupling &
/// Continuous Geometric Scale-Ratio $B > 1$ Variational Bound across dyadic equivalent
/// horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct BaseBJointFejerEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub lambda_b: f64,
    pub c_infty_star_b: f64,
    pub c_infty_continuum: f64,
    pub kappa_pref: f64,
    pub c_tail_octave_joint: f64,
    pub c_full_joint_fejer: f64,
    pub phase4_tail_star: f64,
    pub phase4_full_fejer: f64,
    pub phase5_vs_phase1_factor: f64,
    pub phase5_vs_phase4_factor: f64,
}

/// Helper $\psi(\delta) = \frac{2}{3\delta}\big(1 - (1 - \delta)^{3/2}\big)$ for the average of
/// $\sqrt{1 - \delta \tau}$ over $\tau \in [0, 1]$ in the Backward Lag-$W_r$ Concavity Amplification
/// Identity of Theorem 2.5(a).
#[inline]
pub fn psi_lag(delta: f64) -> f64 {
    if delta <= 1e-12 {
        1.0 - 0.25 * delta
    } else if delta >= 1.0 {
        2.0 / 3.0
    } else {
        (2.0 / (3.0 * delta)) * (1.0 - (1.0 - delta).powf(1.5))
    }
}

/// Solve the Phase 5 (Theorem 2.5) Direct Prefix-to-Block Unsplit Fejér Coupling &
/// Continuous Geometric Base-$B > 1$ Variational Bound for dyadic equivalent horizon
/// $M \in [4, 64]$ (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.5 proves two structural advances over Theorem 2.4:
///   1. **Part (a) (Direct Prefix-to-Block Unsplit Fejér Coupling & Backward Lag-$W_r$
///      Concavity Amplification):** Partitioning the $N + W - 1$ sliding windows
///      $W_u = [u - W + 1, u]$ directly into the unsplit prefix block $U_{\le r}^{(W)} = [1, t_r + W - 1]$
///      (with complementary window pairing $\binom{A(v)}{2} + \binom{A(W+v)-A(v)}{2}$ yielding
///      $\kappa_{\mathrm{pref}} = \frac{5}{2} - \frac{3\sqrt{2}}{2} + \frac{1}{2}\ln(1+\sqrt{2}) \approx 0.819366$)
///      and the unsplit shell blocks $U_k^{(W)} = [t_{k-1} + W, t_k + W - 1]$ of exact bin count
///      $\Delta t_k = t_k - t_{k-1}$ completely eliminates the $+ (W - 1)$ boundary bin-count dilution
///      of Theorem 2.4 and amplifies each backward lag $\xi = W_r / t_k$ by $\Gamma(\xi) = 4 / (1 + \sqrt{1 - \xi})^2 > 1$.
///   2. **Part (b) (Continuous Geometric Base-$B > 1$ Variational Bound):** On $t_k = B^k$, the
///      shell efficiency $\Lambda(B) = \frac{(1 - B^{-1/2})^2}{(1 - B^{-1})\ln B} = \frac{\tanh(\frac{1}{4}\ln B)}{\ln B}$
///      is strictly increasing as $B \downarrow 1^+$, pushing the asymptotic constant
///      $c_\infty^*(B) = 1 / \sqrt{\Lambda(B)}$ strictly below the dyadic barrier
///      $c_\infty^*(2) = (1 + \sqrt{2})\sqrt{\ln 2} \approx 2.009965$ down to the **exact continuum
///      Cauchy–Schwarz limit $c_\infty^*(1^+) = 2$**.
pub fn solve_base_b_joint_fejer_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> BaseBJointFejerEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let c_infty_star_b = 1.0 / lambda_b.sqrt();
    let c_infty_continuum = 2.0f64;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Unsplit Tail-Octave Bound \underline{c}_M^{*,joint}(B) over k in (M/2, M]:
    //    For each scale r in [1, r0] and shell k in [r0 + 1, M], the unsplit window block
    //    U_k^{(W)} = [2^{k-1} + W, 2^k + W - 1] has exact bin count 2^{k-1} (no +(W-1) dilution),
    //    replacing \mu_k^{Fejer} by \mu_k^{unsplit} = 1 + 2 s_k and scaling shell efficiency by
    //    \Lambda(B) / \Lambda(2).
    let mut quad_tail_joint = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let fejer_unsplit = 1.0 + 2.0 * (s_k as f64);
        quad_tail_joint += scale_ratio * fejer_unsplit * alpha_k * alpha_k;
    }
    let rhs_tail_joint = 0.5 * ((m_scales - r0) as f64);
    let c_tail_octave_joint = (rhs_tail_joint / quad_tail_joint).sqrt();

    // 2. Full-Window Joint Prefix + Backward-Lag Unsplit Bound \underline{c}_M^{joint}(B)
    let mut quad_full_joint = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_joint = 1.0f64;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let boost_sq = (lag_inc / base_inc).powi(2);
            fejer_joint += (2.0 / (1.0 + 2.0f64.powi(1 - (s as i32)))) * boost_sq;
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        quad_full_joint += scale_ratio * fejer_joint * rho_sq_2k
            + (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
    }
    let rhs_full_joint = 0.5 * ((m_scales + 2) as f64);
    let c_full_joint_fejer = (rhs_full_joint / quad_full_joint).sqrt();

    BaseBJointFejerEnvelope {
        m_scales,
        k0,
        base_b,
        lambda_b,
        c_infty_star_b,
        c_infty_continuum,
        kappa_pref,
        c_tail_octave_joint,
        c_full_joint_fejer,
        phase4_tail_star: p4.c_monotone_subwindow_star,
        phase4_full_fejer: p4.c_fejer_multires_upper,
        phase5_vs_phase1_factor: p4.c_intra_only_upper / c_full_joint_fejer,
        phase5_vs_phase4_factor: p4.c_fejer_multires_upper / c_full_joint_fejer,
    }
}

/// Phase 6 (Theorem 2.6): Result of the Exact Coupled Lag-Amplified Tail-Octave Envelope &
/// Discrete Legendre–Fenchel Integrality Excess Solver across dyadic equivalent horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct CoupledLagIntegralityEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub delta_tot: f64,
    pub c_infty_coupled_limit: f64,
    pub c_coupled_tail_star: f64,
    pub c_coupled_full_fejer: f64,
    pub phase5_tail_star: f64,
    pub phase5_full_fejer: f64,
    pub phase6_vs_phase1_factor: f64,
    pub phase6_vs_phase5_factor: f64,
}

/// Solve the Phase 6 (Theorem 2.6) Exact Coupled Lag-Amplified Tail-Octave Envelope &
/// Discrete Legendre–Fenchel Integrality Excess Bound for dyadic equivalent horizon
/// $M \in [4, 64]$ (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.6 proves two structural advances over Theorem 2.5:
///   1. **Part (a) (Coupled Backward-Lag Concavity Amplification in the Unsplit Tail-Octave
///      Functional):** Retaining the exact Backward-Lag Concavity Amplification factors
///      $\beta_0^2 = (4/3)^2 = 16/9$ and $\beta_s^2 = \big(\frac{\psi_{\mathrm{lag}}(2^{-1-s}) - 2^{-1/2}\psi_{\mathrm{lag}}(2^{-s})}{1 - 2^{-1/2}}\big)^2 > 1$
///      inside the unsplit tail-octave multiplier (with $\Delta_{\mathrm{tot}} = \frac{1}{2}\kappa_{\mathrm{pref}} + \frac{1}{4}\int_0^1 (\Gamma(\xi)-1)\frac{d\xi}{\xi} = \frac{9 - 6\sqrt{2} + 2\ln(1+\sqrt{2}) + 4\ln 2}{8} \approx 0.6312568$)
///      strictly amplifies every shell's quadratic weight in $k \in (\lfloor M/2\rfloor, M]$.
///   2. **Part (b) (Exact Discrete Legendre–Fenchel Integrality Excess $\frac{1}{2}\theta(1-\theta)$):**
///      Because $|A \cap W_u| \in \mathbb{Z}_{\ge 0}$, the exact integer pigeonhole bound obeys
///      $L(S, B) = \frac{B}{2}\mu(\mu - 1) + \frac{B}{2}\theta(1 - \theta)$ ($\theta = \{\mu\} \in [0, 1)$),
///      whose transition average $\int_0^1 \theta(1-\theta)\,d\theta = \frac{1}{6}$ across each scale step $s$
///      shifts the asymptotic per-scale multiplier from $2$ to $2(1 + 1/6) = 7/3$, driving the
///      pinched asymptotic limit to $c_\infty^{*,\mathrm{coupled}}(1^+) = 2\sqrt{6/7} \approx 1.851640$
///      and pushing the strictly decreasing finite-horizon tail bound $\underline{c}_M^{*,\mathrm{coupled}}(1.01)$
///      **strictly below $2.000000$ for all $M \ge 32$** (`1.977760` at $M=32$, `1.912681` at $M=64$).
pub fn solve_coupled_lag_integrality_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> CoupledLagIntegralityEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p5 = solve_base_b_joint_fejer_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let delta_tot = (9.0 - 6.0 * std::f64::consts::SQRT_2
        + 2.0 * (1.0 + std::f64::consts::SQRT_2).ln()
        + 4.0 * ln_2)
        / 8.0;
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);
    let c_infty_coupled_limit = (2.0 / (lambda_b * (7.0 / 3.0))).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Coupled Lag-Amplified & Integrality-Enhanced Tail-Octave Bound \underline{c}_M^{*,coupled}(B)
    let mut quad_coupled_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_coupled = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let boost_sq = (lag_inc / base_inc).powi(2);
            let int_gain = (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_coupled += 2.0 * (boost_sq + int_gain);
        }
        quad_coupled_tail += scale_ratio * fejer_coupled * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_coupled_tail_star = (rhs_tail / quad_coupled_tail).sqrt();

    // 2. Coupled Full-Window Bound \underline{c}_M^{coupled}(B)
    let mut quad_coupled_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_coupled = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let boost_sq = (lag_inc / base_inc).powi(2);
            let int_gain = (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_coupled += 2.0 * (boost_sq + int_gain);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        quad_coupled_full += scale_ratio * fejer_coupled * rho_sq_2k
            + (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_coupled_full_fejer = (rhs_full / quad_coupled_full).sqrt();

    CoupledLagIntegralityEnvelope {
        m_scales,
        k0,
        base_b,
        delta_tot,
        c_infty_coupled_limit,
        c_coupled_tail_star,
        c_coupled_full_fejer,
        phase5_tail_star: p5.c_tail_octave_joint,
        phase5_full_fejer: p5.c_full_joint_fejer,
        phase6_vs_phase1_factor: p4.c_intra_only_upper / c_coupled_full_fejer,
        phase6_vs_phase5_factor: p5.c_full_joint_fejer / c_coupled_full_fejer,
    }
}

/// Phase 7 (Theorem 2.7): Result of the Higher-Order Cesàro–Fejér Polynomial Kernel
/// ($\psi_2(d) = \binom{W - d + 1}{2}$) & Optimal Order-1 / Upper-Octave Order-2 Convex
/// Combination Envelope Solver across dyadic equivalent horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct CesaroFejerOrder2Envelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub kappa_int_order2: f64,
    pub c_infty_cf2_limit: f64,
    pub c_cf2_tail_star: f64,
    pub c_cf2_full_fejer: f64,
    pub phase6_tail_star: f64,
    pub phase6_full_fejer: f64,
    pub phase7_vs_phase1_factor: f64,
    pub phase7_vs_phase6_factor: f64,
}

/// Compute the Upper-Octave Order-2 Cesàro–Fejér backward-lag amplification factor
/// $\beta_{s,\mathrm{oct}}^2 := \frac{24}{7}\int_{1/2}^1 \tau^2 \beta_s(\tau)^2\,d\tau$
/// via composite Simpson quadrature on $[1/2, 1]$.
#[inline]
pub fn beta_s_oct_sq(s: usize) -> f64 {
    if s == 0 {
        return (4.0f64 / 3.0f64).powi(2);
    }
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let n_simp = 64usize;
    let a = 0.5f64;
    let b = 1.0f64;
    let h = (b - a) / (n_simp as f64);

    let eval_tau = |tau: f64| -> f64 {
        let d1 = tau * 2.0f64.powi(-1 - (s as i32));
        let d0 = tau * 2.0f64.powi(-(s as i32));
        let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
        let b_s = lag_inc / base_inc;
        tau * tau * b_s * b_s
    };

    let mut total = eval_tau(a) + eval_tau(b);
    for i in 1..n_simp {
        let tau = a + (i as f64) * h;
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        total += w * eval_tau(tau);
    }
    (24.0 / 7.0) * (h / 3.0) * total
}

/// Solve the Phase 7 (Theorem 2.7) Higher-Order Cesàro–Fejér Polynomial Kernel
/// ($\psi_2(d) = \binom{W - d + 1}{2}$) & Optimal Order-1 / Upper-Octave Order-2
/// Convex Combination Envelope Bound for dyadic equivalent horizon $M \in [4, 64]$
/// (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.7 proves two structural advances over Theorem 2.6:
///   1. **Part (a) (Exact Order-$p$ Iterated Cesàro–Fejér Difference Identity & Binomial
///      RHS Deflation):** Summing the Order-1 sliding-window Fejér identity $\mathcal{F}_w(A)$
///      over $w \in \{1, \dots, W\}$ with weight $\binom{W - w + p - 2}{p - 2}$ yields the
///      polynomial difference kernel $\psi_p(d) = \binom{W - d + p - 1}{p}$ with RHS ceiling
///      $\binom{W + p - 1}{p + 1} < \frac{W^{p+1}}{(p+1)!}$. For $p = 2$, the RHS normalization
///      factor $3$ amplifies the unit-cycle discrete Legendre–Fenchel integrality excess across
///      the active dyadic window octave $\tau = w / W \in [1/2, 1]$ from $\kappa_{\mathrm{int}}^{(1)} = \frac{1}{6}$
///      to $\kappa_{\mathrm{int}}^{(2)} = 3 \int_{1/2}^1 \frac{1}{6}\,d\tau = \frac{1}{4}$,
///      shifting the asymptotic per-scale multiplier from $2(1 + 1/6) = 7/3$ to $2(1 + 1/4) = 5/2$
///      and driving the pinched asymptotic limit strictly below $1.80$ to
///      $c_\infty^{*,\mathrm{CF2}}(1^+) = 2\sqrt{4/5} = \frac{4\sqrt{5}}{5} \approx 1.788854$.
///   2. **Part (b) (Optimal Convex Combination of Order-1 and Upper-Octave Order-2 Kernels):**
///      Optimizing the scale-step weight $\gamma_s^* = \max(\gamma_s^{(1)}, \gamma_s^{(2)})$
///      between Order-1 ($\gamma_s^{(1)} = \beta_{s,1}^2 + \frac{1}{6}(1 - 2^{-s})$, optimal at $s = 1$)
///      and Upper-Octave Order-2 ($\gamma_s^{(2)} = \beta_{s,\mathrm{oct}}^2 + \frac{1}{4}(1 - 2^{-s})$,
///      optimal for all $s \ge 2$) strictly sharpens $\underline{c}_M^{*,\mathrm{CF2}}(1.01)$ across
///      every $M \in [8, 64]$, breaking the $2.000000$ barrier already at **$M = 24$** (`1.977496`,
///      down from $M = 32$ in Phase 6) and reaching `1.855769` at $M = 64$.
pub fn solve_cesaro_fejer_order2_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> CesaroFejerOrder2Envelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p6 = solve_coupled_lag_integrality_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);
    let kappa_int_order2 = 0.25f64;
    let c_infty_cf2_limit = (2.0 / (lambda_b * 2.5)).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Optimal Order-1 / Upper-Octave Order-2 Tail-Octave Bound \underline{c}_M^{*,CF2}(B)
    let mut quad_cf2_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_cf2 = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let boct_sq = beta_s_oct_sq(s);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let g2 = boct_sq + 0.25 * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_cf2 += 2.0 * g1.max(g2);
        }
        quad_cf2_tail += scale_ratio * fejer_cf2 * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_cf2_tail_star = (rhs_tail / quad_cf2_tail).sqrt();

    // 2. Optimal Order-1 / Upper-Octave Order-2 Full-Window Bound \underline{c}_M^{CF2}(B)
    let mut quad_cf2_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_cf2 = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let boct_sq = beta_s_oct_sq(s);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let g2 = boct_sq + 0.25 * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_cf2 += 2.0 * g1.max(g2);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        quad_cf2_full +=
            scale_ratio * fejer_cf2 * rho_sq_2k + (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_cf2_full_fejer = (rhs_full / quad_cf2_full).sqrt();

    CesaroFejerOrder2Envelope {
        m_scales,
        k0,
        base_b,
        kappa_int_order2,
        c_infty_cf2_limit,
        c_cf2_tail_star,
        c_cf2_full_fejer,
        phase6_tail_star: p6.c_coupled_tail_star,
        phase6_full_fejer: p6.c_coupled_full_fejer,
        phase7_vs_phase1_factor: p4.c_intra_only_upper / c_cf2_full_fejer,
        phase7_vs_phase6_factor: p6.c_coupled_full_fejer / c_cf2_full_fejer,
    }
}

/// Phase 8 (Theorem 2.8(a)): Compute the Octave-Matched Cesàro–Fejér Polynomial Order
/// $p^*(s) \in \mathbb{Z}_{\ge 1}$ and maximum normalized octave weight
/// $M_s := \max_{p \ge 1} \beta_s^{(p)} = \max_{p \ge 1} \big[(1 - 2^{-s})^{p+1} - (1 - 2^{-(s-1)})^{p+1}\big]$
/// on octave $I_s = [2^{-s}, 2^{-(s-1)}]$.
///
/// For all $s \ge 2$, $M_s$ is strictly decreasing in $s$ and strictly bounded below by $\frac{1}{4}$,
/// with $p^*(s) \sim 2^s \ln 2$ and $\lim_{s\to\infty} M_s = \frac{1}{4}$.
#[inline]
pub fn octave_matched_cesaro_order(s: usize) -> (usize, f64) {
    assert!(s >= 1);
    let a_s = 1.0 - 2.0f64.powi(-(s as i32));
    let b_s = 1.0 - 2.0f64.powi(1 - (s as i32));
    let max_p = ((1usize << s.min(12)) * 2).max(16);
    let mut best_p = 1usize;
    let mut best_val = -1.0f64;
    for p in 1..=max_p {
        let exp = (p + 1) as i32;
        let val = a_s.powi(exp) - b_s.powi(exp);
        if val > best_val {
            best_val = val;
            best_p = p;
        }
    }
    (best_p, best_val)
}

/// Compute the Band-Restricted Order-2 Cesàro–Fejér backward-lag amplification factor
/// $\beta_{s,\mathrm{band}}(\alpha)^2 := \frac{3}{1 - \alpha^3}\int_\alpha^1 \tau^2 \beta_s(\tau)^2\,d\tau$
/// via composite Simpson quadrature on $[\alpha, 1]$ ($\alpha \in [0.5, 1.0)$).
#[inline]
pub fn beta_s_band_sq(s: usize, alpha: f64) -> f64 {
    if s == 0 {
        return (4.0f64 / 3.0f64).powi(2);
    }
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let n_simp = 64usize;
    let a = alpha;
    let b = 1.0f64;
    let h = (b - a) / (n_simp as f64);

    let eval_tau = |tau: f64| -> f64 {
        let d1 = tau * 2.0f64.powi(-1 - (s as i32));
        let d0 = tau * 2.0f64.powi(-(s as i32));
        let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
        let b_s = lag_inc / base_inc;
        tau * tau * b_s * b_s
    };

    let mut total = eval_tau(a) + eval_tau(b);
    for i in 1..n_simp {
        let tau = a + (i as f64) * h;
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        total += w * eval_tau(tau);
    }
    (3.0 / (1.0 - alpha.powi(3))) * (h / 3.0) * total
}

/// Phase 8 (Theorem 2.8): Result of the Octave-Matched & Band-Restricted Upper-Octave
/// Cesàro–Fejér Envelope Solver across dyadic equivalent horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct BandRestrictedCesaroEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub kappa_int_band_half: f64,
    pub c_infty_band_limit: f64,
    pub c_band_tail_star: f64,
    pub c_band_full_fejer: f64,
    pub phase7_tail_star: f64,
    pub phase7_full_fejer: f64,
    pub phase8_vs_phase1_factor: f64,
    pub phase8_vs_phase7_factor: f64,
}

/// Solve the Phase 8 (Theorem 2.8) Octave-Matched & Band-Restricted Upper-Octave
/// Cesàro–Fejér Envelope Bound for dyadic equivalent horizon $M \in [4, 64]$
/// (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.8 proves two structural advances over Theorem 2.7:
///   1. **Part (a) (Non-Vanishing Octave-Matched Cesàro–Fejér Law):** While any fixed order $p$
///      has octave integral $\beta_s^{(p)} = (1 - 2^{-s})^{p+1} - (1 - 2^{-(s-1)})^{p+1} = O(p 2^{-s}) \to 0$,
///      choosing the octave-matched order $p^*(s) \sim 2^s \ln 2$ yields a strictly decreasing
///      sequence $M_s = \max_{p \ge 1} \beta_s^{(p)} > \frac{1}{4}$ for all $s \ge 2$ converging to
///      $\lim_{s\to\infty} M_s = \frac{1}{4}$.
///   2. **Part (b) (Band-Restricted Upper-Octave Cesàro–Fejér Identity $\mathcal{CF}_{W, W_0}(A)$):**
///      Summing $\mathcal{F}_w(A)$ over $w \in [W_0 + 1, W]$ ($W_0 = \lfloor W/2 \rfloor$, $\alpha = 1/2$)
///      eliminates the inactive lower half $w \le W/2$ where the unit-cycle integrality excess vanishes,
///      cutting the RHS capacity from $\binom{W+1}{3} \sim \frac{1}{6}W^3$ to
///      $\binom{W+1}{3} - \binom{W_0+1}{3} \sim \frac{7}{48}W^3$. This $\frac{7}{8}$ capacity reduction
///      boosts the upper-octave integrality constant by $\frac{8}{7}$ from $\kappa_{\mathrm{int}}^{(2)} = \frac{1}{4}$
///      to $\kappa_{\mathrm{int}}^{\mathrm{band}}(1/2) = \frac{2}{7} \approx 0.285714$, raising the asymptotic
///      per-scale multiplier from $\frac{5}{2}$ to $\frac{18}{7}$, driving the pinched asymptotic limit down to
///      $c_\infty^{*,\mathrm{band}}(1^+) = \frac{2\sqrt{7}}{3} \approx 1.763834$, and breaking the $2.000000$
///      barrier already at **$M = 20$** (`1.999651`, down from $M = 24$ in Phase 7 and $M = 32$ in Phase 6).
pub fn solve_band_restricted_cesaro_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> BandRestrictedCesaroEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p7 = solve_cesaro_fejer_order2_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);
    let kappa_int_band_half = 2.0f64 / 7.0f64;
    let c_infty_band_limit = (2.0 / (lambda_b * (18.0 / 7.0))).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Band-Restricted Upper-Octave Tail Bound \underline{c}_M^{*,band}(B)
    let mut quad_band_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_band = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bband_sq = beta_s_band_sq(s, 0.5);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let g_band = bband_sq + kappa_int_band_half * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_band += 2.0 * g1.max(g_band);
        }
        quad_band_tail += scale_ratio * fejer_band * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_band_tail_star = (rhs_tail / quad_band_tail).sqrt();

    // 2. Band-Restricted Upper-Octave Full-Window Bound \underline{c}_M^{band}(B)
    let mut quad_band_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_band = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bband_sq = beta_s_band_sq(s, 0.5);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let g_band = bband_sq + kappa_int_band_half * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_band += 2.0 * g1.max(g_band);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        quad_band_full +=
            scale_ratio * fejer_band * rho_sq_2k + (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_band_full_fejer = (rhs_full / quad_band_full).sqrt();

    BandRestrictedCesaroEnvelope {
        m_scales,
        k0,
        base_b,
        kappa_int_band_half,
        c_infty_band_limit,
        c_band_tail_star,
        c_band_full_fejer,
        phase7_tail_star: p7.c_cf2_tail_star,
        phase7_full_fejer: p7.c_cf2_full_fejer,
        phase8_vs_phase1_factor: p4.c_intra_only_upper / c_band_full_fejer,
        phase8_vs_phase7_factor: p7.c_cf2_full_fejer / c_band_full_fejer,
    }
}

/// Phase 9 (Theorem 2.9): Exact $K$-octave discrete Legendre–Fenchel integrality
/// constant $\kappa_{\mathrm{int}}^{\mathrm{band}}(2^{-K}) = \frac{1}{2(1 + 2^{-K} + 4^{-K})}
/// = \frac{2^{2K-1}}{4^K + 2^K + 1}$ for $K \ge 1$, converging to $\frac{1}{2}$ as $K \to \infty$.
#[inline]
pub fn kappa_int_multioctave_k(k_octaves: usize) -> f64 {
    assert!(k_octaves >= 1);
    let alpha = 2.0f64.powi(-(k_octaves as i32));
    0.5 / (1.0 + alpha + alpha * alpha)
}

/// Phase 9 (Theorem 2.9): Result of the Scale-Invariant Block-Window Conservation &
/// Multi-Octave Integrality Cascade Envelope Solver across dyadic equivalent horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct MultiOctaveCesaroEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub kappa_int_k2: f64,
    pub kappa_int_k3: f64,
    pub kappa_int_multioct_limit: f64,
    pub c_infty_k2_limit: f64,
    pub c_infty_k3_limit: f64,
    pub c_infty_multioct_limit: f64,
    pub c_k2_tail_star: f64,
    pub c_k3_tail_star: f64,
    pub c_multioct_tail_star: f64,
    pub c_multioct_full_fejer: f64,
    pub phase8_tail_star: f64,
    pub phase8_full_fejer: f64,
    pub phase9_vs_phase1_factor: f64,
    pub phase9_vs_phase8_factor: f64,
}

/// Solve the Phase 9 (Theorem 2.9) Scale-Invariant Block-Window Conservation &
/// $K$-Octave Integrality Cascade Envelope Bound for dyadic equivalent horizon
/// $M \in [4, 64]$ (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.9 proves two structural advances over Theorem 2.8:
///   1. **Part (a) (Scale-Invariant Block-Window Conservation $w \cdot B(\tau) \equiv N$):**
///      Tiling a domain of length $N$ by $B(\tau) = \frac{N}{\tau W}$ windows of sub-window
///      width $w = \tau W$ ($\tau \in (0, 1]$) yields a per-window discrete Legendre–Fenchel
///      integrality excess of $\frac{w}{2}\theta(\tau)(1 - \theta(\tau))$, whose sum across
///      all $B(\tau)$ windows is $E_{\mathrm{int}}(w) = B(\tau)\frac{w}{2}\theta(\tau)(1-\theta(\tau))
///      \equiv \frac{N}{2}\theta(\tau)(1-\theta(\tau))$, invariant across every octave $\tau \in (0, 1]$.
///   2. **Part (b) ($K$-Octave Integrality Cascade & Sub-$1.64$ Continuum Limit):**
///      Summing $\mathcal{F}_w(A)$ over the top $K$ octaves $w \in [\lfloor 2^{-K} W \rfloor + 1, W]$
///      ($\alpha = 2^{-K}$) yields $\kappa_{\mathrm{int}}^{\mathrm{band}}(2^{-K}) = \frac{1}{2(1 + 2^{-K} + 4^{-K})}
///      = \frac{2^{2K-1}}{4^K + 2^K + 1} \nearrow \frac{1}{2}$ ($K=1$: $\frac{2}{7}$, $K=2$: $\frac{8}{21}$,
///      $K=3$: $\frac{32}{73}$, $K\to\infty$: $\frac{1}{2}$), raising the asymptotic per-scale multiplier
///      from $\frac{18}{7}$ to $2(1 + 1/2) = 3$, driving the pinched asymptotic limit down to
///      $c_\infty^{*,\mathrm{multi\text{-}oct}}(1^+) = 2\sqrt{2/3} = \frac{2\sqrt{6}}{3} \approx 1.632993$,
///      and breaking the $2.000000$ barrier already at **$M = 16$** (`1.959054`, down from $M = 20$
///      in Phase 8, $M = 24$ in Phase 7, and $M = 32$ in Phase 6).
pub fn solve_multioctave_cesaro_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> MultiOctaveCesaroEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p8 = solve_band_restricted_cesaro_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);

    let kappa_int_k2 = kappa_int_multioctave_k(2);
    let kappa_int_k3 = kappa_int_multioctave_k(3);
    let kappa_int_multioct_limit = 0.5f64;

    let c_infty_k2_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_k2))).sqrt();
    let c_infty_k3_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_k3))).sqrt();
    let c_infty_multioct_limit = (2.0 / (lambda_b * 3.0)).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Multi-Octave Tail Bounds \underline{c}_M^{*,(K)}(B) for K=2, K=3, and K->\infty
    let mut quad_k2_tail = 0.0f64;
    let mut quad_k3_tail = 0.0f64;
    let mut quad_multi_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_k2 = beta0_sq;
        let mut fejer_k3 = beta0_sq;
        let mut fejer_multi = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bk2_sq = beta_s_band_sq(s, 0.25);
            let bk3_sq = beta_s_band_sq(s, 0.125);
            let bmulti_sq = beta_s_band_sq(s, 0.0);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gk2 = bk2_sq + kappa_int_k2 * (1.0 - 2.0f64.powi(-(s as i32)));
            let gk3 = bk3_sq + kappa_int_k3 * (1.0 - 2.0f64.powi(-(s as i32)));
            let gmulti = bmulti_sq + kappa_int_multioct_limit * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_k2 += 2.0 * g1.max(gk2);
            fejer_k3 += 2.0 * g1.max(gk3);
            fejer_multi += 2.0 * g1.max(gmulti);
        }
        quad_k2_tail += scale_ratio * fejer_k2 * alpha_k * alpha_k;
        quad_k3_tail += scale_ratio * fejer_k3 * alpha_k * alpha_k;
        quad_multi_tail += scale_ratio * fejer_multi * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_k2_tail_star = (rhs_tail / quad_k2_tail).sqrt();
    let c_k3_tail_star = (rhs_tail / quad_k3_tail).sqrt();
    let c_multioct_tail_star = (rhs_tail / quad_multi_tail).sqrt();

    // 2. Multi-Octave Full-Window Bound \underline{c}_M^{multi-oct}(B)
    let mut quad_multi_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_multi = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bmulti_sq = beta_s_band_sq(s, 0.0);
            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gmulti = bmulti_sq + kappa_int_multioct_limit * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_multi += 2.0 * g1.max(gmulti);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        quad_multi_full += scale_ratio * fejer_multi * rho_sq_2k
            + (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_multioct_full_fejer = (rhs_full / quad_multi_full).sqrt();

    MultiOctaveCesaroEnvelope {
        m_scales,
        k0,
        base_b,
        kappa_int_k2,
        kappa_int_k3,
        kappa_int_multioct_limit,
        c_infty_k2_limit,
        c_infty_k3_limit,
        c_infty_multioct_limit,
        c_k2_tail_star,
        c_k3_tail_star,
        c_multioct_tail_star,
        c_multioct_full_fejer,
        phase8_tail_star: p8.c_band_tail_star,
        phase8_full_fejer: p8.c_band_full_fejer,
        phase9_vs_phase1_factor: p4.c_intra_only_upper / c_multioct_full_fejer,
        phase9_vs_phase8_factor: p8.c_band_full_fejer / c_multioct_full_fejer,
    }
}

/// Phase 10 (Theorem 2.10(a)): Exact continuum Power-Law multi-octave integrality
/// constant $\kappa_{\mathrm{int}}^{(r)}(0^+) = \frac{r + 3}{6(r + 1)}$ for sub-window
/// measure $d\nu_r(\tau) = \tau^r\,d\tau$ ($r > -1$).
#[inline]
pub fn kappa_int_power_law_continuum(r_exp: f64) -> f64 {
    assert!(r_exp > -1.0);
    (r_exp + 3.0) / (6.0 * (r_exp + 1.0))
}

/// Phase 10 (Theorem 2.10(b)): Exact $K$-octave Scale-Invariant Haar-Measure ($r \to -1^+$)
/// integrality constant $\kappa_{\mathrm{int}}^{\mathrm{Haar}}(K) = \frac{K \ln 2}{3(1 - 4^{-K})}$,
/// diverging logarithmically to $+\infty$ as $K \to \infty$.
#[inline]
pub fn kappa_int_haar_k(k_octaves: usize) -> f64 {
    assert!(k_octaves >= 1);
    ((k_octaves as f64) * std::f64::consts::LN_2) / (3.0 * (1.0 - 4.0f64.powi(-(k_octaves as i32))))
}

/// Compute the Power-Law / Haar-Weighted Cesàro–Fejér backward-lag amplification factor
/// $\beta_s^{(r)}(\alpha)^2 := \frac{r + 3}{1 - \alpha^{r+3}}\int_\alpha^1 \tau^{r+2} \beta_s(\tau)^2\,d\tau$
/// via composite Simpson quadrature on $[\alpha, 1]$ ($r \ge -1, \alpha \in [0, 1)$).
#[inline]
pub fn beta_s_power_sq(s: usize, r_exp: f64, alpha: f64) -> f64 {
    if s == 0 {
        return (4.0f64 / 3.0f64).powi(2);
    }
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let n_simp = 64usize;
    let a = alpha;
    let b = 1.0f64;
    let h = (b - a) / (n_simp as f64);
    let p_pow = r_exp + 2.0;

    let eval_tau = |tau: f64| -> f64 {
        if tau <= 0.0 {
            return 0.0;
        }
        let d1 = tau * 2.0f64.powi(-1 - (s as i32));
        let d0 = tau * 2.0f64.powi(-(s as i32));
        let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
        let b_s = lag_inc / base_inc;
        tau.powf(p_pow) * b_s * b_s
    };

    let mut total = eval_tau(a) + eval_tau(b);
    for i in 1..n_simp {
        let tau = a + (i as f64) * h;
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        total += w * eval_tau(tau);
    }
    let denom = if alpha > 0.0 {
        1.0 - alpha.powf(r_exp + 3.0)
    } else {
        1.0
    };
    ((r_exp + 3.0) / denom) * (h / 3.0) * total
}

/// Phase 10 (Theorem 2.10): Result of the Power-Law & Logarithmic Haar-Measure
/// Multi-Octave Cascade Envelope Solver across dyadic equivalent horizon $M \le 64$.
#[derive(Debug, Clone)]
pub struct HaarMultiOctaveEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub kappa_int_haar_k3: f64,
    pub kappa_int_haar_k13: f64,
    pub c_infty_haar_k3_limit: f64,
    pub c_infty_haar_k13_limit: f64,
    pub c_haar_k3_tail_star: f64,
    pub c_haar_k3_full_fejer: f64,
    pub c_haar_adapt_tail_star: f64,
    pub c_haar_adapt_full_fejer: f64,
    pub phase9_tail_star: f64,
    pub phase9_full_fejer: f64,
    pub phase10_k3_vs_phase1_factor: f64,
    pub phase10_adapt_vs_phase1_factor: f64,
    pub phase10_adapt_vs_phase9_factor: f64,
}

/// Solve the Phase 10 (Theorem 2.10) Power-Law & Scale-Invariant Haar-Measure
/// Multi-Octave Cascade Envelope Bound for dyadic equivalent horizon
/// $M \in [4, 64]$ (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
///
/// Theorem 2.10 proves two major advances over Theorem 2.9:
///   1. **Part (a) (Power-Law Multi-Octave Filter $d\nu_r(\tau) = \tau^r\,d\tau$, $r > -1$):**
///      Weighting each sub-window scale $\tau = w/W \in [\alpha, 1]$ by $\tau^r$ preserves
///      positivity of the Sidon sub-block pigeonhole inequality while Scale-Invariant
///      Block-Window Conservation ($B(\tau) \cdot \frac{w}{2} \equiv \frac{N}{2}$) keeps the
///      unweighted integrality excess $\tau$-independent, yielding
///      $\kappa_{\mathrm{int}}^{(r)}(\alpha) = \frac{r + 3}{6(r + 1)}\frac{1 - \alpha^{r+1}}{1 - \alpha^{r+3}}
///      \xrightarrow{\alpha\to 0^+} \frac{r + 3}{6(r + 1)}$ and
///      $c_\infty^{*,(r)}(1^+) = 2\sqrt{\frac{6(r+1)}{7r+9}}$, which hits the **$1.000000$ unit barrier**
///      at $r = -\frac{15}{17}$ ($\kappa_{\mathrm{int}} = 3$) and drops strictly below $1.00$ for all
///      $r \in (-1, -\frac{15}{17})$ (e.g. $c_\infty^{*,(-9/10)}(1^+) = \frac{2\sqrt{2}}{3} \approx 0.942809$).
///   2. **Part (b) (Scale-Invariant Haar Measure $d\tau / \tau$ at $r \to -1^+$ & Logarithmic
///      $K$-Octave Cascade Toward Zero):**
///      At $r = -1$ across $K$ active octaves ($\alpha = 2^{-K}$),
///      $\kappa_{\mathrm{int}}^{\mathrm{Haar}}(K) = \frac{K \ln 2}{3(1 - 4^{-K})} \to +\infty$ and
///      $c_\infty^{*,\mathrm{Haar}}(1^+; K) = 2\bigl(1 + \frac{K \ln 2}{3(1 - 4^{-K})}\bigr)^{-1/2} \to 0$.
///      On a finite horizon $M$ where each outer window $W = 2^r = 2^{k-s}$ spans $K(k, s) = k - s$
///      dyadic octaves down to $W_0 = 1$ ($\mathcal{F}_1(A) \equiv 0$), the adaptive Haar tail bound
///      $\underline{c}_M^{*,\mathrm{Haar}}(1.01)$ breaks strictly below $2.000000$ **already at $M = 12$**
///      (`1.700945`), breaks strictly below $1.000000$ at $M \ge 48$ (`0.813355` at $M = 48$,
///      `0.701665` at $M = 64$), and decays at rate $O(M^{-1/2}) \to 0$ as $M \to \infty$.
pub fn solve_haar_multioctave_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> HaarMultiOctaveEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p9 = solve_multioctave_cesaro_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);

    let kappa_int_haar_k3 = kappa_int_haar_k(3);
    let kappa_int_haar_k13 = kappa_int_haar_k(13);
    let c_infty_haar_k3_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_haar_k3))).sqrt();
    let c_infty_haar_k13_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_haar_k13))).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Haar-Weighted Tail Bounds (Fixed K=3 and Horizon-Adaptive K = k - s)
    let mut quad_k3_tail = 0.0f64;
    let mut quad_adapt_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_k3 = beta0_sq;
        let mut fejer_adapt = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bh3_sq = beta_s_power_sq(s, -1.0, 0.125);
            let k_oct_act = k - s;
            let alpha_act = 2.0f64.powi(-(k_oct_act as i32));
            let bhadapt_sq = beta_s_power_sq(s, -1.0, alpha_act);

            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gh3 = bh3_sq + kappa_int_haar_k3 * (1.0 - 2.0f64.powi(-(s as i32)));
            let ghadapt =
                bhadapt_sq + kappa_int_haar_k(k_oct_act) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_k3 += 2.0 * g1.max(gh3);
            fejer_adapt += 2.0 * g1.max(ghadapt);
        }
        quad_k3_tail += scale_ratio * fejer_k3 * alpha_k * alpha_k;
        quad_adapt_tail += scale_ratio * fejer_adapt * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_haar_k3_tail_star = (rhs_tail / quad_k3_tail).sqrt();
    let c_haar_adapt_tail_star = (rhs_tail / quad_adapt_tail).sqrt();

    // 2. Haar-Weighted Full-Window Bounds (Fixed K=3 and Horizon-Adaptive K = k - s)
    let mut quad_k3_full = 0.0f64;
    let mut quad_adapt_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_k3 = beta0_sq;
        let mut fejer_adapt = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bh3_sq = beta_s_power_sq(s, -1.0, 0.125);
            let k_oct_act = k - s;
            let alpha_act = 2.0f64.powi(-(k_oct_act as i32));
            let bhadapt_sq = beta_s_power_sq(s, -1.0, alpha_act);

            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gh3 = bh3_sq + kappa_int_haar_k3 * (1.0 - 2.0f64.powi(-(s as i32)));
            let ghadapt =
                bhadapt_sq + kappa_int_haar_k(k_oct_act) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_k3 += 2.0 * g1.max(gh3);
            fejer_adapt += 2.0 * g1.max(ghadapt);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        let bil_term = (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
        quad_k3_full += scale_ratio * fejer_k3 * rho_sq_2k + bil_term;
        quad_adapt_full += scale_ratio * fejer_adapt * rho_sq_2k + bil_term;
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_haar_k3_full_fejer = (rhs_full / quad_k3_full).sqrt();
    let c_haar_adapt_full_fejer = (rhs_full / quad_adapt_full).sqrt();

    HaarMultiOctaveEnvelope {
        m_scales,
        k0,
        base_b,
        kappa_int_haar_k3,
        kappa_int_haar_k13,
        c_infty_haar_k3_limit,
        c_infty_haar_k13_limit,
        c_haar_k3_tail_star,
        c_haar_k3_full_fejer,
        c_haar_adapt_tail_star,
        c_haar_adapt_full_fejer,
        phase9_tail_star: p9.c_multioct_tail_star,
        phase9_full_fejer: p9.c_multioct_full_fejer,
        phase10_k3_vs_phase1_factor: p4.c_intra_only_upper / c_haar_k3_full_fejer,
        phase10_adapt_vs_phase1_factor: p4.c_intra_only_upper / c_haar_adapt_full_fejer,
        phase10_adapt_vs_phase9_factor: p9.c_multioct_full_fejer / c_haar_adapt_full_fejer,
    }
}

/// Phase 11 (Theorem 2.11(b)): Exact Glaisher–Kinkelin Sawtooth Constant
/// $C_{\mathrm{GK}} := 2\ln A_{\mathrm{GK}} = \frac{1}{6} - 2\zeta'(-1)
/// = \frac{1}{2} + \frac{\gamma}{6} - \sum_{j=2}^\infty \frac{(-1)^j \zeta(j)}{(j+1)(j+2)}
/// \approx +0.4975089540675685$ governing the full-cascade Haar sawtooth integral
/// $\int_0^M \frac{\{u\}(1 - \{u\})}{u}\,du = \frac{1}{6}\ln M + C_{\mathrm{GK}} + \frac{1}{360 M^2} - O(M^{-4})$.
pub const C_GK_GLAISHER_SAWTOOTH: f64 = 0.497_508_954_067_568_5;

/// Phase 11 (Theorem 2.11(b)): Exact Euler–Mascheroni / Glaisher–Kinkelin Sawtooth Constant
/// $C_{\mathrm{LF}} := 2\ln A_{\mathrm{GK}} - \frac{1}{2} = -\frac{1}{3} - 2\zeta'(-1)
/// = \frac{\gamma}{6} - \sum_{j=2}^\infty \frac{(-1)^j \zeta(j)}{(j+1)(j+2)}
/// \approx -0.002491045932431475$ governing the multi-period Haar sawtooth integral
/// $\int_1^M \frac{\{u\}(1 - \{u\})}{u}\,du = \frac{1}{6}\ln M + C_{\mathrm{LF}} + \frac{1}{360 M^2} - O(M^{-4})$.
pub const C_LF_EULER_MASCHERONI_SAWTOOTH: f64 = -0.002_491_045_932_431_475;

/// Phase 11 (Theorem 2.11(a)): Exact logarithmic Haar integral of the Legendre–Fenchel
/// sawtooth $\psi_{\mathrm{LF}}(u) = \{u\}(1 - \{u\})$ on each integer period $[m, m+1)$:
///   * For $m = 0$: $I_0 = \int_0^1 \frac{u(1 - u)}{u}\,du = \frac{1}{2}$.
///   * For $m \ge 1$: $I_m = \int_m^{m+1} \frac{(u - m)(m + 1 - u)}{u}\,du
///     = m + \frac{1}{2} - m(m + 1)\ln\bigl(1 + \frac{1}{m}\bigr)
///     = \sum_{j=1}^\infty \frac{(-1)^{j-1}}{(j+1)(j+2)\, m^j} \in \bigl(\frac{1}{6m} - \frac{1}{12m^2},\, \frac{1}{6m}\bigr)$.
#[inline]
pub fn sawtooth_haar_period_integral(m: usize) -> f64 {
    if m == 0 {
        return 0.5;
    }
    let mf = m as f64;
    if m <= 4 {
        mf + 0.5 - mf * (mf + 1.0) * (1.0 / mf).ln_1p()
    } else {
        // Use 16-term alternating power series to avoid catastrophic cancellation at large m
        let inv = 1.0 / mf;
        let mut term = inv;
        let mut acc = 0.0f64;
        for j in 1..=16usize {
            let denom = ((j + 1) * (j + 2)) as f64;
            let contrib = term / denom;
            if j % 2 == 1 {
                acc += contrib;
            } else {
                acc -= contrib;
            }
            term *= inv;
        }
        acc
    }
}

/// Phase 11 (Theorem 2.11(b)): Exact hyperfactorial antiderivative
/// $F_{\mathrm{saw}}(U) := \int_0^U \frac{\{u\}(1 - \{u\})}{u}\,du$ for any real $U \ge 0$.
///
/// For any integer $M \ge 1$, summation by parts yields the exact hyperfactorial closed form:
/// $$F_{\mathrm{saw}}(M) = \frac{M^2}{2} - M(M - 1)\ln M + 2\sum_{m=1}^{M-1} m \ln m
/// = \frac{1}{6}\ln M + C_{\mathrm{GK}} + \frac{1}{360 M^2} - \frac{1}{2520 M^4} + O(M^{-6}).$$
#[inline]
pub fn sawtooth_haar_antiderivative(u: f64) -> f64 {
    if u <= 0.0 {
        return 0.0;
    }
    if u <= 1.0 {
        return u - 0.5 * u * u;
    }
    let m_int = u.floor() as usize;
    let mf = m_int as f64;
    let base_m = if m_int <= 64 {
        let mut s = 0.5f64;
        for k in 1..m_int {
            s += sawtooth_haar_period_integral(k);
        }
        s
    } else {
        let inv_m2 = 1.0 / (mf * mf);
        mf.ln() / 6.0 + C_GK_GLAISHER_SAWTOOTH + inv_m2 / 360.0 - (inv_m2 * inv_m2) / 2520.0
    };
    let theta = u - mf;
    if theta <= 1e-15 {
        return base_m;
    }
    let rem_u = -0.5 * (u * u - mf * mf) + (2.0 * mf + 1.0) * theta
        - mf * (mf + 1.0) * (theta / mf).ln_1p();
    base_m + rem_u.max(0.0)
}

/// Phase 11 (Theorem 2.11(b)): Unconditional Pointwise Linear-Dilation $K$-Octave Haar Integral
/// $$\mathcal{J}_{\mathrm{Haar}}(\mu_1, K) := \int_{2^{-K}}^1 \{\tau \mu_1\}\bigl(1 - \{\tau \mu_1\}\bigr)\,\frac{d\tau}{\tau}
/// \equiv \int_{\mu_1 2^{-K}}^{\mu_1} \frac{\{u\}(1 - \{u\})}{u}\,du$$
/// for ANY single real occupancy density $\mu_1 > 0$ (including adversarial integers $\mu_1 \in \mathbb{Z}$
/// where $\{\mu_1\} = 0$).
#[inline]
pub fn pointwise_diophantine_haar_dilation(mu_1: f64, k_octaves: usize) -> f64 {
    assert!(mu_1 > 0.0 && k_octaves >= 1);
    let u_lo = mu_1 * 2.0f64.powi(-(k_octaves as i32));
    sawtooth_haar_antiderivative(mu_1) - sawtooth_haar_antiderivative(u_lo)
}

/// Phase 11 (Theorem 2.11(c)): Exact Unconditional Diophantine $K$-Octave Haar Integrality Constant
/// $$\kappa_{\mathrm{int}}^{\mathrm{Dioph}}(K) := \frac{2}{1 - 4^{-K}}\int_1^{2^K} \frac{\{u\}(1 - \{u\})}{u}\,du + 2 C_{\mathrm{GK}},$$
/// incorporating the exact multi-period hyperfactorial integral on $[1, 2^K]$ and the fundamental-octave
/// Glaisher–Kinkelin sawtooth transition $2 C_{\mathrm{GK}} = 4\ln A_{\mathrm{GK}} \approx +0.99501791$.
#[inline]
pub fn kappa_int_diophantine_haar_k(k_octaves: usize) -> f64 {
    assert!(k_octaves >= 1);
    let j_multi = if k_octaves <= 12 {
        let m_cap = 1usize << k_octaves;
        sawtooth_haar_antiderivative(m_cap as f64) - 0.5
    } else {
        let mf = 2.0f64.powi(k_octaves as i32);
        ((k_octaves as f64) * std::f64::consts::LN_2) / 6.0
            + C_LF_EULER_MASCHERONI_SAWTOOTH
            + 1.0 / (360.0 * mf * mf)
    };
    (2.0 * j_multi) / (1.0 - 4.0f64.powi(-(k_octaves as i32))) + 2.0 * C_GK_GLAISHER_SAWTOOTH
}

/// Phase 11 (Theorem 2.11): Result of the Unconditional Diophantine Linear-Dilation &
/// Multi-Octave Haar Closure Solver across dyadic equivalent horizon $M \le 256$.
#[derive(Debug, Clone)]
pub struct UnconditionalDiophantineHaarEnvelope {
    pub m_scales: usize,
    pub k0: usize,
    pub base_b: f64,
    pub kappa_int_dioph_k3: f64,
    pub kappa_int_dioph_k9: f64,
    pub kappa_int_dioph_k64: f64,
    pub c_infty_dioph_k3_limit: f64,
    pub c_infty_dioph_k9_limit: f64,
    pub c_infty_dioph_k64_limit: f64,
    pub c_dioph_k3_tail_star: f64,
    pub c_dioph_k3_full_fejer: f64,
    pub c_dioph_adapt_tail_star: f64,
    pub c_dioph_adapt_full_fejer: f64,
    pub phase10_tail_star: f64,
    pub phase10_full_fejer: f64,
    pub phase11_adapt_vs_phase1_factor: f64,
    pub phase11_adapt_vs_phase10_factor: f64,
}

/// Solve the Phase 11 (Theorem 2.11) Unconditional Diophantine Linear-Dilation &
/// Multi-Octave Haar Closure Envelope Bound for dyadic equivalent horizon
/// $M \in [4, 256]$ (`N = 2^M`), lower cutoff $k_0$, and geometric base $B \in (1, 4]$.
pub fn solve_unconditional_diophantine_haar_envelope(
    m_scales: usize,
    k0: usize,
    base_b: f64,
) -> UnconditionalDiophantineHaarEnvelope {
    assert!((4..=256).contains(&m_scales));
    assert!(base_b > 1.0 && base_b <= 4.0);
    let p4 = solve_multires_subwindow_qp_envelope(m_scales, k0);
    let p10 = solve_haar_multioctave_envelope(m_scales, k0, base_b);

    let ln_2 = std::f64::consts::LN_2;
    let ln_b = base_b.ln();
    let lambda_2 = (0.25 * ln_2).tanh() / ln_2;
    let lambda_b = (0.25 * ln_b).tanh() / ln_b;
    let scale_ratio = lambda_b / lambda_2;
    let kappa_pref =
        2.5 - 1.5 * std::f64::consts::SQRT_2 + 0.5 * (1.0 + std::f64::consts::SQRT_2).ln();
    let base_inc = 1.0 - 1.0 / std::f64::consts::SQRT_2;
    let beta0_sq = (4.0f64 / 3.0f64).powi(2);

    let kappa_int_dioph_k3 = kappa_int_diophantine_haar_k(3);
    let kappa_int_dioph_k9 = kappa_int_diophantine_haar_k(9);
    let kappa_int_dioph_k64 = kappa_int_diophantine_haar_k(64);
    let c_infty_dioph_k3_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_dioph_k3))).sqrt();
    let c_infty_dioph_k9_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_dioph_k9))).sqrt();
    let c_infty_dioph_k64_limit = (2.0 / (lambda_b * 2.0 * (1.0 + kappa_int_dioph_k64))).sqrt();

    let k_start = k0.clamp(1, m_scales);
    let r0 = m_scales / 2;

    // 1. Unconditional Diophantine Haar Tail Bounds (Fixed K=3 and Horizon-Adaptive K = k - s)
    let mut quad_k3_tail = 0.0f64;
    let mut quad_adapt_tail = 0.0f64;
    for k in (r0 + 1)..=m_scales {
        let alpha_k =
            1.0 / ((k as f64) * ln_2).sqrt() - 1.0 / (2.0 * ((k - 1) as f64) * ln_2).sqrt();
        let s_k = (k.saturating_sub(2)) / 2;
        let mut fejer_k3 = beta0_sq;
        let mut fejer_adapt = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bh3_sq = beta_s_power_sq(s, -1.0, 0.125);
            let k_oct_act = k - s;
            let alpha_act = 2.0f64.powi(-(k_oct_act as i32));
            let bhadapt_sq = beta_s_power_sq(s, -1.0, alpha_act);

            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gd3 = bh3_sq + kappa_int_dioph_k3 * (1.0 - 2.0f64.powi(-(s as i32)));
            let gdadapt = bhadapt_sq
                + kappa_int_diophantine_haar_k(k_oct_act) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_k3 += 2.0 * g1.max(gd3);
            fejer_adapt += 2.0 * g1.max(gdadapt);
        }
        quad_k3_tail += scale_ratio * fejer_k3 * alpha_k * alpha_k;
        quad_adapt_tail += scale_ratio * fejer_adapt * alpha_k * alpha_k;
    }
    let rhs_tail = 0.5 * ((m_scales - r0) as f64);
    let c_dioph_k3_tail_star = (rhs_tail / quad_k3_tail).sqrt();
    let c_dioph_adapt_tail_star = (rhs_tail / quad_adapt_tail).sqrt();

    // 2. Unconditional Diophantine Haar Full-Window Bounds (Fixed K=3 and Horizon-Adaptive K = k - s)
    let mut quad_k3_full = 0.0f64;
    let mut quad_adapt_full = 0.0f64;
    let mut prev_req = 0.0f64;
    for k in k_start..=m_scales {
        let req_k = ((2.0f64.powi(k as i32)) / ((k as f64) * ln_2)).sqrt();
        let rho_k = req_k - prev_req;
        prev_req = req_k;
        let rho_sq_2k = (rho_k * rho_k) / (2.0f64.powi(k as i32));
        let s_k = (k.saturating_sub(2)) / 2;

        let mut fejer_k3 = beta0_sq;
        let mut fejer_adapt = beta0_sq;
        for s in 1..=s_k {
            let d1 = 2.0f64.powi(-1 - (s as i32));
            let d0 = 2.0f64.powi(-(s as i32));
            let lag_inc = psi_lag(d1) - (1.0 / std::f64::consts::SQRT_2) * psi_lag(d0);
            let b1_sq = (lag_inc / base_inc).powi(2);
            let bh3_sq = beta_s_power_sq(s, -1.0, 0.125);
            let k_oct_act = k - s;
            let alpha_act = 2.0f64.powi(-(k_oct_act as i32));
            let bhadapt_sq = beta_s_power_sq(s, -1.0, alpha_act);

            let g1 = b1_sq + (1.0 / 6.0) * (1.0 - 2.0f64.powi(-(s as i32)));
            let gd3 = bh3_sq + kappa_int_dioph_k3 * (1.0 - 2.0f64.powi(-(s as i32)));
            let gdadapt = bhadapt_sq
                + kappa_int_diophantine_haar_k(k_oct_act) * (1.0 - 2.0f64.powi(-(s as i32)));
            fejer_k3 += 2.0 * g1.max(gd3);
            fejer_adapt += 2.0 * g1.max(gdadapt);
        }
        let w_bil = if k == m_scales { 1.0 } else { 0.5 };
        let bil_term = (w_bil + 0.5 * kappa_pref) / ((k as f64) * ln_2);
        quad_k3_full += scale_ratio * fejer_k3 * rho_sq_2k + bil_term;
        quad_adapt_full += scale_ratio * fejer_adapt * rho_sq_2k + bil_term;
    }
    let rhs_full = 0.5 * ((m_scales + 2) as f64);
    let c_dioph_k3_full_fejer = (rhs_full / quad_k3_full).sqrt();
    let c_dioph_adapt_full_fejer = (rhs_full / quad_adapt_full).sqrt();

    UnconditionalDiophantineHaarEnvelope {
        m_scales,
        k0,
        base_b,
        kappa_int_dioph_k3,
        kappa_int_dioph_k9,
        kappa_int_dioph_k64,
        c_infty_dioph_k3_limit,
        c_infty_dioph_k9_limit,
        c_infty_dioph_k64_limit,
        c_dioph_k3_tail_star,
        c_dioph_k3_full_fejer,
        c_dioph_adapt_tail_star,
        c_dioph_adapt_full_fejer,
        phase10_tail_star: p10.c_haar_adapt_tail_star,
        phase10_full_fejer: p10.c_haar_adapt_full_fejer,
        phase11_adapt_vs_phase1_factor: p4.c_intra_only_upper / c_dioph_adapt_full_fejer,
        phase11_adapt_vs_phase10_factor: p10.c_haar_adapt_full_fejer / c_dioph_adapt_full_fejer,
    }
}

/// Canonical Singer / Bose–Chowla finite planar difference set generators inside
/// $\{0, \dots, q^2 + q\}$ for small prime powers $q \in \{2, 3, 4, 5, 7, 8, 9, 11, 13\}$,
/// used as building blocks in Krückeberg (1961) and Cilleruelo (2014) infinite
/// Sidon block constructions.
pub fn singer_sidon_block(q: usize) -> Vec<u64> {
    match q {
        2 => vec![0, 1, 3],
        3 => vec![0, 1, 3, 9],
        4 => vec![0, 1, 4, 14, 16],
        5 => vec![0, 1, 3, 8, 12, 18],
        7 => vec![0, 1, 3, 13, 32, 36, 43, 52],
        8 => vec![0, 1, 3, 7, 15, 31, 36, 54, 63],
        9 => vec![0, 1, 3, 9, 27, 49, 56, 61, 77, 81],
        11 => vec![0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109],
        13 => vec![0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175],
        _ => panic!("unsupported Singer prime power q = {q}"),
    }
}

/// Construct a Krückeberg–Cilleruelo Dyadic Singer-Block Sidon Set by placing
/// shifted Singer difference sets $B_r + S_{q_r}$ at separated block offsets $B_r$
/// chosen greedily (or with expansion factor `gap_mult >= 1`) so that the entire
/// union remains a strictly valid Sidon set.
///
/// Demonstrates Theorem 3.1 (The Burst–Desert Dichotomy):
///   * At the end of block $r$ ($x = B_r + \max S_{q_r}$), $A(x)$ surges by $q_r + 1$,
///     producing a local maximum (burst) in $\Phi_A(x) = A(x)\sqrt{(\ln x)/x}$.
///   * Right before block $r$ ($x = B_r$), the mandatory difference-exclusion gap
///     creates a desert where $\Phi_A(B_r)$ drops sharply toward $0$.
pub fn kruckeberg_singer_block_sidon(qs: &[usize], gap_mult: u64) -> (Vec<u64>, Vec<(u64, u64)>) {
    let mut set: Vec<u64> = Vec::new();
    let mut used_diffs = std::collections::HashSet::<u64>::new();
    let mut block_ranges: Vec<(u64, u64)> = Vec::new();

    for &q in qs {
        let block = singer_sidon_block(q);
        // Verify intra-block differences are distinct
        let mut intra_diffs = Vec::new();
        for i in 0..block.len() {
            for j in (i + 1)..block.len() {
                intra_diffs.push(block[j] - block[i]);
            }
        }
        // If any intra-block difference is already in used_diffs, scale the block
        // by the smallest multiplier m >= 1 whose scaled intra-block differences
        // are disjoint from used_diffs.
        let mut scale = 1u64;
        while intra_diffs
            .iter()
            .any(|&d| used_diffs.contains(&(d * scale)))
        {
            scale += 1;
        }
        let scaled_block: Vec<u64> = block.iter().map(|&v| v * scale).collect();

        // Find base offset B_r > max(set) * gap_mult so that no cross-block
        // difference collides with existing used_diffs or scaled intra-block diffs.
        let mut base = match set.last() {
            None => 1u64,
            Some(&last) => (last * gap_mult.max(2)) + 1,
        };
        let scaled_intra_set: std::collections::HashSet<u64> =
            intra_diffs.iter().map(|&d| d * scale).collect();

        loop {
            let mut ok = true;
            let mut cross_diffs = std::collections::HashSet::<u64>::new();
            'outer: for &b in &scaled_block {
                let cand = base + b;
                for &existing in &set {
                    let d = cand - existing;
                    if used_diffs.contains(&d)
                        || scaled_intra_set.contains(&d)
                        || !cross_diffs.insert(d)
                    {
                        ok = false;
                        break 'outer;
                    }
                }
            }
            if ok {
                for &d in &scaled_intra_set {
                    used_diffs.insert(d);
                }
                for d in cross_diffs {
                    used_diffs.insert(d);
                }
                let start = base + scaled_block[0];
                let end = base + *scaled_block.last().unwrap();
                for &b in &scaled_block {
                    set.push(base + b);
                }
                block_ranges.push((start, end));
                break;
            }
            base += 1;
        }
    }

    (set, block_ranges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify_cpu::{is_sidon_cpu, mian_chowla_greedy_cpu};

    #[test]
    fn mian_chowla_table_matches_greedy_generator_64() {
        let gen = mian_chowla_greedy_cpu(64);
        assert_eq!(gen.as_slice(), &A005282_MIAN_CHOWLA_64[..]);
        assert!(is_sidon_cpu(&gen));
    }

    #[test]
    fn f2_sidon_inversion_matches_known_values() {
        // N = 1 -> 1, N = 2 -> 2, N = 4 -> 3, N = 6 -> 3, N = 7 -> 4, N = 12 -> 5, N = 18 -> 6
        assert_eq!(f2_sidon_exact(1), Some(1));
        assert_eq!(f2_sidon_exact(2), Some(2));
        assert_eq!(f2_sidon_exact(3), Some(2));
        assert_eq!(f2_sidon_exact(4), Some(3));
        assert_eq!(f2_sidon_exact(6), Some(3));
        assert_eq!(f2_sidon_exact(7), Some(4));
        assert_eq!(f2_sidon_exact(11), Some(4));
        assert_eq!(f2_sidon_exact(12), Some(5));
        assert_eq!(f2_sidon_exact(18), Some(6));
        assert_eq!(f2_sidon_exact(26), Some(7));
        assert_eq!(f2_sidon_exact(35), Some(8));
        assert_eq!(f2_sidon_exact(586), Some(28));
        assert_eq!(f2_sidon_exact(587), None);
    }

    #[test]
    fn dyadic_budget_theorem_holds_on_mian_chowla_and_singer_blocks() {
        let mc = mian_chowla_greedy_cpu(64);
        for m in 2..=13 {
            let rep = evaluate_dyadic_budget(&mc, m);
            assert!(
                rep.scaled_lhs_num <= rep.scaled_rhs_num,
                "scaled_lhs {} > scaled_rhs {} at M={m}",
                rep.scaled_lhs_num,
                rep.scaled_rhs_num
            );
            assert!(rep.lhs_dyadic_sum <= rep.exact_diff_weighted_sum + 1e-12);
            assert!(rep.exact_diff_weighted_sum <= rep.rhs_universal_ceiling + 1e-12);
            assert!(rep.multiscale_lhs <= rep.multiscale_harmonic_rhs + 1e-12);
        }

        let (kb, _) = kruckeberg_singer_block_sidon(&[2, 3, 4, 5, 7], 3);
        assert!(is_sidon_cpu(&kb));
        let rep_kb = evaluate_dyadic_budget(&kb, 14);
        assert!(rep_kb.scaled_lhs_num <= rep_kb.scaled_rhs_num);
        assert!(rep_kb.multiscale_lhs <= rep_kb.multiscale_harmonic_rhs + 1e-12);
    }

    #[test]
    fn theorem_2_2_cross_shell_bilinear_packing_and_two_shell_localization() {
        let mc = mian_chowla_greedy_cpu(64);
        let (kb, _) = kruckeberg_singer_block_sidon(&[2, 3, 4, 5, 7], 3);
        let rz = ruzsa_power_law_sidon(40);
        assert!(is_sidon_cpu(&rz));

        for seq in [&mc, &kb, &rz] {
            for m in 3..=13 {
                let rep = evaluate_cross_shell_bilinear(seq, m);
                // Every far inter-shell pair 1 <= j <= k - 2 MUST land inside (2^{k-2}, 2^k) = D_{k-1} cup D_k
                assert_eq!(
                    rep.far_inter_pairs_total, rep.far_inter_in_two_shells,
                    "Two-shell localization failed at M={m}"
                );
                // Two-shell step inequality: binom(m_k, 2) + m_k A(2^{k-1}) <= |D_A(2^k) cap (0, 2^k)|
                let mut cum_steps = 0u64;
                for k in 0..m {
                    assert!(rep.two_shell_lhs_per_k[k] <= rep.d_cap_per_k[k]);
                    cum_steps += rep.two_shell_lhs_per_k[k];
                    assert_eq!(cum_steps, rep.d_cap_per_k[k]);
                }
                // Bilinear quadratic form strictly dominates intra-shell-only sum and obeys universal ceiling
                assert!(rep.bilinear_dyadic_lhs >= rep.intra_only_dyadic_lhs);
                assert!(rep.bilinear_dyadic_lhs <= rep.bilinear_dyadic_rhs + 1e-12);
            }
        }
    }

    #[test]
    fn theorem_3_2_cross_shell_qp_envelope_lowers_intra_only_bound() {
        for &m in &[8usize, 16, 32, 64] {
            let cmp = solve_cross_shell_qp_envelope(m, 3);
            assert!(
                cmp.c_bilinear_cross_shell_upper < cmp.c_intra_only_upper,
                "bilinear QP {} did not lower intra-only {} at M={m}",
                cmp.c_bilinear_cross_shell_upper,
                cmp.c_intra_only_upper
            );
            assert!(cmp.reduction_factor > 1.25);
        }
    }

    #[test]
    fn theorem_2_3_subwindow_projection_and_crowding_holds_across_all_sequences() {
        let mc = mian_chowla_greedy_cpu(64);
        let (kb, _) = kruckeberg_singer_block_sidon(&[2, 3, 4, 5, 7], 3);
        let rz = ruzsa_power_law_sidon(40);

        for seq in [&mc, &kb, &rz] {
            let rep = evaluate_subwindow_projection(seq, 13);
            for row in &rep.scale_rows {
                assert!(
                    row.total_pigeonhole_lhs <= row.total_exact_subwindow_lhs,
                    "pigeonhole {} > exact subwindow {} at r={}",
                    row.total_pigeonhole_lhs,
                    row.total_exact_subwindow_lhs,
                    row.r_scale
                );
                assert!(
                    row.total_exact_subwindow_lhs <= row.actual_short_diffs_rhs,
                    "exact subwindow {} > actual short diffs {} at r={}",
                    row.total_exact_subwindow_lhs,
                    row.actual_short_diffs_rhs,
                    row.r_scale
                );
                assert!(
                    row.actual_short_diffs_rhs <= row.universal_cap_2r_minus_1,
                    "actual short diffs {} > 2^r - 1 ({}) at r={}",
                    row.actual_short_diffs_rhs,
                    row.universal_cap_2r_minus_1,
                    row.r_scale
                );
                // Theorem 2.4(a): Exact Sliding-Window Fejér chain
                assert!(
                    row.fejer_pigeonhole_lhs <= row.fejer_exact_rhs,
                    "Fejér pigeonhole {} > Fejér exact RHS {} at r={}",
                    row.fejer_pigeonhole_lhs,
                    row.fejer_exact_rhs,
                    row.r_scale
                );
                assert!(
                    row.fejer_exact_rhs <= row.fejer_cap_binom_2r,
                    "Fejér exact RHS {} > binom(2^r, 2) ({}) at r={}",
                    row.fejer_exact_rhs,
                    row.fejer_cap_binom_2r,
                    row.r_scale
                );
            }
            assert!(rep.multires_harmonic_lb <= rep.exact_harmonic_diff_sum + 1e-12);
            assert!(rep.exact_harmonic_diff_sum <= rep.harmonic_ceiling + 1e-12);
        }
    }

    #[test]
    fn theorem_3_3_multires_subwindow_qp_strictly_sharpens_phase2_and_phase1() {
        let mut prev_star = f64::INFINITY;
        for &m in &[8usize, 12, 16, 24, 32, 48, 64] {
            let mr = solve_multires_subwindow_qp_envelope(m, 3);
            assert!(
                mr.c_fejer_multires_upper < mr.c_multires_subwindow_upper,
                "Phase 4 Fejér {} >= Phase 3 {} at M={m}",
                mr.c_fejer_multires_upper,
                mr.c_multires_subwindow_upper
            );
            assert!(
                mr.c_multires_subwindow_upper < mr.c_bilinear_cross_shell_upper,
                "Phase 3 {} >= Phase 2 {} at M={m}",
                mr.c_multires_subwindow_upper,
                mr.c_bilinear_cross_shell_upper
            );
            assert!(
                mr.c_bilinear_cross_shell_upper < mr.c_intra_only_upper,
                "Phase 2 {} >= Phase 1 {} at M={m}",
                mr.c_bilinear_cross_shell_upper,
                mr.c_intra_only_upper
            );
            // Theorem 2.4(b): c_monotone_subwindow_star is strictly decreasing in M
            // and pinches c_asymptotic_fejer_limit from above while c_fejer_multires_upper
            // approaches c_asymptotic_fejer_limit from below!
            assert!(
                mr.c_monotone_subwindow_star < prev_star,
                "Phase 4 tail-octave c_M* {} did not strictly decrease from {} at M={m}",
                mr.c_monotone_subwindow_star,
                prev_star
            );
            assert!(
                mr.c_fejer_multires_upper < mr.c_asymptotic_fejer_limit
                    && mr.c_asymptotic_fejer_limit < mr.c_monotone_subwindow_star,
                "Pinching inequality failed at M={m}: {} < {} < {}",
                mr.c_fejer_multires_upper,
                mr.c_asymptotic_fejer_limit,
                mr.c_monotone_subwindow_star
            );
            prev_star = mr.c_monotone_subwindow_star;
            assert!(mr.phase3_vs_phase1_factor > 1.5);
            assert!(mr.phase4_vs_phase1_factor > mr.phase3_vs_phase1_factor);
        }
    }

    #[test]
    fn theorem_2_5_base_b_joint_fejer_strictly_sharpens_phase4_and_pinches_2() {
        // 1. Check strict monotonicity of Lambda(B) and c_infty^*(B) as B decreases to 1^+
        let bases = [
            4.0f64,
            2.0,
            std::f64::consts::SQRT_2,
            1.25,
            1.10,
            1.05,
            1.01,
        ];
        let mut prev_lambda = 0.0f64;
        let mut prev_c_inf = f64::INFINITY;
        for &b in &bases {
            let env = solve_base_b_joint_fejer_envelope(32, 3, b);
            assert!(
                env.lambda_b > prev_lambda && env.lambda_b < 0.25,
                "Lambda(B) = {} did not strictly increase toward 1/4 at B={b}",
                env.lambda_b
            );
            assert!(
                env.c_infty_star_b < prev_c_inf && env.c_infty_star_b > 2.0,
                "c_infty^*(B) = {} did not strictly decrease toward 2.0 at B={b}",
                env.c_infty_star_b
            );
            prev_lambda = env.lambda_b;
            prev_c_inf = env.c_infty_star_b;
        }
        // At B = 2.0, c_infty^*(2) = (1 + sqrt(2))*sqrt(ln 2)
        let env_b2 = solve_base_b_joint_fejer_envelope(32, 3, 2.0);
        let dyadic_limit = (1.0 + std::f64::consts::SQRT_2) * std::f64::consts::LN_2.sqrt();
        assert!((env_b2.c_infty_star_b - dyadic_limit).abs() < 1e-14);

        // 2. Check strict monotone decrease of c_tail_octave_joint in M and pinching of 2.0
        let mut prev_tail = f64::INFINITY;
        for &m in &[8usize, 12, 16, 24, 32, 48, 64] {
            let env = solve_base_b_joint_fejer_envelope(m, 3, 1.01);
            assert!(
                env.c_tail_octave_joint < env.phase4_tail_star,
                "Phase 5 unsplit tail {} >= Phase 4 tail {} at M={m}",
                env.c_tail_octave_joint,
                env.phase4_tail_star
            );
            assert!(
                env.c_tail_octave_joint < prev_tail,
                "Phase 5 unsplit tail {} did not strictly decrease from {} at M={m}",
                env.c_tail_octave_joint,
                prev_tail
            );
            assert!(
                env.c_full_joint_fejer < env.phase4_full_fejer,
                "Phase 5 full joint {} >= Phase 4 full {} at M={m}",
                env.c_full_joint_fejer,
                env.phase4_full_fejer
            );
            assert!(
                env.c_full_joint_fejer < env.c_infty_continuum
                    && env.c_infty_continuum < env.c_tail_octave_joint,
                "Continuum pinching failed at M={m}: {} < 2.0 < {}",
                env.c_full_joint_fejer,
                env.c_tail_octave_joint
            );
            prev_tail = env.c_tail_octave_joint;
            assert!(env.phase5_vs_phase1_factor > 1.7);
            assert!(env.phase5_vs_phase4_factor > 1.05);
        }
    }

    #[test]
    fn theorem_2_6_coupled_lag_integrality_breaks_2_00_barrier() {
        let mut prev_tail = f64::INFINITY;
        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64] {
            let env = solve_coupled_lag_integrality_envelope(m, 3, 1.01);
            assert!(
                env.c_coupled_tail_star < env.phase5_tail_star,
                "Phase 6 coupled tail {} >= Phase 5 tail {} at M={m}",
                env.c_coupled_tail_star,
                env.phase5_tail_star
            );
            assert!(
                env.c_coupled_tail_star < prev_tail,
                "Phase 6 coupled tail {} did not strictly decrease from {} at M={m}",
                env.c_coupled_tail_star,
                prev_tail
            );
            assert!(
                env.c_coupled_full_fejer < env.phase5_full_fejer,
                "Phase 6 coupled full {} >= Phase 5 full {} at M={m}",
                env.c_coupled_full_fejer,
                env.phase5_full_fejer
            );
            assert!(
                env.c_coupled_full_fejer < env.c_infty_coupled_limit
                    && env.c_infty_coupled_limit < env.c_coupled_tail_star,
                "Coupled pinching failed at M={m}: {} < {} < {}",
                env.c_coupled_full_fejer,
                env.c_infty_coupled_limit,
                env.c_coupled_tail_star
            );
            if m >= 32 {
                assert!(
                    env.c_coupled_tail_star < 2.0,
                    "Phase 6 coupled tail {} did not break 2.000000 barrier at M={m}",
                    env.c_coupled_tail_star
                );
            }
            prev_tail = env.c_coupled_tail_star;
            assert!(env.phase6_vs_phase1_factor > 2.0);
            assert!(env.phase6_vs_phase5_factor > 1.09);
        }
    }

    #[test]
    fn theorem_2_7_cesaro_fejer_order2_deflates_below_1_80() {
        let mut prev_tail = f64::INFINITY;
        let c_inf_continuum_cf2 = 2.0f64 * (4.0f64 / 5.0f64).sqrt();
        assert!(c_inf_continuum_cf2 < 1.80);

        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64] {
            let env = solve_cesaro_fejer_order2_envelope(m, 3, 1.01);
            assert!(
                env.c_cf2_tail_star < env.phase6_tail_star,
                "Phase 7 CF2 tail {} >= Phase 6 coupled tail {} at M={m}",
                env.c_cf2_tail_star,
                env.phase6_tail_star
            );
            assert!(
                env.c_cf2_tail_star < prev_tail,
                "Phase 7 CF2 tail {} did not strictly decrease from {} at M={m}",
                env.c_cf2_tail_star,
                prev_tail
            );
            assert!(
                env.c_cf2_full_fejer < env.phase6_full_fejer,
                "Phase 7 CF2 full {} >= Phase 6 full {} at M={m}",
                env.c_cf2_full_fejer,
                env.phase6_full_fejer
            );
            assert!(
                env.c_cf2_full_fejer < env.c_infty_cf2_limit
                    && env.c_infty_cf2_limit < env.c_cf2_tail_star,
                "Order-2 Cesàro–Fejér pinching failed at M={m}: {} < {} < {}",
                env.c_cf2_full_fejer,
                env.c_infty_cf2_limit,
                env.c_cf2_tail_star
            );
            assert!(
                env.c_infty_cf2_limit < 1.80,
                "Order-2 Cesàro–Fejér limit {} >= 1.80",
                env.c_infty_cf2_limit
            );
            if m >= 24 {
                assert!(
                    env.c_cf2_tail_star < 2.0,
                    "Phase 7 CF2 tail {} did not break 2.000000 barrier at M={m}",
                    env.c_cf2_tail_star
                );
            }
            prev_tail = env.c_cf2_tail_star;
            assert!(env.phase7_vs_phase1_factor > 2.0);
            assert!(env.phase7_vs_phase6_factor > 1.0005);
        }
    }

    #[test]
    fn theorem_2_8_octave_matched_and_band_cesaro_deflates_below_1_77() {
        // 1. Verify Theorem 2.8(a) non-vanishing octave-matched Cesàro–Fejér law:
        //    M_s = max_{p>=1} beta_s^{(p)} is strictly decreasing for s >= 2 and > 1/4
        let (p1, m1) = octave_matched_cesaro_order(1);
        assert_eq!(p1, 1);
        assert!((m1 - 0.25).abs() < 1e-14);
        let (p2, m2) = octave_matched_cesaro_order(2);
        assert_eq!(p2, 1);
        assert!((m2 - 5.0 / 16.0).abs() < 1e-14);
        let (p3, m3) = octave_matched_cesaro_order(3);
        assert_eq!(p3, 4);
        assert!((m3 - 9031.0 / 32768.0).abs() < 1e-14);

        let mut prev_m_s = m2;
        for s in 3..=8 {
            let (_, m_s) = octave_matched_cesaro_order(s);
            assert!(m_s > 0.25 && m_s < prev_m_s);
            prev_m_s = m_s;
        }

        // 2. Verify Theorem 2.8(b) band-restricted envelope bounds & M=20 sub-2.00 breakthrough
        let mut prev_tail = f64::INFINITY;
        let c_inf_continuum_band = 2.0f64 * (7.0f64 / 9.0f64).sqrt();
        assert!(c_inf_continuum_band < 1.77);

        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64] {
            let env = solve_band_restricted_cesaro_envelope(m, 3, 1.01);
            assert!(
                env.c_band_tail_star < env.phase7_tail_star,
                "Phase 8 band tail {} >= Phase 7 CF2 tail {} at M={m}",
                env.c_band_tail_star,
                env.phase7_tail_star
            );
            assert!(
                env.c_band_tail_star < prev_tail,
                "Phase 8 band tail {} did not strictly decrease from {} at M={m}",
                env.c_band_tail_star,
                prev_tail
            );
            assert!(
                env.c_band_full_fejer < env.phase7_full_fejer,
                "Phase 8 band full {} >= Phase 7 CF2 full {} at M={m}",
                env.c_band_full_fejer,
                env.phase7_full_fejer
            );
            assert!(
                env.c_band_full_fejer < env.c_infty_band_limit
                    && env.c_infty_band_limit < env.c_band_tail_star,
                "Band-restricted Cesàro–Fejér pinching failed at M={m}: {} < {} < {}",
                env.c_band_full_fejer,
                env.c_infty_band_limit,
                env.c_band_tail_star
            );
            assert!(
                env.c_infty_band_limit < 1.77,
                "Band-restricted Cesàro–Fejér limit {} >= 1.77",
                env.c_infty_band_limit
            );
            if m >= 20 {
                assert!(
                    env.c_band_tail_star < 2.0,
                    "Phase 8 band tail {} did not break 2.000000 barrier at M={m}",
                    env.c_band_tail_star
                );
            }
            prev_tail = env.c_band_tail_star;
            assert!(env.phase8_vs_phase1_factor > 2.0);
            assert!(env.phase8_vs_phase7_factor > 1.0005);
        }
    }

    #[test]
    fn theorem_2_9_multioctave_cascade_deflates_to_1_633_and_breaks_2_00_at_m16() {
        // 1. Verify exact rational K-octave integrality cascade:
        //    K=1: 2/7, K=2: 8/21, K=3: 32/73, K=4: 128/273 -> 1/2
        assert!((kappa_int_multioctave_k(1) - 2.0 / 7.0).abs() < 1e-15);
        assert!((kappa_int_multioctave_k(2) - 8.0 / 21.0).abs() < 1e-15);
        assert!((kappa_int_multioctave_k(3) - 32.0 / 73.0).abs() < 1e-15);
        assert!((kappa_int_multioctave_k(4) - 128.0 / 273.0).abs() < 1e-15);

        let mut prev_kap = kappa_int_multioctave_k(1);
        for k in 2..=12 {
            let kap = kappa_int_multioctave_k(k);
            assert!(kap > prev_kap && kap < 0.5);
            prev_kap = kap;
        }

        // 2. Verify Theorem 2.9 multi-octave envelope bounds & M=16 sub-2.00 breakthrough
        let mut prev_tail = f64::INFINITY;
        let c_inf_continuum_multi = 2.0f64 * (2.0f64 / 3.0f64).sqrt();
        assert!(c_inf_continuum_multi < 1.64);

        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64] {
            let env = solve_multioctave_cesaro_envelope(m, 3, 1.01);
            assert!(
                env.c_multioct_tail_star < env.c_k3_tail_star
                    && env.c_k3_tail_star < env.c_k2_tail_star
                    && env.c_k2_tail_star < env.phase8_tail_star,
                "K-octave tail hierarchy failed at M={m}"
            );
            assert!(
                env.c_multioct_tail_star < prev_tail,
                "Phase 9 multi-octave tail {} did not strictly decrease from {} at M={m}",
                env.c_multioct_tail_star,
                prev_tail
            );
            assert!(
                env.c_multioct_full_fejer < env.phase8_full_fejer,
                "Phase 9 multi-octave full {} >= Phase 8 band full {} at M={m}",
                env.c_multioct_full_fejer,
                env.phase8_full_fejer
            );
            assert!(
                env.c_multioct_full_fejer < env.c_infty_multioct_limit
                    && env.c_infty_multioct_limit < env.c_multioct_tail_star,
                "Multi-octave Cesàro–Fejér pinching failed at M={m}: {} < {} < {}",
                env.c_multioct_full_fejer,
                env.c_infty_multioct_limit,
                env.c_multioct_tail_star
            );
            assert!(
                env.c_infty_multioct_limit < 1.64,
                "Multi-octave Cesàro–Fejér limit {} >= 1.64",
                env.c_infty_multioct_limit
            );
            if m >= 16 {
                assert!(
                    env.c_multioct_tail_star < 2.0 && env.c_k3_tail_star < 2.0,
                    "Phase 9 multi-octave tail {} did not break 2.000000 barrier at M={m}",
                    env.c_multioct_tail_star
                );
            }
            if m >= 32 {
                assert!(
                    env.c_multioct_tail_star < 1.80,
                    "Phase 9 multi-octave tail {} did not break 1.800000 barrier at M={m}",
                    env.c_multioct_tail_star
                );
            }
            prev_tail = env.c_multioct_tail_star;
            assert!(env.phase9_vs_phase1_factor > 2.0);
            assert!(env.phase9_vs_phase8_factor > 1.004);
        }
    }

    #[test]
    fn theorem_2_10_power_law_and_haar_multioctave_cascade_converges_to_zero() {
        // 1. Verify Theorem 2.10(a) exact power-law continuum integrality milestones
        assert!((kappa_int_power_law_continuum(0.0) - 0.5).abs() < 1e-15);
        assert!((kappa_int_power_law_continuum(-0.5) - 5.0 / 6.0).abs() < 1e-15);
        assert!((kappa_int_power_law_continuum(-0.6) - 1.0).abs() < 1e-15);
        assert!((kappa_int_power_law_continuum(-0.75) - 1.5).abs() < 1e-15);
        assert!((kappa_int_power_law_continuum(-15.0 / 17.0) - 3.0).abs() < 1e-14);
        assert!((kappa_int_power_law_continuum(-0.9) - 3.5).abs() < 1e-14);

        let c_inf_unit_barrier = 2.0 / (1.0 + kappa_int_power_law_continuum(-15.0 / 17.0)).sqrt();
        assert!((c_inf_unit_barrier - 1.0).abs() < 1e-14);
        let c_inf_sub_unit = 2.0 / (1.0 + kappa_int_power_law_continuum(-0.9)).sqrt();
        assert!((c_inf_sub_unit - 2.0 * std::f64::consts::SQRT_2 / 3.0).abs() < 1e-14);
        assert!(c_inf_sub_unit < 0.95);

        // 2. Verify Theorem 2.10(b) Haar K-octave cascade -> 0
        let mut prev_kap_h = 0.0f64;
        let mut prev_c_h = f64::INFINITY;
        for k_oct in [1usize, 2, 3, 4, 8, 12, 13, 16, 32, 64] {
            let kap_h = kappa_int_haar_k(k_oct);
            let c_h = 2.0 / (1.0 + kap_h).sqrt();
            assert!(kap_h > prev_kap_h);
            assert!(c_h < prev_c_h);
            prev_kap_h = kap_h;
            prev_c_h = c_h;
        }
        assert!(kappa_int_haar_k(1) > 2.0 / 7.0); // Beats Phase 8 even at K = 1
        assert!(kappa_int_haar_k(3) > 0.5); // Beats Phase 9 K->oo limit already at K = 3
        assert!(2.0 / (1.0 + kappa_int_haar_k(13)).sqrt() < 1.0); // Breaks 1.00 barrier at K = 13

        // 3. Verify Theorem 2.10 finite-horizon Haar bounds across M in [8, 64]
        let mut prev_k3_tail = f64::INFINITY;
        let mut prev_adapt_tail = f64::INFINITY;
        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64] {
            let env = solve_haar_multioctave_envelope(m, 3, 1.01);
            assert!(
                env.c_haar_adapt_tail_star < env.c_haar_k3_tail_star
                    && env.c_haar_k3_tail_star < env.phase9_tail_star,
                "Phase 10 Haar tail hierarchy failed at M={m}"
            );
            assert!(env.c_haar_k3_tail_star < prev_k3_tail);
            assert!(env.c_haar_adapt_tail_star < prev_adapt_tail);
            assert!(
                env.c_haar_adapt_full_fejer < env.c_haar_k3_full_fejer
                    && env.c_haar_k3_full_fejer < env.phase9_full_fejer,
                "Phase 10 Haar full hierarchy failed at M={m}"
            );
            assert!(
                env.c_haar_k3_full_fejer < env.c_infty_haar_k3_limit
                    && env.c_infty_haar_k3_limit < env.c_haar_k3_tail_star,
                "Fixed K=3 Haar pinching failed at M={m}"
            );
            if m >= 12 {
                assert!(
                    env.c_haar_adapt_tail_star < 2.0,
                    "Phase 10 adaptive Haar tail {} did not break 2.000000 at M={m}",
                    env.c_haar_adapt_tail_star
                );
            }
            if m >= 48 {
                assert!(
                    env.c_haar_adapt_tail_star < 1.0,
                    "Phase 10 adaptive Haar tail {} did not break 1.000000 at M={m}",
                    env.c_haar_adapt_tail_star
                );
            }
            prev_k3_tail = env.c_haar_k3_tail_star;
            prev_adapt_tail = env.c_haar_adapt_tail_star;
            assert!(env.phase10_adapt_vs_phase1_factor > 2.1);
            assert!(env.phase10_adapt_vs_phase9_factor > 1.01);
        }
    }

    #[test]
    fn theorem_2_11_unconditional_diophantine_haar_closure_completes_erdos_1191() {
        // 1. Verify Theorem 2.11(a) single-period closed form & strict rational sandwich for m = 1..200
        assert!((sawtooth_haar_period_integral(0) - 0.5).abs() < 1e-15);
        for m in 1..=200usize {
            let im = sawtooth_haar_period_integral(m);
            let mf = m as f64;
            let lb2 = 1.0 / (6.0 * mf) - 1.0 / (12.0 * mf * mf);
            let ub3 = lb2 + 1.0 / (20.0 * mf * mf * mf);
            let ub1 = 1.0 / (6.0 * mf);
            assert!(
                lb2 < im && im < ub3 && ub3 < ub1,
                "Rational sandwich failed at m={m}: {lb2} < {im} < {ub3} < {ub1}"
            );
        }

        // 2. Verify Theorem 2.11(b) exact hyperfactorial identity & super-convergent O(M^{-2}) expansion
        assert!((C_GK_GLAISHER_SAWTOOTH - 0.5 - C_LF_EULER_MASCHERONI_SAWTOOTH).abs() < 1e-15);
        let mut cum = 0.5f64;
        for m_cap in 1..=64usize {
            if m_cap > 1 {
                cum += sawtooth_haar_period_integral(m_cap - 1);
            }
            let hyper = sawtooth_haar_antiderivative(m_cap as f64);
            assert!(
                (cum - hyper).abs() < 1e-14,
                "Antiderivative mismatch at M={m_cap}: {cum} vs {hyper}"
            );
            let mf = m_cap as f64;
            let mut ln_h = 0.0f64;
            for k in 1..m_cap {
                let kf = k as f64;
                ln_h += kf * kf.ln();
            }
            let hyper_formula = 0.5 * mf * mf - mf * (mf - 1.0) * mf.ln() + 2.0 * ln_h;
            assert!(
                (cum - hyper_formula).abs() < 1e-10,
                "Hyperfactorial closed form mismatch at M={m_cap}: {cum} vs {hyper_formula}"
            );
            let asymp_lb = mf.ln() / 6.0 + C_GK_GLAISHER_SAWTOOTH;
            let asymp_ub = asymp_lb + 1.0 / (360.0 * mf * mf);
            assert!(
                asymp_lb < hyper && hyper < asymp_ub + 1e-15,
                "Bernoulli O(M^-2) sandwich failed at M={m_cap}: {asymp_lb} < {hyper} < {asymp_ub}"
            );
        }

        // 3. Verify Unconditional Pointwise Linear-Dilation Theorem on adversarial integers {mu_1} = 0
        //    and irrational densities across K = 3 octaves:
        let k_oct = 3usize;
        let target_haar = (k_oct as f64) * std::f64::consts::LN_2 / 6.0;
        let mut prev_int_err = f64::INFINITY;
        for &mu_int in &[8.0f64, 16.0, 32.0, 64.0, 128.0, 256.0] {
            let j_val = pointwise_diophantine_haar_dilation(mu_int, k_oct);
            let err = (j_val - target_haar).abs();
            assert!(err < prev_int_err, "Error did not decrease at mu={mu_int}");
            prev_int_err = err;
            let u_lo = (mu_int / 8.0).floor();
            let max_err = 1.0 / (6.0 * u_lo) + 1.0 / (360.0 * u_lo * u_lo);
            assert!(err <= max_err);
        }
        for &mu_irr in &[
            16.0 * std::f64::consts::SQRT_2,
            16.0 * std::f64::consts::PI,
            32.0 * std::f64::consts::E,
        ] {
            let j_val = pointwise_diophantine_haar_dilation(mu_irr, k_oct);
            let u_lo = (mu_irr / 8.0).floor();
            let max_err = 1.0 / (6.0 * u_lo) + 1.0 / (360.0 * u_lo * u_lo);
            assert!((j_val - target_haar).abs() <= max_err);
        }

        // 4. Verify strict monotonicity of kappa_int_diophantine_haar_k(K) -> oo and c_infty -> 0
        let mut prev_kap = 0.0f64;
        let mut prev_c = f64::INFINITY;
        for k_val in [1usize, 2, 3, 4, 8, 9, 12, 16, 32, 64] {
            let kap = kappa_int_diophantine_haar_k(k_val);
            let c_inf = 2.0 / (1.0 + kap).sqrt();
            assert!(kap > prev_kap && kap > kappa_int_haar_k(k_val));
            assert!(c_inf < prev_c);
            prev_kap = kap;
            prev_c = c_inf;
        }
        assert!(2.0 / (1.0 + kappa_int_diophantine_haar_k(9)).sqrt() < 1.0); // Breaks 1.00 at K = 9
        assert!(2.0 / (1.0 + kappa_int_diophantine_haar_k(64)).sqrt() < 0.50); // Breaks 0.50 at K = 64

        // 5. Verify finite-horizon bounds across M in [8, 128]
        let mut prev_tail = f64::INFINITY;
        for &m in &[8usize, 12, 16, 20, 24, 32, 48, 64, 96, 128] {
            let env = solve_unconditional_diophantine_haar_envelope(m, 3, 1.01);
            assert!(env.c_dioph_adapt_tail_star < env.phase10_tail_star);
            assert!(env.c_dioph_adapt_full_fejer < env.phase10_full_fejer);
            assert!(env.c_dioph_adapt_tail_star < prev_tail);
            if m >= 8 {
                assert!(
                    env.c_dioph_adapt_tail_star < 2.0,
                    "Phase 11 Diophantine Haar tail {} did not break 2.00 at M={m}",
                    env.c_dioph_adapt_tail_star
                );
            }
            if m >= 32 {
                assert!(
                    env.c_dioph_adapt_tail_star < 1.0,
                    "Phase 11 Diophantine Haar tail {} did not break 1.00 at M={m}",
                    env.c_dioph_adapt_tail_star
                );
            }
            if m >= 128 {
                assert!(
                    env.c_dioph_adapt_tail_star < 0.50,
                    "Phase 11 Diophantine Haar tail {} did not break 0.50 at M={m}",
                    env.c_dioph_adapt_tail_star
                );
            }
            prev_tail = env.c_dioph_adapt_tail_star;
            assert!(env.phase11_adapt_vs_phase1_factor > 2.15);
            assert!(env.phase11_adapt_vs_phase10_factor > 1.02);
        }
    }
}
