//! ============================================================================
//! `gpu.rs` — CUDA driver context and runtime NVRTC compilation for `erdos1191`.
//! ============================================================================
//!
//! Uses `cudarc` with `dynamic-loading`, so `libcuda.so` and `libnvrtc.so` are
//! loaded via `dlopen` only when a GPU path is invoked. CPU-only commands and
//! the default CI test suite run on hosts without an Nvidia GPU.

use anyhow::{Context, Result};
use cudarc::driver::{CudaContext, CudaFunction, CudaModule, CudaStream, LaunchConfig};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use std::sync::Arc;

pub const KERNEL_SRC: &str = include_str!("kernels.cu");

pub struct Gpu {
    pub ctx: Arc<CudaContext>,
    pub stream: Arc<CudaStream>,
    pub module: Arc<CudaModule>,
    pub fn_dfs_maximin_sidon: CudaFunction,
    pub fn_eval_sidon_dyadic_budget: CudaFunction,
    pub fn_profile_sidon_cross_shell: CudaFunction,
    pub fn_eval_sidon_subwindow_spectrum: CudaFunction,
    pub fn_eval_base_b_joint_fejer: CudaFunction,
    pub fn_eval_coupled_lag_integrality_fejer: CudaFunction,
    pub fn_eval_cesaro_fejer_order2: CudaFunction,
    pub fn_eval_band_restricted_cesaro_fejer: CudaFunction,
    pub fn_eval_multioctave_cesaro_fejer: CudaFunction,
    pub fn_eval_haar_multioctave_cesaro_fejer: CudaFunction,
    pub fn_eval_unconditional_diophantine_haar_closure: CudaFunction,
    pub device_name: String,
    pub sm_count: i32,
}

impl Gpu {
    pub fn init(device_ordinal: usize, arch: &str) -> Result<Self> {
        let ctx = CudaContext::new(device_ordinal)
            .with_context(|| format!("failed to initialize CUDA device {device_ordinal}"))?;
        let stream = ctx.default_stream();
        let device_name = ctx
            .name()
            .unwrap_or_else(|_| "Unknown CUDA GPU".to_string());
        let sm_count = ctx
            .attribute(
                cudarc::driver::sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,
            )
            .unwrap_or(82);

        let arch_flag = format!("--gpu-architecture={arch}");
        let opts = CompileOptions {
            options: vec![arch_flag, "--use_fast_math".to_string()],
            ..Default::default()
        };
        let ptx = compile_ptx_with_opts(KERNEL_SRC, opts)
            .context("NVRTC compilation of src/kernels.cu failed")?;
        let module = ctx
            .load_module(ptx)
            .context("failed to load compiled PTX module")?;

        let fn_dfs_maximin_sidon = module
            .load_function("k_dfs_maximin_sidon")
            .context("missing kernel k_dfs_maximin_sidon")?;
        let fn_eval_sidon_dyadic_budget = module
            .load_function("k_eval_sidon_dyadic_budget")
            .context("missing kernel k_eval_sidon_dyadic_budget")?;
        let fn_profile_sidon_cross_shell = module
            .load_function("k_profile_sidon_cross_shell")
            .context("missing kernel k_profile_sidon_cross_shell")?;
        let fn_eval_sidon_subwindow_spectrum = module
            .load_function("k_eval_sidon_subwindow_spectrum")
            .context("missing kernel k_eval_sidon_subwindow_spectrum")?;
        let fn_eval_base_b_joint_fejer = module
            .load_function("k_eval_base_b_joint_fejer")
            .context("missing kernel k_eval_base_b_joint_fejer")?;
        let fn_eval_coupled_lag_integrality_fejer = module
            .load_function("k_eval_coupled_lag_integrality_fejer")
            .context("missing kernel k_eval_coupled_lag_integrality_fejer")?;
        let fn_eval_cesaro_fejer_order2 = module
            .load_function("k_eval_cesaro_fejer_order2")
            .context("missing kernel k_eval_cesaro_fejer_order2")?;
        let fn_eval_band_restricted_cesaro_fejer = module
            .load_function("k_eval_band_restricted_cesaro_fejer")
            .context("missing kernel k_eval_band_restricted_cesaro_fejer")?;
        let fn_eval_multioctave_cesaro_fejer = module
            .load_function("k_eval_multioctave_cesaro_fejer")
            .context("missing kernel k_eval_multioctave_cesaro_fejer")?;
        let fn_eval_haar_multioctave_cesaro_fejer = module
            .load_function("k_eval_haar_multioctave_cesaro_fejer")
            .context("missing kernel k_eval_haar_multioctave_cesaro_fejer")?;
        let fn_eval_unconditional_diophantine_haar_closure = module
            .load_function("k_eval_unconditional_diophantine_haar_closure")
            .context("missing kernel k_eval_unconditional_diophantine_haar_closure")?;

        Ok(Self {
            ctx,
            stream,
            module,
            fn_dfs_maximin_sidon,
            fn_eval_sidon_dyadic_budget,
            fn_profile_sidon_cross_shell,
            fn_eval_sidon_subwindow_spectrum,
            fn_eval_base_b_joint_fejer,
            fn_eval_coupled_lag_integrality_fejer,
            fn_eval_cesaro_fejer_order2,
            fn_eval_band_restricted_cesaro_fejer,
            fn_eval_multioctave_cesaro_fejer,
            fn_eval_haar_multioctave_cesaro_fejer,
            fn_eval_unconditional_diophantine_haar_closure,
            device_name,
            sm_count,
        })
    }

    pub fn launch_cfg(&self, n_items: u32) -> LaunchConfig {
        let block_dim = 256u32;
        let grid_dim = n_items.div_ceil(block_dim).max(1);
        LaunchConfig {
            grid_dim: (grid_dim, 1, 1),
            block_dim: (block_dim, 1, 1),
            shared_mem_bytes: 0,
        }
    }
}
