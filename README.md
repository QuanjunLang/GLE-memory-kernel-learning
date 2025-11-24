# Learning Memory Kernels in Generalized Langevin Equations

This repository contains research code for **learning memory kernels** in generalized Langevin equations (GLEs) from trajectory data. It implements the full pipeline developed in

> Quanjun Lang and Jianfeng Lu, “Learning Memory Kernels in Generalized Langevin Equations,” SIAM Journal on Mathematics of Data Science, 2025.

including nonparametric kernel estimation via a regularized Prony method, Laplace–domain reconstruction, and regression with RKHS regularization.

---

## Features

The current MATLAB code supports:

- **Synthetic GLE data generation**
  - Random Prony–type kernels, power-law and polynomial kernels (via `memoryInfo.kernel_type`)
  - Linear and double-well potentials for the drift force
  - Generation of correlated noise and trajectories using spectral methods

- **Estimation of correlation functions**
  - Velocity autocorrelation and forcing autocorrelation from sparsely observed, possibly noisy trajectories
  - Multiple temporal/spatial averaging strategies (e.g. `temporal_spacial_method = 'M_only'`)

- **Prony-based interpolation**
  - Regularized Prony method to fit correlation functions by exponentials
  - Support for different weightings, polynomial coefficient solvers, and root normalization options
  - Simultaneous estimation of
    - `h(t)`: velocity autocorrelation  
    - `g(t)`: force autocorrelation  
    - derivatives where needed

- **Recovery of the memory kernel**
  - Laplace-domain reconstruction of the memory kernel `γ(t)` from Prony fits to `h` and `g`
  - Combined-loss regression in a weighted `L^2(ρ)` space using spline bases and RKHS regularization to obtain a smooth estimator `θ(t)`
  - Optional Tikhonov–Fourier deconvolution (FFT-based) as an alternative estimator

- **Diagnostics and coercivity analysis**
  - Scripts for computing `L^2(ρ)`-errors between estimated and true kernels
  - Coercivity checks and spectral diagnostics for the learned kernel
  - Publication-quality plots (linear and log scale) for the true vs. estimated memory kernels

---

## Requirements

- MATLAB R2021a or later (earlier versions may work but are not tested)
- Signal Processing Toolbox (for FFT-based routines; can be replaced by core `fft` if needed)

The code assumes a standard MATLAB path setup with all subfolders added. The helper function `addPathAll` recursively adds the project folders to the MATLAB path.

---

## Getting Started

### 1. Clone the repository

    git clone https://github.com/your-username/gle-memory-kernel-learning.git
    cd gle-memory-kernel-learning

### 2. Start MATLAB and add paths

In MATLAB:

    cd /path/to/gle-memory-kernel-learning
    addPathAll    % adds all subfolders of the project to the MATLAB path

### 3. Directory conventions

The main scripts assume that you define the following variables:

    fig_dir  = './figures/';   % directory to save figures
    data_dir = './data/';      % directory to save intermediate data

Create these folders if they do not exist:

    if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end
    if ~exist(data_dir, 'dir'); mkdir(data_dir); end

---

## Main Example Pipeline

A typical driver script (similar to the example in this repository) follows the steps below.

### 1. Generate a true memory kernel

    memoryInfo.kernel_type = 'RandomProny_p_3_q_1';  % or 'PowerLaw', 'Polynomial', etc.
    memoryInfo.plotON      = 0;
    memoryInfo = generate_memory_kernel(memoryInfo);

### 2. Set up the GLE system

    sysInfo.m      = 1;
    sysInfo.beta   = 1;
    sysInfo.mu     = 1;

    sysInfo.F_type = 'linear';   % or 'double_well'
    sysInfo.G_type = '0';

    sysInfo.fig_dir  = fig_dir;
    sysInfo.data_dir = data_dir;

    sysInfo = update_system_info(sysInfo, memoryInfo);  % generates true autocorrelation for linear F

### 3. Generate correlated trajectories

    trajInfo.loadON = 1;
    trajInfo.plotON = 0;
    trajInfo.wu     = 80*pi;      % spectral cutoff
    trajInfo.N      = 8000;       % number of spectral points
    trajInfo.L      = 2000;
    trajInfo.dt     = pi / trajInfo.wu;
    trajInfo.M      = 5000;       % number of trajectories

    random_seed = 20;
    rng(random_seed);

    trajInfo = generate_trajectories(memoryInfo, sysInfo, trajInfo, random_seed);

    sysInfo.dt    = trajInfo.dt;
    sysInfo.tgrid = trajInfo.tgrid;

### 4. Generate sparse/noisy observations

    obsInfo.plotON = 1;
    obsInfo.gap    = 10;
    obsInfo.std    = 0;
    obsInfo.len    = floor(trajInfo.L/obsInfo.gap);

    obsInfo = generate_observed_trajectory(obsInfo, trajInfo);

### 5. Estimate correlation functions

    sysInfo.omega = 2.5;
    sysInfo.rho   = @(x) exp(-2*sysInfo.omega*x);

    temporal_spacial_method = 'M_only';

    [corr_func, corr_func_obs] = get_correlations_true_obs( ...
        sysInfo, trajInfo, obsInfo, temporal_spacial_method);

### 6. Prony interpolation of correlations

    VACF_obs_len = 50;

    % Velocity autocorrelation
    PronyInfo_h.prony_p             = 20;
    PronyInfo_h.prony_N             = VACF_obs_len;
    PronyInfo_h.obs_dx              = obsInfo.dt_obs;
    PronyInfo_h.obs_h_grid          = corr_func_obs.h(1:VACF_obs_len);
    PronyInfo_h.rho                 = sysInfo.rho;
    PronyInfo_h.tgrid_obs           = obsInfo.tgrid_obs;
    PronyInfo_h.polycoef_method     = 'MP';
    PronyInfo_h.weight_method       = 'LS';
    PronyInfo_h.root_normalization  = 1;
    PronyInfo_h.lambda_augmentation = 1;
    PronyInfo_h.drop_0              = 0;

    % Forcing autocorrelation
    PronyInfo_g              = PronyInfo_h;
    PronyInfo_g.weight_method = 'LS';
    PronyInfo_g.obs_h_grid    = corr_func_obs.f(1:VACF_obs_len);
    PronyInfo_g.prony_p       = PronyInfo_h.prony_p;
    PronyInfo_g.drop_0        = 0;

    [EstResult, PronyDetails] = estimate_h_g_prony( ...
        PronyInfo_h, PronyInfo_g, sysInfo, corr_func, corr_func_obs);

    Prony_h  = EstResult.h;
    Prony_dg = EstResult.dg;
    Prony_g  = EstResult.g;

### 7. Recover the memory kernel via Laplace transform

    switch sysInfo.mu
        case 0
            result_gamma = get_gamma_from_prony_h(PronyInfo_h, PronyDetails.result_h);
        otherwise
            result_gamma = get_gamma_from_prony_h_esternal_force( ...
                PronyInfo_h, PronyInfo_g, PronyDetails.result_h, PronyDetails.result_phi);
    end

    theta_Prony = @(x) real(result_gamma.gamma(x));

### 8. Combined-loss regression in L^2(ρ)

    regression_para.lb        = 0;
    regression_para.rb        = 30;
    regression_para.knot_num  = 50;
    regression_para.deg       = 3;
    regression_para.free_bdry = 1;
    regression_para.N         = 2001;
    regression_para.omega     = sysInfo.omega;
    regression_para.reg_method        = 'RKHS';
    regression_para.convolution_method = 'direct';
    regression_para.basis_type         = 'spline';
    regression_para.lam                = result_gamma.lam_seq;

    h_est  = @(x) real(Prony_h(x));
    dg_est = @(x) real(Prony_dg(x));
    g_est  = @(x) real(Prony_g(x));

    [theta_1, result_details_regression_1] = ...
        get_theta_regression_combined(regression_para, h_est, dg_est, g_est, 0);

    [theta_2, result_details_regression_2] = ...
        get_theta_regression_combined(regression_para, h_est, dg_est, g_est, 1);

    [theta_0, result_details_regression_0] = ...
        get_theta_regression_combined(regression_para, h_est, dg_est, g_est, EstResult.comb_coef);

### 9. Optional: Tikhonov–Fourier deconvolution

The example script also includes an FFT-based estimator `theta_F` built from `h_est` and `g_est` via Tikhonov-regularized deconvolution. This provides an additional baseline to compare against the combined-loss regression.

### 10. Plotting and error metrics

The scripts produce

- comparison plots of the true kernel vs. estimated kernels (`θ`, `θ_1`, `θ_2`, `θ_L`, etc.),
- log–log plots of the kernel,
- `L^2(ρ)` errors such as

    gamma_regression_0_error = sqrt(sum((theta_0(tgrid') - R(tgrid')).^2 .* rho(tgrid')) * dt);

---

## Reproducing Figures From the Paper

To reproduce the main experiments and figures:

1. Choose `memoryInfo.kernel_type` (e.g. random Prony kernel, power-law kernel).
2. Set the trajectory parameters (`trajInfo`) and observation parameters (`obsInfo`) to match the desired regime.
3. Run the corresponding driver script in the `scripts/` or `examples/` directory.
4. Figures (e.g. `kernel_est.pdf`) are saved to `fig_dir` as PDF files suitable for inclusion in LaTeX.

---

## Citation

If you use this code in your work, please cite:

    @article{lang2025learningGLE,
      author  = {Lang, Quanjun and Lu, Jianfeng},
      title   = {Learning Memory Kernels in Generalized Langevin Equations},
      journal = {SIAM Journal on Mathematics of Data Science},
      year    = {2025},
    }

---

## License

This project is released under the MIT License. See the `LICENSE` file for details.
