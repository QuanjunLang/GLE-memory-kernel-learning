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
  - Simultaneous estimation of:
    - \( h(t) \): velocity autocorrelation
    - \( g(t) \): force autocorrelation
    - Their derivatives

- **Recovery of the memory kernel**
  - **Laplace-domain reconstruction** of the memory kernel \(\gamma(t)\) from Prony fits to \(h\) and \(g\)
  - **Combined-loss regression** in a weighted \(L^2(\rho)\) space using spline bases and RKHS regularization to obtain a smooth estimator \(\theta(t)\)
  - Optional **Tikhonov–Fourier deconvolution** (FFT-based) as an alternative estimator using
    \[
    \widehat{\gamma}(\omega) = \frac{\overline{\widehat{h}(\omega)}\,\widehat{g}(\omega)}{|\widehat{h}(\omega)|^2 + \lambda}.
    \]

- **Diagnostics and coercivity analysis**
  - Scripts for computing \(L^2(\rho)\)-errors between estimated and true kernels
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

```bash
git clone https://github.com/QuanjunLang/gle-memory-kernel-learning.git
cd gle-memory-kernel-learning
