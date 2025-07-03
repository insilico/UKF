# Unscented Kalman Filter (UKF) for Parameter Estimation in Nonlinear Dynamical Systems

This R package implements the Unscented Kalman Filter (UKF) for parameter estimation and state tracking in nonlinear dynamical systems, with a focus on applications such as fMRI time series analysis and coupled oscillator models.

---

## Overview

The Unscented Kalman Filter (UKF) is a recursive Bayesian filter designed for nonlinear systems. Unlike the Extended Kalman Filter (EKF), which linearizes the system using Jacobians, the UKF uses a deterministic sampling technique (the unscented transform) to more accurately capture the mean and covariance of the system state under nonlinear transformations.

This package provides tools for:
- State and parameter estimation in nonlinear ODE models
- Iterative optimization of model parameters using UKF
- Visualization and smoothing of time series data

---

## Mathematical Background

Given a nonlinear dynamical system:

$$
\begin{aligned}
\mathbf{x}_{k+1} &= f(\mathbf{x}_k, \mathbf{u}_k) + \mathbf{w}_k \\
\mathbf{y}_k &= h(\mathbf{x}_k) + \mathbf{v}_k
\end{aligned}
$$

where:

- $\mathbf{x}_k$ is the state vector at time $k$
- $\mathbf{u}_k$ is the control input (if any)
- $\mathbf{y}_k$ is the observation vector
- $f(\cdot)$ and $h(\cdot)$ are nonlinear functions
- $\mathbf{w}_k \sim \mathcal{N}(0, Q)$ is process noise
- $\mathbf{v}_k \sim \mathcal{N}(0, R)$ is measurement noise

In this package, the state vector is augmented to include both system states and unknown model parameters, allowing for simultaneous parameter and state estimation.

---

## Installation

```r
library(devtools)
install_github("insilico/UKF")
library(UKF)
```

---

## Usage Example

See `inst/model_voxels_2param.R` for a full example.  
Here is a minimal workflow:

```r
# Prepare your data: first column is time, others are observed variables
data <- read.csv("data/your_timeseries.csv")

# Define your ODE model as a function
my_ode_model <- function(t, y, p) {
  # y: state vector, p: parameter vector
  # return dydt as a vector
}

# Set initial parameter guess and UKF settings
param_guess <- c(1, 1) # example for two parameters
N_p <- 2
N_y <- ncol(data) - 1
dt <- 0.1
dT <- 1.0

# Run UKF parameter estimation
result <- iterative_param_optim(
  param_guess, t_dummy = 0, ts_data = data, ode_model = my_ode_model,
  N_p = N_p, N_y = N_y, dt = dt, dT = dT,
  param_tol = 1e-4, MAXSTEPS = 100,
  R_scale = 0.3, Q_scale = 0.01, trace = TRUE, forcePositive = TRUE
)

# Access estimated parameters and diagnostics
result$param_est
result$chisq
plot_ukf_chi_square_loss(result)
```

---

## Tuning Q and R

- `Q_scale`: Controls process noise (model/parameter uncertainty). Higher values allow parameters to adapt more quickly but may introduce noise.
- `R_scale`: Controls measurement noise. Higher values make the filter trust the model more than the data.

**Tip:** Start with `Q_scale = 0.01` and `R_scale = 0.3`, then adjust based on the stability and fit of your results.

---

## Visualization

The package includes plotting functions for:
- Raw and smoothed time series (`plot_voxels_and_smooth_ggplot`)
- UKF state and parameter estimates (`plot_ukf_and_smoothed_ggplot`)
- Chi-square loss over iterations (`plot_ukf_chi_square_loss`)

---

## Dependencies

Some functions require:
```r
install.packages(c('KernSmooth', 'ggplot2', 'pracma'))
```

---

## Academic Use

This package is developed for research in nonlinear time series analysis, parameter estimation, and model validation.  
If you use this package in your research, please cite appropriately and contact the maintainer for collaboration or questions.

---

## Contact

For questions, bug reports, or collaboration, contact:  
[brett.mckinney@gmail.com](mailto:brett.mckinney@gmail.com)