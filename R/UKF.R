###################################
# UKF.R Engine of library
# Unscented Kalman Filter Functions
###################################

#=====================================================#
#' UKF_dT
#' Unscented Kalman Filter (UKF) for time step dT
#' Doesn't use data at this step;
#' UKF_blend fn later blends UKF with data after dT step.
#' This function was known previously as
#' multiple_param_unscented_kalman_filter.
#' @param t_dummy a dummy time variable, because ode models don't have explcit time
#' @param ode_model model with variables y and N_p parameters
#' @param xhat UFK posterior prediction of augmented states (param + ind vars)
#' @param Pxx augmented covariance matrix for sigma points
#' @param y independent variable state vector
#' @param N_p number of unknown model params
#' @param N_y number of ind variables
#' @param R covariance matrix of y state. initialized in UKF_blend and updated with sigma points as Pyy in current function.
#' @param dt smaller time step size within dT
#' @param dT time step size that comes from time series data step
#' @param R_scale related to standard deviation of Gaussian measurement noise of ind variables. A number that must be specified by user.
#' @param Q_scale related to standard deviation of process noise. Noise related to model parameters
#' @param forcePositive logical, if TRUE, ensures all parameters stay positive
#' @return list: xhat (augmented Kalman update), Pxx (augmented covariance at sigma points), and K (Kalman gain matrix)
#' @examples
#' Example
#' @export
UKF_dT <- function(t_dummy, ode_model, xhat, Pxx, y, N_p, N_y, R, dt, dT,
                    R_scale, Q_scale, forcePositive = FALSE) {
  N_x <- N_p + N_y # number of augmented states
  N_sigma <- 2 * N_x # number sigma points
  xsigma <- t(chol(N_x * Pxx, pivot = T))
  # Unscented Sigma Points
  Xa <- xhat %*%
    matrix(rep(1, length = N_sigma), nrow = 1, ncol = N_sigma) +
    cbind(xsigma, -xsigma)

  # X <- kura_ode_model(t,N_p,Xa,dT)
  X <- propagate_model(t_dummy, ode_model, dt, dT, N_p, Xa)
  xtilde <- t(t(rowSums(X)) / N_sigma)

  Pxx <- matrix(rep(0, length = N_x * N_x), nrow = N_x, ncol = N_x)
  for(i in 1:N_sigma){
    Pxx <- Pxx + ((X[, i] - xtilde) %*% t(X[, i] - xtilde)) / N_sigma
  }

  # Add process noise Q at each step (only to parameter block)
  Q <- matrix(0, nrow = N_x, ncol = N_x)
  Q[1:N_p, 1:N_p] <- Q_scale * diag(N_p)
  Pxx <- Pxx + Q  # Add process noise to covariance

  # grab first N_y ind vars
  # skip N_p params from augmented
  # previously used multiple_param_obs_fct
  Y <- X[(N_p + 1):(N_p + N_y),]
  ytilde <- t(t(rowSums(Y)) / N_sigma)

  Pyy <- R
  for(i in 1:N_sigma){
    Pyy <- Pyy + ((Y[, i] - ytilde) %*% t(Y[, i] - ytilde)) / N_sigma
  }

  Pxy <- matrix(rep(0, length = N_x * N_y), nrow = N_x, ncol = N_y)
  for(i in 1:N_sigma){
    Pxy <- Pxy + ((X[, i] - xtilde) %*% t(Y[, i] - ytilde)) / N_sigma
  }
  K <- Pxy %*% solve(Pyy)  # Kalman Gain Matrix

  xhat <- xtilde + K %*% (y - ytilde)
  
  if (forcePositive) {
    param_min <- 1e-8  # Smallest allowed value for parameters
    xhat[1:N_p] <- pmax(param_min, xhat[1:N_p])  # Ensure positive
  }
  
  Pxx <- Pxx - K %*% t(Pxy)

  # xhat is Kalman augmented prediction
  return(list(xhat = xhat, Pxx = Pxx, K = K))
} # ENDFN UKF

#' propagate_model
#' Used inside UKF_dT, previously named kura_ode_model
#' Take the ode model and current augmented state
#' and propagate ind vars to next dT using runge-kutta
#' @param t dummy time variable, because ode models don't have explcit time in our examples
#' @param ode_model model with ind vars y and N_p parameters
#' @param dt smaller time step size within dT
#' @param dT time step size that comes from time series data step
#' @param N_p number of unknown model params
#' @param x augmented state vector
#' @return propagated augmented column vector
#' @examples
#' Example
#' @export
propagate_model <- function(t, ode_model, dt, dT, N_p, x) {
  # dt <- 0.1*dT
  nn <- round(dT / dt)
  p <- x[1:N_p, ]  # grab parameters
  y <- x[(N_p + 1):length(x[, 1]), ]  # grab ind vars
  # Runge-Kutta
  for(n in 1:nn){
    k1 <- dt * ode_model(t, y, p)
    k2 <- dt * ode_model(t, y + k1 / 2, p)
    k3 <- dt * ode_model(t, y + k2 / 2, p)
    k4 <- dt * ode_model(t, y + k3, p)
    y <- y + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
  }
  r <- rbind(x[(1:N_p), ], y) # returns new augmented
  return(r)
} # ENDFN propagate_model

#=====================================================#
#' UKF_blend
#' Apply UKF_dT to all time points in ts_data, which is
#' multi-dimensional time-series, where time is the first
#' column of ts_data. Orignal data file is in rows, but we
#' transpose before the call of UKF_blend.
#' Data should have N_y+1 cols, N_y is num ind vars.
#' Blends UKF_dT with ts_data at each dT.
#' Updates state y and the N_p ode_model parameters.
#' Returns chi-square error output, which can be used
#' as objective function for additional ode_model param
#' optimization. iterative_param_optim and optim_param fns.
#' UKF_blend only steps through time series one time.
#' Iterating UKF_blend with parameters from previous run can
#' improve parameter estimates. Iteration of UKF is what
#' iterative_param_optim does.
#' Note: scalar coefficients for Q and R are hard-coded
#' This function was known previously as myUKFkura.
#' @param t_dummy a dummy time variable, because ode models don't have explcit time
#' @param ts_data panel of time series data. First column holds time values, other N_y columns hold variable values
#' @param ode_model model function with ind variables y and N_p parameters
#' @param N_p number of unknown model params
#' @param N_y number of ind variables
#' @param param_guess vector of length N_p with initial guess for ode_model parameters
#' @param dt smaller time step size within dT for solving ode
#' @param dT time step size that comes from time series data step
#' @param R_scale related to standard deviation of Gaussian measurement noise of ind variables. A number that must be specified by user.
#' @param Q_scale related to standard deviation of process noise. Noise related to model parameters. User choice. Can it be 0?
#' @param forcePositive logical, if TRUE, ensures all parameters stay positive
#' @param seeded logical, if TRUE, sets seed for reproducibility
#' @return list: param_est (N_p model parameter estiamtes after run through time series), xhat (augmented Kalman update), error (error from Pxx augmented covariance at sigma points), and chisq (chi-square goodness of fit of prediction to data)
#' @examples
#' Example
#' @export
UKF_blend <- function(t_dummy, ts_data, ode_model, N_p, N_y,
                      param_guess, dt, dT, R_scale = 0.3, Q_scale = 0.015,
                      forcePositive = FALSE, seeded = FALSE) {
  time_points <- ts_data[, 1]
  num_time <- length(time_points)
  N_x <- N_p + N_y
  # Modified the length of xhat to be N_x*num_time
  xhat <- matrix(rep(0, length = N_x * num_time), nrow = N_x, ncol = num_time)
  # intialize Pxx
  Pxx <- vector(mode = "list", length = num_time)
  for(i in 1:num_time){
    Pxx[[i]] <- matrix(rep(0, length = N_x * N_x), nrow = N_x, ncol = N_x)
  }
  # Modified the length of errors to be N_x*num_time
  errors <- matrix(rep(0, length = N_x * num_time),
                   nrow = N_x, ncol = num_time)
  # intialize Ks
  Ks <- vector(mode = "list", length = num_time)
  for(i in 1:num_time){
    Ks[[i]] <- matrix(rep(0, length = N_x * N_y), nrow = N_x, ncol = N_y)
  }

  # intialize x augmented
  # z stacked on top of y
  # z is param guess repeated for all time points
  z <- t(t(param_guess)) %*%
    matrix(rep(1, length = num_time), nrow = 1, ncol = num_time)
  # y0 is intial y unaugmented state
  # transpose makes each column y values that's num_time wide
  y0 <- t(ts_data[, -1])
  # 2nd and 3rd cols of dat, first col (-1 remove) is time
  # x columns are intial augmented states
  x <- rbind(z, y0) # stack
  xhat[, 1] <- x[, 1]  # right now xhat is all data
  #Q <- 0.015
  #                [grab y rows from x augmented]
  R <- (R_scale)^2 * cov(t(x[(N_p + 1):(N_y + N_p), ])) # only y
  Pxx[[1]] <- pracma::blkdiag(Q_scale * diag(N_p), R)
  #max_steps <- 300
  #steps <- 0
  #convergence <- 1
  #param_array <- list()
  # ?while loop would start here for repeats through time series
  # intial augmented x
  # z is the parameter guesses repeated for all time values
  # y0 is the observed data from the time series
  #x <- rbind(z,y0)
  #xhat[,1] <- x[,1]
  if (seeded) set.seed(1)
  #    grabbing y from x augmented
  y <- x[(N_p + 1):(N_y + N_p), ] + (pracma::sqrtm(R)$B) %*%
    matrix(rnorm(2 * num_time), nrow = 2, ncol = num_time)

  # Create intial dummy UKF_kstep as a default value
  # because cholesky, and hence UKF_dT in for-loop,
  # sometimes fail
  UKF_kstep <- list(xhat = xhat[, 1], Pxx = Pxx[[1]], K = Ks[[1]])
  for(k in 2:num_time){
    # t_k <- time_points[k-1]
    # bam: added tryCatch because by chance chol will fail
    # due to input matrix being non positive definite.
    UKF_kstep <- tryCatch(
      {
        UKF_dT(t_dummy, ode_model, xhat[, k - 1], Pxx[[k - 1]], y[, k],
                N_p, N_y, R, dt, dT,
                R_scale, Q_scale, forcePositive = forcePositive)
      },
      error=function(cond) {
        message("By chance, matrix caused Cholesky to fail.")
        message("Here's the original error message:")
        message(cond)
        message("Skipping this iteration.")
        return(UKF_kstep)
        # if error occurs, return previous out
      }
    ) # end tryCatch
    xhat[, k] <- UKF_kstep$xhat
    Pxx[[k]] <- UKF_kstep$Pxx
    Ks[[k]] <- UKF_kstep$K
    errors[, k] <- sqrt(diag(Pxx[[k]]))

    } # end application of UKF blend to all time points

    #est <- t(xhat[(1:N_p),num_time])
    # parameters at last time point
    param_estimated <- t(t(xhat[(1:N_p), num_time]))
    # bam: Maybe we wnat to divide by the abs average of
    # param_esimated to normalize
    #convergence <- sum(abs(param_estimated - param_guess))
    #param_array[[steps+1]] <- param_estimated
    #param_guess <- param_estimated
    #steps <- steps + 1
    #if(steps >= max_steps){
    #  break
    #}
    # this would end a while loop for repeats through ts
  chisq <- 0
  for(i in 1:N_y){
    chisq <- chisq + (x[(N_p + i), ] - xhat[(N_p + i), ])^2
  }
  chisq <- mean(chisq)
  error <- t(errors[(1:N_p), num_time])

  return(list(param_est = param_estimated,
              xhat = xhat, error = error, chisq = chisq))
} # END FN UKF_blend

#' optim_params
#' Uses pracma optimization methods for either Simulated
#' Annealing or Nelder-Meade. Calls UKF_blend to apply UKF
#' to time series data and model to update parameters and return
#' chi-square goodness of fit for the objective function.
#' @param param_guess vector of length N_p with initial guess for ode_model parameters
#' @param method string to specify optimization method,"L-BFGS-B" or "SANN"
#' @param lower_lim lower bound of solution (L-BFGS only)
#' @param upper_lim upper bound of solution (L-BFGS only)
#' @param maxit maximum number of iterations (SANN only)
#' @param temp annealing temperature
#' @param t_dummy a dummy time variable, because ode models don't have explcit time
#' @param ts_data panel of time series data. First column holds time values, other N_y columns hold variable values
#' @param ode_model model function with ind variables y and N_p parameters
#' @param N_p number of unknown model params
#' @param N_y number of ind variables
#' @param dt smaller time step size within dT for solving ode
#' @param dT time step size that comes from time series data step
#' @param R_scale related to standard deviation of Gaussian measurement noise of ind variables. A number that must be specified by user.
#' @param Q_scale related to standard deviation of process noise. Noise related to model parameters
#' @param forcePositive logical, if TRUE, ensures all parameters stay positive
#' @param seeded logical, if TRUE, sets seed for reproducibility
#' @return list: par (vector of N_p optimized parameters) and value (final chi-square goodness of fit to time series)
#' @examples
#' Example
#' @export
optim_params <- function(param_guess, method = "L-BFGS-B",
                        lower_lim, upper_lim, maxit, temp = 20,
                        t_dummy, ts_data, ode_model, N_p, N_y, dt, dT,
                        R_scale, Q_scale, forcePositive = FALSE, seeded = FALSE) {
  # Optimize the model parameters
  # if method=L-BFGS-B
  #     lower/upper constraints used, no maxit
  # if method=SANN,
  #      temp used, maxit used to stop run, no constraints
  # base stats OPTIM method, Nelder-Meade with
  # param_guess = c(0,0) for example
  # lower_lim=-20 and upper_lim=20 constraints
  # opt <- optim_params(param_guess=c(0,0),
  #                     lower_lim=-20,upper_lim=20,
  #                     t_dummy,smoothed_data,
  #                     kuramoto_ode_model,
  #                     N_y,N_p,N_x,dt,dT)
  # opt$par   # params
  # opt$value # objective function value, chi-square

  # Initialize ukf_obj as NULL in parent environment to store UKF_blend output
  ukf_obj <- NULL

  chisq_objective <- function(par_vec){
    # objective function
    # one full pass through the time series.
    # Assign the result of UKF_blend to ukf_obj in the parent environment
    ukf_obj <<- UKF_blend(t_dummy, ts_data, ode_model,
                          N_p, N_y, par_vec, dt, dT,
                          R_scale, Q_scale, forcePositive = forcePositive,
                          seeded = seeded)
    return(ukf_obj$chisq)
  }
  # from stats base library
  if (method == "SANN") {
    # simulated annealing
    opt <- optim(param_guess, chisq_objective, method = "SANN",
                 control = list(maxit = maxit, temp = temp))
  } else {
    # Broyden-Fletcher-Goldfarb-Shannon
    opt <- optim(param_guess, chisq_objective, method = "L-BFGS-B",
                 lower = lower_lim, upper = upper_lim)
  }
  # Modified the return parameters, param_est and xhat added
  return(list(par = opt$par, value = opt$value,
              param_est = opt$par, xhat = ukf_obj$xhat))
}

#' iterative_param_optim
#' Optimize the model parameters by iteratively running
#' the UKF_blend through the time series data.
#' In other words, repeated calls of UKF_blend to apply UKF
#' to time series data and model to update parameters and return
#' chi-square goodness of fit for the objective function.
#' Convergence of the parameter vector between runs is the
#' stopping criterion
#' @param param_guess vector of length N_p with initial guess for ode_model parameters
#' @param t_dummy a dummy time variable, because ode models don't have explcit time
#' @param ts_data panel of time series data. First column holds time values, other N_y columns hold variable values
#' @param ode_model model function with ind variables y and N_p parameters
#' @param N_p number of unknown model params
#' @param N_y number of ind variables
#' @param dt smaller time step size within dT for solving ode
#' @param dT time step size that comes from time series data step
#' @param param_tol small tolerance for convergence of the parameters to be optimized
#' @param MAXSTEPS max number of iterations through time series
#' @param R_scale related to standard deviation of Gaussian measurement noise of ind variables. A number that must be specified by user.
#' @param Q_scale related to standard deviation of process noise. Noise related to model parameters
#' @param trace logical, if TRUE, returns the best step and chi-square value
#' @param forcePositive logical, if TRUE, ensures all parameters stay positive
#' @param seeded logical, if TRUE, sets seed for reproducibility
#' @return list: par (vector of N_p optimized parameters) and value (final chi-square goodness of fit to time series)
#' @examples
#' Example
#' @export
iterative_param_optim <- function(param_guess, t_dummy, ts_data, ode_model,
                                  N_p, N_y, dt, dT, param_tol = 0.01, MAXSTEPS = 30,
                                  R_scale = 1, Q_scale = 1, trace = FALSE,
                                  forcePositive = FALSE, seeded = FALSE) {
  done <- FALSE
  steps <- 0
  chisq_history <- numeric()  # Store chi-square values for each iteration

  # Initialize best parameter tracking
  best_param <- param_guess
  best_chisq <- Inf
  best_xhat <- NULL

  while (!done) {
    # One run through whole time series
    ukf_run <- UKF_blend(t_dummy, ts_data, ode_model,
                         N_p, N_y, param_guess, dt, dT,
                         R_scale, Q_scale, forcePositive = forcePositive,
                         seeded = seeded)

    param_new <- ukf_run$param_est
    chisq_history <- c(chisq_history, ukf_run$chisq)  # Append chi-square

    # Track best parameters
    if (ukf_run$chisq < best_chisq) {
      best_chisq <- ukf_run$chisq
      best_param <- param_new
      best_xhat <- ukf_run$xhat
    }

    steps <- steps + 1
    param_norm <- abs(sum(param_new - param_guess))
    converged <- param_norm < param_tol
    done <- converged | steps >= MAXSTEPS

    param_guess <- param_new  # Update guess for next iteration
  }
  # If trace is TRUE, return the best step and chi-square value
  if (trace) {
    return(list(
      par = best_param,
      value = best_chisq,
      param_norm = param_norm,
      steps = steps,
      param_est = best_param,
      xhat = best_xhat,
      chisq = chisq_history  # chi-square over iterations
    ))
  } else {
    # Return includes full chi-square trace
    return(list(
      par = param_new,
      value = ukf_run$chisq,
      param_norm = param_norm,
      steps = steps,
      param_est = param_new,
      xhat = ukf_run$xhat,
      chisq = chisq_history  # NEW: chi-square over iterations
    ))
  }
}
