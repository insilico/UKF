######################################
# model_voxels_2param.R
# example script for applying UKF
# to a pair of voxel time series
# and optimizing coupling parameters
# for an ode model of synchronization
######################################

# parameters to optimize, p1 p2, in coupled_osc_model
param_guess <- c(4,4)  # intial guess, oscillator coupling
# script runs multiple optimization methods

########################
# load fMRI voxel data
# smooth and make plots
########################

# first row is time, then voxel A and voxel B rows
# UKF_blend assumes column data, so transpose
vox.A.B.data <- read.delim("data/voxel_A&B_data.txt",sep="",header=F)
# make column-wise
vox.A.B.data <- t(vox.A.B.data)
# show panel plot of each voxel with kernel smooth
# return kernel smoothed data
# plot_and_smooth.R
smoothed_data <- plot_voxels_and_smooth(vox.A.B.data)

############################
# Set up Model Dimensions  #
############################

t_vec <- vox.A.B.data[,1] # time is first column
dT <- t_vec[2]-t_vec[1]   # assume uniform time steps dT=1
# smaller steps size for propagating model between dT steps
dt <- 0.1*dT
#nn <- round(dT/dt)
# num observed ind vars, first col is time, so -1
N_y <- ncol(vox.A.B.data)-1  # N_y=2
# number of unknown model parameters to be estimated
N_p <- 2
# size of augmented state vector
N_x <- N_p + N_y

###########################
# Specify Model Structure #
###########################
coupled_osc_model <- function(t,x,p){
  # Coupled Oscillator for Synchronization
  # Coupled pendulum approximation.
  # How related to Kuramoto?
  # p are unknown parameters that we will estimate.
  # These represent the coupling between voxels.
  # The other parameters could also be optimized or
  # determined by other means.
  # There is no explicit t in model but keep anyway
  w1 <- 0 #1
  w2 <- 0 #0.53
  a1 <-  1 #0.8
  a4 <-  1 #1.2
  r <- rbind((w1 + a1*sin(x[1,]) + p[1,]*sin(x[2,])),
             w2 + p[2,]*sin(x[1,]) + a4*sin(x[2,]))
  return(r)
} # ENDFN coupled_osc_model

###################################
# Run UKF with model and RAW data
# One pass through the time series
###################################
ukf_out <- UKF_blend(t_dummy,vox.A.B.data,coupled_osc_model,
                     N_p,N_y,param_guess,dt,dT)
ukf_out$param_est
ukf_out$chisq
# function in plot_and_smooth.R
plot_ukf_and_smoothed(ukf_out, smoothed_data,
                      top_title='One UKF Step, Raw Data')

#########################################
# Run UKF with model and SMOOTHED data as input
# One pass through the time series
#########################################
ukf_from_smooth <- UKF_blend(t_dummy,smoothed_data,
                             coupled_osc_model,
                     N_p,N_y,param_guess,dt,dT)
ukf_from_smooth$param_est
ukf_from_smooth$chisq
plot_ukf_and_smoothed(ukf_from_smooth, smoothed_data,
                    top_title = 'One UKF Step, Smooth Data')

#########################################
# Simulated Annealing with smoothed data as input
# to optimize model parameters, chisquare goodness of fit
#########################################
### Might want to make a separate UKF_blend that does not
### change the parameters because the params get changed a lot
### and yet the chisquare is the same.
### Rugged fitness landscape?
# Simulated Annealing, SANN ignores lower/upper limits
opt <- optim_params(param_guess,method="SANN",
                    lower_lim=-20,upper_lim=20,
                    maxit=500,temp=20,
                    t_dummy,smoothed_data,coupled_osc_model,
                    N_p,N_y,dt,dT)
opt$par   # params
opt$value # objective function value, chi-square
# plug optim parameter back into UKF and plot
# TODO: Make UKF_blend that just updates the y state,
# not the parameters. A run through UKF_blend can drastically
# change the guess parameters that were optimized, which
# makes you wonder about the fitness landscape.
ukf_optimized <- UKF_blend(t_dummy,smoothed_data,
                           coupled_osc_model,
                           N_p,N_y,opt$par,dt,dT)
ukf_optimized$param_est
ukf_optimized$chisq
plot_ukf_and_smoothed(ukf_optimized, smoothed_data,
                      top_title = 'Simulated Annealing')

#########################################
# Nelder-Mead with smoothed data as input
# to optimize the parameter using UKF chisquare
#########################################
# L-BFGS-B, Nelder-Meade ignores maxit and temp
opt2 <- optim_params(param_guess,method="L-BFGS-B",
                    lower_lim=-20,upper_lim=20,
                    maxit=500,temp=20,
                    t_dummy,smoothed_data,coupled_osc_model,
                    N_p,N_y,dt,dT)
opt2$par   # params
opt2$value # objective function value, chi-square
# plug optim parameter back into UKF and plot
ukf_optimized <- UKF_blend(t_dummy,smoothed_data,
                           coupled_osc_model,
                           N_p,N_y,opt2$par,dt,dT)
ukf_optimized$param_est
ukf_optimized$chisq
plot_ukf_and_smoothed(ukf_optimized, smoothed_data,
                      top_title = 'Nelder-Meade')

#########################################
# Use the iterative approach to optimize params.
# Convergence of parameters between iterations.
# Run UKF model through smoothed data iteratively.
#########################################
iter_opt <- iterative_param_optim(param_guess,
                      t_dummy, smoothed_data,
                      coupled_osc_model,
                      N_p,N_y,dt,dT,
                      param_tol=.01,MAXSTEPS=30)
iter_opt$par   # params
iter_opt$value # chi-square
iter_opt$steps
iter_opt$param_norm
# plug iterative optim parameters back into UKF and plot
ukf_iter_opt <- UKF_blend(t_dummy,smoothed_data,
                           coupled_osc_model,
                           N_p,N_y,iter_opt$par,dt,dT)
ukf_iter_opt$param_est
ukf_iter_opt$chisq
plot_ukf_and_smoothed(ukf_iter_opt, smoothed_data,
                      top_title = 'Iterative Optimization')
