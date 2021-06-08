#########################################
# Simulate data and fit osc params
#########################################

#########################################
# 1. Simulate time series data
#    so we can fit params p1 and p2
#########################################
# Parameters for the coupled pendulum equations
# These give linear trend for one state,
# which could be useful as an example in fMRI data
#(parms <- c(w1 = 1, w2 = 0.53, a1 = 0.8, a4 = 1.2,
#            p1 = 0.1, p2 = 0.1))
# p1 and p2 are parameters we want to fit (coupling).
parms <- c(p1 = 2*pi/3, p2 = 2*pi/3)
# We are later going to try to find p1 and p2
## vector of timesteps

###########################
# Specify Model Structure #
###########################
coupled_osc_model <- function(t,x,p){
  # Coupled Oscillator for Synchronization
  # There is no explicit t in model but keep anyway
  w1 <- 0; w2 <- 0; a1 <- pi/2; a4 <- pi/2
  # x and p are matrices where columns are all time points
  r <- rbind((w1 + a1*sin(x[1,]) + p[1,]*sin(x[2,])),
             w2 + p[2,]*sin(x[1,]) + a4*sin(x[2,]))
  return(r)
} # ENDFN coupled_osc_model

## intial conditions
#(ystart <- c(y1=0.01,y2=0.05))
yinit <- c(y1=0,y2=.2)
### Using our own RK solver so we can sample
### data at specific time points more easily
### dT: uniform sampling time points of simulation
dT=1  
times <- seq(0, 176, by=dT)
## Solve ode to create data
dt <- 0.01*dT  # dt are intermediate values between dT steps
# we only want to predict values at dT
nn <- round(dT/dt)
# propagate_model uses the augmented state, so as a workaround, 
# make a sort of artificial augmented state (x). 
# function will get y and parms p1 and p2 from x.
# x is a matrix where each column is the augmented state
# and columns are each time point. Intialize for all time:
# 2 ode parameters
dq <- 2
dy <- length(yinit)
# for dt time point columns
# intialize x augmented
x<-matrix(rep(c(parms,yinit),nn),nrow=dq+dy,ncol=nn, byrow=F)
# intialize y solution at times vector
osc.sol.mat <- 
  matrix(rep(yinit,length(times)),nrow=dy,ncol=length(times))
for (i in seq(2,length(times))){
  # step solution forward with dt steps with current dT
  temp_dt <-
    propagate_model(t,coupled_osc_model,dt,dT,dq,x)
  # update augmented x with last dt value for next dT step
  # for this parms doesn't change, nn columns for dt
  x <- matrix(rep(temp_dt[,nn],nn),nrow=dq+dy,ncol=nn)
  #x<-matrix(rep(c(parms,temp_dt[,nn]),nn),nrow=dq+dy,ncol=nn)
  osc.sol.mat[,i] <- temp_dt[(dq+1):(dq+dy),nn]
}

## this is for plotting
osc.sol.df <- data.frame(time=times, y1=osc.sol.mat[1,],
                         y2=osc.sol.mat[2,])
library(reshape2)
library(ggplot2)
osc.sol.melted <- melt(osc.sol.df, id = "time")
g <- ggplot(data = osc.sol.melted, aes(x = time, y = value, color = variable)) 
#g <- g + geom_point(size=2) 
g <- g + geom_line() 
g <- g + xlab("time") +  ylab("Amplitude") 
g <- g + ggtitle("Coupled Oscillators")  
g <- g + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
                  legend.position = c(.8,.2),
                  legend.key = element_rect(fill="white")) 
show(g)

#########################################
# 2. Fit osc params from simulated data
#########################################
#########################################
# Simulated Annealing with smoothed data as input  
# to optimize model parameters, chisquare goodness of fit
#########################################
### Might want to make a separate UKF_blend that does not 
### change the parameters because the params get changed a lot
### and yet the chisquare is the same.
### Rugged fitness landscape?
# Simulated Annealing, SANN ignores lower/upper limits
opt <- optim_params(param_guess=c(0,0),method="SANN",
                    maxit=500,temp=20, 
                    t_dummy,osc_sol.df,coupled_osc_model,
                    dy,dq,dx,dt,dT)
opt$par   # params
opt$value # objective function value, chi-square
# plug optim parameter back into UKF and plot
# TODO: Make UKF_blend that just updates the y state,
# not the parameters. A run through UKF_blend can drastically
# change the guess parameters that were optimized, which 
# makes you wonder about the fitness landscape. 
ukf_optimized <- UKF_blend(t_dummy,smoothed_data,
                           coupled_osc_model,
                           dy,dq,dx,opt$par,dt,dT)
ukf_optimized$param_est
ukf_optimized$chisq
plot_ukf_and_smoothed(ukf_optimized, smoothed_data,
                      top_title = 'Simulated Annealing')

# Not using prama ode45 because not sure how to
# control time sampling

# Model for coupled pendulum written as
# a system of first order for ode45
coupled_osc2 <- function(t, y, parms){
  # fixed parameters
  w1 <- 0; w2 <- 0; a1 <- pi/2; a4 <- pi/2
  # p1 and p2 are the coupling paramters
  # they come from parms input
  y1 <- y[1]
  y2 <- y[2]
  with(as.list(c(parms, y)), {
    dy1 <- w1 + a1*sin(y1) + p1*sin(y2)
    dy2 <- w2 + p2*sin(y1) + a4*sin(y2)
    dy <- matrix(data = c(dy1, dy2), nrow = 2, ncol = 1)
  })
}

library(pracma)
osc.sol <- ode45(f=function(t,y){coupled_osc2(t,y,parms)},
                 y=yinit, t0=min(times), tfinal=max(times))

osc.sol.df <- data.frame(time=osc.sol$t, y1=osc.sol$y[,1],
                         y2=osc.sol$y[,2])

osc.sol.melted <- melt(osc.sol.df, id = "time")
g <- ggplot(data = osc.sol.melted, aes(x = time, y = value, color = variable)) 
#g <- g + geom_point(size=2) 
g <- g + geom_line() 
g <- g + xlab("time") +  ylab("Amplitude") 
g <- g + ggtitle("Coupled Oscillators")  
g <- g + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(.8,.2),
        legend.key = element_rect(fill="white")) 
show(g)
