######################################
# plot_and_smooth.R contains utility
# functions for smoothing fMRI
# time series and making some plots
# todo: translate to ggplot
######################################
#install.packages("KernSmooth")
library(KernSmooth)  # dpill function
# ksmooth from stats library

#' plot_voxels_and_smooth
#' Plot a pair of voxel time series and kernel smoothed data.
#' Panel for each voxel with raw and smoothed.
#' Also returns the smoothed data for later use.
#' @param voxel_input three columns: column 1 = time, column 2 = voxel 1 signal, column 3 = voxel 2 signal
#' @return matrix: same shape as input data but kernel smoothed
#' @examples
#' Example
#' @export
plot_voxels_and_smooth <- function(voxel_input){
  voxel_Time <- voxel_input[,1]
  voxel_data <- cbind(voxel_Time,voxel_input[,2],
                        voxel_input[,3])
  voxel.df <- data.frame(times=voxel_data[,1],
                             VoxA=voxel_data[,2],
                             VoxB=voxel_data[,3])
  y1_obs <- voxel.df[,2]
  y2_obs <- voxel.df[,3]

  n <- length(t(voxel_Time))
  m_points <- seq(1,n,by=1)

  #                 mar=c(bottom,left,top,right)
  par(mfrow=c(2,1), mar=c(4.2, 4.5, 1.3, 0.8))

  #plot(m_points, y1_obs)
  h <- tryCatch( # kernel bandwidth using KernSmooth
    {
      dpill(m_points, y1_obs)
    },
    error=function(cond) {
      message("Install/load KernSmooth for better smoothing.")
      message(cond)
      message("Using default kernel bandwidth h=.5.")
      return(0.5)
    }
  ) # end tryCatch
  #h <- dpill(m_points, y1_obs)
  fit1 <- ksmooth(m_points, y1_obs, kernel="normal",
                  bandwidth = h, n.points=176)
  #lines(fit1,col='red',lwd=2)

  #plot(m_points, y2_obs)
  h <- dpill(m_points, y2_obs)
  fit2 <- ksmooth(m_points, y2_obs, kernel="normal",
                  bandwidth = h, n.points=176)
  #lines(fit2,col='red',lwd=2)

  y1_obs.processed <- fit1$y
  y2_obs.processed <- fit2$y

  plot(m_points,y1_obs,pch=3,
       main="real data and kernel smooth",xlab="",
       ylab="Voxel A")
  lines(m_points,y1_obs.processed,col='red',lwd=2)
  plot(m_points,y2_obs,pch=3,xlab="Times",ylab="Voxel B")
  lines(m_points,y2_obs.processed,col='red',lwd=2)

  smoothed.data <- cbind(t(t(m_points)),
                          t(t(y1_obs.processed)),
                          t(t(y2_obs.processed)))
  par(mfrow=c(1,1))
  return(smoothed.data)
}

#' plot_ukf_and_smoothed
#' Plot a pair of voxel time series and UKF blended data.
#' Panel for each voxel with input data and UKF estimate.
#' The legend says Kernel for input data, but could input raw.
#' So, might want to make legend text an input
#' @param ukf_out output object from UKF_blend
#' @param smoothed_data time series panel used for UKF_blend
#' @param top_title string title typically to describe how parameters were optimized
#' @return Just a plot
#' @examples
#' Example
#' @export
plot_ukf_and_smoothed <- function(ukf_out, smoothed_data,
                top_title="UKF estimate and kernel smooth"){
  par(mfrow=c(2,1), mar=c(4.2, 4.5, 1.3, 0.8))
  plot(ukf_out$xhat[3,], xlab='Time', ylab='Voxel A')
  lines(smoothed_data[,2],col='red')
  title(top_title)
  legend('topright', legend=c("UKF", "Kernel"), bty='n',
         col=c("black", "red"), pch=c('o', NA), lty=c(NA,1))

  plot(ukf_out$xhat[4,], ylab='Voxel B', xlab='Time')
  lines(smoothed_data[,3],col='red')
  title(paste('p1, p2 = ',
            round(ukf_out$param_est[1],digits=2),
            round(ukf_out$param_est[2],digits=2)))
  legend('topright', legend=c("UKF", "Kernel"), bty='n',
         col=c("black", "red"), pch=c('o', NA), lty=c(NA,1))
}
