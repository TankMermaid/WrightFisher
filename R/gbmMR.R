#' Mean-Reverting geometric brownian motion
#' 
#' @export
#' @param nspecies number of companies to simulate
#' @param mu Number or vector of length \code{nspecies} for the drift
#' @param sigma Number of vector of length \code{nspecies} for the volatility
#' @param Tmax End of time-window for simulation
#' @param X0 initial abundances
#' @param nsamples Number of samples to collect between [0,Tmax]
#' @param dt size of timesteps - an alternative to nsamples for specifying numerical integration. If input, the output will be of length Tmax/dt
#' @examples 
#' library(WrightFisher)
#' library(plotrix)
#' 
#' x <- gbmMR(10)
#' par(mfrow=c(1,2))
#' stackpoly(x$Market,stack=T,main='Market Dynamics')
#' stackpoly(x$Shares,stack=T,main='Market Share Dynamics')
gbmMR <- function(nspecies,mu=1,sigma=1,rho=0,Tmax=1,X0=NULL,nsamples=1000,dt=NULL){
  
  #simulates brownian motions & then exponentiates
  
  
  X <- matrix(0,nrow=nsamples,ncol=nspecies)
  X[1,] <- rnorm(nspecies)
  x <- X[1,]
  
  if (is.null(dt)){
    dt <- Tmax/nsamples
  }
  times <- seq(0,Tmax,by=dt)
  
  
  
  for (nn in 2:nsamples){
    dw <- rnorm(nspecies,sd = sqrt(dt))
    X[nn,] <- X[nn-1,]+mu*(rho-X[nn-1,])*dt+sigma*dw
  }
  
  output <- list(times,exp(X),exp(X)/rowSums(exp(X)))
  names(output) <- c('time','Market','Shares')
  return(output)
}