#' Simulate volatility-stabilized market model
#' 
#' @param n Number of companies in market
#' @param delta Growth rate of market
#' @param lambda Alternative parameterization via underlying Wright-Fisher Process
#' @param rho Alternative parameterization via underlying Wright-Fisher Process
#' @param dt Size of timestep in numerical integration
#' @param Tmax Time window of numerical integration [0,Tmax]
#' @param X0 initial community size, default is rlnorm(n)
#' @param tol Tolerance - reflecting wall close to zero to prevent negative values
#' @return Two member list. First element is time vector. Second element is matrix whose columns are companies in the market and whose rows are timepoints
#' @export
#' @examples 
#' library(plotrix)
#' set.seed(1)
#' X0=rmultinom(1,size=200,prob=rep(1,5))
#' X <- vsm(5,Tmax=1,nsamples=5e3,dt=1e-4,X0=X0)
#' tms <- X$time
#' X <- X$Market
#' R <- X/rowSums(X)
#' par(mfrow=c(2,1))
#' stackpoly(X,stack=T,main='Market Dynamics',ylim=c(0,350))
#' stackpoly(R,stack=T,main='Market Shares')
#' 
#' 
#' ### now let's see how the variance scales with the mean
#' n=50
#' lambda=1
#' fs <- function(S) S
#' X0 <- rlnorm(50)
#' #the step below takes ~10s
#' X <- vsm(n,lambda=lambda,Tmax=1,nsamples=5e3,X0=X0,dt=1e-4,fs=fs)
#' tms <- X$time
#' X <- X$Market
#' par(mfrow=c(2,2))
#' stackpoly(X,stack=T,main='Market Dynamics')
#' R <- X/rowSums(X)
#' stackpoly(R,stack=T,main='Relative Abundances')

#' ms <- colMeans(X) 
#' vs <- apply(X,MARGIN=2,var)
#' model <- glm(log(vs)~log(ms))
#' model$coefficients[2]
#' prediction <- exp(predict(model))
#' 

#' plot(ms,vs,xlab='mean',ylab='variance',main='Variance Scaling of VSM',log='xy')
#' lines(ms,prediction,lwd=2,col='blue')
#' legend(min(ms),.3*max(vs),legend=paste('slope=',substr(toString(model$coefficients[2]),1,4),sep=''),col='blue',lty=1,lwd=2)
#' 
#' ### Finally, let's perform a neutral test
#' Pvals <- CVTest(500,R)
#' plot(ecdf(Pvals),main='CVT-Neutrality Test',xlim=c(0,1),ylim=c(0,1))
#' lines(c(0,1),c(0,1),lwd=2,col='red')
#' ksP <- ks.test(Pvals,runif(length(Pvals)))$p.value
#' ksP <- round(1000*ksP)/1000
#' legend(0,.9,legend=c('Neut. Expectation',paste('P=',toString(ksP))),lty=c(1,NA),col=c('red',NA),lwd=c(2,NA))
vsm <- function(n=2,delta=NULL,lambda=NULL,rho=NULL,fs=NULL,dt=1e-5,Tmax=1,X0=NULL,tol=1e-16,nsamples=1000){
  
  if (is.null(lambda)==F){
    if (is.null(rho)){
      rho=rep(1/n,n)
    } else { 
      if (length(rho)!=n){stop('If input lambda, must also input compositional vector rho of length n')}
    }
    delta=rho*lambda
  } else {
    delta = rep(2/n,n)
  }
  
  if (is.null(fs)){
    fs=function(S) S
  }

  timesteps = Tmax/dt
  times <- seq(0,Tmax,length.out=nsamples)
  X <- matrix(0,nrow=nsamples,ncol=n)
  samples <- round(seq(1,timesteps,length.out = nsamples))
  if (is.null(X0)){
    X[1,] <- rlnorm(n)   #initial abundances
  } else {
    X[1,] <- X0
  }
  
  colnames(X) <- unlist(lapply(1:n,FUN = function(x,y) paste(y,x),y='species'))

  x <- X[1,]
  ns <- 1
  for (tt in 2:timesteps){
    S <- sum(x)
    dw <- rnorm(n,sd = sqrt(dt*abs(fs(S))*x))
    x <- x+0.5*delta*fs(S)*dt+dw
    x[x<tol]=tol
    if (tt %in% samples){
      ns <- ns+1
      X[ns,] <- x
    }
  }
  output <- list(times,X)
  names(output) <- c('time','Market')
  return(output)
}
