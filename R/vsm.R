#' Simulate volatility-stabilized market model
#' 
#' @param n Number of companies in market
#' @param delta Growth rate of market
#' @param dt Size of timestep in numerical integration
#' @param Tmax Time window of numerical integration [0,Tmax]
#' @param X0 initial community size, default is rlnorm(n)
#' @param tol Tolerance - reflecting wall close to zero to prevent negative values
#' @param rho Mean to which market model reverts. Default is 0
#' @return Two member list. First element is time vector. Second element is matrix whose columns are companies in the market and whose rows are timepoints
#' @export
#' @examples 
#' set.seed(1)
#' X0=rmultinom(1,size=200,prob=rep(1,5))
#' X <- vsm(5,Tmax=1,rho=1e3,samples=5e3,dt=1e-4,X0=X0)
#' tms <- X$time
#' X <- X$Market
#' R <- X/rowSums(X)
#' par(mfrow=c(1,2))
#' stackpoly(X,stack=T,main='Market Dynamics',ylim=c(0,350))
#' stackpoly(R,stack=T,main='Market Shares')
#' 
#' 
#' ### now let's see how the variance scales with the mean
#' set.seed(2)
#' n=50
#' X0 <- rmultinom(1,size=250,prob=rep(1,n))
#' #the step below takes ~10s
#' X <- vsm(50,Tmax=10,delta=-3e-2,samples=5e3,X0=X0,dt=1e-4)
#' tms <- X$time
#' X <- X$Market
#' par(mfrow=c(1,1))
#' plot(rowSums(X),type='l',main='Market Size')
#' R <- X/rowSums(X)


#' ms <- colMeans(X) 
#' vs <- apply(X,MARGIN=2,var)
#' model <- glm(log(vs)~log(ms)-1)
#' model$coefficients
#' prediction <- exp(predict(model))
#' 
#' par(mfrow=c(1,2))
#' plot(ms,vs,xlab='mean',ylab='variance',main='Variance Scaling of VSM',log='xy')
#' lines(ms,prediction,lwd=2,col='blue')
#' 
#' ### Finally, let's perform a neutral test
#' Pvals <- CVTest(1000,R)
#' plot(ecdf(Pvals),main='Washburne Test',xlim=c(0,1),ylim=c(0,1))
#' lines(c(0,1),c(0,1),lwd=2,col='red')

vsm <- function(n=2,delta=0.1,dt=1e-5,Tmax=1,X0=NULL,tol=1e-16,rho=0,samples=1000){

  timesteps = Tmax/dt
  times <- seq(0,Tmax,length.out=samples)
  X <- matrix(0,nrow=n,ncol=samples)
  samples <- round(seq(1,timesteps,length.out = samples))
  if (is.null(X0)){
    X[,1] <- rlnorm(n)   #initial abundances
  } else {
    X[,1] <- X0
  }
  
  rownames(X) <- unlist(lapply(1:n,FUN = function(x,y) paste(y,x),y='species'))

  x <- X[,1]
  ns <- 1
  for (tt in 2:timesteps){
    S <- sum(x)
    dw <- rnorm(n,sd = dt*S*x)
    x <- x+delta*(rho-S)*dt+dw
    x[x<tol]=tol
    if (tt %in% samples){
      ns <- ns+1
      X[,ns] <- x
    }
  }
  output <- list(times,t(X))
  names(output) <- c('time','Market')
  return(output)
}
