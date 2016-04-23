#' Simulate Wright-Fisher Process
#' 
#' @param nSpecies One way to parameterize WFP by number of species. \code{rho} will be rep(1/nSpecies,nSpecies)
#' @param rho Compositional vector signifying the mean of the community
#' @param lambda Rate of migration or mean-reversion
#' @param X0 optional initial community state - must be a compositional vector of same size as rho
#' @param dt size of timesteps for numerical integration
#' @param Tmax Size of time window of numerical integration
#' @param sample Number of timepoints between 0 and T to sample. Output will be a matrix with nrow=length(rho) and ncol=sample
#' @export
#' @return Outputs a two-element list. The first element is the time vector for the time of samples, and the second element is the matrix of the fluctuating community.
#' @examples 
#' library(plotrix)
#' X <- WFP(rep(.2,5),dt=1e-6,T=0.25)
#' 
#' time <- X[[1]]
#' X <- X[[2]]
#' stackpoly(X,stack=T)
WFP <- function(nSpecies=NULL,rho=rep(1/nSpecies,nSpecies),lambda=15,X0=rho,dt=1e-5,Tmax=1,sample=1000){
  #Simulate Wright-Fisher Process using log-ratio transformation and subsequent inversion as described in Washburne's (2015) thesis
  
  p <- rho
  l <- lambda
  n <- round(Tmax/dt)
  d <- length(p)
  m <- length(X0)
  if (d!=m){stop('dimensions of rho and X0 are not equal')}
  x <- matrix(0,m,n+1)
  x[,1] <- log(X0/(1-X0))
  
  dW <- sqrt(2*dt)*matrix(rnorm(m*(m-1)*n/2),m*(m-1)/2,n)
  
  P <- matrix(0,m,m*(m-1)/2)  
  j <- 1
  i1 <- NULL
  i2 <- i1
  
  for (i in 1:(m-1)){
    P[i:m,(j:(j+m-i-1))] <- rbind(rep(1,m-i),-diag(m-i))
    i1 <- c(i1,i*rep(1,m-i))
    i2 <- c(i2,(i+1):m)
    j <- j+m-i
  }
  
  Dm <- diag(m)
  Di <- diag(length(i2))
           
  
  for (t in 1:n){
    x[,t+1] <- x[,t]  +  4*cosh(x[,t]/2)^2*(l*p+(1-d)/d+(1-l)/(1+exp(-x[,t])))*dt  +  (4*cosh(x[,t]/2)^2)*Dm %*% P %*% (sqrt((1/(1+exp(-x[i1,t])))*(1/(1+exp(-x[i2,t]))))*Di)%*%dW[,t]
  }            
  
  time <- seq(0,T,length=sample)
  ind <- round(seq(1,n+1,length=sample))
  x <- 1/(1+exp(-x[,ind]))
  x <- t(x)
  return(list(time,x))
}