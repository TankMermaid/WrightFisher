#' Performs Washburne et al.'s (2016) Neutrality Test with corrected ks-test
#' 
#' @export
#' @param R Compositional matrix whose columns are species and whose rows are timepoints
#' @param tm Optional time-points. Use if time intervals between samples are not uniform.
#' @param ntests Number of constant-volatility transforms to use in test. Must be less than or equal to 2^dim(R)[2]
#' @param regress String, either 'f', or 'tm', based on whether heteroskedasticity should regress against the state, f, or the time, tm. Default is 'f'
#' @param ncores Number of cores for parallelization of constant-volatility tests.
#' @param formula Formula for regression of drift. Default DF~f+I(f^2)
#' @param varformula Formula for regression of residuals in Breush-Pagan test. Default DF~f+I(f^2)
#' @param standard Logical indicating whether to use a standard, reproducible choice of CVTs by setting seed within NeutralCovTest
#' @param ... optional input arguments for bptest
#' @examples
#' library(WrightFisher)
#' library(plotrix)
#' 
#' set.seed(10)
#' # For a simple test, we can produce a dataset:
#' R <- wfp(nspecies=12,dt=1e-3)$Community
#' Rg <- gbmMR(nspecies=12,dt=1e-3,sigma=3,mu=2)$Shares
#' 
#' # Visualizing Market/Community Dynamics
#' par(mfrow=c(2,2))
#' stackpoly(R,stack=T,main='Wright-Fisher Process',ylab='Relative Abundance',xlab='time')
#' stackpoly(Rg,stack=T,main='Mean-Reverting GBM',ylab='Relative Abundance',xlab='time') 
#' 
#' # Neutrality Covariance Testing
#' P_wfp <- NeutralCovTest(R)
#' P_gbm <- NeutralCovTest(Rg)
#' P_wfp
#' P_gbm
#' 
#' # In the script below, we'll simulate 'reps' neutral communities through simulation of a Volatility-stabilized market, and then
#' # we'll perform a NeutralCovTest on each, using 1000 constant-volatility transformations. If performed on a computer with 4 or more
#' # cores, the parallelized version will be run for a speed-up.  
#'  
#' reps <- 20
#' Pvals_vsm <- rep(NA,reps)
#' Pvals_gbm <- rep(NA,reps)
#' parallelize=F
#' 
#' for (nn in 1:reps){
#' 
#' ## Simulating Communities
#' R <- vsm(nspecies=12,lambda=25,dt=1e-4,Tmax=1,nsamples=125)$Shares
#' R2 <- gbmMR(12,Tmax=1,nsamples=125,sigma=4,mu=15)$Shares
#' 
#' ## Neutrality Testing
#' if (detectCores()>=3 && parallelize==T){
#' # For species-rich datasets, NeutralCovTest may take a while and can be sped-up through paralellization by inputting ncores
#'  ncores=detectCores()-1
#'  Pvals_vsm[nn] <- NeutralCovTest(R,ncores=ncores)
#'  Pvals_gbm[nn] <- NeutralCovTest(R2,ncores=ncores)
#' } else {
#'  Pvals_vsm[nn] <- NeutralCovTest(R)
#'  Pvals_gbm[nn] <- NeutralCovTest(R2)
#' }
#' 
#' }
#' 
#' plot(ecdf(Pvals_vsm),xlab='NeutralCovTest Pvalue',ylab='F(P)',main='Pvalues: NeutralCovTest on VSMs',xlim=c(0,1),ylim=c(0,1))
#' lines(c(0,1),c(0,1),col='blue',lwd=2)
#' legend(.6,.4,legend=c('Simulations','Uniform'),pch=c(16,NA),lwd=c(1,2),col=c('black','blue'))
#' plot(ecdf(Pvals_gbm),xlab='Pvalue',ylab='F(P)',main='Pvalues: NeutralCovTests on mean-reverting GBMs',xlim=c(0,1),ylim=c(0,1))
#' lines(c(0,1),c(0,1),col='blue',lwd=2)
#' legend(.6,.4,legend=c('Simulations','Uniform'),pch=c(16,NA),lwd=c(1,2),col=c('black','blue'))


NeutralCovTest <- function(R,tm=NULL,ntests=NULL,regress='f',ncores=NULL,formula=DF~f+I(f^2),varformula=NULL,standard=T,...){
  
  
  # Pull out
  S <- dim(R)[2]
  m <- dim(R)[1]
  
  
  data("bounds")
  upperbound <- bounds$upperbound
  
  if (is.null(ntests)){
    ntests <- min(16000,floor((upperbound*S)^3))
  } else {
    if (ntests>(upperbound*S)^3){
      warning('ntests input is large relative to number of species, leading to high false-positive rates for P<0.05')
    }
  }
  
  
  if (standard){
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
  }
   
  if (is.null(ncores)){
    if (standard){
      set.seed(10)
    }
    Pvals <- CVTest(ntests,R,regress,formula,varformula,tm,...) 
  } else{
    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl,library(WrightFisher))
    parallel::clusterEvalQ(cl,library(lmtest))
    if (standard){
      seedList <- as.list(101^seq(1,length(cl)))
      parallel::clusterApply(cl,x=seedList,fun=set.seed)
    }
    
    testsEach <- round(ntests/ncores)
    
    Pvals <- unlist(parallel::parLapply(cl,X=as.list(rep(testsEach,ncores)),fun = CVTest,R=R,regress=regress,formula=formula,varformula=varformula))
  
    stopCluster(cl)
    gc()
  }
  
  ntests <- length(Pvals)
  D <- suppressWarnings(ks.test(Pvals,runif(1e5))$statistic)
  
  data('WashburneKSn')
  Q <- predict(WashburneKSn,newdata=data.frame('f'=ntests,'s'=S))
  nstar_est <- ntests/(1+exp(-Q))
  P <- 1-kolmim::pkolmim(D,nstar_est)
  
  
  return(P)
}