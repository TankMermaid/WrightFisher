#' Constant Volatility Test for WFP
#' 
#' @export
#' @param ntests Positive integer - number of randomly selected constant-volatility transformations (CVTs) for heteroskedasticity testing
#' @param R Compositional matrix - columns must be species and rows timepoints.
#' @param regress Either "f" or "t" - whether regrssion should be performed on the state of the CVT or on the time.
#' @param formula Formula for regression of the jumps in CVTs vs time or state of cvts. Should be in form "DF ~ f" or "DF ~ tm + I(tm^2)" or  etc. Default is quadratic regression
#' @param varformla Input argument for \code{bptest}. Default is quadratic regression.
#' @param tm Input vector of timpoints to be used if \code{regress}='time'
#' @return vector of p-values from Bruesh-Pagan heteroskedasticity tests.
#' @references http://biorxiv.org/content/early/2016/03/18/044495
#' @examples
#' library(plotrix)
#' library(WrightFisher)
#' 
#' set.seed(2)
#' 
#' #this may take ~20-30 seconds
#' X <- wfp(nSpecies=15,lambda=25)
#' 
#' t <- X$time
#' R <- X$Community
#' 
#' par(mfrow=c(1,3))
#' stackpoly(R,stack=T,main='Community Dynamics',xlab='time',ylab='Abundance')
#' 
#' Pvalsf <- CVTest(100,R)
#' Pvalst <- CVTest(100,R,regress='time')
#' 
#' plot(ecdf(Pvalsf),xlab='P',ylab='#',main='P-values from regress="f"')
#' lines(c(0,1),c(0,1),lwd=4,col='blue')
#' 
#' plot(ecdf(Pvalst),xlab='P',ylab='#',main='P-values from regress="time"')
#' lines(c(0,1),c(0,1),lwd=4,col='blue')
#' 
#' #Uncorrected KS test - the P-value from the KS-test gives a comparable estimate of the incompatibility of the simulation with a WFP
#' ks <- ks.test(Pvalsf,runif(length(Pvalsf)))
#' ks$p.value
#'
#' ks.test(Pvalst,runif(length(Pvalst)))$p.value
#' 
CVTest <- function(ntests,R,regress='f',formula=NULL,varformula=NULL,tm=NULL,...){
  
  if (!(regress %in% c('f','time'))){stop('unknown inputt "regress" - must be either "f" or "time"')}
  if (is.null(formula)){
    if (regress=='f'){
      formula=DF~f+I(f^2)
    } else 
      formula=DF~tm
  }
  if (is.null(varformula)){
    if (regress=='f'){
      varformula=DF~f+I(f^2)
    } else {
      varformula=DF~tm+I(tm^2)
      }
  }
  
  
  
  p <- rep(0,ntests)
  b <- vector('list',ntests)
  
  if (regress=='time'){
    if (is.null(tm)){
      tm <- 1:dim(R)[1]
      DT <- diff(tm)
      tm <- tm[1:(length(tm)-1)]
    }
  }
  
  for (nn in 1:ntests){
    
    cvt <- drawCVT(R)
    f <- cvt$f
    DF <- cvt$DF

    if (regress=='f'){
      
      dataset=model.frame(formula,environment())
      model <- lm(formula,dataset)
      p[nn] <- lmtest::bptest(model,varformula,...)$p.value
      DF <- model$residuals^2
      b[[nn]] <- glm(varformula)$coefficients
      
    } else {
      
      DF <- DF/sqrt(DT)
      
      dataset=model.frame(formula,environment())
      model <- lm(formula,dataset)
      p[nn] <- lmtest::bptest(model,varformula,...)$p.value
      b[[nn]] <- glm(varformula)$coefficients
    }
    
  }
  
  return(p)
}