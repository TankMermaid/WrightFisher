#' Performs Washburne et al.'s (2016) Neutrality Test with corrected ks-test
#' 
#' 

WashburneNTest <- function(R,ntests=min(16000,2^dim(X)[2]),regress='f',ncores=NULL,formula=DF~f+I(f^2),varformula=~f+I(f^2),...){
  
  n <- dim(R)[2]
  m <- dim(R)[1]
   
  if (is.null(ncores)){
  Pvals <- CVTest(ntests,R,regress,formula,varformula,...) 
  } else{
    cl <- makeCluster(ncores)
    clusterEvalQ(cl,library(WrightFisher))
    clusterEvalQ(cl,library(lmtest))
    
    testsEach <- round(ntests/ncores)
    
    Pvals <- unlist(parLapply(cl,X=as.list(rep(testsEach,ncores)),fun = CVTest,R=R,regress=regress,formula=formula,varformula=varformula))
  }
  
}