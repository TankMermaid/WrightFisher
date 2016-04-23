#' Draws a WFP constant-volatility transformation
#' 
#' @export
#' @param R Compositional matrix. See \code{\link{CVTest}}
#' @param n number of species
#' @param m number of timepoints
#' @return two member list, "f" - the CVT evaluated at all points in the input matrix - and "DF", the jumps in the CVT
drawCVT <- function(R,n=dim(R)[2],m=dim(R)[1]){
  dum=T
  while(dum==T){
    a <- sign(rnorm(n))
    if (abs(sum(a)<n)){dum=F}
  }
  
  vec <- R %*% a
  vec[vec < -1]=-1
  vec[vec > 1]=1
  output <- vector('list',2)
  output[[1]] <- asin(vec)
  output[[2]] <- diff(output[[1]])
  output[[1]] <- output[[1]][1:(m-1)]
  names(output) <- c('f','DF')
  return(output)
}
