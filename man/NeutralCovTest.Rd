% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeutralCovTest.R
\name{NeutralCovTest}
\alias{NeutralCovTest}
\title{Performs Washburne et al.'s (2016) Neutrality Test with corrected ks-test}
\usage{
NeutralCovTest(R, tm = NULL, ntests = NULL, regress = "f",
  ncores = NULL, formula = DF ~ f + I(f^2), varformula = NULL,
  standard = T, method = "logitnorm", ...)
}
\arguments{
\item{R}{Compositional matrix whose columns are species and whose rows are timepoints}

\item{tm}{Optional time-points. Use if time intervals between samples are not uniform.}

\item{ntests}{Number of constant-volatility transforms to use in test. Must be less than or equal to 2^dim(R)[2]}

\item{regress}{String, either 'f', or 'tm', based on whether heteroskedasticity should regress against the state, f, or the time, tm. Default is 'f'}

\item{ncores}{Number of cores for parallelization of constant-volatility tests.}

\item{formula}{Formula for regression of drift. Default DF~f+I(f^2)}

\item{varformula}{Formula for regression of residuals in Breush-Pagan test. Default DF~f+I(f^2)}

\item{standard}{Logical indicating whether to use a standard, reproducible choice of CVTs by setting seed within NeutralCovTest}

\item{method}{String, either "Kolmogorov" or "logitnorm", for approximation of KS-statistic from P-value distribution of CVTests against uniform. Default is logitnorm.}

\item{...}{optional input arguments for bptest}
}
\description{
Performs Washburne et al.'s (2016) Neutrality Test with corrected ks-test
}
\examples{
library(WrightFisher)
library(plotrix)

set.seed(10)
# For a simple test, we can produce a dataset:
R <- wfp(nspecies=12,dt=1e-3)$Community
Rg <- gbmMR(nspecies=12,dt=1e-3,sigma=3,mu=2)$Shares

# Visualizing Market/Community Dynamics
par(mfrow=c(2,2))
stackpoly(R,stack=T,main='Wright-Fisher Process',ylab='Relative Abundance',xlab='time')
stackpoly(Rg,stack=T,main='Mean-Reverting GBM',ylab='Relative Abundance',xlab='time') 

# Neutrality Covariance Testing
P_wfp <- NeutralCovTest(R)
P_gbm <- NeutralCovTest(Rg)
P_wfp
P_gbm

# In the script below, we'll simulate 'reps' neutral communities through simulation of a Volatility-stabilized market, and then
# we'll perform a NeutralCovTest on each, using 1000 constant-volatility transformations. If performed on a computer with 4 or more
# cores, the parallelized version will be run for a speed-up.  
 
reps <- 20
nspecies=50
Pvals_vsm <- rep(NA,reps)
Pvals_gbm <- rep(NA,reps)
parallelize=T

for (nn in 1:reps){

## Simulating Communities
R <- vsm(nspecies=nspecies,lambda=25,dt=1e-4,Tmax=1,nsamples=125)$Shares
R2 <- gbmMR(nspecies=nspecies,Tmax=1,nsamples=125,sigma=4,mu=15)$Shares

## Neutrality Testing
if (detectCores()>=3 && parallelize==T){
# For species-rich datasets, NeutralCovTest may take a while and can be sped-up through paralellization by inputting ncores
 ncores=detectCores()-1
 Pvals_vsm[nn] <- NeutralCovTest(R,ncores=ncores)
 Pvals_gbm[nn] <- NeutralCovTest(R2,ncores=ncores)
} else {
 Pvals_vsm[nn] <- NeutralCovTest(R)
 Pvals_gbm[nn] <- NeutralCovTest(R2)
}

}

plot(ecdf(Pvals_vsm),xlab='NeutralCovTest Pvalue',ylab='F(P)',main='Pvalues: NeutralCovTest on VSMs',xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col='blue',lwd=2)
legend(.6,.4,legend=c('Simulations','Uniform'),pch=c(16,NA),lwd=c(1,2),col=c('black','blue'))
plot(ecdf(Pvals_gbm),xlab='Pvalue',ylab='F(P)',main='Pvalues: NeutralCovTests on mean-reverting GBMs',xlim=c(0,1),ylim=c(0,1))
lines(c(0,1),c(0,1),col='blue',lwd=2)
legend(.6,.4,legend=c('Simulations','Uniform'),pch=c(16,NA),lwd=c(1,2),col=c('black','blue'))
}

