## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)
library(doSNOW)

## Functions

# computing Tyler's shape matrix
TylerSigma = function(X, tol=1e-5, maxit=100, depth=F){
  n = nrow(X); p = ncol(X)
  iSig = diag(rep(1,p))
  
  # whether to use depth weights
  if(depth){
    mult = EPQD(X,X)[,p+1]
    mult = max(mult) - mult
  }
  else{
    mult = rep(1,n)
  }
  
  for(i in 1:maxit){
    iiSig = matrix(0,p,p)
    inv.iSig = solve(iSig)
    for(j in 1:n){
      xj = as.matrix(X[j,])
      iiSig = iiSig + mult[j]^2 * (xj %*% t(xj))/ as.numeric(t(xj) %*% inv.iSig %*% xj)
    }
    
    iiSig = iiSig/det(iiSig)^(1/p)
    if(norm(iSig - iiSig, type="F") < tol){
      break
    }
    else{
      iSig = iiSig
    }
  }
  
  iSig
}

## setup 1: Sigma = diag(1,1)
set.seed(02032015)

n = 1000
Sig = diag(c(2,1))
simiter = 100
err = matrix(0, nrow=simiter, ncol=2)

for(i in 1:simiter){
  iX = matrix(rnorm(2*n), ncol=2) %*% sqrt(Sig)
  
  err[i,1] = norm(Sig - TylerSigma(iX), type="F")
  err[i,2] = norm(Sig - TylerSigma(iX, depth=T), type="F")
}

colSums(err)/simiter
