## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)
library(doSNOW)

## Functions

## setup 1: Sigma = diag(1,2)
FSE.norm = function(n, p, iter=1e3){
  set.seed(12182014)
  v = c(rep(0,p-1), 1)
  lam = 1:p
  Sigma = diag(lam)
  
  MSE.mat = matrix(0, nrow=length(n), ncol=2)
  for(i in 1:length(n)){
    
    # function to compute stuff 1000 times for a given n
    loopfun = function(j){
      #source('misc_functions.R')
      iv = rep(0,3)
      
      # get sample and construct sign matrix
      iX = my.mvrnorm(n, mu=c(0,p), Sig=Sigma)
      iXnorm = sqrt(iX^2 %*% rep(1,p))
      iS = iX / (iXnorm %*% rep(1,p))
      
      # PCA on original sample
      iP = princomp(iX)
      iv[1] = abs(sum(v * iP$loadings[,1]))
      
      # PCA on SCM
      iPsign = princomp(iS)
      iv[2] = abs(sum(v * iPsign$loadings[,1]))
      
      # PCA on depth-CM
      idep = EPQD(iX,iX)[,p+1]
      iXd = iS * (max(idep) - idep)
      iPdepth = princomp(iXd)
      iv[3] = abs(sum(v * iPdepth$loadings[,1]))
      
      iv
    }
    
    # parallel code: compute MSE elements 1000 times
    cl = makeCluster(detectCores())
    registerDoSNOW(cl)
    system.time(eff.v <- foreach(j=1:iter) %dopar% loopfun(j))
    stopCluster(cl)
    
    # get MSE and return
    eff.v = matrix(unlist(eff.v), ncol=3, byrow=T)
    (MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
    MSE.mat[i,] = MSE.vec[2:3]/MSE.vec[1]
  }
  
  MSE.mat
}

FSE.norm(c(20,50,100,300), 2, 100)

