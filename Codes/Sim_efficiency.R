## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

library(parallel)
library(doSNOW)

## Functions
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## function to calculate weighted projection quantile depth
EPQD = function(X, grid, nu=1e3){
  
  p = ncol(X)
  b = apply(X, 2, median)
  X0 = X - ones(nrow(X),1) %*% b
  grid0 = grid - ones(nrow(grid),1) %*% b
  
  ## get matrix of weighted PQDs for all points
  npt = dim(grid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(p)); u = u/sqrt(sum(u^2))
    uecdf = ecdf(X0%*%u)
    Fuxu.mat[,iu] = uecdf(grid0%*%u)
  }
  EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  
  return(cbind(grid,EPQD.vec))
  
}


## setup 1: Sigma = diag(2,1)
require(mvtnorm)
set.seed(12182014)
v = c(1,0)
lam = c(2,1)
Sigma = diag(lam)
n = c(20,50,100,300)
iter = 1e3

loopfun = function(i){
  require(mvtnorm)
  iv = rep(0,3)
  
  # get sample and construct sign matrix
  #iX = my.mvrnorm(n[3], mu=c(0,0), Sig=Sigma)
  iX = rmvt(n[2], sigma=Sigma, df=5)
  iXnorm = sqrt(iX^2 %*% rep(1,2))
  iS = iX / (iXnorm %*% rep(1,2))
  
  # PCA on original sample
  iP = princomp(iX)
  iv[1] = abs(sum(v * iP$loadings[,1]))
  
  # PCA on SCM
  iPsign = princomp(iS)
  iv[2] = abs(sum(v * iPsign$loadings[,1]))
  
  # PCA on depth-CM
  idep = EPQD(iX,iX)[,3]
  iXd = iS * (max(idep) - idep)
  iPdepth = princomp(iXd)
  iv[3] = abs(sum(v * iPdepth$loadings[,1]))
  
  iv
}

set.seed(12182014)
cl = makeCluster(detectCores())
registerDoSNOW(cl)
system.time(eff.v <- foreach(i=1:iter) %dopar% loopfun(i))
stopCluster(cl)

eff.v = matrix(unlist(eff.v), ncol=3, byrow=T)
(MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
(eff.vec = MSE.vec[1]/MSE.vec[2:3])
