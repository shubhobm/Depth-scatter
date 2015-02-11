## misc_functions: miscellaneous functions for ease of use
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

# function giving matrix of ones
ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## calculate weighted projection quantile depth
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
  #EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  EPQD.vec = apply(1-Fuxu.mat, 1, min)
  
  return(cbind(grid,EPQD.vec))
  
}

# compute Tyler's shape matrix, optionally with depth weights
TylerSig = function(X, tol=1e-5, maxit=100, weight=NULL){
  n = nrow(X); p = ncol(X)
  iSig = diag(rep(1,p))
  
  # whether to use depth weights
  if(is.null(weight)){
    weight = rep(1,n)
  }
  
  for(i in 1:maxit){
    iiSig = matrix(0,p,p)
    inv.iSig = solve(iSig)
    for(j in 1:n){
      xj = as.matrix(X[j,])
      iiSig = iiSig + weight[j]^2 * (xj %*% t(xj))/ as.numeric(t(xj) %*% inv.iSig %*% xj)
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