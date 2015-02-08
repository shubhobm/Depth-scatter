## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)
library(doSNOW)

## Functions

## setup 1: Sigma = diag(2,1)
set.seed(02052015)
lam = c(2,1)
Sig = diag(lam)
Z = matrix(rnorm(1e3),ncol=2)
X = Z %*% sqrt(Sig)

# make grid of points
pts = seq(-2, 2, by=.1)
lengrid = length(pts)
xcoord = rep(pts, rep(lengrid,lengrid))
ycoord = rep(pts, lengrid)
xygrid = cbind(xcoord,ycoord)
rm(xcoord,ycoord)

# Sample covariance matrix
r = sqrt(rowSums(xygrid^2))
Ugrid = xygrid / r
IFnorm.S = r^2*abs(Ugrid[,1] * Ugrid[,2] * sqrt(lam[1]*lam[2])/(lam[1] - lam[2]))
persp(pts, pts, matrix(IFnorm.S, nrow=lengrid, byrow=T),
      xlab="x1", ylab="x2", zlab="IF(x0)", ticktype="detailed",
      theta=45, phi=45)

# Influence fn plot for DCM
# get eigenvalues of DCM
lam.Zsq = (Z * Z) %*% Sig
sum.lam.Zsq = rowSums(lam.Zsq)

DZ = EPQD(Z, Z)[,3]
DZ = max(Z) - Z
lamDS = colMeans((DZ^2/sum.lam.Zsq) * lam.Zsq)

# get htped at grid points
Dgrid = EPQD(X, xygrid)[,3]
Dgrid = max(Dgrid) - Dgrid

# get norms of influence fns for eigenvectors
mult = sqrt(lam[1]*lam[2])/(lamDS[1] - lamDS[2])
IFnorm = sqrt(abs(mult * xygrid[,1] * xygrid[,2] * Dgrid^2 / diag(xygrid %*% Sig %*% t(xygrid))))

# plot result
persp(pts, pts, matrix(IFnorm, nrow=lengrid, byrow=T),
      xlab="x1", ylab="x2", zlab="IF(x0)", ticktype="detailed",
      theta=45, phi=45)

# Influence fn plot for Tyler's scatter matrix
Ugrid = xygrid / sqrt(rowSums(xygrid^2))
IFnorm.tyler = 4*sqrt(abs(Ugrid[,1] * Ugrid[,2] * sqrt(lam[1]*lam[2])/(lam[1] - lam[2])))
persp(pts, pts, matrix(IFnorm.tyler, nrow=lengrid, byrow=T),
      xlab="x1", ylab="x2", zlab="IF(x0)", ticktype="detailed",
      theta=45, phi=45)
