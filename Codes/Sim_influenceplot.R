## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)
library(doSNOW)

## Scatterplot of data and D-rank
n = 1e3
set.seed(120214)
Gamma = matrix(c(1,-1,1,1), nrow=2)/sqrt(2)
sig = Gamma %*% diag(c(9,1)) %*% t(Gamma)
X = my.mvrnorm(n, mu=c(0,0), Sig=2*sig)
uX = X / sqrt(rowSums(X^2))
dX = EPQD(X, X)[,3]
Xrank = uX * (max(dX) - dX)

par(mfrow=c(1,2))
plot(X, pch=19, cex=.5)
plot(Xrank, pch=19, cex=.5)
par(mfrow=c(1,1))


## Influence function plots
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

default = par()
par(mfrow=c(2,2), mai=rep(.5,4))
# Sample covariance matrix
r = sqrt(rowSums(xygrid^2))
Ugrid = xygrid / r
IFnorm = abs(Ugrid[,1] * Ugrid[,2] * sqrt(lam[1]*lam[2])/(lam[1] - lam[2]))

IFnorm.Sig = r^2*IFnorm
persp(pts, pts, matrix(IFnorm.Sig, nrow=lengrid, byrow=T),
      main="(a)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for SCM
# get eigenvalues of SCM
lam.Zsq = (Z * Z) %*% Sig
sum.lam.Zsq = rowSums(lam.Zsq)
lamS = colMeans(lam.Zsq / sum.lam.Zsq)

# get norms of influence fns for eigenvectors
multS = sqrt(lam[1]*lam[2])/(lamS[1] - lamS[2])
IFnorm.S = sqrt(abs(multS * Ugrid[,1] * Ugrid[,2]))

# plot result
persp(pts, pts, matrix(IFnorm.S, nrow=lengrid, byrow=T),
      main="(b)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for Tyler's scatter matrix
IFnorm.tyler = 4*IFnorm
persp(pts, pts, matrix(IFnorm.tyler, nrow=lengrid, byrow=T),
      main="(c)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for DCM
# get eigenvalues of DCM
DZ = EPQD(Z, Z)[,3]
DZ = max(Z) - Z
lamDS = colMeans((DZ^2/sum.lam.Zsq) * lam.Zsq)

# get htped at grid points
Dgrid = EPQD(X, xygrid)[,3]
Dgrid = max(Dgrid) - Dgrid

# get norms of influence fns for eigenvectors
mult = sqrt(lam[1]*lam[2])/(lamDS[1] - lamDS[2])
IFnorm.D = abs(mult * xygrid[,1] * xygrid[,2] * Dgrid^2 / diag(xygrid %*% Sig %*% t(xygrid)))

# plot result
persp(pts, pts, matrix(IFnorm.D, nrow=lengrid, byrow=T),
      main="(d)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

par(default)
