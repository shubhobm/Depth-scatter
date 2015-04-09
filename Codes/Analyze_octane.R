## Analyze_octane: analysis of octane data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")
source('misc_functions.R')

## Octane data
require(ellipse)
require(ChemoSpec)
require(rrcov)
Octane <- read.csv("../Data/Octane.csv")
Octane.X = Octane[,-c(1,2)]
n = nrow(Octane.X)
p = ncol(Octane.X)

distanceplot = function(data.X, npc, ...){
  # vanila PCA
  pcamod = PcaClassic(Octane.X, k=2)
  
  sd = pcamod@sd
  od = pcamod@od
  csd = pcamod@cutoff.sd
  cod = pcamod@cutoff.od
  indices = 1:n
  which.ind = which(sd > csd | od > cod)
  
  ## distance-distance plots
  par(mfrow=c(1,2))
  
  plot(sd, od,
       main="Classical PCA", xlab="Score distance", ylab="Orthogonal distance", ...)
  abline(v=csd, col="red")
  abline(h=cod, col="red")
  if(length(which.ind>0)){
    text(sd[which.ind], od[which.ind], indices[which.ind], pos=1, cex=.7)
  }
  
  # rank PCA
  pcarank = PcaRank(Octane.X, k=2)
  
  sdrank = pcarank@sd
  odrank = pcarank@od
  csdrank = pcarank@cutoff.sd
  codrank = pcarank@cutoff.od
  indices = 1:n
  which.ind = which(sdrank > csdrank | odrank > codrank)
  
  ## distance-distance plots
  plot(sdrank, odrank,
       main="Depth PCA", xlab="Score distance", ylab="Orthogonal distance", ...)
  abline(v=csdrank, col="red")
  abline(h=codrank, col="red")
  if(length(which.ind>0)){
    text(sdrank[which.ind], odrank[which.ind], indices[which.ind], pos=1, cex=.7)
  }  
  par(mfrow=c(1,1))  
}

scoreplot = function(data.X, npc){
  svd.oct = svd(data.X)
  P.oct = svd.oct$v[,1:npc]
  scores.oct = as.matrix(data.X) %*% P.oct
  
  ## score plot
  par(mfrow=c(1,2))
  
  plot(scores.oct, pch=19, cex=.7)
  lines(ellipse(cov(scores.oct)), lwd=2)
  
  # get ranks
  norms = sqrt(data.X^2 %*% rep(1,p))
  signs = data.X / (norms %*% rep(1,p))
  
  # calculate depth
  require(fda.usc)
  depths = mdepth.RP(data.X, data.X)$dep
  depths = max(depths) - depths
  data.rank = scale(signs * depths)
  
  svd.oct.rank = svd(data.rank)
  P.oct.rank = svd.oct.rank$v[,1:npc]
  scores.oct.rank = as.matrix(data.rank) %*% P.oct.rank
  
  ## score plot
  plot(scores.oct.rank, pch=19, cex=.7)
  lines(ellipse(cov(scores.oct.rank)), lwd=2)
  
  par(mfrow=c(1,1))
  
}


distanceplot(Octane.X, 2, xlim=c(0,10), ylim=c(0,1.3), pch=19, col="blue")
scoreplot(Octane.X, 2)
