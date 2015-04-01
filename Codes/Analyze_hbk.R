## Analyze_octane: analysis of octane data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")

## Octane data
library(ellipse)
library(cepp)
# data(Colon)
# data.X = scale(Colon$X[which(Colon$Y==2),])

library(robustbase)
hbk.X = hbk[,-4]
data.X = scale(stack.x)
n = nrow(data.X)
p = ncol(data.X)

distanceplot = function(data.X, npc, ...){
  svd.oct = svd(data.X)
  P.oct = svd.oct$v[,1:npc]
  scores.oct = as.matrix(data.X) %*% P.oct
  eig.oct = svd.oct$d[1:npc]
  
  # calculate score distance and orthogonal distance
  SD = sqrt(rowSums(scores.oct^2 / matrix(eig.oct, nrow=n, ncol=npc, byrow=T)))
  OD = sqrt(rowSums((data.X - t(P.oct %*% t(scores.oct)))^2))
  
  # get cutoffs
  SD.cutoff = sqrt(qchisq(0.975, 2))
  ODt = OD^(2/3)
  t = mean(ODt); s = sd(ODt)
  OD.cutoff = (t+s*qnorm(.975))^(3/2)
  
  indices = 1:n
  which.ind = which(SD > SD.cutoff | OD > OD.cutoff)
  
  ## distance-distance plots
  par(mfrow=c(1,2))
  
  plot(SD, OD, ...)
  abline(v=SD.cutoff, col="red")
  abline(h=OD.cutoff, col="red")
  if(length(which.ind>0)){
    text(SD[which.ind], OD[which.ind], indices[which.ind], pos=1)
  }
  
  
  # get ranks
  norms = sqrt(data.X^2 %*% rep(1,p))
  signs = data.X / (norms %*% rep(1,p))
  
  # calculate depth
  require(fda.usc)
  depths = mdepth.RP(data.X, data.X, proj=2000)$dep
  depths = max(depths) - depths
  data.rank = scale(signs * depths)
  
  svd.oct.rank = svd(data.rank)
  P.oct.rank = svd.oct.rank$v[,1:npc]
  scores.oct.rank = as.matrix(data.rank) %*% P.oct.rank
  eig.oct.rank = svd.oct.rank$d[1:npc]
  
  # calculate score distance and orthogonal distance
  SDrank = sqrt(rowSums(scores.oct.rank^2 / matrix(eig.oct.rank, nrow=n, ncol=npc, byrow=T)))
  ODrank = sqrt(rowSums((data.rank - t(P.oct.rank %*% t(scores.oct.rank)))^2))
  
  # get cutoffs
  SDrank.cutoff = sqrt(qchisq(0.975, 2))
  ODtrank = ODrank^(2/3)
  trank = mean(ODtrank); srank = sd(ODtrank)
  ODrank.cutoff = (trank+srank*qnorm(.975))^(3/2)
  
  indices = 1:n
  which.ind = which(SDrank > SDrank.cutoff | ODrank > ODrank.cutoff)
  
  
  plot(SDrank, ODrank, ...)
  abline(v=SDrank.cutoff, col="red")
  abline(h=ODrank.cutoff, col="red")
  if(length(which.ind>0)){
    text(SDrank[which.ind], ODrank[which.ind], indices[which.ind], pos=1)
  }
  
  par(mfrow=c(1,1))
  
}

scoreplot = function(data.X, npc, ...){
  svd.oct = svd(data.X)
  P.oct = svd.oct$v[,1:npc]
  scores.oct = as.matrix(data.X) %*% P.oct
  
  ## score plot
  par(mfrow=c(1,2))
  
  plot(scores.oct, ...)
  lines(ellipse(cov(scores.oct), level=.975), lwd=2)
  
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
  plot(scores.oct.rank, ...)
  lines(ellipse(cov(scores.oct.rank), level=.975), lwd=2)
  
  par(mfrow=c(1,1))  
}


distanceplot(data.X, 2, pch=19, col="blue")
scoreplot(data.X, 2, pch=19, cex=.7, col="blue")
