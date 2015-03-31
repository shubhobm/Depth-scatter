## Comm_analyze: analysis of communities data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")
commdata = read.csv("../Data/communities.data.txt", header=F)
NAcols = c(102:118, 122:125, 127) # columns with most NA entries
commdata.X = commdata[,-c(1:5, NAcols)]

# convert all columns into numeric
for(i in 1:ncol(commdata.X)){
  if(class(commdata.X[,i]) != "numeric"){
    commdata.X[,i] = as.numeric(paste(commdata.X[,i]))
  }
}

Y = commdata.X[,101]+1
commdata.X = data.frame(scale(commdata.X[,-101]))
Y = Y[complete.cases(commdata.X)]
commdata.X = commdata.X[complete.cases(commdata.X),]

## wil do vanilla and rank PCA, then
## principal components regression on increasing number of PCs

## vanilla PCA
pcamod = princomp(commdata.X)
scores = pcamod$scores

## get ranks
n = nrow(commdata.X)
p = ncol(commdata.X)
norms = sqrt(commdata.X^2 %*% rep(1,p))
signs = commdata.X / (norms %*% rep(1,p))

# calculate depth
require(fda.usc)
depths = mdepth.RP(commdata.X, commdata.X)$dep
depths = max(depths) - depths
commdata.rank = scale(signs * depths)

# rank PCA
pca.rank = princomp(commdata.rank)
scores.rank = pca.rank$scores

rsq.mat = matrix(0, nrow=p, ncol=2)
for(npc in 1:p){
  # PC regression on vanilla PCA scores
  lm.PC = lm(1/Y^3~scores[,1:npc]) # inverse cube transformation chosen by boxcox
  
  # PC-regression on rank PCA scores
  lm.PCrank = lm(1/Y^3~scores.rank[,1:npc])
  
  rsq.mat[npc,] = c(summary(lm.PC)$r.squared, summary(lm.PCrank)$r.squared)
}

plot(1:p, rsq.mat[,1], type='l', ylim=c(0.2,.8), lwd=2)
lines(1:p, rsq.mat[,2], type='l', lty=2, lwd=2)

