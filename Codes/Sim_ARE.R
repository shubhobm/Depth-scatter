## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

library(parallel)
library(doSNOW)

# Functions
ARE.norm = function(rho, depth, ns=1e4){
  require(fda.usc)
  Z = matrix(rnorm(2*ns), ncol=2)
  
  # get depth values
  if(depth=='HS'){
    DZ = mdepth.HS(Z,Z)$dep
  }
  else if(depth=='MhD'){
    DZ = mdepth.MhD(Z,Z)$dep
  }
  else {
    DZ = mdepth.RP(Z,Z)$dep
  }
  DZ = mdepth.RP(Z,Z)$dep
  DZ = (max(DZ) - DZ)
  
  
  DZ=1
  
  # calculate required quantities
  Z2 = Z*Z
  Z.Sig.Z = Z2[,1] + rho*Z2[,2]
  (E1 = mean(DZ^4 * Z2[,1] * Z2[,2] / Z.Sig.Z^2))
  (E2 = mean(DZ^2 * (Z2[,1] - rho*Z2[,2]) / Z.Sig.Z)^2)
  ARE = E2 / (E1*(1-rho)^2)
  
  ARE  
}

ARE.norm(.5, 'HS', ns=1e3)
ARE.norm(.5, 'MhD', ns=1e4)
ARE.norm(.5, 'RP', ns=1e4)

# Z = matrix(rnorm(2*1e4), ncol=2)
# r = sqrt(rowSums(Z*Z))
# mean(r^4)

lamS1 = mean(DZ^2*Z2[,1]/Z.Sig.Z)
lamS2 = rho*mean(DZ^2*Z2[,2]/Z.Sig.Z)
(AV1 = rho/(lamS1-lamS2)^2 * mean(DZ^4*Z2[,1]*Z2[,2]/Z.Sig.Z^2))
(AV2 = rho/(1-rho)^2)
AV2/AV1
