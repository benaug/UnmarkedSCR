NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  p0 ~ dunif(0,1)
  # sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y.true2D[i,1:J] ~ dBinVector(pd=pd[i,1:J],K1D[1:J],z=z[i]) #vectorized obs mod
  }
  capcounts <- Getcapcounts(y.true2D=y.true2D[1:M,1:J],y.true3D=y.true3D[1:M,1:J,1:K])
  N <- sum(z[1:M])
})#model