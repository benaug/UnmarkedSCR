NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi~dunif(0,1)
  p0~dunif(0,1)
  # sigma~dunif(0,100)
  #--------------------------------------------------------------
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
  }
  for(j in 1:J){
    pd.j[j] <- 1-(prod(1-pd[1:M,j]*z[1:M]))
    y.j[j] ~ dbinom(pd.j[j],K1D[j])
  }
  N <- sum(z[1:M])
})#model