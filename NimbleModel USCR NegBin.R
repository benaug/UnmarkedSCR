NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  theta ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  # sigma ~ dunif(0,100)
  sigma ~ dgamma(25,scale=0.02) #informative prior with mean 0.5
  #--------------------------------------------------------------
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    p[i,1:J] <- theta/(theta+lam[i,1:J])
    y.true[i,1:J] ~ dNBVector(p=p[i,1:J],theta=theta*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of NB RVs is NB with theta=N*theta
  }
  #calculate number of inds captured
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  n <- Getncap(capcounts=capcounts[1:M])
  N <- sum(z[1:M])
})#model