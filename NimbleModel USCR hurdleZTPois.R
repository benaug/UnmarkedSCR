NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  p0 ~ dunif(0,1)
  lambda ~ dunif(0,10)
  sigma ~ dunif(0,100) #uninformative prior
  # sigma ~ dgamma(25,scale=0.02) #informative prior with mean 0.5
  #--------------------------------------------------------------
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y.true[i,1:J,1:K] ~ dHurdleZTPoisVector(pd=pd[i,1:J],K2D=K2D[1:J,1:K],z=z[i],lambda=lambda) #vectorized obs mod
  }
  #calculate number of inds captured
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  n <- Getncap(capcounts=capcounts[1:M])
  N <- sum(z[1:M])
})