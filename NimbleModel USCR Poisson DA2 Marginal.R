NimModel <- nimbleCode({
  # priors
  lambda.N ~ dunif(0,1000) #expected abundance
  lam0 ~ dunif(0,10) #baseline detection rate
  # sigma ~ dunif(0,10) #detection spatial scale parameter, uninformative
  sigma ~ dgamma(20,40) #informative around 0.5
  N ~ dpois(lambda.N) #realized abundance
  #data augmentation "under the hood", jointly update N/z, 
  #no distribution induced on z, just turns obsmod on/off, used in y.true/ID update
  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
  }
  #this speeds up lam0/sigma updates a bit by skipping z=0 inds
  bigLam[1:J] <- GetbigLam(lam=lam[1:M,1:J],z=z[1:M])
  lam.noID[1:J] <- bigLam[1:J]*K1D[1:J]
  y.noID[1:J] ~ dPoissonVector(lam.noID[1:J],z=1) #plug in z=1 to reuse dPoissonVector
})#model