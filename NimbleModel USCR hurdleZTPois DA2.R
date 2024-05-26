NimModel <- nimbleCode({
  # priors
  lambda.N ~ dunif(0,1000) #expected abundance
  logit(p0) ~ dlogis(0,1) #baseline detection probability on logit scale
  # sigma ~ dunif(0,10) #detection spatial scale parameter, uninformative
  sigma ~ dgamma(20,40) #informative around 0.5
  lambda ~ dunif(0,10) #parameter for counts|detection
  N ~ dpois(lambda.N) #realized abundance
  #data augmentation "under the hood", jointly update N/z, 
  #no distribution induced on z, just turns obsmod on/off, used in y.true/ID update
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y.true[i,1:J,1:K] ~ dHurdleZTPoisMatrix(pd=pd[i,1:J],K2D=K2D[1:J,1:K],z=z[i],lambda=lambda) #vectorized obs mod
  }
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J,1:K])
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples])
})#model
