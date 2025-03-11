#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dNBVector <- nimbleFunction(
  run = function(x = double(1), p = double(1), theta = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dnbinom(x, p = p, size = theta, log = TRUE))
      return(logProb)
    }
  }
)

#dummy random vector generator to make nimble happy
rNBVector <- nimbleFunction(
  run = function(n = integer(0),p = double(1),theta = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(p)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(ID=double(1),M=double(0)){
    returnType(double(1))
    n.samples <- nimDim(ID)[1]
    capcounts <- numeric(M, value = 0)
    for(l in 1:n.samples){
      capcounts[ID[l]] <- capcounts[ID[l]] + 1
    }
    return(capcounts)
  }
)

Getncap <- nimbleFunction(
  run = function(capcounts=double(1)){
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

## sampler to update y[1:M,1:J]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j<-control$this.j
    M<-control$M
    J <- control$J
    n.samples <- control$n.samples
    calcNodes <- model$getDependencies(c("y.true","ID"))
  },
  
  run = function() {
    z <- model$z
    y.true <- model$y.true
    ID.curr <- model$ID
    
    #precalculate log likelihoods at individual by trap level
    ll.y <- matrix(0,nrow=M,ncol=J)
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          #if theta is a function of individual or trap covariates, need to fix the line below, e.g., size=model$theta[i,j]*model$K1D[j]
          ll.y[i,j] <-  dnbinom(y.true[i,j],size=model$theta[1]*model$K1D[j],prob=model$p[i,j], log = TRUE)
        }
      }
    }
    ll.y.cand <- ll.y
    ID.cand <- ID.curr
    y.true.cand <- y.true
    for(l in 1:n.samples){#loop over samples
      #proposal distribution
      propprobs <- model$lam[1:M,this.j[l]]*z
      sumpropprobs <- sum(propprobs)
      if(sumpropprobs>0){ #abort here to avoid dividing by 0. I don't think this will actually every happen, though. Can always propose to yourself.
        propprobs <- propprobs/sumpropprobs
        ID.cand[l] <- rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l]] <- y.true[ID.curr[l],this.j[l]]-1
          y.true.cand[ID.cand[l],this.j[l]] <- y.true[ID.cand[l],this.j[l]]+1
          #if theta is a function of individual or trap covariates, need to fix these 2 lines below, e.g., size=model$theta[swapped[1],this.j[l]]*model$K1D[this.j[l]]
          ll.y.cand[swapped[1],this.j[l]] <- dnbinom(y.true.cand[swapped[1],this.j[l]],size=model$theta[1]*model$K1D[this.j[l]],prob=model$p[swapped[1],this.j[l]],log=TRUE)
          ll.y.cand[swapped[2],this.j[l]] <- dnbinom(y.true.cand[swapped[2],this.j[l]],size=model$theta[1]*model$K1D[this.j[l]],prob=model$p[swapped[2],this.j[l]],log=TRUE)
          #select sample to move proposal probabilities
          focalprob <- y.true[ID.curr[l],this.j[l],this.k[l]]/n.samples
          focalbackprob <- y.true.cand[ID.cand[l],this.j[l],this.k[l]]/n.samples

          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped,this.j[l]])
          lp_proposed <- sum(ll.y.cand[swapped,this.j[l]])
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l]] <- y.true.cand[swapped[1],this.j[l]]
            y.true[swapped[2],this.j[l]] <- y.true.cand[swapped[2],this.j[l]]
            ll.y[swapped[1],this.j[l]] <- ll.y.cand[swapped[1],this.j[l]]
            ll.y[swapped[2],this.j[l]] <- ll.y.cand[swapped[2],this.j[l]]
            ID.curr[l] <- ID.cand[l]
          }else{
            y.true.cand[swapped[1],this.j[l]] <- y.true[swapped[1],this.j[l]]
            y.true.cand[swapped[2],this.j[l]] <- y.true[swapped[2],this.j[l]]
            ID.cand[l] <- ID.curr[l]
          }
        }
      }
    }
    #put everything back into the model$stuff
    model$y.true <<- y.true
    model$ID <<- ID.curr
    model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)