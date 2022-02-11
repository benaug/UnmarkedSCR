#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBinVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D=double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, prob = pd, size = K1D, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1),K1D=double(1), z = double(0)) {
    returnType(double(1))
    J=nimDim(pd)[1]
    out=numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
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
    this.k<-control$this.k
    n.samples <- control$n.samples
    M<-control$M
    J <- control$J
    calcNodes <- model$getDependencies(target)
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
          ll.y[i,j] <-  dbinom(y.true[i,j],size=model$K1D[j],prob=model$pd[i,j], log = TRUE)
        }
      }
    }
    ll.y.cand <- ll.y
    
    for(l in 1:n.samples){#loop over samples
      ID.cand=ID.curr
      y.true.cand=y.true
      #proposal distribution
      propprobs=model$pd[1:M,this.j[l]]*z[1:M]
      backprobs=propprobs #store for faster calculation below
      #remove matches that are illegal under bernoulli obsmod. will remove self matches
      rem.idx=which(this.j==this.j[l]&this.k==this.k[l])
      for(i in 1:length(rem.idx)){
        propprobs[ID.curr[rem.idx[i]]]=0
      }
      sumpropprobs=sum(propprobs)
      if(sumpropprobs>0){ #abort here to avoid dividing by 0. I don't think this will actually every happen, though. Can always propose to yourself.
        propprobs <- propprobs/sumpropprobs
        ID.cand[l]=rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){ #bernoulli exclusions will prevent selecting current ID.. leaving in anyway
          swapped=c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob=propprobs[swapped[2]]
          #calculate backwards propprobs
          for(i in 1:length(rem.idx)){
            backprobs[ID.cand[rem.idx[i]]]=0
          }
          backprobs=backprobs/sum(backprobs)
          backprob=backprobs[swapped[1]]
          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l]]=y.true[ID.curr[l],this.j[l]]-1
          y.true.cand[ID.cand[l],this.j[l]]=y.true[ID.cand[l],this.j[l]]+1
          
          #if theta is a function of individual or trap covariates, need to fix these 2 lines below, e.g., size=model$theta[swapped[1],this.j[l]]*model$K1D[this.j[l]]
          ll.y.cand[swapped[1],this.j[l]]=dbinom(y.true.cand[swapped[1],this.j[l]],size=model$K1D[this.j[l]],prob=model$pd[swapped[1],this.j[l]],log=TRUE)
          ll.y.cand[swapped[2],this.j[l]]=dbinom(y.true.cand[swapped[2],this.j[l]],size=model$K1D[this.j[l]],prob=model$pd[swapped[2],this.j[l]],log=TRUE)
          #select sample to move proposal probabilities
          #P(select a sample for this ID)*P(select this j|this ID)
          #n.samples cancels out in MH ratio. Including for clarity
          focalprob=(sum(ID.curr==swapped[1])/n.samples)*(y.true[swapped[1],this.j[l]])/sum(y.true[swapped[1],1:J])
          focalbackprob=(sum(ID.cand==swapped[2])/n.samples)*(y.true.cand[swapped[2],this.j[l]])/sum(y.true.cand[swapped[2],1:J])
          
          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped,this.j[l]])
          lp_proposed <- sum(ll.y.cand[swapped,this.j[l]])
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l]]=y.true.cand[swapped[1],this.j[l]]
            y.true[swapped[2],this.j[l]]=y.true.cand[swapped[2],this.j[l]]
            ll.y[swapped[1],this.j[l]]=ll.y.cand[swapped[1],this.j[l]]
            ll.y[swapped[2],this.j[l]]=ll.y.cand[swapped[2],this.j[l]]
            ID.curr[l]=ID.cand[l]
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