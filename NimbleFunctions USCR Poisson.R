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

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J=nimDim(lam)[1]
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
    trapup <- control$trapup
    M<-control$M
    J <- control$J
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    z <- model$z
    for(j in 1:length(trapup)){ #only iterate through traps with samples
      lam.curr <- model$lam[,trapup[j]] #individual by trap expected counts
      lam.curr[z==0] <- 0 #can zero out z==0 here, will already be zeroed out using supplied nimble functions
      lam.use <- lam.curr
      sum.lam.use=sum(lam.use)
      if(sum.lam.use>0){ #abort if 0. Should never happen?
        #full conditional for identity update at this j
        fullcond <- lam.use/sum.lam.use
        idx <- which(this.j==trapup[j]) #which samples were recorded at this trap?
        for(l in 1:length(idx)){ #update the individual identities of these samples one by one
          #update derived parameter ID
          ID.curr <- model$ID[idx[l]]
          ID.prop <- rcat(1,fullcond)
          model$ID[idx[l]] <<- ID.prop
          #update y
          model$y.true[ID.curr,trapup[j]] <<- model$y.true[ID.curr,trapup[j]]-1 #subtract out old ID obs
          model$y.true[ID.prop,trapup[j]] <<- model$y.true[ID.prop,trapup[j]]+1 #add in new ID obs
        }
      }
    }
    #update lp
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)