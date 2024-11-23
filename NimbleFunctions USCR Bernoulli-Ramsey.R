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
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)
Getcapcounts <- nimbleFunction(
  run = function(y.true2D=double(2),y.true3D=double(3)){ #dummy y.true3D so i can use in ID update
    returnType(double(0))
    capcounts <- sum(y.true2D)
    return(capcounts)
  }
)

## sampler to update y[1:M,1:J]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    n.samples <- control$n.samples
    M <- control$M
    J <- control$J
    K <- control$K
    this.j <- control$this.j
    this.k <- control$this.k
    y.jk <- control$y.jk
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    #do y update with 3D data, then reduce to 2D ot feed back to nimble
    z <- model$z
    y.true3D <- model$y.true3D
    n.idx <- length(this.j)

    # precalculate log likelihoods at individual by trap by occasion level
    ll.y <- array(0,dim=c(M,J,K))
    for(i in 1:M){
      if(z[i]==1){
        for(l in 1:n.idx){
          ll.y[i,this.j[l],this.k[l]] <-  dbinom(y.true3D[i,this.j[l],this.k[l]],size=1,prob=model$pd[i,this.j[l]], log = TRUE)
        }
      }
    }
    
    ll.y.cand <- ll.y
    y.true3D.cand <- y.true3D
    for(l in 1:n.idx){#for each j-k with a detection
      propprobs <- model$pd[1:M,this.j[l]]*z[1:M]
      #propose full i x j x k candidate
      y.true3D.cand[,this.j[l],this.k[l]] <- rbinom(M,1,propprobs)
      if(sum(y.true3D.cand[,this.j[l],this.k[l]])>=y.jk[this.j[l],this.k[l]]){ #always reject if 3D candidate inconsistent with observed data
        #proposal probs and likelihood use 3D data
        log.prop.for <- sum(dbinom(y.true3D.cand[,this.j[l],this.k[l]],1,propprobs,log=TRUE))
        log.prop.back <- sum(dbinom(y.true3D[,this.j[l],this.k[l]],1,propprobs,log=TRUE))
        ll.y.cand[,this.j[l],this.k[l]] <- dbinom(y.true3D.cand[,this.j[l],this.k[l]],1,prob=model$pd[,this.j[l]], log = TRUE)
        #sum log likelihoods and do MH step
        lp_initial <- sum(ll.y[,this.j[l],this.k[l]])
        lp_proposed <- sum(ll.y.cand[,this.j[l],this.k[l]])
        log_MH_ratio <- (lp_proposed+log.prop.back) - (lp_initial+log.prop.for)
        accept <- decide(log_MH_ratio)
        #proposal is full conditional. except when proposal invalid (propose too few detections)... switch to gibbs?
        if(accept){
          y.true3D[,this.j[l],this.k[l]] <- y.true3D.cand[,this.j[l],this.k[l]]
          ll.y[,this.j[l],this.k[l]] <- ll.y.cand[,this.j[l],this.k[l]]
        }
      }
    }
    y.true2D <- matrix(0,M,J)
    for(i in 1:M){
      for(j in 1:J){
        y.true2D[i,j] <- sum(y.true3D[i,j,])
      }
    }
    #put everything back into the model$stuff
    model$y.true2D <<- y.true2D
    model$y.true3D <<- y.true3D
    model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
