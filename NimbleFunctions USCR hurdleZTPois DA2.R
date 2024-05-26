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

dHurdleZTPoisMatrix <- nimbleFunction(
  run = function(x = double(2), pd = double(1),K2D = double(2), z = double(0), lambda = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(x)[1]
      K <- nimDim(x)[2]
      logProb <- 0
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]==1){
            if(x[j,k]==0){
              logProb <- logProb + log(1-pd[j])
            }else{
              logProb <- logProb + log(pd[j]) + log(dpois(x[j,k],lambda=lambda)/(1-exp(-lambda)))
            }
          }
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rHurdleZTPoisMatrix <- nimbleFunction(
  run = function(n = integer(0), pd = double(1),K2D = double(2), z = double(0), lambda = double(0)) {
    returnType(double(2))
    J <- nimDim(pd)[1]
    K <- nimDim(pd)[2]
    out <- matrix(0,J,K)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(3)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    K <- nimDim(y.true)[3]
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:J,1:K])
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

# sampler to update y[1:M,1:J,1:K]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j <- control$this.j
    this.k <- control$this.k
    n.samples <- control$n.samples
    M <- control$M
    J <- control$J
    K <- control$K
    IDups <- control$IDups
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    z <- model$z
    y.true <- model$y.true
    ID.curr <- model$ID
    lambda <- model$lambda[1]
    
    #precalculate likelihoods at i x j x k level
    ll.y <- array(0,dim=c(M,J,K))
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(model$K2D[j,k]==1){
              if(y.true[i,j,k]==0){
                ll.y[i,j,k] <- log(1-model$pd[i,j])
              }else{
                #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                ll.y[i,j,k] <- log(model$pd[i,j])
                ll.y[i,j,k] <- ll.y[i,j,k] + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
              }
            }
          }
        }
      }
    }
    ll.y.cand <- ll.y
    ID.cand <- ID.curr
    y.true.cand <- y.true
    
    for(up in 1:IDups){
      for(l in 1:n.samples){ #loop over samples
        #get proposal distribution for sample l
        lp.y.prop <- rep(0,M)
        for(i in 1:M){
          if(z[i]==1){
            if(i!=ID.curr[l]){ #new state
              y.tmp1 <- y.true[i,this.j[l],this.k[l]] + 1 #if we add sample here
              y.tmp2 <- y.true[ID.curr[l],this.j[l],this.k[l]] - 1 #if we add sample here
              if(y.tmp1==0){ #if not captured
                lp.y.prop[i] <- dbinom(0,size=1,prob=model$pd[i,this.j[l]],log=TRUE)
              }else{ #if captured
                lp.y.prop[i] <- dbinom(1,size=1,prob=model$pd[i,this.j[l]],log=TRUE) + 
                  log(dpois(y.tmp1,lambda=lambda)/(1-dpois(0,lambda=lambda)))
              }
              if(y.tmp2==0){ #if not captured
                lp.y.prop[i] <- lp.y.prop[i] + dbinom(0,size=1,prob=model$pd[i,this.j[l]],log=TRUE)
              }else{ #if captured
                lp.y.prop[i] <- lp.y.prop[i] + dbinom(1,size=1,prob=model$pd[i,this.j[l]],log=TRUE) + 
                  log(dpois(y.tmp2,lambda=lambda)/(1-dpois(0,lambda=lambda)))
              }
            }else{ #can't propose this guy if z==0
              lp.y.prop[i] <- -Inf
            }
          }else{
            lp.y.prop[i] <- -Inf
          }
        }
        prop.probs <- exp(lp.y.prop)
        prop.probs <- prop.probs/sum(prop.probs)
        ID.cand[l] <- rcat(1,prob=prop.probs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])

          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l],this.k[l]] <- y.true[ID.curr[l],this.j[l],this.k[l]] - 1
          y.true.cand[ID.cand[l],this.j[l],this.k[l]] <- y.true[ID.cand[l],this.j[l],this.k[l]] + 1
          
          #update ll.y
          if(y.true.cand[swapped[1],this.j[l],this.k[l]]==0){
            ll.y.cand[swapped[1],this.j[l],this.k[l]] <- log(1-model$pd[swapped[1],this.j[l]])
          }else{
            #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
            ll.y.cand[swapped[1],this.j[l],this.k[l]] <- log(model$pd[swapped[1],this.j[l]])
            ll.y.cand[swapped[1],this.j[l],this.k[l]] <- ll.y.cand[swapped[1],this.j[l],this.k[l]]+
              log(dpois(y.true.cand[swapped[1],this.j[l],this.k[l]],lambda=lambda)/(1-exp(-lambda)))
          }
          if(y.true.cand[swapped[2],this.j[l],this.k[l]]==0){
            ll.y.cand[swapped[2],this.j[l],this.k[l]] <- log(1-model$pd[swapped[2],this.j[l]])
          }else{
            #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
            ll.y.cand[swapped[2],this.j[l],this.k[l]] <- log(model$pd[swapped[2],this.j[l]])
            ll.y.cand[swapped[2],this.j[l],this.k[l]] <- ll.y.cand[swapped[2],this.j[l],this.k[l]]+
              log(dpois(y.true.cand[swapped[2],this.j[l],this.k[l]],lambda=lambda)/(1-exp(-lambda)))
          }
          
          #get backwards proposal probs (not symmetric)
          lp.y.prop.back <- rep(0,M)
          for(i in 1:M){
            if(z[i]==1){
              if(i!=ID.cand[l]){#new state
                y.tmp1 <- y.true.cand[i,this.j[l],this.k[l]] + 1 #if we add sample here
                y.tmp2 <- y.true.cand[ID.cand[l],this.j[l],this.k[l]] - 1 #if we add sample here
                if(y.tmp1==0){ #if not captured
                  lp.y.prop.back[i] <- dbinom(0,size=1,prob=model$pd[i,this.j[l]],log=TRUE)
                }else{ #if captured
                  lp.y.prop.back[i] <- dbinom(1,size=1,prob=model$pd[i,this.j[l]],log=TRUE) + 
                    log(dpois(y.tmp1,lambda=lambda)/(1-dpois(0,lambda=lambda)))
                }
                if(y.tmp2==0){ #if not captured
                  lp.y.prop.back[i] <- lp.y.prop.back[i] + dbinom(0,size=1,prob=model$pd[i,this.j[l]],log=TRUE)
                }else{ #if captured
                  lp.y.prop.back[i] <- lp.y.prop.back[i] + dbinom(1,size=1,prob=model$pd[i,this.j[l]],log=TRUE) + 
                    log(dpois(y.tmp2,lambda=lambda)/(1-dpois(0,lambda=lambda)))
                }
              }else{ #can't propose this guy if z==0
                lp.y.prop.back[i] <- -Inf
              }
            }else{
              lp.y.prop.back[i] <- -Inf
            }
          }
          prop.probs.back <- exp(lp.y.prop.back)
          prop.probs.back <- prop.probs.back/sum(prop.probs.back)
          
          prop.prob.for <- prop.probs[swapped[2]]
          prop.prob.back <- prop.probs.back[swapped[1]]
          
          #probability we select this y[i,j,k] to update by selecting a sample ID at random
          select.prob.for <- sum(ID.curr==ID.curr[l]&this.j==this.j[l]&this.k==this.k[l])/n.samples
          select.prob.back <- sum(ID.cand==ID.cand[l]&this.j==this.j[l]&this.k==this.k[l])/n.samples
          
          #sum log likelihoods and do MH step
          lp_initial <- ll.y[swapped[1],this.j[l],this.k[l]] + ll.y[swapped[2],this.j[l],this.k[l]]
          lp_proposed <- ll.y.cand[swapped[1],this.j[l],this.k[l]] + ll.y.cand[swapped[2],this.j[l],this.k[l]]
          log_MH_ratio <- (lp_proposed+log(prop.prob.back)+log(select.prob.back)) - 
            (lp_initial+log(prop.prob.for)+log(select.prob.for))
          
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l],this.k[l]] <- y.true.cand[swapped[1],this.j[l],this.k[l]]
            y.true[swapped[2],this.j[l],this.k[l]] <- y.true.cand[swapped[2],this.j[l],this.k[l]]
            ll.y[swapped[1],this.j[l],this.k[l]] <- ll.y.cand[swapped[1],this.j[l],this.k[l]]
            ll.y[swapped[2],this.j[l],this.k[l]] <- ll.y.cand[swapped[2],this.j[l],this.k[l]]
            ID.curr[l] <- ID.cand[l]
          }else{ #set back
            y.true.cand[swapped[1],this.j[l],this.k[l]] <- y.true[swapped[1],this.j[l],this.k[l]]
            y.true.cand[swapped[2],this.j[l],this.k[l]] <- y.true[swapped[2],this.j[l],this.k[l]]
            ll.y.cand[swapped[1],this.j[l],this.k[l]] <- ll.y[swapped[1],this.j[l],this.k[l]]
            ll.y.cand[swapped[2],this.j[l],this.k[l]] <- ll.y[swapped[2],this.j[l],this.k[l]]
            ID.cand[l] <- ID.curr[l]
          }
        }else{ #set back
          ID.cand[l] <- ID.curr[l]
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

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)