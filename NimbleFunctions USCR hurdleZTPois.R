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

dHurdleZTPoisVector <- nimbleFunction(
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
      J=nimDim(x)[1]
      K=nimDim(x)[2]
      logProb = 0
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]==1){
            if(x[j,k]==0){
              logProb = logProb + log(1-pd[j])
            }else{
              logProb = logProb + log(pd[j]) + log(dpois(x[j,k],lambda=lambda)/(1-exp(-lambda)))
            }
          }
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rHurdleZTPoisVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1),K2D = double(2), z = double(0), lambda = double(0)) {
    returnType(double(2))
    J=nimDim(pd)[1]
    K=nimDim(pd)[2]
    out=matrix(0,J,K)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(3)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    K <- nimDim(y.true)[3]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J,1:K])
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


## sampler to update y[1:M,1:J,1:K]
#this sampler updates all samples at a j-k index at once. The version commented out
#below updates sample IDs one at a time. For this model, this algorithm (one at a time) is not sufficient.
#use this one. It is a little confusing because I am still keeping up with sample IDs.
IDSampler_jk <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j<-control$this.j
    this.k<-control$this.k
    M<-control$M
    J <- control$J
    K <- control$K
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    z <- model$z
    y.true <- model$y.true
    ID.curr <- model$ID
    lambda <- model$lambda
    #precalculate log likelihoods at individual by trap level
    ll.y <- array(0,dim=c(M,J,K))
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(model$K2D[j,k]==1){
              if(y.true[i,j,k]==0){
                ll.y[i,j,k] = log(1-model$pd[i,j])
              }else{
                #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                ll.y[i,j,k] = log(model$pd[i,j])
                ll.y[i,j,k] = ll.y[i,j,k] + log(dpois(y.true[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
              }
            }
          }
        }
      }
    }
    ll.y.cand <- ll.y
    
    ID.cand=ID.curr
    y.true.cand=y.true
    for(j in 1:J){
      if(any(this.j==j)){
        #same proposal distribution for all traps
        propprobs=model$pd[,j]*z
        if(sum(propprobs)>0){ #abort if propprobs sum to 0 so we don't divide by 0
          propprobs=propprobs/sum(propprobs)
          for(k in 1:K){
            if(any(this.k==k&this.j==j)){
              these.samps=which(this.j==j&this.k==k)
              n.these.samps=length(these.samps)
              #sample the .jk capture history
              prop.jk=rmulti(1,n.these.samps,propprobs)
              #reformat proposal in terms of sample IDs
              #scramble "these.samps" so posterior ID probs are correct. bc no "sample" function...
              if(n.these.samps>1){
                tmp=rep(1,n.these.samps)
                scramble.idx=rep(0,n.these.samps)
                for(xx in 1:(n.these.samps)){
                  probs=tmp/sum(tmp)
                  choose=rcat(1,probs)
                  scramble.idx[xx]=choose
                  tmp[choose]=0
                }
                these.samps=these.samps[scramble.idx]
              }
              
              idx=1
              for(i in 1:M){
                if(prop.jk[i]>0){
                  for(i2 in 1:prop.jk[i]){
                    ID.cand[these.samps[idx]]=i
                    idx=idx+1
                  }
                }
              }
              y.true.cand[,j,k]=prop.jk
              #update likelihood
              for(i in 1:M){ #can be made more efficient here...
                if(z[i]>0){
                  if(y.true.cand[i,j,k]==0){
                    ll.y.cand[i,j,k] = log(1-model$pd[i,j])
                  }else{
                    #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                    ll.y.cand[i,j,k] = log(model$pd[i,j])
                    ll.y.cand[i,j,k] = ll.y.cand[i,j,k] + log(dpois(y.true.cand[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
                  }
                }
              }
              
              #get proposal probabilities
              forprob=dmulti(y.true.cand[,j,k],n.these.samps,propprobs)
              backprob=dmulti(y.true[,j,k],n.these.samps,propprobs)
              
              #sum log likelihoods and do MH step
              lp_initial <- sum(ll.y[,j,k])
              lp_proposed <- sum(ll.y.cand[,j,k])
              log_MH_ratio <- (lp_proposed+log(backprob)) - (lp_initial+log(forprob))
              accept <- decide(log_MH_ratio)
              if(accept){
                y.true[,j,k]=y.true.cand[,j,k]
                ll.y[,j,k]=ll.y.cand[,j,k]
                ID.curr[these.samps]=ID.cand[these.samps]
              }
            }
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


## sampler to update y[1:M,1:J,1:K]
#This version updates 1 sample ID at a time. This algorithm is not sufficient. Use the one above.
# IDSampler <- nimbleFunction(
#   contains = sampler_BASE,
#   setup = function(model, mvSaved, target, control) {
#     this.j<-control$this.j
#     this.k<-control$this.k
#     n.samples <- control$n.samples
#     M<-control$M
#     J <- control$J
#     K <- control$K
#     calcNodes <- model$getDependencies(target)
#   },
#   
#   run = function() {
#     z <- model$z
#     y.true <- model$y.true
#     ID.curr <- model$ID
#     lambda <- model$lambda
#     #precalculate log likelihoods at individual by trap level
#     ll.y <- array(0,dim=c(M,J,K))
#     for(i in 1:M){
#       if(z[i]==1){
#         for(j in 1:J){
#           for(k in 1:K){
#             if(model$K2D[j,k]==1){
#               if(y.true[i,j,k]==0){
#                 ll.y[i,j,k] = log(1-model$pd[i,j])
#               }else{
#                 #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
#                 ll.y[i,j,k] = log(model$pd[i,j])
#                 ll.y[i,j,k] = ll.y[i,j,k] + log(dpois(y.true[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
#               }
#             }
#           }
#         }
#       }
#     }
#     ll.y.cand <- ll.y
#     
#     for(l in 1:n.samples){#loop over samples
#       ID.cand=ID.curr
#       y.true.cand=y.true
#       #proposal distribution
#       propprobs=model$pd[1:M,this.j[l]]*z
#       sumpropprobs=sum(propprobs)
#       if(sumpropprobs>0){ #abort here to avoid dividing by 0. I don't think this will actually every happen, though. Can always propose to yourself.
#         propprobs <- propprobs/sumpropprobs
#         ID.cand[l]=rcat(1,prob=propprobs)
#         if(ID.cand[l]!=ID.curr[l]){
#           swapped=c(ID.curr[l],ID.cand[l])
#           #new sample proposal probabilities
#           forprob=propprobs[swapped[2]]
#           backprob=propprobs[swapped[1]]
#           #new y.true's - move sample from ID to ID.cand
#           y.true.cand[ID.curr[l],this.j[l],this.k[l]]=y.true[ID.curr[l],this.j[l],this.k[l]]-1
#           y.true.cand[ID.cand[l],this.j[l],this.k[l]]=y.true[ID.cand[l],this.j[l],this.k[l]]+1
#           #update ll.y
#           if(y.true.cand[swapped[1],this.j[l],this.k[l]]==0){
#             ll.y.cand[swapped[1],this.j[l],this.k[l]]=log(1-model$pd[swapped[1],this.j[l]])
#           }else{
#             #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
#             ll.y.cand[swapped[1],this.j[l],this.k[l]]=log(model$pd[swapped[1],this.j[l]])
#             ll.y.cand[swapped[1],this.j[l],this.k[l]]=ll.y.cand[swapped[1],this.j[l],this.k[l]]+
#               log(dpois(y.true.cand[swapped[1],this.j[l],this.k[l]],lambda=lambda[1])/(1-exp(-lambda[1])))
#           }
#           if(y.true.cand[swapped[2],this.j[l],this.k[l]]==0){
#             ll.y.cand[swapped[2],this.j[l],this.k[l]]=log(1-model$pd[swapped[2],this.j[l]])
#           }else{
#             #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
#             ll.y.cand[swapped[2],this.j[l],this.k[l]]=log(model$pd[swapped[2],this.j[l]])
#             ll.y.cand[swapped[2],this.j[l],this.k[l]]=ll.y.cand[swapped[2],this.j[l],this.k[l]]+
#               log(dpois(y.true.cand[swapped[2],this.j[l],this.k[l]],lambda=lambda[1])/(1-exp(-lambda[1])))
#           }
#           #select sample to move proposal probabilities
#           #P(select a sample for this ID)*P(select this j and k|this ID)
#           #n.samples cancels out in MH ratio. Including for clarity
#           focalprob=(sum(ID.curr==swapped[1])/n.samples)*
#             (y.true[swapped[1],this.j[l],this.k[l]]/sum(y.true[swapped[1],,]))
#           
#           focalbackprob=(sum(ID.cand==swapped[2])/n.samples)*
#             (y.true.cand[swapped[2],this.j[l],this.k[l]]/sum(y.true.cand[swapped[2],,]))
#           
#           #sum log likelihoods and do MH step
#           lp_initial <- ll.y[swapped[1],this.j[l],this.k[l]]+ll.y[swapped[2],this.j[l],this.k[l]]
#           lp_proposed <- ll.y.cand[swapped[1],this.j[l],this.k[l]]+ll.y.cand[swapped[2],this.j[l],this.k[l]]
#           log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
#           accept <- decide(log_MH_ratio)
#           if(accept){
#             y.true[swapped[1],this.j[l],this.k[l]]=y.true.cand[swapped[1],this.j[l],this.k[l]]
#             y.true[swapped[2],this.j[l],this.k[l]]=y.true.cand[swapped[2],this.j[l],this.k[l]]
#             ll.y[swapped[1],this.j[l],this.k[l]]=ll.y.cand[swapped[1],this.j[l],this.k[l]]
#             ll.y[swapped[2],this.j[l],this.k[l]]=ll.y.cand[swapped[2],this.j[l],this.k[l]]
#             ID.curr[l]=ID.cand[l]
#           }
#         }
#       }
#     }
#     #put everything back into the model$stuff
#     model$y.true <<- y.true
#     model$ID <<- ID.curr
#     model$calculate(calcNodes) #update logprob
#     copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
#   },
#   methods = list( reset = function () {} )
# )
