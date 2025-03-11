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
rHurdleZTPoisVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1),K2D = double(2), z = double(0), lambda = double(0)) {
    returnType(double(2))
    J <- nimDim(pd)[1]
    K <- nimDim(pd)[2]
    out <- matrix(0,J,K)
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


## sampler to update y[1:M,1:J,1:K]
#this sampler updates all samples at a j-k index at once. The version commented out
#below updates sample IDs one at a time. For this model, this algorithm (one at a time) is not sufficient.
#use this one. It is a little confusing because I am still keeping up with sample IDs.
IDSampler_jk <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j <- control$this.j
    this.k <- control$this.k
    M <- control$M
    J <- control$J
    K <- control$K
    j.indicator <- control$j.indicator
    jk.indicator <- control$jk.indicator
    K2D <- control$K2D
    cluster.ups <- control$cluster.ups
    local.eval <- control$local.eval
    swap.rad.multiplier <- control$swap.rad.multiplier
    calcNodes <- model$getDependencies(c("y.true","z","ID"))
  },
  
  run = function() {
    z <- model$z
    y.true <- model$y.true
    ID.curr <- model$ID
    lambda <- model$lambda
    pd <- model$pd
    #precalculate log likelihoods at individual by trap by occasion level
    ll.y <- array(0,dim=c(M,J,K))
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(K2D[j,k]==1){
              if(y.true[i,j,k]==0){
                ll.y[i,j,k] <- log(1-pd[i,j])
              }else{
                #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                ll.y[i,j,k] <- log(pd[i,j])
                ll.y[i,j,k] <- ll.y[i,j,k] + log(dpois(y.true[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
              }
            }
          }
        }
      }
    }
    ll.y.cand <- ll.y
    
    ID.cand <- ID.curr
    y.true.cand <- y.true
    for(j in 1:J){
      if(any(this.j==j)){
        #same proposal distribution for all samples at this trap
        propprobs <- pd[,j]*z
        if(sum(propprobs)>0){ #abort if propprobs sum to 0 so we don't divide by 0
          propprobs <- propprobs/sum(propprobs)
          for(k in 1:K){
            if(any(this.k==k&this.j==j)){
              these.samps <- which(this.j==j&this.k==k)
              n.these.samps <- length(these.samps)
              #sample the .jk capture history
              prop.jk <- rmulti(1,n.these.samps,propprobs)
              #reformat proposal in terms of sample IDs
              #scramble "these.samps" so posterior ID probs are correct. bc no "sample" function...
              if(n.these.samps>1){
                tmp <- rep(1,n.these.samps)
                scramble.idx <- rep(0,n.these.samps)
                for(xx in 1:(n.these.samps)){
                  probs <- tmp/sum(tmp)
                  choose <- rcat(1,probs)
                  scramble.idx[xx] <- choose
                  tmp[choose] <- 0
                }
                these.samps <- these.samps[scramble.idx]
              }
              
              idx <- 1
              for(i in 1:M){
                if(prop.jk[i]>0){
                  for(i2 in 1:prop.jk[i]){
                    ID.cand[these.samps[idx]] <- i
                    idx <- idx+1
                  }
                }
              }
              y.true.cand[,j,k] <- prop.jk
              #update likelihood
              for(i in 1:M){ #can be made more efficient here...
                if(z[i]>0){
                  if(y.true.cand[i,j,k]==0){
                    ll.y.cand[i,j,k] <- log(1-pd[i,j])
                  }else{
                    #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                    ll.y.cand[i,j,k] <- log(pd[i,j])
                    ll.y.cand[i,j,k] <- ll.y.cand[i,j,k] + log(dpois(y.true.cand[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
                  }
                }
              }
              
              #get proposal probabilities
              forprob <- dmulti(y.true.cand[,j,k],n.these.samps,propprobs)
              backprob <- dmulti(y.true[,j,k],n.these.samps,propprobs)
              
              #sum log likelihoods and do MH step
              lp_initial <- sum(ll.y[,j,k])
              lp_proposed <- sum(ll.y.cand[,j,k])
              log_MH_ratio <- (lp_proposed+log(backprob)) - (lp_initial+log(forprob))
              accept <- decide(log_MH_ratio)
              if(accept){
                y.true[,j,k] <- y.true.cand[,j,k]
                ll.y[,j,k] <- ll.y.cand[,j,k]
                ID.curr[these.samps] <- ID.cand[these.samps]
              }
            }
          }
        }
      }
    }
    #put everything back into the model$stuff
    model$y.true <<- y.true
    model$ID <<- ID.curr
    
    #Now, we do a joint z-ID update
    if(cluster.ups>0){ #skip if cluster.ups=0
      y.true <- model$y.true
      pd <- model$pd
      z <- model$z
      #precalculate log likelihoods at individual by trap by occasion level
      ll.y <- array(0,dim=c(M,J,K))
      for(i in 1:M){
        if(z[i]==1){
          for(j in 1:J){
            for(k in 1:K){
              if(K2D[j,k]==1){
                if(y.true[i,j,k]==0){
                  ll.y[i,j,k] <- log(1-pd[i,j])
                }else{
                  ll.y[i,j,k] <- log(pd[i,j])
                  ll.y[i,j,k] <- ll.y[i,j,k] + log(dpois(y.true[i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
                }
              }
            }
          }
        }
      }
      ID.curr2 <- model$ID #can't reuse object with same name but different dimensions, adding "2" to make nimble happy
      swap.rad <- swap.rad.multiplier*model$sigma[1] #radius for traps to update sample IDs around a focal individual
      for(up in 1:cluster.ups){ #how many times to do this update per iteration?
        # select z==1 to turn off
        z.cand <- z
        ID.cand2 <- ID.curr2
        pd.cand <- pd
        y.cand <- y.true
        ll.y.cand <- ll.y
        z.on <- which(z==1)
        n.z.on <- length(z.on)
        if(n.z.on>1){ #Cannot turn off anyone unless there are at least 2 guys on. samples must belong to someone!
          if(n.z.on>1){
            pick <- rcat(1,rep(1/n.z.on,n.z.on))
            focal.i <- z.on[pick]
          }else{
            focal.i <- z.on[1]
          }
          z.cand[focal.i] <- 0
          pd.cand[focal.i,] <- 0
          p.select.z.for <- 1/n.z.on
          if(local.eval==TRUE){# find local traps with samples
            dists <- sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
            focal.traps <- which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps <- which(j.indicator) #j.indicator removes traps with 0 samples
          }
          total.log.j.probs.for <- 0
          total.log.j.probs.back <- 0
          n.focal.traps <- length(focal.traps)
          abort <- FALSE #abort if any propprobs so small we get underflow. Would be rejected if there were no underflow.
          #almost never happens...
          if(n.focal.traps>0){
            # repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              propprobs.for <- pd.cand[,focal.traps[j]]*z.cand
              propprobs.back <- pd[,focal.traps[j]]*z
              sum.propprobs.for <- sum(propprobs.for)
              if(sum.propprobs.for==0){
                abort <- TRUE
              }
              propprobs.for <- propprobs.for/sum.propprobs.for
              propprobs.back <- propprobs.back/sum(propprobs.back)
              for(k in 1:K){
                if(jk.indicator[focal.traps[j],k]){ #if samples at this j-k
                  these.samps <- which(this.j==focal.traps[j]&this.k==k)
                  n.these.samps <- length(these.samps)
                  for(i in 1:M){
                    y.cand[i,focal.traps[j],k] <- 0
                  }
                  for(l in 1:n.these.samps){
                    pick <- rcat(1,prob=propprobs.for)
                    ID.cand2[these.samps[l]] <- pick
                    y.cand[ID.cand2[these.samps[l]],focal.traps[j],k] <- 
                      y.cand[ID.cand2[these.samps[l]],focal.traps[j],k]+1
                  }
                  total.log.j.probs.for <- total.log.j.probs.for+dmulti(y.cand[,focal.traps[j],k],n.these.samps,prob=propprobs.for,log=TRUE)
                  total.log.j.probs.back <- total.log.j.probs.back+dmulti(y.true[,focal.traps[j],k],n.these.samps,prob=propprobs.back,log=TRUE)
                  
                  #update ll.y.cand - only focal traps with samples here
                  for(i in 1:M){
                    if(z.cand[i]==1){
                      if(K2D[focal.traps[j],k]==1){
                        if(y.cand[i,focal.traps[j],k]==0){
                          ll.y.cand[i,focal.traps[j],k] <- log(1-pd.cand[i,focal.traps[j]])
                        }else{
                          ll.y.cand[i,focal.traps[j],k] <- log(pd.cand[i,focal.traps[j]])
                          ll.y.cand[i,focal.traps[j],k] <- ll.y.cand[i,focal.traps[j],k] +
                            log(dpois(y.cand[i,focal.traps[j],k],lambda=lambda[1])/(1-exp(-lambda[1])))
                        }
                      }
                    }
                  }
                }
              }
            }
            #update ll.y.cand for focal.i
            for(j in 1:J){
              for(k in 1:K){
                ll.y.cand[focal.i,j,k] <- 0
              }
            }
          }else{#if we only turn off a z and no local samples to reallocate
            for(j in 1:J){
              for(k in 1:K){
                ll.y.cand[focal.i,j,k] <- 0
              }
            }
          }
          if(!abort){#if propprobs didn't have underflow
            ll.z.curr <- dbinom(z[focal.i],1,model$psi[1],log=TRUE)
            ll.z.cand <- dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
            
            z.off <- which(z.cand==0)
            p.select.z.back <- 1/length(z.off)
            
            logforprob <- log(p.select.z.for)+total.log.j.probs.for
            logbackprob <- log(p.select.z.back)+total.log.j.probs.back
            
            if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
              ll.total.curr <- ll.z.curr #just summing full y likelihood for ease
              ll.total.cand <- ll.z.cand
              for(i in 1:M){
                for(j in 1:J){
                  for(k in 1:K){
                    ll.total.curr <- ll.total.curr+ll.y[i,j,k]
                    ll.total.cand <- ll.total.cand+ll.y.cand[i,j,k]
                  }
                }
              }
            }else{#y likelihood for focal ind only, all traps
              ll.total.curr <- ll.z.curr #just summing full y likelihood for ease
              ll.total.cand <- ll.z.cand
              for(j in 1:J){
                for(k in 1:K){
                  ll.total.curr <- ll.total.curr+ll.y[focal.i,j,k]
                  ll.total.cand <- ll.total.cand+ll.y.cand[focal.i,j,k]
                }
              }
            }
            log_MH_ratio <- (ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
            accept <- decide(log_MH_ratio)
            if(accept){
              if(n.focal.traps>0){
                for(i in 1:M){
                  for(j in 1:n.focal.traps){
                    for(k in 1:K){
                      y.true[i,focal.traps[j],k] <- y.cand[i,focal.traps[j],k]
                      ll.y[i,focal.traps[j],k] <- ll.y.cand[i,focal.traps[j],k]
                    }
                  }
                }
                ID.curr2 <- ID.cand2
              }
              for(j in 1:J){
                for(k in 1:K){
                  ll.y[focal.i,j,k] <- ll.y.cand[focal.i,j,k]
                }
              }
              z[focal.i] <- z.cand[focal.i]
              pd[focal.i,] <- pd.cand[focal.i,]
            }
          }
        }
        
        #select z==0 to turn on
        z.cand <- z
        ID.cand2 <- ID.curr2
        y.cand <- y.true
        pd.cand <- pd
        ll.y.cand <- ll.y
        z.off <- which(z==0)
        n.z.off <- length(z.off)
        if(n.z.off>0){
          if(n.z.off>1){
            pick <- rcat(1,rep(1/n.z.off,n.z.off))
            focal.i <- z.off[pick]
          }else{
            focal.i <- z.off[1]
          }
          z.cand[focal.i] <- 1
          
          p.select.z.for <- 1/length(z.off)
          #find local traps
          dists <- sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
          dists <- sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
          if(local.eval==TRUE){# find local traps with samples
            focal.traps <- which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps <- which(j.indicator) #j.indicator removes traps with 0 samples
          }
          pd.cand[focal.i,] <- model$p0[1]*exp(-dists^2/(2*model$sigma[1]^2))
          total.log.j.probs.for <- 0
          total.log.j.probs.back <- 0
          n.focal.traps=length(focal.traps)
          if(n.focal.traps>0){
            # repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              propprobs.for <- pd.cand[,focal.traps[j]]*z.cand
              propprobs.back <- pd[,focal.traps[j]]*z
              propprobs.for <- propprobs.for/sum(propprobs.for)
              propprobs.back <- propprobs.back/sum(propprobs.back)
              for(k in 1:K){
                if(jk.indicator[focal.traps[j],k]){ #if samples at this j-k
                  these.samps <- which(this.j==focal.traps[j]&this.k==k)
                  n.these.samps <- length(these.samps)
                  for(i in 1:M){
                    y.cand[i,focal.traps[j],k] <- 0
                  }
                  for(l in 1:n.these.samps){
                    pick <- rcat(1,prob=propprobs.for)
                    ID.cand2[these.samps[l]] <- pick
                    y.cand[ID.cand2[these.samps[l]],focal.traps[j],k] <- 
                      y.cand[ID.cand2[these.samps[l]],focal.traps[j],k]+1
                  }
                  total.log.j.probs.for <- total.log.j.probs.for+dmulti(y.cand[,focal.traps[j],k],n.these.samps,prob=propprobs.for,log=TRUE)
                  total.log.j.probs.back <- total.log.j.probs.back+dmulti(y.true[,focal.traps[j],k],n.these.samps,prob=propprobs.back,log=TRUE)
                  
                  #update ll.y.cand - only focal traps with samples here
                  for(i in 1:M){
                    if(z.cand[i]==1){
                      if(K2D[focal.traps[j],k]==1){
                        if(y.cand[i,focal.traps[j],k]==0){
                          ll.y.cand[i,focal.traps[j],k] <- log(1-pd.cand[i,focal.traps[j]])
                        }else{
                          #breaking into two because nimble "can't do math with arrays of more than 2 dimensions"
                          ll.y.cand[i,focal.traps[j],k] <- log(pd.cand[i,focal.traps[j]])
                          ll.y.cand[i,focal.traps[j],k] <- ll.y.cand[i,focal.traps[j],k] +
                            log(dpois(y.cand[i,focal.traps[j],k],lambda=lambda[1])/(1-exp(-lambda[1])))
                        }
                      }
                    }
                  }
                }
              }
            }
            #update focal.i likelihood for all j-k
            for(j in 1:J){
              for(k in 1:K){
                if(K2D[j,k]==1){
                  if(y.cand[focal.i,j,k]==0){
                    ll.y.cand[focal.i,j,k] <- log(1-pd.cand[focal.i,j])
                  }else{
                    ll.y.cand[focal.i,j,k] <- log(pd.cand[focal.i,j])
                    ll.y.cand[focal.i,j,k] <- ll.y.cand[focal.i,j,k] +
                      log(dpois(y.cand[focal.i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
                  }
                }
              }
            }
          }else{#if we only turn on a z and no local samples to reallocate
            #update focal.i likelihood for all j-k
            for(j in 1:J){
              for(k in 1:K){
                if(K2D[j,k]==1){
                  if(y.cand[focal.i,j,k]==0){
                    ll.y.cand[focal.i,j,k] <- log(1-pd.cand[focal.i,j])
                  }else{
                    ll.y.cand[focal.i,j,k] <- log(pd.cand[focal.i,j])
                    ll.y.cand[focal.i,j,k] <- ll.y.cand[focal.i,j,k] +
                      log(dpois(y.cand[focal.i,j,k],lambda=lambda[1])/(1-exp(-lambda[1])))
                  }
                }
              }
            }
          }
          
          ll.z.curr <- dbinom(z[focal.i],1,model$psi[1],log=TRUE)
          ll.z.cand <- dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
          
          z.on <- which(z.cand==1)
          p.select.z.back <- 1/length(z.on)
          
          logforprob <- log(p.select.z.for)+total.log.j.probs.for
          logbackprob <- log(p.select.z.back)+total.log.j.probs.back
          
          if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
            ll.total.curr <- ll.z.curr #just summing full y likelihood for ease
            ll.total.cand <- ll.z.cand
            for(i in 1:M){
              for(j in 1:J){
                for(k in 1:K){
                  ll.total.curr <- ll.total.curr+ll.y[i,j,k]
                  ll.total.cand <- ll.total.cand+ll.y.cand[i,j,k]
                }
              }
            }
          }else{#y likelihood for focal ind only, all traps
            ll.total.curr <- ll.z.curr #just summing full y likelihood for ease
            ll.total.cand <- ll.z.cand
            for(j in 1:J){
              for(k in 1:K){
                ll.total.curr <- ll.total.curr+ll.y[focal.i,j,k]
                ll.total.cand <- ll.total.cand+ll.y.cand[focal.i,j,k]
              }
            }
          }
          
          log_MH_ratio <- (ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
          accept <- decide(log_MH_ratio)
          if(accept){
            if(n.focal.traps>0){
              for(i in 1:M){
                for(j in 1:n.focal.traps){
                  for(k in 1:K){
                    y.true[i,focal.traps[j],k] <- y.cand[i,focal.traps[j],k]
                    ll.y[i,focal.traps[j],k] <- ll.y.cand[i,focal.traps[j],k]
                  }
                }
              }
              ID.curr2 <- ID.cand2
            }
            for(j in 1:J){
              for(k in 1:K){
                ll.y[focal.i,j,k] <- ll.y.cand[focal.i,j,k]
              }
            }
            z[focal.i] <- z.cand[focal.i]
            pd[focal.i,] <- pd.cand[focal.i,]
          }
        }
      }
      
      #update model$stuff
      model$y.true <<- y.true
      model$ID <<- ID.curr2
      model$z <<- z #pd updated due to z dependencies in calcNodes below
    }#end joint z-ID update
    
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
#     this.j <- control$this.j
#     this.k <- control$this.k
#     n.samples <- control$n.samples
#     M <- control$M
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
