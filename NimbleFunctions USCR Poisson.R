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
    j.indicator <- control$j.indicator
    M<-control$M
    J <- control$J
    cluster.ups <- control$cluster.ups
    local.eval <- control$local.eval
    swap.rad.multiplier <- control$swap.rad.multiplier
    calcNodes <- model$getDependencies(c("y.true","z"))
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
          ID.cand <- rcat(1,fullcond)
          model$ID[idx[l]] <<- ID.cand
          #update y
          model$y.true[ID.curr,trapup[j]] <<- model$y.true[ID.curr,trapup[j]]-1 #subtract out old ID obs
          model$y.true[ID.cand,trapup[j]] <<- model$y.true[ID.cand,trapup[j]]+1 #add in new ID obs
        }
      }
    }
    
    #Now, we do a joint z-ID update
    if(cluster.ups>0){ #skip if cluster.ups=0
      y.true <- model$y.true
      lam <- model$lam
      K1D <- model$K1D
      z <- model$z
      #precalculate ll.y
      ll.y <- matrix(0,M,J)
      for(i in 1:M){
        if(z[i]==1){
          ll.y[i,] = dpois(y.true[i,],K1D*lam[i,],log=TRUE)
        }
      }
      ID.curr2 <- model$ID #can't reuse object with same name but different dimensions, adding "2" to make nimble happy
      swap.rad=swap.rad.multiplier*model$sigma[1] #radius for traps to update sample IDs around a focal individual
      for(up in 1:cluster.ups){ #how many times to do this update per iteration?
        # select z==1 to turn off
        z.cand=z
        ID.cand2=ID.curr2
        y.cand=y.true
        ll.y.cand=ll.y
        lam.cand=lam
        z.on=which(z==1)
        n.z.on=length(z.on)
        if(n.z.on>1){ #Cannot turn off anyone unless there are at least 2 guys on. samples must belong to someone!
          if(n.z.on>1){
            pick=rcat(1,rep(1/n.z.on,n.z.on))
            focal.i=z.on[pick]
          }else{
            focal.i=z.on[1]
          }
          z.cand[focal.i]=0
          lam.cand[focal.i,]=0
          p.select.z.for=1/n.z.on
          if(local.eval==TRUE){# find local traps with samples
            dists=sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
            focal.traps=which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps=which(j.indicator) #j.indicator removes traps with 0 samples
          }
          total.log.j.probs.for=0
          total.log.j.probs.back=0
          n.focal.traps=length(focal.traps)
          abort=FALSE #abort if any propprobs so small we get underflow. Would be rejected if there were no underflow.
          #almost never happens...
          if(n.focal.traps>0){
            # repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              these.samps=which(this.j==focal.traps[j])
              n.these.samps=length(these.samps)
              propprobs.for=lam.cand[,focal.traps[j]]*z.cand
              propprobs.back=lam[,focal.traps[j]]*z
              sum.propprobs.for=sum(propprobs.for)
              if(sum.propprobs.for==0){
                abort=TRUE
              }
              propprobs.for=propprobs.for/sum.propprobs.for
              propprobs.back=propprobs.back/sum(propprobs.back)
              y.cand[,focal.traps[j]]=0
              for(l in 1:n.these.samps){
                pick = rcat(1,prob=propprobs.for)
                ID.cand2[these.samps[l]]=pick
                y.cand[ID.cand2[these.samps[l]],focal.traps[j]]=y.cand[ID.cand2[these.samps[l]],focal.traps[j]]+1
              }
              total.log.j.probs.for=total.log.j.probs.for+dmulti(y.cand[,focal.traps[j]],n.these.samps,prob=propprobs.for,log=TRUE)
              total.log.j.probs.back=total.log.j.probs.back+dmulti(y.true[,focal.traps[j]],n.these.samps,prob=propprobs.back,log=TRUE)
              
              #update ll.y.cand - only focal traps with samples here
              for(i in 1:M){
                if(z.cand[i]==1){
                  ll.y.cand[i,focal.traps[j]]=dpois(y.cand[i,focal.traps[j]],
                                                    K1D[focal.traps[j]]*lam.cand[i,focal.traps[j]],log=TRUE)
                }
              }
            }
            #update ll.y.cand for focal.i
            ll.y.cand[focal.i,]=0
          }else{#if we only turn off a z and no local samples to reallocate
            ll.y.cand[focal.i,]=0
          }
          if(!abort){#if propprobs didn't have underflow
            ll.z.curr=dbinom(z[focal.i],1,model$psi[1],log=TRUE)
            ll.z.cand=dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
            
            z.off=which(z.cand==0)
            p.select.z.back=1/length(z.off)
            
            logforprob=log(p.select.z.for)+total.log.j.probs.for
            logbackprob=log(p.select.z.back)+total.log.j.probs.back
            
            if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
              ll.total.curr=sum(ll.y)+ll.z.curr #just summing full y likelihood for ease
              ll.total.cand=sum(ll.y.cand)+ll.z.cand
            }else{#y likelihood for focal ind only, all traps
              ll.total.curr=sum(ll.y[focal.i,])+ll.z.curr
              ll.total.cand=sum(ll.y.cand[focal.i,])+ll.z.cand
            }
            log_MH_ratio=(ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
            accept <- decide(log_MH_ratio)
            if(accept){
              if(n.focal.traps>0){
                y.true[,focal.traps]=y.cand[,focal.traps]
                ll.y[,focal.traps]=ll.y.cand[,focal.traps]
                ID.curr2=ID.cand2
              }
              lam[focal.i,]=lam.cand[focal.i,]
              ll.y[focal.i,]=ll.y.cand[focal.i,]
              z[focal.i]=z.cand[focal.i]
            }
          }
        }
        
        #select z==0 to turn on
        z.cand=z
        ID.cand2=ID.curr2
        y.cand=y.true
        ll.y.cand=ll.y
        lam.cand=lam
        z.off=which(z==0)
        n.z.off=length(z.off)
        if(n.z.off>0){
          if(n.z.off>1){
            pick=rcat(1,rep(1/n.z.off,n.z.off))
            focal.i=z.off[pick]
          }else{
            focal.i=z.off[1]
          }
          z.cand[focal.i]=1

          p.select.z.for=1/length(z.off)
          
          dists=sqrt((model$s[focal.i,1]-model$X[,1])^2+(model$s[focal.i,2]-model$X[,2])^2)
          if(local.eval==TRUE){# find local traps with samples
            focal.traps=which(dists<swap.rad&j.indicator) #j.indicator removes traps with 0 samples
          }else{
            focal.traps=which(j.indicator) #j.indicator removes traps with 0 samples
          }
          lam.cand[focal.i,]=model$lam0[1]*exp(-dists^2/(2*model$sigma[1]^2))
          total.log.j.probs.for=0
          total.log.j.probs.back=0
          n.focal.traps=length(focal.traps)
          if(n.focal.traps>0){
            #repropose all samples at focal.traps
            for(j in 1:n.focal.traps){
              these.samps=which(this.j==focal.traps[j])
              n.these.samps=length(these.samps)
              propprobs.for=lam.cand[,focal.traps[j]]*z.cand
              propprobs.back=lam[,focal.traps[j]]*z
              propprobs.for=propprobs.for/sum(propprobs.for)
              propprobs.back=propprobs.back/sum(propprobs.back)
              y.cand[,focal.traps[j]]=0
              for(l in 1:n.these.samps){
                pick = rcat(1,prob=propprobs.for)
                ID.cand2[these.samps[l]]=pick
                y.cand[ID.cand2[these.samps[l]],focal.traps[j]]=y.cand[ID.cand2[these.samps[l]],focal.traps[j]]+1
              }
              total.log.j.probs.for=total.log.j.probs.for+dmulti(y.cand[,focal.traps[j]],n.these.samps,prob=propprobs.for,log=TRUE)
              total.log.j.probs.back=total.log.j.probs.back+dmulti(y.true[,focal.traps[j]],n.these.samps,prob=propprobs.back,log=TRUE)
              #update ll.y.cand - only focal traps with samples here
              for(i in 1:M){
                if(z.cand[i]==1){
                  ll.y.cand[i,focal.traps[j]]=dpois(y.cand[i,focal.traps[j]],
                                                    K1D[focal.traps[j]]*lam.cand[i,focal.traps[j]],log=TRUE)
                }
              }
            }
            ll.y.cand[focal.i,]=dpois(y.cand[focal.i,],K1D*lam.cand[focal.i,],log=TRUE)
          }else{#if we only turn on a z and no local samples to reallocate
            ll.y.cand[focal.i,]=dpois(y.cand[focal.i,],K1D*lam.cand[focal.i,],log=TRUE)
          }
          ll.z.curr=dbinom(z[focal.i],1,model$psi[1],log=TRUE)
          ll.z.cand=dbinom(z.cand[focal.i],1,model$psi[1],log=TRUE)
          
          z.on=which(z.cand==1)
          p.select.z.back=1/length(z.on)

          logforprob=log(p.select.z.for)+total.log.j.probs.for
          logbackprob=log(p.select.z.back)+total.log.j.probs.back

          if(n.focal.traps>0){#y likelihood of everyone at focal traps and all traps for focal individual
            ll.total.curr=sum(ll.y)+ll.z.curr #just summing full likelihood for ease
            ll.total.cand=sum(ll.y.cand)+ll.z.cand
          }else{#y likelihood for focal ind only, all traps
            ll.total.curr=sum(ll.y[focal.i,])+ll.z.curr
            ll.total.cand=sum(ll.y.cand[focal.i,])+ll.z.cand
          }

          log_MH_ratio=(ll.total.cand+logbackprob)-(ll.total.curr+logforprob)
          accept <- decide(log_MH_ratio)
          if(accept){
            if(n.focal.traps>0){
              y.true[,focal.traps]=y.cand[,focal.traps]
              ID.curr2=ID.cand2
              ll.y[,focal.traps]=ll.y.cand[,focal.traps]
            }
            ll.y[focal.i,]=ll.y.cand[focal.i,]
            z[focal.i]=z.cand[focal.i]
            lam[focal.i,]=lam.cand[focal.i,]
          }
        }
      }
      
      #update model$stuff
      model$y.true <<- y.true
      model$ID <<- ID.curr2
      model$z <<- z #lam will be updated with calculate below
    }
    #update lp
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)