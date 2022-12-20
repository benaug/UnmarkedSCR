#both base sampler and cluster buster appear to work n=1. test. I think. Compare 1 run each again.

#Multisession unmarked SCR. See single session test scripts for more
#comments about the single session models

library(nimble)
library(coda)
source("simUSCR.multisession.R")
source("simUSCR.R")
source("init.data.USCR.multisession.R")
source("init.data.USCR.R")
source("NimbleModel USCR Multisession Poisson.R")
source("NimbleFunctions USCR Multisession Poisson.R")
source("sSampler Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session=3
D = rep(0.2,N.session) #expected density in units of sigma and X
lam0=rep(0.25,N.session)
sigma=rep(0.5,N.session)
K=c(5,6,7) #number of occasions
buff=rep(3,N.session) #state space buffer
#make trapping arrays
X1=expand.grid(3:11,3:11)
X2=expand.grid(3:12,3:12)
X3=expand.grid(3:13,3:13)
X=list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area=getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda=D*area
lambda #expected N in each session

obstype="poisson"

#Simulate some data
data=simUSCR.multisession(N.session=N.session,lambda=lambda,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member in each session
head(t(data$this.j))
head(t(data$this.k))#not used in this sampler for 2D data, but required if using 3D data

#Data augmentation level
M=c(200,250,300)
maxM=max(M)

#trap operation matrix - plugging in perfect operation here
J=unlist(lapply(X,nrow))
maxJ=max(J)
K1D=matrix(NA,nrow=N.session,ncol=maxJ)
for(g in 1:N.session){
  K1D[g,1:J[g]]=rep(K[g],J[g])
}

inits=list(lam0=lam0,sigma=sigma)#initial values for lam0, sigma to build data. in practice, don't use simulated values as inits.
nimbuild=init.data.USCR.multisession(data=data,M=M,inits=inits)

#inits for nimble
N.init=rowSums(nimbuild$z,na.rm=TRUE) #N.init must be consistent with z.init!

Niminits <- list(N=N.init,z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=apply(nimbuild$y.true2D,c(1,2),sum,na.rm=TRUE),
                 y.true=nimbuild$y.true2D,lam0.fixed=0.5,sigma.fixed=0.5,D=0.5)

#constants for Nimble
constants<-list(N.session=N.session,M=M,J=J,K1D=K1D,n.samples=nimbuild$n.samples,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,area=area)

#supply data to nimble
Nimdata<-list(y.true=array(NA,dim=c(N.session,maxM,maxJ)),
              ID=matrix(NA,N.session,max(nimbuild$n.samples)),
              z=matrix(NA,N.session,maxM),X=nimbuild$X,capcounts=matrix(NA,N.session,maxM))

# set parameters to monitor
parameters<-c('lam0.fixed','sigma.fixed','N','n','D','lambda')

#can also monitor a different set of parameters with a different thinning rate
# parameters2 <- c("s","z") #use this if you want to check s acceptance rates when z=1 below
parameters2 <- c("ID")
nt=1 #thinning rate
nt2=1#thin more

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, monitors2=parameters2,thin2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Three *required* sampler replacements

##z/N update
z.ups=round(M*0.25) # how many z proposals per iteration per session?
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[",g,",","1:",M[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,y.nodes,lam.nodes)
  
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".

#Some notes on the joint z-ID update. First, you may not need it. It slows down computation.
#But mixing will be improved substantially in certain situations, particularly when overall detection 
#probability is high. 
#"cluster.ups" controls how many joint z-ID updates to do per iteration.
#So cluster.ups=0 does the default algorithm. (1 z addition proposal and 1 z subtraction proposal).
#Perhaps just cluster.ups=1 is a good compromise between computation speed and mixing for data sets that benefit from it.
#"local.eval" determines if we only consider the local traps during the update. This makes the update faster.
#If local.eval=TRUE, "swap.rad.multiplier" determines the local neighborhood of traps around the 
#focal z to consider. The traps considered are a function of sigma, a radius of 
#sigma*swap.rad.multiplier of the focal z. If you set "swap.rad.multipler" too low, the 
#proposal may no longer be reversible and will not work correctly. To be safe, set local.eval=FALSE.
#Otherwise, 6 seems to work fine in almost all cases (sometimes can fail for very sparse data sets).
conf$removeSampler("y.true")
cluster.ups=0 #how many joint z-ID updates per iteration? 
for(g in 1:N.session){
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[",g,",","1:",M[g],",1:",J[g],"]"))
  calcNodes=c(y.nodes,paste("capcounts[",g,",1:",M[g],"]"),paste("n[",g,"]"))
  if(cluster.ups>0){
    lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
    calcNodes=c(calcNodes,lam.nodes,paste("N[",g,"]")) #need N and lam nodes for joint update
  }
  trapup=unique(data$this.j[g,1:nimbuild$n.samples[g]])
  j.indicator=1:J[g]%in%data$this.j[g,1:nimbuild$n.samples[g]]
  conf$addSampler(target = y.nodes,
                  type = 'IDSampler',control = list(g=g,M=M[g],J=J[g],this.j=data$this.j[g,1:nimbuild$n.samples[g]],trapup=trapup,
                                                    j.indicator=j.indicator,calcNodes=calcNodes,z.nodes=NA,
                                                    n.samples=nimbuild$n.samples[g],
                                                    cluster.ups=cluster.ups,local.eval=FALSE,swap.rad.multiplier=6),
                  silent = TRUE)
}


# ###Two *optional* sampler replacements:
#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#use block update for lam0 and sigma. bc correlated posteriors.
conf$removeSampler(c("lam0.fixed","sigma.fixed"))
conf$addSampler(target = c(paste("lam0.fixed"),paste("sigma.fixed")),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



# Run the model.
start.time2<-Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$N #true realized N
data$n #true number of inds captured


aa=mcmc.list(mcmc(mvSamples[2:7501,]),mcmc(mvSamples2[2:7501,]))
plot(aa)
