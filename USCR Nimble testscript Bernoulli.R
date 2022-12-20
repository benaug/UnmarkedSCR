library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR Bernoulli.R")
source("NimbleFunctions USCR Bernoulli.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
N=38
p0=0.25 #simulator treats p0 as lam0 for bernoulli obsmod
sigma=0.50
K=5
buff=3 #state space buffer. Should be at least 3 sigma.
X<- expand.grid(3:11,3:11)
obstype="bernoulli"

#Simulate some data
data=simUSCR(N=N,p0=p0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)#not used in this sampler for 2D data, but required if using 3D data

#Data augmentation level
M=300

#trap operation matrix
J=nrow(X)
K1D=rep(K,J)

inits=list(p0=0.8,sigma=0.4,psi=0.5)#initial values for lam0, sigma, and psi to build data
nimbuild=init.data.USCR(data=data,M=M,inits=inits,obstype="bernoulli")

#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true2D),
                 y.true=nimbuild$y.true2D,p0=inits$p0,sigma=inits$sigma)

#constants for Nimble
constants<-list(M=M,J=J,K=data$K,K1D=K1D,n.samples=nimbuild$n.samples,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata<-list(y.true=matrix(NA,nrow=M,ncol=J),
              ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters<-c('psi','p0','sigma','N','n')

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

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,this.j=data$this.j,this.k=data$this.k,
                                                  n.samples=nimbuild$n.samples),
                silent = TRUE)


# ###Two *optional* sampler replacements:
#replace default activity center sampler that updates x and y locations separately with a joint update
#should be a little more efficient. Could use AFslice or block random walk.
#BUT! I suggest using "sSampler", which is a RW block update for the x and y locs with no covariance,
#AND only tuned for when z=1. When z=0, it draws from the prior, assumed to be uniform. 
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
  # block RW option
  # do not adapt covariance bc samples not deterministically linked to individuals
  # longer adapt interval to average over more data configurations for each s_i
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}




#use block update for p0 and sigma. bc correlated posteriors.
conf$removeSampler(c("p0","sigma"))
conf$addSampler(target = c(paste("p0"),paste("sigma")),
                type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



# Run the model.
start.time2<-Sys.time()
Cmcmc$run(10000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

#n is number of individuals captured. True value:
data$n


####Look an posterior pairwise sample match probs
#assuming you monitored ID in 2nd monitor
library(MCMCglmm)
mvSamples2 = as.matrix(Cmcmc$mvSamples2)
burnin=1000
IDpost=posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),]))
#For simulated data sets, comparing posterior mode ID to truth.
#Numbers will not be the same, but all samples with same true ID will have
#same ID in posterior mode when posterior mode is exactly correct. Numbers just don't match up.
cbind(data$ID,round(IDpost))

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
n.samples=length(data$this.j)
n.iter=nrow(mvSamples2)
pair.probs=matrix(NA,n.samples,n.samples)
for(i in 1:n.samples){
  for(j in 1:n.samples){
    count=0
    for(iter in burnin:n.iter){
      count=count+1*(mvSamples2[iter,j]==mvSamples2[iter,i])
    }
    pair.probs[i,j]=count/(n.iter-burnin+1)
  }
}

this.samp=1 #sample number to look at
round(pair.probs[this.samp,],3) #probability this sample is from same individual as all other samples
round(pair.probs[this.samp,data$ID==data$ID[this.samp]],3) #for simulated data, these are the other samples truly from same individual


### If using s update that only tunes when z=1...
#Crude check of acceptance rates when z=1 if you monitor s and z.
#crude bc it is overestimate counting 1 acceptance every time  a z is turned off and back on.
#can get exact with more work
burn=1000
mvSamples2 = as.matrix(Cmcmc$mvSamples2)
z.idx=grep("z",colnames(mvSamples2))
z.post=mvSamples2[burn:nrow(mvSamples2),z.idx]
s.idx=grep("s",colnames(mvSamples2))
s.post=mvSamples2[burn:nrow(mvSamples2),s.idx]

checkID=1
1-rejectionRate(mcmc(s.post[,checkID][z.post[,checkID]==1]))
