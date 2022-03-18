#This script uses a "site use - zero-truncated Poisson observation model"
#This is a model where 1) detection is a function of distance from activity center
#and 2) given detection, counts follow a zero-truncated Poisson with a fixed lambda
#parameter that is not a function of distance from the activity center.
#We assume y.det[i,j,k] ~ Bernoulli(p[i,j]) and
# y.count[i,j,k] ~ ZTPois(lambda*y.det[i,j,k]) (lambda is zeroed out if no detection)
#To fit the model, we marginalize over y.det.

#This model provides worse estimates than regular USCR. I believe
#disconnecting the counts from the activity centers and detection probability is one cause of this.
#Counts > 1 do not inform the activity center locations or detection probability. 
#Then, adding more counts per detection adds more uncertainty about how many individuals
#contributed those counts without adding any more information about the detection process.

#Because this model sucks more than regular SCR, this test script is set up with an 
#informative prior for sigma centered around the current simulation value of 0.5.
#Switch back to an uninformative prior if you want to see how the model performs in that
#case, or change the informative prior if you change the simulation value or use your
#own data set.

#A note about the MCMC. Unlike the other USCR MCMC samplers, I have chosen to update
#all samples at a trap-occasion at once instead of one by one. You cannot update just
#1 ID at a time for this model, particularly with larger lambda values--there will be 
#convergence problems and the posterior won't be fully explored. An implication of this
#is that my approach for handling partial IDs that requires updating 1 sample ID at a time
#won't work. Wah wah.

#I expect this model to also be subject to "trap saturation" where if p0 is large enough, there
#is no spatial correlation in detections and the model parameters will no longer be identifiable
#See Ramsey et al. (2015).

library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR siteUseZTPois.R")
source("NimbleFunctions USCR siteUseZTPois.R")
source("sSampler.R")

#make sure to run this line!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 

#simulate some data
N=38
p0=0.125 #baseline detection probability
lambda=0.25 #count parameter given detection
sigma=0.50
K=10
buff=3 #state space buffer. Should be at least 3 sigma.
X<- expand.grid(3:11,3:11)
obstype="siteUseZTPois"

#Simulate some data
data=simUSCR(N=N,p0=p0,lambda=lambda,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)

#Data augmentation level
M=200

#trap operation matrix
J=nrow(X)
K2D=matrix(1,J,K) #2 dimensional. Must be either 0 or 1.

inits=list(p0=0.05,lambda=1,sigma=0.75,psi=0.5)#initial values for lam0, sigma, and psi to build data
nimbuild=init.data.USCR(data=data,M=M,inits=inits,obstype="siteUseZTPois")

#inits for nimble - using 3D data here
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true3D),
                 y.true=nimbuild$y.true3D,p0=inits$p0,sigma=inits$sigma,lambda=inits$lambda)

#constants for Nimble
constants<-list(M=M,J=J,K=data$K,K2D=K2D,n.samples=nimbuild$n.samples,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata<-list(y.true=array(NA,dim=c(M,J,K)),
              ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters<-c('psi','p0','lambda','sigma','N','n')

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
# conf$removeSampler("y.true")
# conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
#                 type = 'IDSampler',control = list(M=M,J=J,K=K,this.j=data$this.j,this.k=data$this.k,
#                                                   n.samples=length(data$this.j)),
#                 silent = TRUE)

conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDSampler_jk',control = list(M=M,J=J,K=K,this.j=data$this.j,this.k=data$this.k),
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




#use block update for lam0 and sigma. bc correlated posteriors.
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c(paste("lam0"),paste("sigma")),
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
burnin=1000
mvSamples2 = as.matrix(Cmcmc$mvSamples2)
z.idx=grep("z",colnames(mvSamples2))
z.post=mvSamples2[burnin:nrow(mvSamples2),z.idx]
s.idx=grep("s",colnames(mvSamples2))
s.post=mvSamples2[burnin:nrow(mvSamples2),s.idx]

checkID=1
1-rejectionRate(mcmc(s.post[,checkID][z.post[,checkID]==1]))
