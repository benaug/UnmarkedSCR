#CDL and marginal samplers for Ramsey et al. (2015) model. Binary i x j x k detections reduced to j x k presence/absence
#sigma fixed to simulated value in model file. AC sampler tuning parameter set to fixed value.

#CDL version first. By this I mean we update the latent capture histories instead of marginalize over ID
library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR Bernoulli-Ramsey.R")
source("NimbleFunctions USCR Bernoulli-Ramsey.R")
source("sSampler.R")

#make sure to run this line!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 

#simulate some data
N=38
p0=0.15 #simulator treats p0 as lam0 for bernoulli obsmod
sigma=0.50
K=10
buff=3 #state space buffer. Should be at least 3 sigma.
X<- expand.grid(3:12,3:12)
obstype="bernoulli" #simulate from bernoulli obstype

#Simulate some data
data=simUSCR(N=N,p0=p0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)#not used in this sampler for 2D data, but required if using 3D data

#What if we observe even less information? Say only j x k binary data? Like Ramsey et al. (2015)
y.jk=apply(data$y.true,c(2,3),sum)
y.jk[y.jk>1]=1

#add to data
data$y.jk=y.jk

#Data augmentation level
M=150

#trap operation matrix
J=nrow(X)
K1D=rep(K,J)

inits=list(p0=0.25,sigma=0.5,psi=0.5)#initial values for lam0, sigma, and psi to build data
nimbuild=init.data.USCR(data=data,M=M,inits=inits,obstype="ramsey") #build with ramsey obstype

#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,y.true2D=nimbuild$y.true2D,p0=inits$p0,sigma=inits$sigma,y.true3D=nimbuild$y.true3D)

#constants for Nimble
constants<-list(M=M,J=J,K=data$K,K1D=K1D,xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata<-list(z=rep(NA,M),X=as.matrix(data$X))

# set parameters to monitor
parameters<-c('psi','p0','sigma','N','capcounts')
nt=1 #thinning rate


# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
conf$removeSampler("y.true2D")
conf$addSampler(target = paste0("y.true2D[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K=K,this.j=nimbuild$jk.idx[,1],this.k=nimbuild$jk.idx[,2],
                                                  y.jk=data$y.jk),
                silent = TRUE)


# ###Two *optional* sampler replacements:
#sSampler is a RW block update for the x and y locs with no covariance,
#AND only tuned for when z=1. When z=0, it draws from the prior, assumed to be uniform. 
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=0.125),silent = TRUE,adaptive=FALSE)
}

#scale parameter set to 1/4 simulated sigma, unless you changed that...
#set adaptive=TRUE to let nimble tune, but often tuned poorly



# #use block update for p0 and sigma. bc correlated posteriors.
# conf$removeSampler(c("p0","sigma"))
# conf$addSampler(target = c(paste("p0"),paste("sigma")),
#                 type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

#true number of capture events
sum(data$y.true)


#compare to marginal obsmod
source("NimbleModel USCR Bernoulli-Ramsey Marginal.R")
#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,p0=inits$p0,sigma=inits$sigma)

#constants for Nimble
constants<-list(M=M,J=J,K=data$K,K1D=K1D,xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata<-list(y.j=rowSums(data$y.jk),z=rep(NA,M),X=as.matrix(data$X))

# set parameters to monitor
parameters<-c('psi','p0','sigma','N')
nt=1 #thinning rate


# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE) 

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

mvSamples2 = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples2[2:nrow(mvSamples2),]))

