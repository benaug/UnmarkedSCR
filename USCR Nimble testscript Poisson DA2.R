#this version uses an alternative data augmentation approach that runs faster and allows a poisson
#prior on N. Mixing sucks in these models, may need to run 100K or more iterations.

library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR Poisson DA2.R")
source("NimbleFunctions USCR Poisson DA2.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
N <- 38
lam0 <- 0.25
sigma <- 0.5 #change prior if you change sigma, set up with informative prior around 0.5
K <- 10
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11)
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
diff(xlim)*diff(ylim) #state space area
obstype <- "poisson"

#Simulate some data
data <- simUSCR(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)#not used in this sampler for 2D data, but required if using 3D data

#Data augmentation level
M <- 300

#trap operation matrix
J <- nrow(X)
K1D <- rep(K,J)

inits <- list(lam0=1,sigma=1)#ballpark initial values for lam0 and sigma to build data
nimbuild <- init.data.USCR(data=data,M=M,inits=inits)

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #must initialize N to be the sum of z init
                 lambda.N=sum(nimbuild$z), #initializing lambda.N to be consistent with N.init
                 s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true2D),
                 y.true=nimbuild$y.true2D,lam0=inits$lam0,sigma=inits$sigma)

#constants for Nimble
constants <- list(M=M,J=J,K1D=K1D,n.samples=nimbuild$n.samples,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=J),
              ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters <- c('lambda.N','lam0','sigma','N','n')
nt <- 1 #thinning rate

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2 <- 25 #thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c("lambda.N","lam0","sigma")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,
                      monitors2=parameters2,thin2=nt2,useConjugacy = FALSE) 

# conf$removeSampler("y.true")
IDups <- 2 #how many times to propose updates for each sample ID per iteration. No idea what is optimal in specific scenarios.
trapup <- unique(data$this.j)
j.indicator <- 1:J%in%data$this.j
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,this.j=data$this.j,trapup=trapup,
                                                  j.indicator=j.indicator,
                                                  IDups=IDups),
                silent = TRUE)


z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 50% of M here.
#nodes used for update
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 y.nodes=y.nodes,lam.nodes=lam.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

#"sSampler", which is a RW block update for the x and y locs with no covariance,
#and only tuned for when z=1. When z=0, it draws from the prior, assumed to be uniform. 
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#use block update for  correlated posteriors. Can use "tries" to control how many times per iteration
conf$addSampler(target = c("lam0","sigma","lambda.N"),
                type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

cor(mcmc(mvSamples[250:nrow(mvSamples),]))
#n is number of individuals captured. True value:
data$n


####Look an posterior pairwise sample match probs
#assuming you monitored ID in 2nd monitor
library(MCMCglmm)
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
burnin <- 5
IDpost <- posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),]))
#For simulated data sets, comparing posterior mode ID to truth.
#Numbers will not be the same, but all samples with same true ID will have
#same ID in posterior mode when posterior mode is exactly correct. Numbers just don't match up.
cbind(data$ID,round(IDpost))

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
n.samples <- length(data$this.j)
n.iter <- nrow(mvSamples2)
pair.probs <- matrix(NA,n.samples,n.samples)
for(i in 1:n.samples){
  for(j in 1:n.samples){
    count <- 0
    for(iter in burnin:n.iter){
      count <- count+1*(mvSamples2[iter,j]==mvSamples2[iter,i])
    }
    pair.probs[i,j] <- count/(n.iter-burnin+1)
  }
}

this.samp <- 1 #sample number to look at
round(pair.probs[this.samp,],3) #probability this sample is from same individual as all other samples
round(pair.probs[this.samp,data$ID==data$ID[this.samp]],3) #for simulated data, these are the other samples truly from same individual
