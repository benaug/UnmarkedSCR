#This script uses a hurdle zero-truncated Poisson observation model
#This is a model where 1) detection is a function of distance from activity center
#and 2) given detection, counts follow a zero-truncated Poisson with a fixed lambda
#parameter that is not a function of distance from the activity center.
#We assume y.det[i,j,k] ~ Bernoulli(pd[i,j]) and
# y.count[i,j,k] ~ ZTPois(lambda*y.det[i,j,k]) (lambda is zeroed out if no detection)
#To fit the model, we marginalize over y.det.
#Mixing sucks in these models, may need to run 100K or more iterations.

library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR hurdleZTPois DA2.R")
source("NimbleFunctions USCR hurdleZTPois DA2.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
N <- 38
p0 <- 0.25
lambda <- 0.25
sigma <- 0.5 #change prior if you change sigma, set up with informative prior around 0.5
K <- 10
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11)
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
diff(xlim)*diff(ylim) #state space area
obstype <- "hurdleZTPois"

#Simulate some data
data <- simUSCR(N=N,p0=p0,lambda=lambda,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)

#Data augmentation level
M <- 300

#trap operation matrix
J <- nrow(X)
K2D <- matrix(1,J,K) #2 dimensional. Must be either 0 or 1.

inits <- list(p0=0.05,lambda=1,sigma=1)#initial values for p0, sigma, and psi to build data
nimbuild <- init.data.USCR(data=data,M=M,inits=inits,obstype="hurdleZTPois")

#inits for nimble - using 3D data here
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #must initialize N to be sum(z)
                 lambda.N=sum(nimbuild$z), #initializing lambda.N to be consistent with N.init
                 s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true3D),
                 y.true=nimbuild$y.true3D,
                 logit_p0=qlogis(inits$p0),sigma=inits$sigma,lambda=inits$lambda)

#constants for Nimble
constants <- list(M=M,J=J,K=data$K,K2D=K2D,n.samples=nimbuild$n.samples,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata <- list(y.true=array(NA,dim=c(M,J,K)),
              ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters <- c('lambda.N','p0','lambda','sigma','N','n')

#can also monitor a different set of parameters with a different thinning rate
# parameters2 <- c("s","z") #use this if you want to check s acceptance rates when z=1 below
parameters2 <- c("ID")
nt <- 5 #thinning rate
nt2 <- 20 #thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c("lambda.N","logit_p0","sigma","lambda")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, monitors2=parameters2,thin2=nt2,nodes=config.nodes,
                      useConjugacy = FALSE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("y.true")
IDups <- 2 #how many times to propose updates for each sample ID per iteration. No idea what is optimal in specific scenarios.
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K=K,n.samples=nimbuild$n.samples,
                                                  this.j=data$this.j,this.k=data$this.k,
                                                     K2D=K2D,IDups=IDups),
                silent = TRUE)

z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
#nodes used for update
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,",1:",K,"]"))
pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,pd.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,
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
# conf$addSampler(target = c("logit_p0","log_sigma","log_lambda.N"),
#                 type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)
conf$addSampler(target = c("logit_p0","sigma","lambda.N"),
                type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
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
burnin <- 10
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
