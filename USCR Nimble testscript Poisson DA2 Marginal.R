#this version uses an alternative data augmentation approach that runs faster and allows a poisson
#prior on N. Also uses observation model marginalized over individuals with results from Herliansyah et al.
#(2024) https://link.springer.com/article/10.1007/s13253-023-00598-3
#Marginal only possible with Poisson detections

library(nimble)
library(coda)
source("simUSCR.R")
source("init.data.USCR.R")
source("NimbleModel USCR Poisson DA2 Marginal.R")
source("NimbleFunctions USCR Poisson DA2 Marginal.R")
source("sSampler Poisson Marginal.R")

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
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #z and N inits must be consistent
                 lambda.N=sum(nimbuild$z), #converges faster if you set N and lambda.N at similar values
                 s=nimbuild$s,sigma=inits$sigma,lam0=inits$lam0)

#constants for Nimble
J <- nrow(data$X)
constants <- list(M=M,J=J,xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
y.noID <- tabulate(nimbuild$this.j,J) #number of unidentified counts by trap
Nimdata <- list(X=as.matrix(X),K1D=K1D,y.noID=y.noID)


# set parameters to monitor
parameters <- c('lambda.N','lam0','sigma','N')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c("lambda.N","lam0","sigma")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE) 

###*required* sampler replacement
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
#nodes used for update
y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
bigLam.nodes <- Rmodel$getDependencies("bigLam") #only need this in calcNodes
lam.noID.nodes <- Rmodel$getDependencies("lam.noID")
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                 y.noID.nodes=y.noID.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

#must use this activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,J=J,scale=1),silent = TRUE)
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
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

cor(mcmc(mvSamples[250:nrow(mvSamples),]))