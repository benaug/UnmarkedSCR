e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.data.USCR=function(data=NA,M=NA,inits=inits,obstype="poisson"){
  library(abind)
  this.j=data$this.j
  this.k=data$this.k
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- data$K
  n.samples=length(this.j)
  
  #state space extent
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  psi <- inits$psi
  sigma <- inits$sigma
  #If bernoulli data, add constraints that prevent y.true[i,j,k]>1
  constraints=matrix(1,nrow=n.samples,ncol=n.samples)
  if(obstype=="bernoulli"){
    # idx=which(y.obs>0,arr.ind=TRUE)
    for(i in 1:n.samples){
      for(j in 1:n.samples){
        if(i!=j){
          # if(all(idx[i,2:3]==idx[j,2:3])){
          if(!(this.j[i]==this.j[j]&this.k[i]==this.k[j])){
            constraints[i,j]=0 #can't combine samples from same trap and occasion in binomial model
            constraints[j,i]=0
          }
        }
      }
    }
  }
  
  
  #assign random activity centers
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  z=rbinom(M,1,psi)
  D=e2dist(s, X)
  
  if(obstype=="poisson"){
    lam0<- inits$lam0
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    #Build y.true
    y.true=array(0,dim=c(M,J,K))
    ID=rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist=z*lamd[,this.j[l]]
      propdist=propdist/sum(propdist)
      ID[l]=sample(1:M,1,replace=FALSE,prob=propdist)
      y.true[ID[l],this.j[l],this.k[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
    }
  }else{
    p0<- inits$p0
    pd<- p0*exp(-D*D/(2*sigma*sigma))
    #Build y.true
    y.true=array(0,dim=c(M,J,K))
    ID=rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist=z*pd[,this.j[l]]
      propdist=propdist*(1-y.true[,this.j[l],this.k[l]]) #zero out ID's that already have a sample at this j-k
      propdist=propdist/sum(propdist)
      ID[l]=sample(1:M,1,replace=FALSE,prob=propdist)
      y.true[ID[l],this.j[l],this.k[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
    }
  }
  
  
    
  y.true2D=apply(y.true,c(1,2),sum)
    
  #Optimize s after assigning samples
  idx=which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  sigma<- inits$sigma
  D=e2dist(s, X)
  if(obstype=="bernoulli"){
    p0<- inits$p0
    pd<- p0*exp(-D*D/(2*sigma*sigma))
    ll.y=dbinom(y.true2D,K,pd*z,log=TRUE)
  }else if(obstype=="poisson"){
    lam0<- inits$lam0
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y=dpois(y.true2D,K*lamd*z,log=TRUE)
  }
  if(!is.finite(sum(ll.y)))stop("Starting obs model likelihood is not finite")
  
  return(list(y.true2D=y.true2D,y.true3D=y.true,s=s,z=z,
              ID=ID,n.samples=n.samples,xlim=xlim,ylim=ylim))
}