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
  
  #assign random activity centers
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  z=rbinom(M,1,psi)
  D=e2dist(s, X)
  
  if(obstype%in%c("poisson","negbin")){
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
  }else if(obstype=="bernoulli"){
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
  }else if(obstype=="siteUseZTPois"){
    p0<- inits$p0
    pd<- p0*exp(-D*D/(2*sigma*sigma))
    y.true=array(0,dim=c(M,J,K))
    ID=rep(NA,n.samples)
    for(l in 1:n.samples){
      propdist=z*pd[,this.j[l]]
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
  }else if(obstype=="negbin"){
    lam0<- inits$lam0
    theta<- inits$theta
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y=y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,]=dnbinom(y.true2D[i,],mu=lamd[i,],size=theta*K,log=TRUE)
      }
    }
  }else if(obstype=="siteUseZTPois"){
    p0<- inits$p0
    lambda<- inits$lambda
    pd<- p0*exp(-D*D/(2*sigma*sigma))
    ll.y=y.true*0
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(y.true[i,j,k]==0){
              ll.y[i,j,k]=log(1-pd[j])
            }else{
              ll.y[i,j,k]=log(pd[j]) + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
            }
          }
        }
      }
    }
  }else{
    stop("obstype not recognized")
  }
  if(!is.finite(sum(ll.y)))stop("Starting obs model likelihood is not finite")
  
  return(list(y.true2D=y.true2D,y.true3D=y.true,s=s,z=z,
                ID=ID,n.samples=n.samples,xlim=xlim,ylim=ylim))
}