e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getArea <- function(X=X,buff=buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(a in 1:N.session){
    xlim=c(min(X[[a]][,1]),max(X[[a]][,1]))+c(-buff[[a]],buff[[a]])
    ylim=c(min(X[[a]][,2]),max(X[[a]][,2]))+c(-buff[[a]],buff[[a]])
    area[a]=diff(xlim)*diff(ylim)
  }
  return(area)
}


simUSCR.multisession <-
  function(N.session=NA,lambda=NA,lam0=NA,p0=NA,sigma=NA,theta=NA,
           K=NA,X=NA,buff=NA,obstype="poisson"){
    if(length(lambda)!=N.session)stop("lambda must be of length N.session")
    if(obstype%in%c("poisson","negbin")){
      if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    }
    if(obstype=="bernoulli"){
      if(length(p0)!=N.session)stop("lam0 must be of length N.session")
    }
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
    if(length(buff)!=N.session)stop("buff must be of length N.session")
    if(obstype=="negbin"){
      if(length(theta)!=N.session)stop("theta must be of length N.session")
    }else{
      theta=rep(theta,N.session)
    }

    #realized N
    N=rpois(N.session,lambda)
    
    library(abind)
    xlim=ylim=matrix(NA,N.session,2)
    s=D=vector("list",N.session)
    J=rep(NA,N.session)
    
    for(g in 1:N.session){
      X[[g]]=as.matrix(X[[g]])
      xlim[g,]=c(min(X[[g]][,1]),max(X[[g]][,1]))+c(-buff[g],buff[g])
      ylim[g,]=c(min(X[[g]][,2]),max(X[[g]][,2]))+c(-buff[g],buff[g])
      s[[g]]<- cbind(runif(N[g], xlim[g,1],xlim[g,2]), runif(N[g],ylim[g,1],ylim[g,2]))
      D[[g]]<- e2dist(s[[g]],X[[g]])
      J[g]=nrow(X[[g]])
    }
    
    #simulate sessions one at a time
    data=vector("list",N.session)
    for(g in 1:N.session){
      data[[g]]=simUSCR(N=N[g],lam0=lam0[g],p0=p0[g],sigma=sigma[g],K=K[g],X=X[[g]],buff=buff[g],
                        obstype=obstype,theta=theta[g])
    }
    
    #combine session data
    n.samples=rep(NA,N.session)
    for(g in 1:N.session){
      n.samples[g]=length(data[[g]]$this.j)
    }
    n.samples.max=max(n.samples)
    n=rep(NA,N.session)
    this.j=this.k=ID=matrix(NA,N.session,n.samples.max)
    y=s=vector("list",N.session)
    for(g in 1:N.session){
      n[g]=data[[g]]$n
      this.j[g,1:n.samples[g]]=data[[g]]$this.j
      this.k[g,1:n.samples[g]]=data[[g]]$this.k
      ID[g,1:n.samples[g]]=data[[g]]$ID
      y[[g]]=data[[g]]$y
      s[[g]]=data[[g]]$s
    }
    
    out<-list(this.j=this.j,this.k=this.k, #observed data
              y=y,s=s, ID=ID,N=N,n=n,#true data
              X=X,K=K,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }