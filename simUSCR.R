e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
simUSCR <-
  function(N=120,lam0=NA,p0=NA,sigma=0.50,theta=NA,lambda=NA,
           K=10,X=X,buff=3,obstype="poisson"){
    # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    J<- nrow(X)
    
    # Capture individuals
    y=array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(p0))stop("must provide p0 for bernoulli obstype")
      pd<- p0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rbinom(1,1,pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd<- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta))stop("Must provide theta for negbin obstype")
      lamd<- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rnbinom(1,mu=lamd[i,j],size=theta)
          }
        }
      } 
    }else if(obstype=="hurdleZTPois"){
      if(is.na(p0))stop("must provide p0 for siteUse_ZTpois obstype")
      if(is.na(lambda))stop("must provide lambda for siteUse_ZTpois obstype")
      library(VGAM)
      pd<- p0*exp(-D*D/(2*sigma*sigma))
      y.det=array(0,dim=c(N,J,K))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.det[i,j,k]=rbinom(1,1,pd[i,j])
            if(y.det[i,j,k]==1){
              y[i,j,k]=rzapois(1,lambda,pobs0=0)
            }
          }
        }
      }
    }else{
      stop("obstype not recognized")
    }

    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught=which(apply(y,c(1),sum)>0)
    y.true=y
    y=y[caught,,]
    if(K==1){
      y=array(y,dim=c(dim(y),1))
    }
    n=length(caught)
    n.samples=sum(y)
    
    ID=this.j=this.k=rep(NA,n.samples)
    idx=1
    for(i in 1:length(caught)){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              ID[idx]=i
              this.j[idx]=j
              this.k[idx]=k
              idx=idx+1
            }
          }
        }
      }
    }
    
    #reconstruct y to make sure this algorithm works
    y.check=y*0
    for(l in 1:n.samples){
      y.check[ID[l],this.j[l],this.k[l]]=y.check[ID[l],this.j[l],this.k[l]]+1
    }
    if(!all(y==y.check))stop("Error rebuilding data. Report bug.")
    
    out<-list(this.j=this.j,this.k=this.k,#observed data
              ID=ID,y=y,y.true=y.true,X=X,K=K,buff=buff,obstype=obstype,s=s,n=nrow(y))
    return(out)
  }