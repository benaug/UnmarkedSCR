e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.data.USCR.multisession=function(data=NA,M=NA,inits=inits,obstype="poisson"){
  N.session=nrow(data$this.j)
  n.samples=rowSums(!is.na(data$this.j))
  init.session=vector("list",N.session)
  
  #split inits by session
  inits.use=vector("list",N.session)
  inits$psi=rep(0.75,N.session) #single session uses psi to init. Giving it one per session here.
  parms=names(inits)
  for(g in 1:N.session){
    inits.use[[g]]=vector("list",length(parms))
    names(inits.use[[g]])=parms
    for(i in 1:length(parms)){
      inits.use[[g]][[i]]=inits[[i]][g]
    }
  }
  
  #initialize sessions one by one
  for(g in 1:N.session){
    data.use=list(this.j=data$this.j[g,1:n.samples[g]],this.k=data$this.k[g,1:n.samples[g]],
                  X=data$X[[g]],buff=data$buff[g],K=data$K[g],xlim=data$xlim[g,],ylim=data$ylim[g,])
    init.session[[g]]=init.data.USCR(data.use,inits.use[[g]],M=M[g],obstype=obstype)
  }
  J=unlist(lapply(data$X,nrow))
  K=data$K
  maxM=max(M)
  s=array(NA,dim=c(N.session,maxM,2))
  z=matrix(NA,N.session,maxM)
  ID=matrix(NA,N.session,max(n.samples))
  y.true2D=array(NA,dim=c(N.session,maxM,max(J)))
  y.true3D=array(NA,dim=c(N.session,maxM,max(J),max(K)))
  
  for(g in 1:N.session){
    s[g,1:M[g],]=init.session[[g]]$s
    z[g,1:M[g]]=init.session[[g]]$z
    ID[g,1:n.samples[g]]=init.session[[g]]$ID
    y.true2D[g,1:M[g],1:J[g]]=init.session[[g]]$y.true2D
    y.true3D[g,1:M[g],1:J[g],1:K[g]]=init.session[[g]]$y.true3D
  }
  #put X in ragged array
  X.new=array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],]=data$X[[g]]
  }
  
  return(list(s=s,z=z,ID=ID,y.true2D=y.true2D,y.true3D=y.true3D,K=K,J=J,X=X.new,
              n.samples=n.samples,this.j=data$this.j,
              xlim=data$xlim,ylim=data$ylim))
}
