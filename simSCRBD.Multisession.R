e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getArea = function (X, buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(g in 1:N.session){
    area[g]=diff(range(X[[g]][,1])+c(-buff[g],buff[g]))*diff(range(X[[g]][,2])+c(-buff[g],buff[g]))
  }
  return(area)
}

simSCRBD.Multisession<-
  function(N.session=3,lambda=NA,p0=NA,lam0=NA,sigma=NA,theta=NA,K=NA,K2D=NA,X=NA,
           omega=NA,t.mu=NA,t.sd=NA,buff=NA,obstype="bernoulli",plot=TRUE){
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
    
    #realized B
    B=rpois(N.session,lambda)
    
    library(abind)
    xlim=ylim=matrix(NA,N.session,2)
    s=D=vector("list",N.session)
    J=rep(NA,N.session)
    
    for(g in 1:N.session){
      X[[g]]=as.matrix(X[[g]])
      xlim[g,]=c(min(X[[g]][,1]),max(X[[g]][,1]))+c(-buff[g],buff[g])
      ylim[g,]=c(min(X[[g]][,2]),max(X[[g]][,2]))+c(-buff[g],buff[g])
      J[g]=nrow(X[[g]])
    }
    
    #simulate sessions one at a time
    data=vector("list",N.session)
    for(g in 1:N.session){
      data[[g]]=simSCRBD(B=B[g],lam0=lam0[g],p0=p0[g],sigma=sigma[g],K=K[g],X=X[[g]],buff=buff[g],
                        obstype=obstype,theta=theta[g],omega=omega[g],t.mu=t.mu[g],t.sd=t.sd[g])
    }
    
    #combine session data
    n=rep(NA,N.session)
    y=s=bday=dday=lifetime=K2D=vector("list",N.session)
    for(g in 1:N.session){
      y[[g]]=data[[g]]$y
      s[[g]]=data[[g]]$s
      bday[[g]]=data[[g]]$bday
      dday[[g]]=data[[g]]$dday
      lifetime[[g]]=data[[g]]$lifetime
      K2D[[g]]=data[[g]]$K2D
      n[g]=data[[g]]$n
    }
    
    out<-list(y=y,s=s,B=B,n=n,bday=bday,dday=dday,lifetime=lifetime,
              X=X,K=K,K2D=K2D,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }
