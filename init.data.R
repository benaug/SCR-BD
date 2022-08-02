init.data=function(data=NA,inits=NA,M=NA,plot=TRUE){
  library(abind)
  buff<- data$buff
  X=data$X
  J=nrow(X)
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  y=data$y
  n=nrow(y)
  K=data$K
  K2D=data$K2D
  bday.range=data$bday.range
  
  t.mu=inits$t.mu
  t.sd=inits$t.sd
  omega=inits$omega
  
  #augment
  y=abind(y,array(0,dim=c(M-n,J,K)),along=1)
  y2D=apply(y,c(1,2),sum)
  #intialize s
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[which(y2D[i,]>0),],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  #initialize birth/death stuff
  y2D2=apply(y,c(1,3),sum)
  first.seen=apply(y2D2[1:n,],1,function(x){which(x>0)[1]})
  last.seen=K-apply(y2D2[1:n,K:1],1,function(x){which(x>0)[1]})+1
  bday=round(rnorm(M,t.mu,t.sd))
  lifetime=round(rexp(M,1/omega))
  dday=bday+lifetime
  #update for detected guys
  bday[1:n]=first.seen-10 #initialize to 1 day before first seen
  dday[1:n]=last.seen+10 #1 day after last seen
  #enforce bday.range if available
  #initialize to mean of bday range
  if(!all(is.na(bday.range))){
    for(i in 1:n){
      if(!is.na(bday.range[i,1])){ #assuming you see both min/max or nothing
        bday[i]=mean(bday.range[i,])
      }else{
        bday.range[i,]=c(-Inf,Inf) #if you dont' see anyting, set to -/+ infinity
      }
    }
    #need to augment bday range. set limits to -/+ infinity so no constraints
    bday.range=rbind(bday.range,matrix(0,M-n,2))
    bday.range[(n+1):M,1]=-Inf
    bday.range[(n+1):M,2]=Inf
  }else{
    bday.range=NA
  }
  
  lifetime[1:n]=dday[1:n]-bday[1:n]
  b=rbinom(M,1,inits$psi)
  b[1:n]=1
  
  #exposure to capture by occasion
  z=matrix(0,M,K) #exposed from the day of birth to day of death
  for(i in 1:M){
    first=bday[i]
    last=dday[i]
    if(first<1&last<1)next
    if(first>K&last>K)next
    if(first<1){
      first=1
    }
    if(last>K){
      last=K
    }
    z[i,first:last]=1
  }
  
  if(plot==TRUE){
    ## Temporal view of bdays and lifetimes
    plot(0, type="n", xlim=c(0,K), ylim=c(1,M),
         xlab="Time", ylab="Individual",main="Initialized Lifetimes")
    for(i in 1:M) {
      segments(bday[i], i, dday[i], i)
    }
    points(bday, 1:M, pch=16, col=3)
    points(dday, 1:M, pch=17, col=4)
    abline(v=1,lty=2)
    abline(v=K,lty=2)
    abline(h=n,lty=2)
  }
  #augment first and last seen
  first.seen=c(first.seen,rep(Inf,M-n))
  last.seen=c(last.seen,rep(-Inf,M-n))
  
  return(list(s=s,z=z,K2D=K2D,b=b,y=y,J=J,K=K,bday=bday,dday=dday,first.seen=first.seen,last.seen=last.seen,
              lifetime=lifetime,xlim=xlim,ylim=ylim,bday.range=bday.range))
  
  
}