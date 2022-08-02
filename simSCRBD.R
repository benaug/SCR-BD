e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
simSCRBD<-
  function(B=NA,p0=NA,lam0=NA,sigma=0.50,theta=NA,X=X,buff=3,obstype="bernoulli",
           omega=NA,t.mu=NA,t.sd=NA,K=NA,K2D=NA,bday.delta=NA,plot=TRUE){
    # simulate a population of activity centers
    xlim=range(X[,1])+c(-buff,buff)
    ylim=range(X[,2])+c(-buff,buff)
    s<- cbind(runif(B, xlim[1],xlim[2]), runif(B,ylim[1],ylim[2]))
    
    #add births, deaths
    bday <- rnorm(B, t.mu,t.sd) # Times of birth events
    lifetime <- rexp(B, 1/omega) # Lifetimes
    dday <- bday+lifetime # Times of mortality events
    
    D<- e2dist(s,X)
    J<- nrow(X)
    
    if(all(is.na(K2D))){
      K2D=matrix(1,J,K)
    }else{
      if(nrow(K2D)!=J|ncol(K2D)!=K)stop("K2D must be a J x K matrix")
    }
    
    #exposure to capture by occasion
    z=matrix(0,B,K) #exposed from the day of birth to the day of death. Entire day. For Poisson, could consider partial exposure on those days.
    for(i in 1:B){
      first=round(max(bday[i], 1)-0.49999)
      last=round(min(K, dday[i])-0.49999)
      if(last<1)next
      if(first>K)next
      if(first<1){
        first=1
      }
      if(last>K){
        last=K
      }
      z[i,first:last]=1
    }
    
    # Capture individuals
    y=array(0,dim=c(B,J,K))
    if(obstype=="bernoulli"){
      pd<- p0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:B){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rbinom(1,z[i,k]*K2D[j,k],pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      lamd<- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:B){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rpois(1,lamd[i,j]*z[i,k]*K2D[j,k])
          }
        }
      } 
    }else if(obstype=="negbin"){
      for(i in 1:B){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rnbinom(1,mu=lamd[i,j],size=theta*z[i,k]*K2D[j,k])
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    
    if(plot==TRUE){
      par(mfrow=c(1,1),ask=FALSE)
      ## Temporal view of bdays and lifetimes
      plot(0, type="n", xlim=c(0, K), ylim=c(1,B),
           xlab="Time", ylab="Individual",main="Simulated Lifetimes")
      for(i in 1:B) {
        segments(bday[i], i, dday[i], i)
      }
      points(bday, 1:B, pch=16, col=3)
      points(dday, 1:B, pch=17, col=4)
      abline(v=1,lty=2)
      abline(v=K,lty=2)
    }

    #discard uncaptured inds
    caught=which(apply(y,c(1),sum)>0)
    y=y[caught,,]
    n=length(caught)
    bday=bday[caught]
    dday=dday[caught]
    lifetime=lifetime[caught]
    
    #make bday windows
    bday.range=matrix(NA,n,2)
    bday.discrete=floor(bday)
    if(!is.na(bday.delta)){
      for(i in 1:n){
        bday.range[i,1]=sample((bday.discrete[i]-bday.delta):(bday.discrete[i]),1)
        bday.range[i,2]=bday.range[i,1]+bday.delta
      }
    }
    
    out<-list(y=y,X=X,K=K,buff=buff,obstype=obstype,s=s,n=n,K=K,z=z,K2D=K2D,
              xlim=xlim,ylim=ylim,bday=bday,dday=dday,lifetime=lifetime,bday.range=bday.range)
    return(out)
  }
