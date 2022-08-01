init.data.Multisession<-
  function(data=data,M=M,inits=inits,plot=TRUE){
    N.session=length(data$y)
    if(length(M)!=N.session)stop("M and data must be of length 'N.session'")
    M.max=max(M)
    J=unlist(lapply(data$X,function(x){nrow(x)}))
    K=data$K
    
    init.session=vector("list",N.session)
    #split inits by session
    inits.use=vector("list",N.session)
    inits$psi=rep(0.5,N.session) #single session uses psi to init. Giving it one per session here.
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
      data.use=list(y=data$y[[g]],X=data$X[[g]],xlim=data$xlim[g,],ylim=data$ylim[g,],
                    J=J[g],K=K[g],K2D=data$K2D[[g]],n=data$n[g],buff=data$buff[g])
      init.session[[g]]=init.data(data.use,inits.use[[g]],M=M[g],plot=plot)
    }
    
    #reformat
    maxM=max(M)
    s=array(NA,dim=c(N.session,maxM,2))
    b=matrix(NA,N.session,maxM)
    bday=dday=lifetime=first.seen=last.seen=matrix(NA,N.session,maxM)
    y=array(NA,dim=c(N.session,maxM,max(J),max(K)))
    xlim=ylim=matrix(NA,N.session,2)
    n=rep(NA,N.session)
    K2D=array(NA,dim=c(N.session,max(J),max(K)))
    for(g in 1:N.session){
      s[g,1:M[g],]=init.session[[g]]$s
      b[g,1:M[g]]=init.session[[g]]$b
      y[g,1:M[g],1:J[g],1:K[g]]=init.session[[g]]$y
      bday[g,1:M[g]]=init.session[[g]]$bday
      dday[g,1:M[g]]=init.session[[g]]$dday
      lifetime[g,1:M[g]]=init.session[[g]]$lifetime
      first.seen[g,1:M[g]]=init.session[[g]]$first.seen
      last.seen[g,1:M[g]]=init.session[[g]]$last.seen
      xlim[g,]=init.session[[g]]$xlim
      ylim[g,]=init.session[[g]]$ylim
      n[g]=data$n[g]
      K2D[g,1:J[g],1:K[g]]=init.session[[g]]$K2D
    }
    #put X in ragged array
    X.new=array(NA,dim=c(N.session,max(J),2))
    for(g in 1:N.session){
      X.new[g,1:J[g],]=data$X[[g]]
    }
    
    return(list(y=y,X=X.new,xlim=xlim,ylim=ylim,K=K,J=J,n=n,K2D=K2D,s=s,b=b,
                bday=bday,dday=dday,first.seen=first.seen,last.seen=last.seen,
                bday=bday,dday=dday,lifetime=lifetime))
  }
