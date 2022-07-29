GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), b=double(0)){ 
    returnType(double(1))
    if(b==0) return(rep(0,J))
    if(b==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBernoulliMatrix <- nimbleFunction(
  run = function(x = double(2), pd = double(1),first = double(0),last = double(0), 
                 first.seen = double(0),last.seen = double(0), b = double(0), K = double(0), K2D=double(2),
                 log = integer(0)) {
    returnType(double(0))
    if(b==0){
      return(0)
    }else if(first>first.seen|last<last.seen){ #test here so we can sum from first:last below
      return(-Inf)
    }else{
      J=nimDim(x)[1]
      K=nimDim(x)[2]
      logProb = 0
      if(last>0&first<=K&last>=first){
        for(j in 1:J){
          for(k in first:last){
            logProb = logProb + dbinom(x[j,k],size=1,prob=pd[j]*K2D[j,k],log = TRUE)
          }
        }
      } #else not exposed to capture
      return(logProb)
    }
  }
)
#make dummy random vector generator to make nimble happy
rBernoulliMatrix <- nimbleFunction(
  run = function(n = integer(0),pd = double(1),first = double(0),last = double(0), first.seen = double(0),last.seen = double(0),
                 b = double(0),K = double(0), K2D = double(2)) {
    returnType(double(2))
    J=nimDim(pd)[1]
    out=matrix(0,J,K)
    return(out)
  }
)