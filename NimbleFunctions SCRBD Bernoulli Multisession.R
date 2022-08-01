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


#Required custom update for number of calls
bSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    J <- control$J
    M <- control$M
    inds.detected <- control$inds.detected
    b.ups <- control$b.ups
    xlim <- control$xlim
    ylim <- control$ylim
    y.nodes <- control$y.nodes
    pd.nodes <- control$pd.nodes
    B.node <- control$B.node
    b.nodes <- control$b.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:b.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown=rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject=FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all b's currently on
        b.on=which(model$b[g,1:M]==1)
        n.b.on=length(b.on)
        pick=rcat(1,rep(1/n.b.on,n.b.on)) #select one of these individuals
        pick=b.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject=TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          
          #get initial logprobs for B and y
          lp.initial.B <- model$getLogProb(B.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new B/b
          model$B[g] <<-  model$B[g] - 1
          model$b[g,pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for B and y
          lp.proposed.B <- model$calculate(B.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.B + lp.proposed.y) - (lp.initial.B + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["B",1][g] <<- model[["B"]][g]
            for(j in 1:J){
              mvSaved["pd",1][g,pick,j] <<- model[["pd"]][g,pick,j]
            }
            mvSaved["b",1][g,pick] <<- model[["b"]][g,pick]
          }else{
            model[["B"]][g] <<- mvSaved["B",1][g]
            for(j in 1:J){
              model[["pd"]][g,pick,j] <<- mvSaved["pd",1][g,pick,j]
            }
            model[["b"]][g,pick] <<- mvSaved["b",1][g,pick]
            model$calculate(y.nodes[pick])
            model$calculate(B.node)
          }
        }
      }else{#add
        if(model$B[g] < M){ #cannot update if b maxed out. Need to raise M
          b.off=which(model$b[g,1:M]==0)
          n.b.off=length(b.off)
          pick=rcat(1,rep(1/n.b.off,n.b.off)) #select one of these individuals
          pick=b.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.B <- model$getLogProb(B.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/b
          model$B[g] <<-  model$B[g] + 1
          model$b[g,pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.B <- model$calculate(B.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.B + lp.proposed.y) - (lp.initial.B + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["B",1][g] <<- model[["B"]][g]
            for(j in 1:J){
              mvSaved["pd",1][g,pick,j] <<- model[["pd"]][g,pick,j]
            }
            mvSaved["b",1][g,pick] <<- model[["b"]][g,pick]
          }else{
            model[["B"]][g] <<- mvSaved["B",1][g]
            for(j in 1:J){
              model[["pd"]][g,pick,j] <<- mvSaved["pd",1][g,pick,j]
            }
            model[["b"]][g,pick] <<- mvSaved["b",1][g,pick]
            model$calculate(y.nodes[pick])
            model$calculate(B.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)