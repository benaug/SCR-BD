#multisession version using RJMCMC

library(nimble)
source("simSCRBD.Multisession.R")
source("simSCRBD.R")
source("init.data.Multisession.R")
source("init.data.R")
source("NimbleModel SCRBD Bernoulli Multisession.R")
source("NimbleFunctions SCRBD Bernoulli Multisession.R")
source("sSampler Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE) 

N.session=3
D=0.3 # Density of births
t.mu=rep(50,N.session) # mean birth date
t.sd=rep(25,N.session) # sd of birth dates
omega=rep(50,N.session) # mean lifetime under exponential lifetime model
K=rep(200,N.session)

#simulate some data
p0=rep(0.05,N.session)
sigma=rep(0.50,N.session)
buff=rep(2,N.session) #state space buffer. Should be at least 3 sigma.
#make trapping arrays
X1=expand.grid(3:11,3:11)
X2=expand.grid(3:12,3:12)
X3=expand.grid(3:13,3:13)
X=list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area=getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda=D*area
lambda #expected births (B) in each session
J=nrow(X)

data=simSCRBD.Multisession(lambda=lambda,p0=p0,sigma=sigma,X=X,buff=buff,obstype="bernoulli",
             omega=omega,t.mu=t.mu,t.sd=t.sd,K=K,K2D=K2D)

M=c(100,110,120)
inits=list(t.mu=t.mu,t.sd=t.sd,omega=omega)

nimbuild=init.data.Multisession(data=data,inits=inits,M=M)

Niminits <- list(s=nimbuild$s,b=nimbuild$b,bday=nimbuild$bday,dday=nimbuild$dday,
                 lifetime=nimbuild$lifetime, p0=0.25,sigma=0.25,omega=50,t.mu=50,t.sd=50,D=0.5,
                 B=rowSums(nimbuild$b,na.rm=TRUE))

J=nrow(data$X)
constants<-list(N.session=N.session,J=nimbuild$J,K=nimbuild$K,M=M,xlim=nimbuild$xlim,ylim=nimbuild$ylim,X=nimbuild$X,K2D=nimbuild$K2D,
                area=area)

b.data=matrix(NA,N.session,max(M))
for(g in 1:N.session){
  b.data[g,1:nimbuild$n[g]]=1
}

Nimdata<-list(y=nimbuild$y,b=b.data,first.seen=nimbuild$first.seen,last.seen=nimbuild$last.seen)

# set parameters to monitor
parameters=c("B","D","lambda","t.mu","t.sd","omega","p0","sigma")
parameters2=c("bday","dday","lifetime")
# parameters2="s"
# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, monitors2=parameters2,thin2=25,
                      useConjugacy = TRUE)


##*required* sampler replacement
b.ups=c(25,25,25) # how many z proposals per iteration per session?
J=nimbuild$J
conf$removeSampler("B")
for(g in 1:N.session){
  #nodes used for update
  y.nodes <- Rmodel$expandNodeNames(paste("y[",g,",","1:",M[g],",1:",J[g],",1:",K[g],"]"))
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[",g,",","1:",M[g],",1:",J[g],"]"))
  B.node <- Rmodel$expandNodeNames(paste("B[",g,"]"))
  b.nodes <- Rmodel$expandNodeNames(paste("b[",g,",","1:",M[g],"]"))
  calcNodes <- c(B.node,y.nodes,pd.nodes)
  
  conf$addSampler(target = c("B"),
                  type = 'bSampler',control = list(inds.detected=1:nimbuild$n[g],b.ups=b.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,pd.nodes=pd.nodes,B.node=B.node,b.nodes=b.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

for(g in 1:N.session){
  conf$removeSampler(paste("s[",g,",1:",M[g],", 1:2]", sep=""))
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1,
                                                   adaptive=TRUE),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned, unless adaptive=FALSE.
    
  }
}


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(15,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
library(coda)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$B #true number of births (B) in each session

