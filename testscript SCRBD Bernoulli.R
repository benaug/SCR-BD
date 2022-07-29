library(nimble)
source("simSCRBD.R")
source("init.data.R")
source("NimbleModel SCRBD Bernoulli.R")
source("NimbleFunctions SCRBD Bernoulli.R")
source("sSampler.R")

B=50 ## nBorn
t.mu=50 ## mean birth date
t.sd=25 ## sd of birth dates
omega=50 ## mean lifetime under exponential lifetime model
K=200

#simulate some data
N=50
p0=0.5
sigma=0.50
buff=3 #state space buffer. Should be at least 3 sigma.
X<- as.matrix(expand.grid(3:11,3:11))
J=nrow(X)
K2D=matrix(1,J,K)

data=simSCRBD(B=B,p0=p0,sigma=sigma,X=X,buff=3,obstype="bernoulli",
             t.mu=t.mu,t.sd=t.sd,K=K,K2D=K2D)

M=100
inits=list(t.mu=t.mu,t.sd=t.sd,omega=omega,psi=0.5)
nimbuild=init.data(data=data,inits=inits,M=M)

Niminits <- list(s=nimbuild$s,b=nimbuild$b,bday=nimbuild$bday,dday=nimbuild$dday,
                 lifetime=nimbuild$lifetime,
                 p0=0.25,sigma=0.25,omega=omega,t.mu=t.mu,t.sd=t.sd)

J=nrow(data$X)
constants<-list(J=J,K=data$K,M=M,xlim=data$xlim,ylim=data$ylim,X=as.matrix(data$X),K2D=data$K2D)

b.data=rep(NA,M)
b.data[1:data$n]=1
Nimdata<-list(y=nimbuild$y,b=b.data,first.seen=nimbuild$first.seen,last.seen=nimbuild$last.seen)

# set parameters to monitor
parameters=c("B","t.mu","t.sd","omega","p0","sigma")
parameters2=c("bday","dday","lifetime")
# parameters2="s"
# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, monitors2=parameters2,thin2=25,
                      useConjugacy = TRUE)


conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=data$xlim,ylim=data$ylim,
                                                 scale=1,adaptive=TRUE),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned, unless adaptive=FALSE.
}


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-compilation run time

mvSamples = as.matrix(Cmcmc$mvSamples)
library(coda)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

