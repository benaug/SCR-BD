NimModel <- nimbleCode({
  D ~ dunif(0,100) #expected birth density
  t.mu ~ dunif(-100,500) # Mean birth date
  t.sd ~ dunif(0,500) # Variance of birth dates
  omega ~ dunif(0,1000) # Mean lifetime
  p0 ~ dunif(0,1) # baseline detection prob
  sigma ~ dunif(0,10) # df scale
  for(g in 1:N.session){
    lambda[g] <- D*area[g] #expected births in session g
    B[g] ~ dpois(lambda[g]) #realized births in session g
    for(i in 1:M[g]) {
      s[g,i,1] ~ dunif(xlim[g,1], xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1], ylim[g,2])
      bday[g,i] ~ dnorm(t.mu, sd=t.sd)
      lifetime[g,i] ~ dexp(1/omega)
      dday[g,i] <- bday[g,i]+lifetime[g,i]
      first[g,i] <- round(max(bday[g,i], 1)-0.49999)
      last[g,i] <- round(min(K[g], dday[g,i])-0.49999)
      pd[g,i,1:J[g]] <- GetDetectionProb(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma, p0=p0, b=b[g,i])
      y[g,i,1:J[g],1:K[g]] ~ dBernoulliMatrix(pd=pd[g,i,1:J[g]],first=first[g,i],last=last[g,i],
                                      first.seen=first.seen[g,i],last.seen=last.seen[g,i],
                                      K=K[g],K2D=K2D[g,1:J[g],1:K[g]],b=b[g,i])
    }
  }
})