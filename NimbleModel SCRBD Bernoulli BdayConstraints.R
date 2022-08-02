NimModel <- nimbleCode({
  t.mu ~ dunif(-100,500) # Mean birth date
  t.sd ~ dunif(0,500) # Variance of birth dates
  omega ~ dunif(0,1000) # Mean lifetime
  psi ~ dunif(0,1) # Data augmentation parameter, E(nBirths)=M*psi
  p0 ~ dunif(0,1) # baseline detection prob
  sigma ~ dunif(0,10) # df scale
  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    b[i] ~ dbern(psi)
    bday[i] ~ dnorm(t.mu, sd=t.sd)
    bday.constraint.data[i] ~ dconstraint(bday[i] >= bday.range[i,1] & bday[i] <= bday.range[i,2])
    lifetime[i] ~ dexp(1/omega)
    dday[i] <- bday[i]+lifetime[i]
    first[i] <- round(max(bday[i], 1)-0.49999)
    last[i] <- round(min(K, dday[i])-0.49999)
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, b=b[i])
    y[i,1:J,1:K] ~ dBernoulliMatrix(pd=pd[i,1:J],first=first[i],last=last[i],
                                  first.seen=first.seen[i],last.seen=last.seen[i],b=b[i],
                                  K=K,K2D=K2D[1:J,1:K])
  }
  B <- sum(b[1:M])
})