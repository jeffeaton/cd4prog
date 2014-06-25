source("likelihood.R")

###################
####  Example  ####
###################

pi0 <- c(0.58, 0.23, 0.16, 0.03)
lambda <- 1/c(6.37, 2.86, 3.54)
mu <- c(0.004, 0.010, 0.026, 0.435)
Lam <- transmat(lambda, mu)

theta.true <- c(pi0[-1], log(lambda), log(mu))

sim.dat <- function(n.i, pi0, lambda, mu, cens.rate=NULL){
  
  s0 <- sample(1:4, n.i, TRUE, pi0)
  prog.t <- t(matrix(rexp(n.i*(NS-1), lambda), NS-1))
  mort.t <- t(matrix(rexp(n.i*NS, mu), NS))
  
  ts <- matrix(NA, n.i, NS)
  t <- rep(NA, n.i)
  d <- rep(NA, n.i)
  
  for(s in 1:(NS-1)){
    ts[s0==s,s] <- 0
    ts[,s+1] <- ts[,s] + prog.t[,s]
    which.mort <- is.na(t) & mort.t[,s] < prog.t[,s]
    ts[which.mort,s+1] <- NA
    t[which.mort] <- ts[which.mort,s] + mort.t[which.mort,s]
    d[which.mort] <- 1
  }
  ts[s0==NS,NS] <- 0
  which.mort <- is.na(t)
  t[which.mort] <- ts[which.mort,NS] + mort.t[which.mort,NS]
  d[which.mort] <- 1

  if(!is.null(cens.rate)){
    cens.t <- rexp(n.i, cens.rate)
    for(s in 1:NS)
      ts[ts[,s] > cens.t,s] <- NA
    d[t > cens.t] <- 0
    t[t > cens.t] <- cens.t[t > cens.t]
  }

  return(list(s0=s0, ts=ts, t=t, d=d))
}

## censoring the last CD4 observation if death
update.progonly <- function(dat){
  dat$t[dat$d == 1] <- apply(dat$ts[dat$d == 1,], 1, max, na.rm=TRUE)
  dat$d[dat$d == 1] <- 0
  return(dat)
}

update.survonly <- function(dat){return(dat)}



### 1. perfectly observed data ###

dat1 <- sim.dat(1000, pi0, lambda, mu)
ll.cd4dat(pi0, lambda, mu, dat1)
mod1.fit <- optim(c(pi0[-1], log(lambda), log(mu)), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat1, hessian=TRUE)
mod1.fit.bfgs <- optim(c(pi0[-1], log(lambda), log(mu)), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat1, hessian=TRUE, method="BFGS")

mod1.fit.alt <- optim(c(0.25, 0.25, 0.25, log(c(0.3, 0.3, 0.3)), log(c(0.1, 0.1, 0.1, 0.1))), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat1, hessian=TRUE)
mod1.fit.alt <- optim(c(0.25, 0.25, 0.25, log(c(0.3, 0.3, 0.3)), log(c(0.1, 0.1, 0.1, 0.1))), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat1, hessian=TRUE, method="BFGS")
        
solve(-mod1.fit$hessian)
mean(cd4dat$t)


### 2. perfectly observed cd4data with random censoring ###

dat2 <- sim.dat(1000, pi0, lambda, mu, 1/5)
ll.cd4dat(pi0, lambda, mu, dat2)
mod2.fit <- optim(c(pi0[-1], log(lambda), log(mu)), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat2, hessian=TRUE)
mod2.fit.bfgs <- optim(c(pi0[-1], log(lambda), log(mu)), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat2, hessian=TRUE, method="BFGS")

mod2.fit.alt <- optim(c(0.25, 0.25, 0.25, log(c(0.3, 0.3, 0.3)), log(c(0.1, 0.1, 0.1, 0.1))), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat2, hessian=TRUE, method="BFGS")


### 3. perfectly observed CD4 transition and mortality, separate ###

dat3.cd4 <- sim.dat(1000, pi0, lambda, mu)
dat3.surv <- sim.dat(1000, pi0, lambda, mu)
dat3.cd4 <- update.progonly(dat3.cd4)
dat3.surv <- update.survonly(dat3.surv)

ll(theta.true, dat3.cd4)
ll(theta.true, dat3.cd4, dat3.surv)
mod3.fit <- optim(theta.true, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat3.cd4, survdat=dat3.surv, hessian=TRUE)
mod3.fit.bfgs <- optim(theta.true, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat3.cd4, survdat=dat3.surv, hessian=TRUE, method="BFGS")

mod3.fit.alt <- optim(c(0.25, 0.25, 0.25, log(c(0.3, 0.3, 0.3)), log(c(0.1, 0.1, 0.1, 0.1))), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat3.cd4, survdat=dat3.surv, hessian=TRUE, method="BFGS")

cbind(theta.true, mod3.fit.bfgs$par + sqrt(diag(solve(-mod3.fit.bfgs$hessian))) %o% c(0, -1,1)*qnorm(0.975))
## poor estimate of mortality for higher CD4 stages, may overestimate mortality from last stage (needs replication on lots of simulated data)


### 4. separately observed CD4 transition and mortality, with censoring

dat4.cd4 <- sim.dat(1000, pi0, lambda, mu, 1/5)
dat4.surv <- sim.dat(1000, pi0, lambda, mu, 1/10)
dat4.cd4 <- update.progonly(dat4.cd4)
dat4.surv <- update.survonly(dat4.surv)

ll(theta.true, dat4.cd4)
ll(theta.true, dat4.cd4, dat4.surv)
mod4.fit <- optim(theta.true, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat4.cd4, survdat=dat4.surv, hessian=TRUE)
mod4.fit.bfgs <- optim(theta.true, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat4.cd4, survdat=dat4.surv, hessian=TRUE, method="BFGS")

mod4.fit.alt <- optim(c(0.25, 0.25, 0.25, log(c(0.3, 0.3, 0.3)), log(c(0.1, 0.1, 0.1, 0.1))), ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=dat4.cd4, survdat=dat4.surv, hessian=TRUE, method="BFGS")

cbind(theta.true, mod4.fit.bfgs$par + sqrt(diag(solve(-mod4.fit.bfgs$hessian))) %o% c(0, -1,1)*qnorm(0.975))
## poor estimate of mortality for higher CD4 stages, may underestimate mortality from last stage (needs replication on lots of simulated data)


### 5. some amount of overlap ###


## Note: correlation between initial CD4 and rate of decline could be modelled be replacing pi0 by pi0_i and modelling correlation with lambda_i and pi0_i


  ## Note: knowing that somebody starts treatment at > CD4 might give some survival information
