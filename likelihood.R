library(expm)

transmat <- function(lambda, mu){
  Lam <- -diag(c(lambda, 0)+mu)
  diag(Lam[,-1]) <- lambda
  return(Lam)
}


## CDF for convoluation of exponential distribution with unequal rates
pexpsum <- Vectorize(function(x, rates){
  if(!length(rates))
    return(1)
  return(1-sum(sapply(seq_along(rates), function(i) prod(rates[-i])/prod(rates[-i] - rates[i])*exp(-x*rates[i]))))
  }, "x")

## calculates expm(t*Lam), but appears slower than expm(t*Lam)
progmat <- function(t, lambda, mu){

  NS <- length(mu)
  rate <- c(lambda, 0) + mu
  psurv <- c(lambda, 0)/rate

  pmat <- matrix(0, NS, NS)

  for(i in 1:NS)
    for(j in i:NS){
      ii <- (i-1)+seq_len(j-i+1)
      iip <- (i-1)+seq_len(j-i)
      pmat[i,j] <- (pexpsum(t, rate[iip]) - pexpsum(t, rate[ii]))*prod(psurv[iip])
    }

  return(pmat)
}
      

## probability for an individual only with observed survival information
ll.surv.ind <- function(t, d, t0=0, pi0, lambda, mu, Lam = NULL){
  ## t: time of last observation
  ## d: 1 = death, 0 = censored
  ## Lam: transition matrix

  if(is.null(Lam))
    Lam <- transmat(lambda, mu)

  pi.t <- pi0 %*% expm(t*Lam) # probability in stage s at time t;; !! THIS IS SLOW
  ## pi.t <- expAtv(t(Lam), pi0, t)$eAtv # probability in stage s at time t;; appears slower for NS = 4
  
  psurv.t0 <- ifelse(t0 > 0, sum(pi0 %*% expm(t0*Lam)), 1) # probability of surviving to t0

  if(d)
    return(log(sum(pi.t * mu)) - log(psurv.t0))
  else
    return(log(sum(pi.t)))
}

## probability for an individual with perfectly observed CD4 transition times (and censoring or death)
ll.cd4.ind <- function(ts, t, d, t0=0, pi0, lambda, mu, Lam = NULL){
  ## ts: time of observed CD4 transitions, ts[initial stage] = 0.
  ## t: time of last observation
  ## d: 1 = death, 0 = censored

  NS <- length(pi0)
  
  if(is.null(Lam))
    Lam <- transmat(lambda, mu)

  s.t0 <- min(which(!is.na(ts)))
  pi.t0 <- if(t0 > 0) as.vector(pi0 %*% expm(t0*Lam)) else pi0
  log.p <- log(pi.t0[s.t0])  # probability of s.t0 is initial stage

  if(s.t0 < NS){
    for(s in s.t0:(NS-1)){  # this could probably be programmed more elegantly...
      if(!is.na(ts[s+1])){
        dts <- ts[s+1] - ts[s]  # duration in stage s
        log.p <- log.p - dts*(lambda[s]+mu[s]) + log(lambda[s]) # P(trans after dts) = exp(-dts*(lambda[s]+mu[s])) * lambda[s]
      } else {
        ## censoring / death interval
        dt <- t - ts[s]
        if(d)
          log.p <- log.p - dt*(lambda[s]+mu[s])  + log(mu[s]) # P(die at t | ts[s]) = exp(-dt*(lambda[s]+mu[s])) * mu[s]
        else
          log.p <- log.p - dt*(lambda[s]+mu[s])          # P(surv ts[s] to t) = exp(-dt*(lambda[s]+mu[s]))
        return(log.p)
      }
    }
  }
  
  dt <- t - ts[NS]
  if(d)
    log.p <- log.p - dt*mu[NS] + log(mu[NS])
  else
    log.p <- log.p - dt*mu[NS]

  return(log.p)

  ## !! note: still informative ART censoring issue that I'm not sure how to deal with...
  ## !!       This likelihood assumes random censoring.
  ## !! note: assumes died from current CD4 stage
}

ll.cd4dat <- function(pi0, lambda, mu, cd4dat){
  Lam <- transmat(lambda, mu)
  if(is.null(cd4dat$t0)) cd4dat$t0 <- rep(0, length(cd4dat$t))
  sum(sapply(1:length(cd4dat$t), function(i) ll.cd4.ind(cd4dat$ts[i,], cd4dat$t[i], cd4dat$d[i], cd4dat$t0[i], pi0, lambda, mu, Lam)))
}

ll.survdat <- function(pi0, lambda, mu, survdat){
  Lam <- transmat(lambda, mu)
  if(is.null(survdat$t0)) survdat$t0 <- rep(0, length(survdat$t))
  return(sum(sapply(1:length(survdat$t), function(i) ll.surv.ind(survdat$t[i], survdat$d[i], survdat$t0[i], pi0, lambda, mu, Lam))))
}

ll <- function(theta, cd4dat, survdat=NULL, NS=4){
  pi0 <- c(1-sum(theta[1:(NS-1)]), theta[1:(NS-1)])  
  lambda <- exp(theta[(NS-1) + 1:(NS-1)])
  mu <- exp(theta[2*(NS-1) + 1:NS])

  if(any(pi0 <= 0) || any(pi0 > 1))
    return(-Inf)

  ll.cd4 <- ll.cd4dat(pi0, lambda, mu, cd4dat)

  if(!is.null(survdat))
    ll.surv <- ll.survdat(pi0, lambda, mu, survdat)
  else
    ll.surv <- 0

  return(ll.cd4+ll.surv)
}
