library(expm)

NS <- 4 # number of stages

transmat <- function(lambda, mu){
  Lam <- -diag(c(lambda, 0)+mu)
  diag(Lam[,-1]) <- lambda
  return(Lam)
}

## probability for an individual only with observed survival information
ll.surv.ind <- function(t, d, pi0, lambda, mu, Lam = NULL){
  ## t: time of last observation
  ## d: 1 = death, 0 = censored
  ## Lam: transition matrix

  if(is.null(Lam))
    Lam <- transmat(lambda, mu)

  pi.t <- pi0 %*% expm(t*Lam) # probability in stage s at time t;; !! THIS IS SLOW
  if(d)
    return(log(sum(pi.t * mu)))
  else
    return(log(sum(pi.t)))
}

## probability for an individual with perfectly observed CD4 transition times (and censoring or death)
ll.cd4.ind <- function(ts, t, d, pi0, lambda, mu){
  ## ts: time of observed CD4 transitions, ts[initial stage] = 0.
  ## t: time of last observation
  ## d: 1 = death, 0 = censored

  s0 <- which(ts==0)
  log.p <- log(pi0[s0])  # probability that s is initial stage

  if(s0 < NS){
    for(s in s0:(NS-1)){  # this could probably be programmed more elegantly...
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
  sum(sapply(1:length(cd4dat$t), function(i) ll.cd4.ind(cd4dat$ts[i,], cd4dat$t[i], cd4dat$d[i], pi0, lambda, mu)))
}

ll.survdat <- function(pi0, lambda, mu, survdat){
  Lam <- transmat(lambda, mu)
  return(sum(sapply(1:length(survdat$t), function(i) ll.surv.ind(survdat$t[i], survdat$d[i], pi0, lambda, mu, Lam))))
}

ll <- function(theta, cd4dat, survdat=NULL){
  pi0 <- c(1-sum(theta[1:(NS-1)]), theta[1:(NS-1)])
  lambda <- exp(theta[(NS-1) + 1:(NS-1)])
  mu <- exp(theta[2*(NS-1) + 1:NS])

  ll.cd4 <- ll.cd4dat(pi0, lambda, mu, cd4dat)

  if(!is.null(survdat))
    ll.surv <- ll.survdat(pi0, lambda, mu, survdat)
  else
    ll.surv <- 0

  return(ll.cd4+ll.surv)
}
