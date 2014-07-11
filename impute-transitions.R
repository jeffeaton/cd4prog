library(mvtnorm)

impute.cd4trans <- function(pat.gam, cd4.thresh=c(500, 350, 200), rand=TRUE){
  ## pat.gam: GAM model output from mcgv package with spline output for individual patient CD4 trajectory
  ## rand: whether to sample random set of coefficients from covariance matrix (TRUE), or use coefficient point estimates (FALSE)

  if(rand)
    coef.sim <- c(rmvnorm(1, coef(pat.gam), vcov(pat.gam)))
  else
    coef.sim <- coef(pat.gam)

  t0 <- min(pat.gam$smooth[[1]]$xp)
  t.final <- max(pat.gam$smooth[[1]]$xp)
  pred.sim <- function(x) sum(predict(pat.gam, list(cd4.t=x), type="lpmatrix") * coef.sim)
  
  cd4.range <- c(pred.sim(t0), pred.sim(t.final))  # !! THIS ASSUMES MONOTONIC DECLINE

  ts <- numeric(length(cd4.thresh)+1)
  ts[1] <- t0
  for(i in 1:length(cd4.thresh))
    ts[i+1] <- ifelse(cd4.range[1] < cd4.thresh[i], t0,          # already below threshold at t0
                      ifelse(cd4.range[2] > cd4.thresh[i], NA,   # never reached CD4 threshold
                             optimize(function(x) (pred.sim(x)-cd4.thresh[i])^2, lower=t0, upper=t.final, tol=0.01)$minimum))

  ts[1:length(ts) < max(which(ts==t0))] <- NA

  return(ts)
}
