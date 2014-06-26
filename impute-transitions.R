library(mvtnorm)

sample.cd4trans <- function(pat.gam, rand=TRUE){
  ## pat.gam: GAM model output from mcgv package with spline output for individual patient CD4 trajectory
  ## rand: whether to sample random set of coefficients from covariance matrix (TRUE), or use coefficient point estimates (FALSE)

  if(rand)
    coef.sim <- c(rmvnorm(1, coef(pat.gam), vcov(pat.gam)))
  else
    coef.sim <- coef(pat.gam)

  t.final <- max(pat.gam$smooth[[1]]$xp)
  cd4.range <- c(predict(pat.gam, list(cd4.t=c(0, t.final)), type="lpmatrix") %*% coef.sim)  # !! THIS ASSUMES MONOTONIC DECLINE

  pred.sim <- function(x) sum(predict(pat.gam, list(cd4.t=x), type="lpmatrix") * coef.sim)
  
  if(cd4.range[1]>500 & cd4.range[2]<500) {
    cd4.500<-optimize(function(x) (pred.sim(x)-500)^2, lower=0, upper=t.final, tol=0.01)
    t.cd4.500<-cd4.500$minimum  # time leaving CD4 500 compartment
  } else if(cd4.range[2]>500) t.cd4.500=-1 else t.cd4.500=0
  
  if(cd4.range[1]>350&cd4.range[2]<350) {
    cd4.350<-optimize(function(x) (pred.sim(x)-350)^2,lower=0,upper=t.final,tol=0.01)
    t.cd4.350<-cd4.350$minimum
  } else if(cd4.range[2]>350) t.cd4.350=-1 else t.cd4.350=0
  
  if(cd4.range[1]>200&cd4.range[2]<200) {
    cd4.200<-optimize(function(x) (pred.sim(x)-200)^2,lower=0,upper=t.final,tol=0.01)
    t.cd4.200<-cd4.200$minimum
  } else if(cd4.range[2]>200) t.cd4.200=-1 else t.cd4.200=0

  ts <- c(0, t.cd4.500, t.cd4.350, t.cd4.200)
  ts[1:4 < max(which(ts==0))] <- NA
  ts[ts==-1] <- NA
  
  return(ts)
}
