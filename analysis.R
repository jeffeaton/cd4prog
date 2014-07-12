#########################
####  Data cleaning  ####
#########################


library(foreign)


####  CD4 data  ####

bas <- read.dta("/Volumes/CD4/Data/pooled data/CASCADE (from Nikos)/received 12 February 2014/MyBas0.dta")
cd4 <- read.dta("/Volumes/CD4/Data/pooled data/CASCADE (from Nikos)/received 12 February 2014/MyCD4_0.dta")

nrow(bas)
nrow(cd4)

## exclude serohow "SC illness" and "other"
bas <- subset(bas, !serohow %in% c("SC illness", "other") & !is.na(serohow)); nrow(bas)

### seroconversion window
bas$seroco.wind <- bas$firstpos - bas$lastneg

with(bas, table(serohow, seroco.wind < 0, useNA="always"))
with(bas, table(floor(seroco.wind/365.25), serohow, useNA="always"))

## limit seroconversion window to < 2 years
bas <- subset(bas, !(serohow == "midpoint" & (seroco.wind < 0 | seroco.wind/365.25 > 2))); nrow(bas)

bas$t <- as.numeric((bas$cens_d - bas$serodate)/365.25)
bas$d <- bas$cens_r == "death"

## exclude missing CENS_D (for now !!! CHECK THIS LATER)
bas <- subset(bas, !is.na(t)); nrow(bas)

## limit to enrolment < 12 months
bas <- subset(bas, enrolltm < 12); nrow(bas)


#### merge BAS and CD4 data
dat <- merge(bas[,c("patient", "serodate", "t", "d")], cd4); nrow(dat)
dat$cd4.t <- as.numeric((dat$cd4date - dat$serodate)/365.25)

table(floor(dat$cd4.t))  # number of CD4s by years post SC
table(dat$cd4.t > dat$t) # CD4s after censor date
dat <- subset(dat, cd4.t >=0 & cd4.t <= t); nrow(dat)
dat <- dat[order(dat$patient, dat$cd4.t),]

dat.spl <- split(dat, dat$patient)

table(sapply(dat.spl, nrow)) # number of CD4 counts per patient

## limit to minimum 5 CD4 measurements
dat.spl <- dat.spl[sapply(dat.spl, nrow) >= 5]; length(dat.spl)
table(sapply(dat.spl, nrow))

sum(sapply(lapply(dat.spl, "[[", "d"), any)) # number of deaths observed

## limit to first CD4 count within 2 years !!! (RELAX THIS LATER) --> implement left censoring
## dat.spl <- dat.spl[sapply(lapply(dat.spl, "[[", "cd4.t"), min) < 2]

dat <- subset(dat, patient %in% names(dat.spl)) # dat <- do.call(rbind, dat.spl)
table(dat$cd4 < 100); sum(cd4$cd4 < 100) # lost more than half of CD4 < 100 measurements


####  Survival data  ####

alpha <- read.dta("/Volumes/CD4/Data/pooled data/Todd 2007 ALPHA survival analysis/alpha3_final.dta")

nrow(alpha)

table(!is.na(alpha$lastneg_date), !is.na(alpha$frstpos_date))
table(!is.na(alpha$seroconv_date), alpha$site)

## limit to seroconverters
alpha <- subset(alpha, !is.na(seroconv_date)); nrow(alpha) 

alpha$seroco.wind <- alpha$frstpos_date - alpha$lastneg_date
table(floor(alpha$seroco.wind/365.25))
table(floor(alpha$seroco.wind/365.25), alpha$site)

alpha <- subset(alpha, seroco.wind/365.25 < 4); nrow(alpha)  # limit to window < 4 years
alpha <- subset(alpha, alpha$exit_date > alpha$seroconv_date); nrow(alpha)  # limit to exit_date after seroconv_date

alpha <- subset(alpha, !is.na(exit_type)) # exclude missing exit_type

alpha <- subset(alpha, alpha$frstpos_date >= alpha$seroconv_date) # exclude first positive date < seroconversion date (4 cases)

alpha$t <- as.numeric((alpha$exit_date - alpha$seroconv_date)/365.25)
alpha$d <- alpha$exit_type == "death"
alpha$t0 <- as.numeric((alpha$frstpos_date - alpha$seroconv_date)/365.25)
  
survdat <- alpha[,c("t", "d", "t0")]


#####################################################
####  Impute CD4 transition times using splines  ####
#####################################################

library(mgcv)

with(dat.spl[[4]], plot(cd4.t, cd4, ylim=c(0, 900)))
with(list(tpred=seq(0, 2.5, 0.01)), lines(tpred, predict(gam(cd4 ~ s(cd4.t, bs="cr", k=3), data=dat.spl[[4]]), list(cd4.t=tpred))))
with(list(tpred=seq(0, 2.5, 0.01)), lines(tpred, predict(gam(cd4 ~ s(cd4.t, bs="cr", fx=TRUE, k=3), data=dat.spl[[4]]), list(cd4.t=tpred)), col=2))
with(list(tpred=seq(0, 2.5, 0.01)), lines(tpred, predict(smooth.spline(dat.spl[[4]]$cd4.t, dat.spl[[4]]$cd4, df=3), tpred)$y, col=3))

gam.fit <- lapply(dat.spl, function(dat) gam(cd4 ~ s(cd4.t, bs="cr", k=3), data=dat))
## gam.fit <- lapply(dat.spl, function(dat) gam(cd4 ~ s(cd4.t, bs="cr", fx=TRUE, k=3), data=dat))

source("impute-transitions.R")

ts.fixed <- t(sapply(gam.fit, impute.cd4trans, cd4.thresh=c(500, 350, 200, 100), rand=FALSE))

s.t0.fixed <- apply(ts.fixed, 1, which.min) # stage at first observation

cd40.fixed <- sapply(gam.fit, predict, list(cd4.t=0)) # CD4 projecting spline to 0
table(cut(cd40.fixed, c(0, 100, 200, 350, 500, Inf), 5:1))

t.dat <- sapply(lapply(dat.spl, "[[", "t"), "[", 1)
d.dat <- sapply(lapply(dat.spl, "[[", "d"), "[", 1)
t0.dat <- sapply(lapply(dat.spl, "[[", "cd4.t"), min)

cd4dat.fixed <- list(ts=ts.fixed, t=t.dat, d=d.dat, t0=t0.dat)


############################
####  Fit Markov model  ####
############################

source("likelihood.R")

pi0.start <- c(0.58, 0.23, 0.16, 0.03)
lambda.start <- 1/c(6.37, 2.86, 3.54)
mu.start <- c(0.004, 0.010, 0.026, 0.435)
theta.start <- c(pi0.start[-1], log(lambda.start), log(mu.start))


####  1: Only CD4 seroconverter data, 4 stages  ####
cd4dat.fixed$ts <- cd4dat.fixed$ts[,-5]
cd4dat.fixed$ts[s.t0.fixed==5, 4] <- ts.fixed[s.t0.fixed==5, 5]

ll.cd4dat(pi0.start, lambda.start, mu.start, cd4dat.fixed)

mod1.fit.bfgs <- optim(theta.start, ll, cd4dat=cd4dat.fixed, method="BFGS", control=list(fnscale=-1, trace=4, maxit=10000), hessian=FALSE)
mod1.fit.bfgs$hessian <- optimHess(mod1.fit.bfgs$par, ll, cd4dat=cd4dat.fixed, control=list(fnscale=-1, trace=4, maxit=10000))

round(cov2cor(solve(-mod1.fit.bfgs$hessian)), 3)

mod1.Lam <- transmat(exp(mod1.fit.bfgs$par[4:6]), exp(mod1.fit.bfgs$par[7:10]))
mod1.pi0 <- c(1-sum(mod1.fit.bfgs$par[1:3]), mod1.fit.bfgs$par[1:3])

integrate(Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod1.pi0 %*% expm(t*mod1.Lam)) - 0.5)^2, c(0, 25)) # median survival 


####  2: Include survival data, 4 stages  ####

system.time(print(ll.survdat(pi0.start, lambda.start, mu.start, survdat)))
system.time(print(ll(theta.start, cd4dat.fixed, survdat)))

mod2.fit.bfgs <- optim(theta.start, ll, cd4dat=cd4dat.fixed, survdat=survdat, method="BFGS", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1),  hessian=FALSE)
mod2.fit.bfgs$hessian <- optimHess(mod2.fit.bfgs$par, ll, cd4dat=cd4dat.fixed, survdat=survdat, control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1))

rbind(mod2.fit.bfgs$par, mod1.fit.bfgs$par)
round(cov2cor(solve(-mod2.fit.bfgs$hessian)), 3)

exp(mod2.fit.bfgs$par[4:6])
exp(tail(mod2.fit.bfgs$par, 4))

mod2.Lam <- transmat(exp(mod2.fit.bfgs$par[4:6]), exp(mod2.fit.bfgs$par[7:10]))
mod2.pi0 <- c(1-sum(mod2.fit.bfgs$par[1:3]), mod2.fit.bfgs$par[1:3])

integrate(Vectorize(function(t) sum(mod2.pi0 %*% expm(t*mod2.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod2.pi0 %*% expm(t*mod2.Lam)) - 0.5)^2, c(0, 25)) # median survival 

library(survival)

print(survfit(Surv(survdat$t, survdat$d)~1))

plot(survfit(Surv(survdat$t, survdat$d)~1), mark.time=FALSE, col="blue", xlab="Years after SC")
lines(survfit(Surv(cd4dat.fixed$t, cd4dat.fixed$d)~1), mark.time=FALSE, col="green")
lines(survfit(Surv(c(survdat$t, cd4dat.fixed$t), c(survdat$d, cd4dat.fixed$d))~1), mark.time=FALSE, col="red")

lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam)))(seq(0, 17, 0.1)), col="darkgreen", lwd=2)
lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod2.pi0 %*% expm(t*mod2.Lam)))(seq(0, 17, 0.1)), col="darkred", lwd=2)

legend("bottomleft", c("CD4 data", "Surv data", "All data", "Mod 1 (CD4 only, 4 stages)", "Mod 2 (CD4+surv, 4 stages)"), lwd=c(1, 1, 1, 2, 2), col=c("green", "blue", "red", "darkgreen", "darkred"))


####  3: add CD4 < 100 stage, CD4 data only  ####

cd4dat.fixed$ts <- ts.fixed # data for all 5 stages

## pi0.start <- c(0.58, 0.23, 0.16, 0.02, 0.01)
## lambda.start <- 1/c(6.37, 2.86, 2.54, 1.00)
## mu.start <- c(0.004, 0.010, 0.026, 0.1, 0.5)
pi0.start <- c(0.6818071, 0.212106518, 0.082567728, 0.014531260, 0.008987408)
lambda.start <- exp(c(-1.547370403, -1.134108219, -1.460319857, -0.824318568))
mu.start <- exp(c(-6.199900226, -5.383513573, -4.238293157, -2.738070149, -0.822751797))
theta.start <- c(pi0.start[-1], log(lambda.start), log(mu.start))

ll.cd4dat(pi0.start, lambda.start, mu.start, cd4dat.fixed)
ll(theta.start, cd4dat.fixed, NS=5)

## mod3.fit.bfgs <- optim(theta.start, ll, cd4dat=cd4dat.fixed, NS=5, method="BFGS", control=list(fnscale=-1, trace=4, maxit=10000), hessian=FALSE) # ERROR: non-finite finite-difference value [4]
mod3.fit.bfgs.not0 <- optim(theta.start, ll, cd4dat=cd4dat.fixed[c("ts", "t", "d")], NS=5, method="BFGS", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE) # fit ignoring left censoring of CD4 data (t0=0) --> will mess up pi0 estimate
mod3.fit.bfgs.not0$hessian <- optimHess(mod3.fit.bfgs.not0$par, ll, cd4dat=cd4dat.fixed[c("ts", "t", "d")], NS=5, control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1))
mod3.fit.nm <- optim(mod3.fit.bfgs.not0$par, ll, cd4dat=cd4dat.fixed, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)
mod3.fit.nm <- optim(mod3.fit.nm$par, ll, cd4dat=cd4dat.fixed, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)
mod3.fit.nm <- optim(mod3.fit.nm$par, ll, cd4dat=cd4dat.fixed, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE) 

## pi0.alt <- c(0.25, 0.25, 0.25, 0.15, 0.10)
## lambda.alt <- 1/c(4, 4, 4, 4)
## mu.alt <- c(0.1, 0.1, 0.1, 0.1, 0.1)o
## theta.alt <- c(pi0.alt[-1], log(lambda.alt), log(mu.alt))
## ll(theta.alt, cd4dat.fixed)

## mod3alt.fit.bfgs <- optim(theta.alt, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, hessian=TRUE, method="BFGS")


####  4: add CD4 < 100 stage, with survival data  ####

system.time(print(ll.survdat(pi0.start, lambda.start, mu.start, survdat)))
system.time(print(ll(theta.start, cd4dat.fixed, survdat, NS=5)))

mod4.fit.bfgs.not0 <- optim(theta.start, ll, cd4dat=cd4dat.fixed[c("ts", "t", "d")], survdat=survdat, NS=5, method="BFGS", control=list(fnscale=-1, trace=4, maxit=10000), hessian=FALSE)
## mod4.fit.bfgs.not0$hessian <- optimHess(mod4.fit.bfgs.not0$par, ll, cd4dat=cd4dat.fixed[c("ts", "t", "d")], survdat=survdat, NS=5, control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1)) ## ERROR: non-finite finite-difference value [4]
mod4.fit.nm <- optim(mod4.fit.bfgs.not0$par, ll, cd4dat=cd4dat.fixed, survdat=survdat, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)
mod4.fit.nm <- optim(mod4.fit.nm$par, ll, cd4dat=cd4dat.fixed, survdat=survdat, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)
mod4.fit.nm <- optim(mod4.fit.nm$par, ll, cd4dat=cd4dat.fixed, survdat=survdat, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)
mod4.fit.nm <- optim(mod4.fit.nm$par, ll, cd4dat=cd4dat.fixed, survdat=survdat, NS=5, method="Nelder-Mead", control=list(fnscale=-1, trace=4, maxit=10000, REPORT=1), hessian=FALSE)


## plot

mod1.pi0 <- c(1-sum(mod1.fit.bfgs$par[1:3]), mod1.fit.bfgs$par[1:3])
mod1.lambda <- exp(mod1.fit.bfgs$par[4:6])
mod1.mu <- exp(mod1.fit.bfgs$par[7:10])
mod1.Lam <- transmat(mod1.lambda, mod1.mu)

mod2.pi0 <- c(1-sum(mod2.fit.bfgs$par[1:3]), mod2.fit.bfgs$par[1:3])
mod2.lambda <- exp(mod2.fit.bfgs$par[4:6])
mod2.mu <- exp(mod2.fit.bfgs$par[7:10])
mod2.Lam <- transmat(mod2.lambda, mod2.mu)

mod3.pi0 <- c(1-sum(mod3.fit.nm$par[1:4]), mod3.fit.nm$par[1:4])
mod3.lambda <- exp(mod3.fit.nm$par[5:8])
mod3.mu <- exp(mod3.fit.nm$par[9:13])
mod3.Lam <- transmat(mod3.lambda, mod3.mu)

mod4.pi0 <- c(1-sum(mod4.fit.nm$par[1:4]), mod4.fit.nm$par[1:4])
mod4.lambda <- exp(mod4.fit.nm$par[5:8])
mod4.mu <- exp(mod4.fit.nm$par[9:13])
mod4.Lam <- transmat(mod4.lambda, mod4.mu)

integrate(Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod1.pi0 %*% expm(t*mod1.Lam)) - 0.5)^2, c(0, 25)) # median survival

integrate(Vectorize(function(t) sum(mod2.pi0 %*% expm(t*mod2.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod2.pi0 %*% expm(t*mod2.Lam)) - 0.5)^2, c(0, 25)) # median survival

integrate(Vectorize(function(t) sum(mod3.pi0 %*% expm(t*mod3.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod3.pi0 %*% expm(t*mod3.Lam)) - 0.5)^2, c(0, 25)) # median survival

integrate(Vectorize(function(t) sum(mod4.pi0 %*% expm(t*mod4.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod4.pi0 %*% expm(t*mod4.Lam)) - 0.5)^2, c(0, 25)) # median survival 

library(survival)

print(survfit(Surv(survdat$t, survdat$d)~1))

plot(survfit(Surv(survdat$t0, survdat$t, survdat$d)~1), mark.time=FALSE, col="blue", xlab="Years after SC")
lines(survfit(Surv(cd4dat.fixed$t0, cd4dat.fixed$t, cd4dat.fixed$d)~1), mark.time=FALSE, col="green")
lines(survfit(Surv(c(survdat$t, cd4dat.fixed$t), c(survdat$d, cd4dat.fixed$d))~1), mark.time=FALSE, col="red")

lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam)))(seq(0, 17, 0.1)), col="darkgreen", lwd=2)
lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod2.pi0 %*% expm(t*mod2.Lam)))(seq(0, 17, 0.1)), col="darkred", lwd=2)
lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod3.pi0 %*% expm(t*mod3.Lam)))(seq(0, 17, 0.1)), col="darkgreen", lwd=2, lty=2)
lines(seq(0, 17, 0.1), Vectorize(function(t) sum(mod4.pi0 %*% expm(t*mod4.Lam)))(seq(0, 17, 0.1)), col="darkred", lwd=2, lty=2)
legend("bottomleft", c("CD4 data", "Surv data", "All data", "Mod 1 (CD4 only, 4 stages)", "Mod 2 (CD4+surv, 4 stages)", "Mod 3 (CD4 only, 5 stages)", "Mod 4 (CD4+surv, 5 stages)"), lwd=c(1, 1, 1, 2, 2, 2, 2), lty=c(1, 1, 1, 1, 1, 2, 2), col=c("green", "blue", "red", "darkgreen", "darkred", "darkgreen", "darkred"))



## summary table of estimated parameters and mean & median survival
rbind(cbind(mod1.pi0, mod2.pi0, mod3.pi0, mod4.pi0),
      cbind(mod1.lambda, mod2.lambda, mod3.lambda, mod4.lambda),
      cbind(mod1.mu, mod2.mu, mod3.mu, mod4.mu),
      c(integrate(Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam))), 0, Inf)$val,
        integrate(Vectorize(function(t) sum(mod2.pi0 %*% expm(t*mod2.Lam))), 0, Inf)$val,
        integrate(Vectorize(function(t) sum(mod3.pi0 %*% expm(t*mod3.Lam))), 0, Inf)$val,
        integrate(Vectorize(function(t) sum(mod4.pi0 %*% expm(t*mod4.Lam))), 0, Inf)$val),
      c(optimize(function(t) (sum(mod1.pi0 %*% expm(t*mod1.Lam)) - 0.5)^2, c(0, 25))$min,
        optimize(function(t) (sum(mod2.pi0 %*% expm(t*mod2.Lam)) - 0.5)^2, c(0, 25))$min,
        optimize(function(t) (sum(mod3.pi0 %*% expm(t*mod3.Lam)) - 0.5)^2, c(0, 25))$min,
        optimize(function(t) (sum(mod4.pi0 %*% expm(t*mod4.Lam)) - 0.5)^2, c(0, 25))$min))

save(mod1.fit.bfgs, mod2.fit.bfgs, mod3.fit.nm, mod4.fit.nm, file="fit-cd4-with-leftcens_2014-07-11.RData")
