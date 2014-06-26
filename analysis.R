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

## seroconversion window
bas$seroco.wind <- bas$firstpos - bas$lastneg

with(bas, table(serohow, seroco.wind < 0, useNA="always"))
with(bas, table(floor(seroco.wind/365.25), serohow, useNA="always"))

## limit seroconversion window to < 2 years
bas <- subset(bas, !(serohow == "midpoint" & (seroco.wind < 0 | seroco.wind/365.25 > 2))); nrow(bas)

bas$t <- as.numeric((bas$cens_d - bas$serodate)/365.25)
bas$d <- bas$cens_r == "death"

## exclude missing CENS_D (for now !!! CHECK THIS LATER)
bas <- subset(bas, !is.na(t)); nrow(bas)

dat <- merge(bas[,c("patient", "serodate", "t", "d")], cd4); nrow(dat)
dat$cd4.t <- as.numeric((dat$cd4date - dat$serodate)/365.25)

table(floor(dat$cd4.t))  # number of CD4s by years post SC
table(dat$cd4.t > dat$t) # CD4s after censor date
dat <- subset(dat, cd4.t >=0 & cd4.t <= t); nrow(dat)


dat.spl <- split(dat, dat$patient)

table(sapply(dat.spl, nrow)) # number of CD4 counts per patient

## limit to minimum 5 CD4 measurements
dat.spl <- dat.spl[sapply(dat.spl, nrow) >= 5]; length(dat.spl)

table(sapply(dat.spl, nrow))

sum(sapply(lapply(dat.spl, "[[", "d"), any)) # number of deaths observed



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

alpha$t <- as.numeric((alpha$exit_date - alpha$seroconv_date)/365.25)
alpha$d <- alpha$exit_type == "death"

survdat <- alpha[,c("t", "d")]


###########################
####  Test spline fit  ####
###########################

library(mgcv)

with(dat.spl[[4]], plot(cd4.t, cd4, ylim=c(0, 900)))
with(list(tpred=seq(0, 2.5, 0.01)), lines(tpred, predict(gam(cd4 ~ s(cd4.t, bs="cr", k=3), data=dat.spl[[4]]), list(cd4.t=tpred))))


gam.fit <- lapply(dat.spl, function(dat) gam(cd4 ~ s(cd4.t, bs="cr", k=3), data=dat))

source("impute-transitions.R")

ts.fixed <- t(sapply(gam.fit, sample.cd4trans, rand=FALSE))

s0.fixed <- apply(ts.fixed == 0, 1, which)

t.dat <- sapply(lapply(dat.spl, "[[", "t"), "[", 1)
d.dat <- sapply(lapply(dat.spl, "[[", "d"), "[", 1)

cd4dat.fixed <- list(ts=ts.fixed, t=t.dat, d=d.dat)


###########################
####  Test likelihood  ####
###########################

source("likelihood.R")

pi0.start <- c(0.58, 0.23, 0.16, 0.03)
lambda.start <- 1/c(6.37, 2.86, 3.54)
mu.start <- c(0.004, 0.010, 0.026, 0.435)
theta.start <- c(pi0.start[-1], log(lambda.start), log(mu.start))


####  1: Only CD4 seroconverter data  ####

ll.cd4dat(pi0.start, lambda.start, mu.start, cd4dat.fixed)

mod1.fit <- optim(theta.start, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, hessian=TRUE)
mod1.fit.bfgs <- optim(theta.start, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, hessian=TRUE, method="BFGS")

mod1.fit.bfgs2 <- optim(test.fit$par, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, hessian=TRUE, method="BFGS")

round(cov2cor(solve(-test.fit.bfgs2$hessian)), 3)

mod1.Lam <- transmat(exp(mod1.fit.bfgs$par[4:6]), exp(mod1.fit.bfgs$par[7:10]))
mod1.pi0 <- c(1-sum(mod1.fit.bfgs$par[1:3]), mod1.fit.bfgs$par[1:3])

integrate(Vectorize(function(t) sum(mod1.pi0 %*% expm(t*mod1.Lam))), 0, Inf) # mean survival
optimize(function(t) (sum(mod1.pi0 %*% expm(t*mod1.Lam)) - 0.5)^2, c(0, 25)) # median survival 




####  2: include survival data  ####

system.time(print(ll.survdat(pi0.start, lambda.start, mu.start, survdat)))
system.time(print(ll(theta.start, cd4dat.fixed, survdat)))

mod2.fit <- optim(theta.start, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, survdat=survdat, hessian=TRUE)
mod2.fit.bfgs <- optim(theta.start, ll, control=list(fnscale=-1, trace=4, maxit=10000), cd4dat=cd4dat.fixed, survdat=survdat, hessian=TRUE, method="BFGS")

rbind(mod2.fit.bfgs$par, test.fit.bfgs$par)
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

legend("bottomleft", c("CD4 data", "Surv data", "All data", "Mod 1 (CD4 only)", "Mod 2 (CD4+surv)"), lwd=c(1, 1, 1, 2, 2), col=c("green", "blue", "red", "darkgreen", "darkred"))


####  3: exponential mortality model  ####
