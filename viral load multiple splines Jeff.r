rm(list=ls())

require(lattice)
require(gam)
require(mgcv)
require(survival)
require(abind)

setwd("//fi--san02/homes/tmangal/")
sero <- read.csv("HIV\\CD4_data\\SPVL_CD4 analysis\\sero.csv", header=T, sep=",")

sero <- sero[!is.na(sero$SPVL),]
sero$ID<-factor(sero$ID)
length(sero$ID) 
sero$LASTNEG_D <- as.Date(sero$LASTNEG_D, "%d/%m/%Y")
sero$FIRSTPOS_D <- as.Date(sero$FIRSTPOS_D, "%d/%m/%Y")
sero$SERO_INT <- as.numeric(sero$FIRSTPOS_D - sero$LASTNEG_D) 
summary(sero$SERO_INT)  # seroconversion interval, days

# select seroconversion interval < 1year	
cd4 <- sero[sero$SERO_INT <=365,]
cd4 <- cd4[cd4$SERO_INT >=0,]
cd4$ID<-factor(cd4$ID)
cd4.patients <- levels(cd4$ID)
n.cd4.patients<-length(cd4.patients)
print(paste("N = ",n.cd4.patients))

# choose cutoff for # of CD4 counts - 5
cd4.length <- aggregate(cd4$CD4_V ~ cd4$ID, FUN=length) # no. viral counts
names(cd4.length)<- c("ID", "N.cd4")
cd4 <- merge(cd4, cd4.length, by="ID")
cd4 <- cd4[cd4$N.cd4>4,]

cd4$PatSPVL <- paste(round(cd4$SPVL, 3), cd4$ID, sep=" ") # SPVL+ID
cd4$PatSPVL <- factor(cd4$PatSPVL)
cd4$DATE <- as.Date(cd4$DATE, "%d/%m/%Y")
cd4 <- cd4[with(cd4, order(PatSPVL, DATE)),] # sort

cd4$ID <- factor(cd4$ID)
cd4.patients <- levels(cd4$ID)
print(paste("N = ", n.cd4.patients <- length(cd4.patients))) 
length(cd4[,1])
length(levels(cd4$ID))


########### Fitting splines to CD4 counts
cd4.times<-data.frame("ID"=character(), "SPVL"=numeric(),
	"t.cd4.500"=numeric(), "t.cd4.350"=numeric(), "t.cd4.200"=numeric(),
	"t.final"=numeric(), stringsAsFactors=FALSE)
	
cd4_1 <- cd4[!is.na(cd4$CD4_V),] # remove rows for censored obs, for spline fitting	
cd4_1$TimePos <- cd4_1$TIMESERO_YR

# ops<-par(mfrow=c(2,2))
for(i in 1:n.cd4.patients) {
	cd4.temp<-cd4_1[cd4_1$ID==cd4.patients[i],]
	main=cd4.temp$ID[1]
	cd4.spl <- smooth.spline(cd4.temp$TimePos, cd4.temp$CD4_V, df=3) # ASSUMPTION TO BE TESTED
	xx <- unique(sort(c(seq(0, max(cd4.temp$TimePos), by = .2), kn <- unique(cd4.temp$TimePos))))
	# plot(cd4.temp$TimePos, cd4.temp$CD4_V, xlim = range(xx)+0.5, ylim=c(0, 1500),
		# xlab="Years since seroconversion", ylab="CD4 count", cex.lab=0.9, cex.axis=0.9, cex.main=0.9,
		# main=main)
	# lines(pp <- predict(cd4.spl, xx), col = "red")
	#abline(500,0)
	#abline(350,0)
	#abline(200,0)

	#i.kn <- match(kn, xx)   # indices of knots within xx
	#points(kn, pp$y[i.kn], pch = 3, col = "dark red") 
	xx<-xx[xx>min(cd4.temp$TimePos)]
	xx<-xx[xx<max(cd4.temp$TimePos)]
	t.final<-max(cd4.temp$TimePos)
	cd4.range<-c(predict(cd4.spl,0)$y, min(predict(cd4.spl,xx)$y))

	if(cd4.range[1]>500 & cd4.range[2]<500) {
	#	cd4.500<-uniroot(function(x) predict(cd4.spl,x)$y-500,lower=min(xx),upper=max(xx))
	#	t.cd4.500<-cd4.500$root
		cd4.500<-optimize(function(x) (predict(cd4.spl, x)$y-500)^2, lower=min(xx),
			upper=max(xx), tol=0.01)
		t.cd4.500<-cd4.500$minimum  # time leaving CD4 500 compartment
	#	points(t.cd4.500, predict(cd4.spl,t.cd4.500)$y, pch = 5, col = "blue") 
	#	abline(v=t.cd4.500)
	} else if(cd4.range[2]>500) t.cd4.500=-1 else t.cd4.500=0

	# if(t.cd4.500>0) rect(0, 500, t.cd4.500, 1200, col=rgb(0,0,1,alpha=0.2), lty=0)
	# if (cd4.range[1]>500&t.cd4.500<0) rect(0,500,t.final,1200,col=rgb(0,0,1,alpha=0.2),lty=0)

	if(cd4.range[1]>350&cd4.range[2]<350) {
	#	cd4.350<-uniroot(function(x) predict(cd4.spl,x)$y-350,lower=min(xx),upper=max(xx))
	#	t.cd4.350<-cd4.350$root
		cd4.350<-optimize(function(x) (predict(cd4.spl,x)$y-350)^2,lower=min(xx),upper=max(xx),tol=0.01)
		t.cd4.350<-cd4.350$minimum
	#	points(t.cd4.350, predict(cd4.spl,t.cd4.350)$y, pch = 5, col = "green") 
	#	abline(v=t.cd4.350)
	} else if(cd4.range[2]>350) t.cd4.350=-1 else t.cd4.350=0

	# if(t.cd4.350>0) rect(max(0,t.cd4.500),350,t.cd4.350,500,col=rgb(0,1,0,alpha=0.2),lty=0)
	# if (cd4.range[1]<500&cd4.range[1]>350&t.cd4.350<0) rect(max(0,t.cd4.500),350,t.final,500,col=rgb(0,1,0,alpha=0.2),lty=0)
	# if(t.cd4.500==0&t.cd4.350<0) rect(0,350,t.final,500,col=rgb(0,1,0,alpha=0.2),lty=0)
	# if(t.cd4.500>0&t.cd4.350<0) rect(t.cd4.500,350,t.final,500,col=rgb(0,1,0,alpha=0.2),lty=0)

	if(cd4.range[1]>200&cd4.range[2]<200) {
	#	cd4.200<-uniroot(function(x) predict(cd4.spl,x)$y-200,lower=min(xx),upper=max(xx))
	#	t.cd4.200<-cd4.200$root
		cd4.200<-optimize(function(x) (predict(cd4.spl,x)$y-200)^2,lower=min(xx),upper=max(xx),tol=0.01)
		t.cd4.200<-cd4.200$minimum
	#	points(t.cd4.200, predict(cd4.spl,t.cd4.200)$y, pch = 5, col = "green") 
	#	abline(v=t.cd4.200)
	} else if(cd4.range[2]>200) t.cd4.200=-1 else t.cd4.200=0

	# if(t.cd4.200>0) rect(max(0,t.cd4.350),200,t.cd4.200,350,col=rgb(1,0,0,alpha=0.2),lty=0)
	# if (t.cd4.350>0&t.cd4.200<0) rect(t.cd4.350,200,t.final,350,col=rgb(1,0,0,alpha=0.2),lty=0)
	# if (cd4.range[1]<350&t.cd4.200<0) rect(0,200,t.final,350,col=rgb(1,0,0,alpha=0.2),lty=0)

	# if(t.cd4.350==0&t.cd4.200<0&!t.cd4.500==0) rect(0,200,t.final,350,col=rgb(1,0,0,alpha=0.2),lty=0)

	# if((t.final>t.cd4.200&t.cd4.200>0)|t.cd4.200==0) rect(t.cd4.200,0,t.final,200,col=rgb(1,0,1,alpha=0.2),lty=0)
	# if(t.cd4.200==0&!t.cd4.350==0) rect(0,0,t.final,200,col=rgb(1,0,1,alpha=0.2),lty=0)

	cd4.times[i,]<-c(cd4.patients[i],cd4.temp[1,]$SPVL,t.cd4.500,t.cd4.350,t.cd4.200,t.final)
}

## Predicting CD4 compartment compared with observed compartment
cd4.times$ID <- factor(cd4.times$ID)
cd4.times$SPVL <- as.numeric(cd4.times$SPVL)
cd4.times$t.cd4.500 <- as.numeric(cd4.times$t.cd4.500) # time leaving CD4=500+ compartment
cd4.times$t.cd4.350 <- as.numeric(cd4.times$t.cd4.350)
cd4.times$t.cd4.200 <- as.numeric(cd4.times$t.cd4.200)
cd4.times$t.final <- as.numeric(cd4.times$t.final) # time of final obs
head(cd4.times)

cd4.cat <- merge(cd4, cd4.times,by="ID")
head(cd4.cat)
cd4.cat <- cd4.cat[!is.na(cd4.cat$CD4_V),]
cd4.cat <- cd4.cat[ order(cd4.cat$ID, cd4.cat$DATE), ]
cd4.cat$CD4obs <- sapply(cd4.cat$CD4_V, function(x) if(x>=500) 1 else(if (x>=350) 2 
	else (if (x>=200) 3 else 4 )))

# 500+
cd4.cat$CD4temp1 <- as.numeric(cd4.cat$t.cd4.500 > 0 & cd4.cat$TIMESERO_YR < cd4.cat$t.cd4.500)
# 350-500
cd4.cat$CD4temp2 <- as.numeric(cd4.cat$t.cd4.350 > 0 & cd4.cat$TIMESERO_YR >= cd4.cat$t.cd4.500 &
	cd4.cat$TIMESERO_YR < cd4.cat$t.cd4.350)
# 200-350	
cd4.cat$CD4temp3 <- as.numeric(cd4.cat$t.cd4.200 > 0 & cd4.cat$TIMESERO_YR >= cd4.cat$t.cd4.350 &
	cd4.cat$TIMESERO_YR < cd4.cat$t.cd4.200)
# 0-200	
cd4.cat$CD4temp4 <- as.numeric(cd4.cat$t.cd4.200 > 0 & cd4.cat$TIMESERO_YR >= cd4.cat$t.cd4.200)

cd4.cat$CD4Pred <- 1*cd4.cat$CD4temp1 + 2*cd4.cat$CD4temp2 + 3*cd4.cat$CD4temp3 + 
	4*cd4.cat$CD4temp4 # creates predicted value of 1-4
	
table(cd4.cat$CD4obs)
table(cd4.cat$CD4Pred)

cd4.cat2 <- cd4.cat[!cd4.cat$CD4Pred==0,]
table(cd4.cat2$CD4obs)
table(cd4.cat2$CD4Pred)

(contingency<-table(cd4.cat2$CD4obs, cd4.cat2$CD4Pred))
prop.table(contingency, 1)


############### Kaplan-Meier ###########################################################
head(cd4.times)

q.spvl<-quantile(cd4.times$SPVL)
cd4.times$VLQ<-sapply(cd4.times$SPVL, function(x) if(x<q.spvl[2]) return(1) 
	else (if(x<q.spvl[3])  return(2) else ( if (x<q.spvl[4]) return(3) else return(4))))
# <25% =1, <50% =2, <75% =3, else =4
cd4.times$VLQ<-factor(cd4.times$VLQ)
nm<-function(x) round(unname(x),3) 
categories<-c(paste("SPVL<=", nm(q.spvl[2])), paste(nm(q.spvl[2]),"<SPVL<=",nm(q.spvl[3])),
	paste(nm(q.spvl[3]), "<SPVL<=",nm(q.spvl[4])), paste(nm(q.spvl[4]),"<SPVL"))
	
print(data.frame(Group=categories, N=table(cd4.times$VLQ)))
(nq<-table(cd4.times$VLQ)[1])

# function for the Kaplan-Meier analyses:
change.2.1000<-function(x) {return(sapply(x,function(y) if (y==-1) return(1000) 
	else return(y)))}
cd4.times$t.cd4.500<-change.2.1000(cd4.times$t.cd4.500)
cd4.times$t.cd4.350<-change.2.1000(cd4.times$t.cd4.350)
cd4.times$t.cd4.200<-change.2.1000(cd4.times$t.cd4.200)

	
########## KM function for multiple transition estimates ##################################	
head(cd4.times)

options(survfit.rmean = "common") 
options(survfit.print.rmean = "common") 

kaplan2 <- function(incl, t.cd4, t.final, VLQ, plot.lines=T) {
	t.final <- t.final[incl]
	VLQ <- VLQ[incl]
	t.cd4 <- t.cd4[incl]
	o.km <- sapply(t.cd4, function(x) if(x==1000) return(0) else return(1))
	# o.km -> 0=right-censored, i.e. no transition, 1=event occured at time t.cd4
	t.km <- pmin(t.cd4, t.final)
	km.surv <- Surv(t.km, o.km) # KM object
	km.curve <- survfit(Surv(t.km, o.km) ~ VLQ) # create survival curve from KM
	km.fit <- summary(km.curve)$table[,5]
	#print(km.curve, print.rmean=TRUE)
	return(list(curve=km.curve, fit=km.fit))
}	

#  covariance matrix
t_sample <- data.frame(cd4.times$t.cd4.500, cd4.times$t.cd4.350, cd4.times$t.cd4.200)
t_sample[t_sample==1000] <- NA # remove for covariance estimate
t_sample[t_sample==0] <- NA
head(t_sample)

# remove individuals with no transitions between any states 
t_sample <- t_sample[!(is.na(t_sample$cd4.times.t.cd4.500)) | !(is.na(t_sample$cd4.times.t.cd4.350))
	| !(is.na(t_sample$cd4.times.t.cd4.200)),]
dim(t_sample)	
head(t_sample)

(Sigma <- cov(t_sample, use="pairwise.complete.obs"))

N <- dim(cd4.times)[1]
Iter <- 100

# dataset for sampling
t_sample2 <- data.frame(cd4.times$VLQ, cd4.times$t.final, cd4.times$t.cd4.500, 
	cd4.times$t.cd4.350, cd4.times$t.cd4.200)	
names(t_sample2) <- c("VLQ", "t.final", "t.cd4.500", "t.cd4.350", "t.cd4.200")
head(t_sample2)
t_sample2[t_sample2==0] <- NA

arr1 <- array(NA, dim=c(Iter, 5, N))
for (i in 1:N){ # for each set of patients
	coefs <- as.numeric(t_sample2[i, 3:5])
	arr1[1:Iter, 1:5, i] <- c(rep(t_sample2$VLQ[i], Iter), rep(t_sample2$t.final[i], Iter),
		mvrnorm(n=Iter, mu=coefs, Sigma, empirical = TRUE)) # sample times
}

arr1[arr1>100] <- 1000 # reset censor values 
arr1[is.na(arr1)] <- 0
colnames(arr1)=colnames(t_sample2) 

# time from seroco to CD4<500
f1 <- matrix(NA, ncol=4, nrow=Iter)	
colnames(f1) <- c("VLQ1", "VLQ2", "VLQ3", "VLQ4")
for (j in 1:Iter){
	arr_sample <- data.frame(t(arr1[j, ,]))
	km.500 <- kaplan2(!(arr_sample$t.cd4.500==0), arr_sample$t.cd4.500, arr_sample$t.final,
	arr_sample$VLQ)	
	
	f1[j,] <- km.500$fit # extract rmean for each set
}

# time from CD4=500 to CD4<350
f2 <- matrix(NA, ncol=4, nrow=Iter)	
colnames(f2) <- c("VLQ1", "VLQ2", "VLQ3", "VLQ4")
for (j in 1:Iter){
	arr_sample <- data.frame(t(arr1[j, ,]))
	
	km.500.to.350 <- kaplan2((!(arr_sample$t.cd4.500==1000)) & (!(arr_sample$t.cd4.350==0)),
	arr_sample$t.cd4.350 - arr_sample$t.cd4.500, arr_sample$t.final - arr_sample$t.cd4.500,
	cd4.times$VLQ)
	
	f2[j,] <- km.500.to.350$fit 
}

# time from CD4=350 to CD4<200
f3 <- matrix(NA, ncol=4, nrow=Iter)	
colnames(f3) <- c("VLQ1", "VLQ2", "VLQ3", "VLQ4")
for (j in 1:Iter){
	arr_sample <- data.frame(t(arr1[j, ,]))
	
	km.350.to.200<-kaplan2((!(arr_sample$t.cd4.350==1000)) & (!(arr_sample$t.cd4.200==0)),
	arr_sample$t.cd4.200 - arr_sample$t.cd4.350, arr_sample$t.final - arr_sample$t.cd4.350,
	arr_sample$VLQ)

	f3[j,] <- km.350.to.200$fit 
}

# extract median + 95% CI for sampled sets
f1MED <- apply(f1[,c("VLQ1","VLQ2", "VLQ3", "VLQ4")], 2, mean)
f2MED <- apply(f2[,c("VLQ1","VLQ2", "VLQ3", "VLQ4")], 2, mean)
f3MED <- apply(f3[,c("VLQ1","VLQ2", "VLQ3", "VLQ4")], 2, mean)
















