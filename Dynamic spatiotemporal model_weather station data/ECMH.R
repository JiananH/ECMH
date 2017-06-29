library(ECMH)

restraint <- read.csv("~/Dropbox/MCMH/restraint on each components/restraint_expand.csv")
#restraint <- read.csv("~/weatherstation/restraint.csv")
restraint <- restraint[,-1]
rest <- as.numeric(restraint[1,])
probability=0.5

###################################################################################################
## Not run:
data("NETemp.dat")
ne.temp <- NETemp.dat

#set.seed(1)

##take a chunk of New England
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3000000,]

##subset first 2 years (Jan 2000 - Dec. 2001)
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months


##scale to km
coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))

##set starting and priors
p <- 2 #number of regression parameters in each month

#starting <- list("beta"=rep(runif(1),N.t*p), "phi"=rep(runif(1,2.5,3.5)/(0.5*max.d), N.t),
#                   "sigma.sq"=rep(runif(1,1.5,2.5),N.t), "tau.sq"=rep(runif(1,0.5,1.5), N.t),
#                   "sigma.eta"=diag(rep(runif(1,0.005,0.015), p)))
starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(5, N.t))

priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))

##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)

n.samples <- 2

m.1 <- spDynLMmodgone(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
                  starting=starting, tuning=tuning, priors=priors, get.fitted =FALSE,
                  cov.model="exponential", n.samples=n.samples, n.report=0.05*n.samples,
		  radiusbeta0=c(rest[1],rest[2]), radiusbeta=c(rest[3],rest[4]), radiustausq = rest[6],
		  radiussigmasq = rest[5], radiusphi = rest[7], radiussigmaEta = 171.67*2, prob=probability)

# m.2 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
#                starting=starting, tuning=tuning, priors=priors, get.fitted =FALSE,
#                cov.model="exponential", n.samples=n.samples, n.report=0.1*n.samples)

esjd <- function(parameters,itr){
  sum(diff(parameters)^2)/(itr-1)
}

esjd.beta0.m1=esjd(m.1$p.beta.0.samples,n.samples)
# esjd.beta0.m2=esjd(m.2$p.beta.0.samples,n.samples)
esjd.beta.m1=esjd(m.1$p.beta.samples,n.samples)
# esjd.beta.m2=esjd(m.2$p.beta.samples,n.samples)
esjd.theta.m1=esjd(m.1$p.theta.samples,n.samples)
# esjd.theta.m2=esjd(m.2$p.theta.samples,n.samples)
esjd.sigmaeta.m1=esjd(m.1$p.sigma.eta.samples,n.samples)
# esjd.sigmaeta.m2=esjd(m.2$p.sigma.eta.samples,n.samples)

m.1.samples=cbind(m.1$p.beta.0.samples,n.samples,m.1$p.beta.samples,m.1$p.theta.samples,m.1$p.sigma.eta.samples)
# m.2.samples=cbind(m.2$p.beta.0.samples,n.samples,m.2$p.beta.samples,m.2$p.theta.samples,m.2$p.sigma.eta.samples)

esjdcombine.m1=esjd(m.1.samples,n.samples)
# esjdcombine.m2=esjd(m.2.samples,n.samples)


#comparing results
burn.in <- 0.75*n.samples

#compute parameter estimates from the long run

#beta.0.samples
beta0.mean.m1=apply(m.1$p.beta.0.samples[burn.in:n.samples,],2,mean)
# beta0.mean.m2=apply(m.2$p.beta.0.samples[burn.in:n.samples,],2,mean)

#beta.samples
beta.mean.m1=apply(m.1$p.beta.samples[burn.in:n.samples,],2,mean)
# beta.mean.m2=apply(m.2$p.beta.samples[burn.in:n.samples,],2,mean)

#theta.samples
theta.mean.m1=apply(m.1$p.theta.samples[burn.in:n.samples,],2,mean)
# theta.mean.m2=apply(m.2$p.theta.samples[burn.in:n.samples,],2,mean)

#sigma.eta.samples
sigma.eta.mean.m1=apply(m.1$p.sigma.eta.samples[burn.in:n.samples,],2,mean)
# sigma.eta.mean.m2=apply(m.2$p.sigma.eta.samples[burn.in:n.samples,],2,mean)

m1.est=c(beta0.mean.m1,beta.mean.m1,theta.mean.m1,sigma.eta.mean.m1)
# m2.est=c(beta0.mean.m2,beta.mean.m2,theta.mean.m2,sigma.eta.mean.m2)

#trueest <- read.csv("/Users/Jianan/Dropbox/MCMH/long run results/trueest100chains.csv")
#trueest=as.numeric(trueest$x)
trueest <- read.csv("~/weatherstation/trueest.csv")
trueest=as.numeric(trueest$x)

# mse.m1=sqrt(sum((m1.est-trueest)^2))
# mse.m2=sqrt(sum((m2.est-trueest)^2))
beta0.mse.m1=mean((beta0.mean.m1-trueest[1:2])^2)
# beta0.mse.m2=mean((beta0.mean.m2-trueest[1:2])^2)

beta.mse.m1=mean((beta.mean.m1-trueest[3:26])^2)
# beta.mse.m2=mean((beta.mean.m2-trueest[3:50])^2)

theta.mse.m1=mean((theta.mean.m1-trueest[27:62])^2)
# theta.mse.m2=mean((theta.mean.m2-trueest[51:122])^2)

sigma.eta.mse.m1=mean((sigma.eta.mean.m1-trueest[63:66])^2)
# sigma.eta.mse.m2=mean((sigma.eta.mean.m2-trueest[123:126])^2)

mse.m1=mean((m1.est-trueest)^2)
# mse.m2=mean((m2.est-trueest)^2)
# mse=mse.m1/mse.m2

res=c(esjd.beta0.m1,esjd.beta.m1,esjd.theta.m1,esjd.sigmaeta.m1,esjdcombine.m1,beta0.mse.m1,beta.mse.m1,theta.mse.m1,sigma.eta.mse.m1,mse.m1)
names(res)=c('esjd.beta0.m1','esjd.beta.m1','esjd.theta.m1','esjd.sigmaeta.m1','esjdcombine.m1','beta0.mse.m1','beta.mse.m1','theta.mse.m1','sigma.eta.mse.m1','mse.m1')

write.csv(res,paste("5percent",jianan,".csv",sep=""))
