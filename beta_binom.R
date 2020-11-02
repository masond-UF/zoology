
####### Today's code:


p<- 0.8
N <- 400
Y <- rbinom(n=1, size=400, prob=p)
Y
Y/N
p.values <- seq(0.001,0.999, by=0.001)
prob.obs.data <- dbinom(x=Y,size=400,prob=p.values)
plot(p.values, prob.obs.data, type="l", lwd=2, col="red")
plot(p.values, prob.obs.data, type="l", lwd=2, col="red", ylab="Likelihood of p", xlim=c(0.6,1))
prand <- rbeta(n=10000, shape1=2,shape2=5)
hist(prand)
mean(prand)
2/(2+5)
abline(v=2/(2+5))
plot(p.values, prob.obs.data/max(prob.obs.data), type="l", lwd=2, col="red", ylab="Likelihood of p", xlim=c(0.6,1))
plot(p.values, prob.obs.data/max(prob.obs.data), type="l", lwd=2, col="red", ylab="Relative Likelihood of p", xlim=c(0.6,1))
abline(h=0.1469)
qchisq(p=0.95,df=1)
exp(-qchisq(p=0.95,df=1)/2)
prands <- rbeta(n=N, shape1=2,shape2=5)
quartz;hist(prands)
newYs <- rbinom(n=N, size=1, prob=prands)
newYs
nsurviv <- sum(newYs)
nsurviv
#assume a binomial probability model
p.mle <- nsurviv/N
p.mle
2/(2+5)
prands
library("extraDistr")
library("MASS")
n <- 100
n
a <- 2
b <- 5
ntrials <- 20
n <- 100
a <- 2
b <- 5
ntrials <- 20
P <- rbeta(n=n, shape1=a, shape2=b)
P
n <- 100 # 100 turtles
a <- 2
b <- 5
ntrials <- 1
P <- rbeta(n=n, shape1=a, shape2=b)
P
YgP <- rbinom(n=n, size=ntrials, prob=P)
YgP
bbnegllike <- function(guess, Yvec, ntrials){
pos.guess <- exp(guess)
a <- pos.guess[1]
b <- pos.guess[2]
bbnegll.vec <- dbbinom(x=Yvec, size=ntrials, alpha=a, beta=b, log=TRUE)
return(-sum(bbnegll.vec))
}
bbnegllike <- function(guess, Yvec, ntrials){
pos.guess <- exp(guess)
a <- pos.guess[1]
b <- pos.guess[2]
bbnegll.vec <- dbbinom(x=Yvec, size=ntrials, alpha=a, beta=b, log=TRUE)
return(-sum(bbnegll.vec))
}
myguess <- log(c(1.5,4))
opt.samp <- optim(par=myguess, fn=bbnegllike, method="Nelder-Mead", Yvec=YgP, ntrials=ntrials)
a.hat <- exp(opt.samp$par[1])
b.hat <- exp(opt.samp$par[2])
print(a.hat)
print(b.hat)
pvals <- seq(.001,0.999, by=0.001)
surv.prob.dist.hat <- dbeta(x=pvals, shape1=a.hat, shape2=b.hat)
quartz();
plot(pvals,surv.prob.dist.hat, type="l", lwd=2, col="red")
abline(v=a.hat/(a.hat+b.hat))
