# Program to compute the likelihood ratio test of the poisson fields ex.
# Data are the number of trees observed inside a quadrat.  

y1.vec <- c(2,7,3,2,4,6)
n1 <- length(y1.vec)

y2.vec <- c(1,1,1,2,2,0,1,1)
n2 <- length(y2.vec)

y3.vec <- c(2,4,3,0,2,4,1)
n3 <- length(y3.vec)

all.quadrats <- c(y1.vec,y2.vec,y3.vec)
n <- n1+n2+n3

mu1.hat <- mean(y1.vec)
mu2.hat <- mean(y2.vec)
mu3.hat <- mean(y3.vec)
mu0.hat <- mean(all.quadrats) 

ln.Lo.hat <- sum(dpois(x=all.quadrats,lambda=mu0.hat,log=TRUE))

ln.L1.hat <- sum(dpois(x=y1.vec,lambda=mu1.hat,log=TRUE)) +

					 sum(dpois(x=y2.vec,lambda=mu2.hat,log=TRUE)) + 
					 
					 sum(dpois(x=y3.vec,lambda=mu3.hat,log=TRUE))

Gsq.obs <- -2*(ln.Lo.hat-ln.L1.hat)

dfs <- 3-1 # 3 parameters estimated under H1 - 1 parameter estimated under Ho
alpha <- 0.05
Gsq.crit <- qchisq(p=1-alpha, df=dfs)
pvalue <- 1-pchisq(q=Gsq.obs,df=dfs)

Chisq.samp <- rchisq(n=100000, df=dfs)

hist(Chisq.samp)
abline(v=Gsq.obs, col="blue")
abline(v=Gsq.crit, col="red")








