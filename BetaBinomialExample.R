bbnegllike <- function(guess, Yvec, ntrials = 20){
	
	pos.guess <- exp(guess)
	a <- pos.guess[1]
	b <- pos.guess[2]
	
	bbnegll.vec <- dbbinom(x= Yvec, size = ntrials, 
												 alpha = a, beta = b, log = TRUE)
	
}

# function for likelihood of observing the data (beta binomial distribution
# that averages out the randomness in the values of p)

n <- 100
a <- 2
b <- 5
ntrials <- 20 # 20 turtles checking 100 times how many survive
 
P <- rbeta(n=n, shape1 = a, shape2 = b) # for every set of 20 turtles we draw a survival probability
YgP <- rbinom(n = n, size = ntrials, prob = P) # 


# optimization routine which tries different values of parameters of interest
# to maximize probability or minimize negative log
myguess <- log(c(1.5,4))
opt.samp <- optim(par = myguess, fn=bbnegllike, method = "Nelder-Mead", 
									Yvec = YgP, ntrials = ntrials)
a.hat <- exp(opt.samp$par[1])
b.hat <- exp(opt, samp$par[2])
print(a.hat)
print(b.hat)
# we were just estimating survival probability, but now we have a + b. This
# allows us to estimate something else.

pvals <- seq(0.001, 0.999, by = 0.001)
surv.prob.dist.hat <- dbeta(x = pvals, shape1 = a.hat, shape2 = b.hat)
quartz();
plot(pvals, surviv.prob.dist.hat, type = "l", lwd = 2, col = "red")
abline(v=a.hat/(a.hat+b.hat)) # mean with

