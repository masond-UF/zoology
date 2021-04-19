# Number of observations is 21 for each trap
# The number of trials is the number of seeds arriving, which is dependent
# upon the MLE of lambda and the proportion of seeds establishing. The proportion
# of seeds establishing is from the establishment trials. 

#find the probability of 7 successes during 20 trials where the probability of
#success on each trial is 0.5

dbinom(x=7, size=20, prob=.5)

# number of successes is establishing individuals?
# number of trials is arriving and that stays the same (MLE)
# 


# Simulating establishment with different seed arrival values ####
# Percent changes from the seed arrival MLE for burn
change.inc <- c(0.25, 0.50, 0.75)
change.dec <- c(-0.75,-0.50,-0.25)

burn.samp.inc <- c()
burn.samp.dec <- c()

for(i in 1:3){
	burn.samp.inc[i] <- abs((burn.lamb.mle*change.inc[i])+burn.lamb.mle)
	burn.samp.dec[i] <- abs((burn.lamb.mle*change.dec[i])+burn.lamb.mle)
}

burn.samp <- c(burn.samp.dec,burn.lamb.mle,burn.samp.inc)

## Percent changes from the seed arrival MLE for control

con.samp.inc <- c()
con.samp.dec <- c()

for(i in 1:3){
	con.samp.inc[i] <- abs((control.lamb.mle*change.inc[i])+control.lamb.mle)
	con.samp.dec[i] <- abs((control.lamb.mle*change.dec[i])+control.lamb.mle)
}

con.samp <- c(con.samp.dec,control.lamb.mle,con.samp.inc)

### For loop estimating plant establishment with varying 
### seed arrival estimates (trials)

boot <- 5000
both <- c(burn.est,cont.est)
Gsq.list <- list()
phat.Ho <- sum(both)/(N*length(both)); # MLE under H0 = sum(xis)/sum(nis)


for (i in 1:7){
	# Estimate the number of establishing plants
	burn.est <- rbinom(n=burn.reps, size=round(burn.samp[i]), prob = burn.p.hat)
	cont.est <- rbinom(n=control.reps, size=round(con.samp[i]), prob = cont.p.hat)
	
	# Calculate likelihood
	lnLo.hat <- lnLo(burn=burn.est, control=cont.est, Ntrials=N)
	lnL1.hat <- lnL1(burn=burn.est, control=cont.est, Ntrials=N)

	Gsq <- -2*(lnLo.hat-lnL1.hat)
	alpha <- 0.001
	Gsq.crit <- qchisq(p=1-alpha, df=3-1)
	pvalue <- 1-pchisq(q=Gsq, df=3-1)
	
	Gsq.list <- list()
	
		for(j in 1:boot){
		
		# Simulate data like the one observed but under the Null hypothesis
		burn.boot <- rbinom(n=burn, size=burn.samp[i], prob=phat.Ho)
		cont.boot <- rbinom(n=control, size=con.samp[i], prob=phat.Ho)
	
		lnLo.boot <- lnLo(burn=burn.boot, control=cont.boot, Ntrials=N)	
		lnL1.boot <- lnL1(burn=burn.boot, control=cont.boot, Ntrials=N)		
	
		Gsq.vec[j] <- -2*(lnLo.boot-lnL1.boot)
		}
	boot.pval <- sum(Gsq.vec > Gsq)/2000 # Proportion of boostrap Gsq's that are bigger than the observed Gsq
	print(boot.pval) # compare to Chisquare pvalue
	print(pvalue)

	hist(Gsq.vec)
	abline(v=Gsq, col="red")
}


	