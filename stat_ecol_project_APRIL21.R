##### David Mason ##### 14 April 2021 ##### STAT ECOL PROJ ###
# Set-up #####
library(tidyverse)
seeds <- read.csv("seed_trap.csv", header = T)
seeds_long <- pivot_longer(seeds, cols = 4:9,names_to = "SPECIES", 
													 values_to = "SEEDS")

burn <- seeds_long %>% 
filter(TREATMENT == "BURN")
control <- seeds_long %>% 
filter(TREATMENT == "CONTROL")
# Calculate seed arrival MLE ####
treatment <- function(burn,control ,plot.it=FALSE){

	# Computing Likelihood for Observed Data ####
	llh_poisson <- function(lambda, y){
		# log(likelihood) by summing 
		llh <- sum(dpois(y, lambda, log=TRUE))
  		return(llh)
	}

	lambdas <- seq(0,100, by=0.01)

	# compute log-likelihood for all lambda values
	burn_ll <- sapply(lambdas,function(x){llh_poisson(x,burn$SEEDS)})
	max.burn.ll <- max(burn_ll)
	burn.mle <- lambdas[burn_ll==max.burn.ll]
	burn.rel.l <- exp(burn_ll - max.burn.ll)
	burn.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	burn.cut.m.rell <- burn.ci.cuttoff - burn.rel.l
	burn.ci.indices <- which(burn.cut.m.rell<0, arr.ind=TRUE)
	burn.len.ciind <- length(burn.ci.indices)
	burn.left.CL <- lambdas[burn.ci.indices[1]]
	burn.right.CL <-lambdas[burn.ci.indices[burn.len.ciind]]
	n.burn <- length(burn$SEEDS)
	
	control_ll <- sapply(lambdas,function(x){llh_poisson(x,control$SEEDS)})
	max.control.ll <- max(control_ll)
	control.mle <- lambdas[control_ll==max.control.ll]
	control.rel.l <- exp(control_ll - max.control.ll)
	control.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	control.cut.m.rell <- control.ci.cuttoff - control.rel.l
	control.ci.indices <- which(control.cut.m.rell<0, arr.ind=TRUE)
	control.len.ciind <- length(control.ci.indices)
	control.left.CL <- lambdas[control.ci.indices[1]]
	control.right.CL <-lambdas[control.ci.indices[control.len.ciind]]
	n.control <- length(control$SEEDS)	

	#### Null model log-likelihood: no difference between burn and control?
	all.data  <- c(as.vector(burn$SEEDS), as.vector(control$SEEDS))
	null_ll <- sapply(lambdas,function(x){llh_poisson(x,all.data)})
	max.null.ll <- max(null_ll)
	null.mle <- lambdas[null_ll==max.null.ll]

	bic.sep.model <- -2*(max.burn.ll+max.control.ll) + 2*log(n.burn+n.control)
	bic.null.model <- -2*max.null.ll + 1*log(n.burn+n.control)
	if(bic.sep.model<bic.null.model){outcome <- 
		"Separate model is best"}else{outcome<- 
			"There's no difference between control and burn"}
	# Save the lambdas and log-likelihoods in a data frame
	burn_df <- data.frame(ll=burn_ll, lambda=lambdas)
	burn_df$TREATMENT <- c("BURN")
	control_df <- data.frame(ll=control_ll, lambda=lambdas)
	control_df$TREATMENT <- c("CONTROL")
	profile.ll <- rbind(burn_df, control_df)
	
	# Save the MLES and the confidence intervals	
	burn.cis <- c(burn.left.CL, burn.mle, burn.right.CL)
	ctl.cis  <- c(control.left.CL, control.mle, control.right.CL)
	out <- rbind(burn.cis, ctl.cis)	
	colnames(out) <- c("LCL", "MLE", "UPL")
	rownames(out) <- c("Burn treat.", "Control")
	
	
	
	if(plot.it==TRUE){
		
		burn.ci <- seq(from=burn.left.CL, to=burn.right.CL, by=0.01)
		ctl.ci <- seq(from=control.left.CL, to=control.right.CL, by=0.01)
		ci.height <- burn.ci.cuttoff
		
		xlims.burn <- c(burn.left.CL-1, burn.right.CL+1)
		xlims.ctl  <- c(control.left.CL-1, control.right.CL+1) 
		par(mfrow=c(1,2), oma=c(1,2,1,1), mar=c(4,5,4,2))		
		plot(lambdas, burn.rel.l, type="l", col="red", lwd=2, xlim=xlims.burn, xlab=expression(lambda), 
		ylab="Relative Profile Likelihood", cex.lab=1.5, bty="l", main="Burn treatment")
		polygon(x=c(burn.ci,rev(burn.ci)), y=c(rep(0,length(burn.ci)),rep(ci.height,length(burn.ci))), col="grey")
		plot(lambdas, control.rel.l, type="l", col="red", lwd=2, xlim=xlims.ctl,xlab=expression(lambda), 
		ylab="Relative Profile Likelihood", cex.lab=1.5, bty="l", main="Control")
		polygon(x=c(ctl.ci,rev(ctl.ci)), y=c(rep(0,length(ctl.ci)),rep(ci.height,length(ctl.ci))), col="grey")		
		
	}
	
	
	return(list(CIS.lambdas=out,bic.sep.model=bic.sep.model, bic.null.model=bic.null.model, outcome=outcome))
}

output <- treatment(burn,control,plot.it = TRUE)
out.mle <- output[[1]]

burn.lamb.mle <- out.mle[1,3]
control.lamb.mle <- out.mle[2,3]
# Calculate proportion of establishment MLE ####
# Number of trials 
N <- 10
burn.reps <- 21
control.reps <- 21

# Simulate the data
burn.samp <- rbinom(n=burn.reps, size=N, prob=0.6)
cont.samp <- rbinom(n=burn.reps, size=N, prob=0.3)

# Calculate phat
burn.p.hat <- sum(burn.samp)/(N*length(burn.samp))
cont.p.hat <- sum(cont.samp)/(N*length(cont.samp))

# Estimate seed establishment for treatments ####
# Number of trials (arriving seeds) is from lambda
burn.est <- rbinom(n=burn.reps, size=round(burn.lamb.mle), prob = burn.p.hat)
cont.est <- rbinom(n=control.reps, size=round(control.lamb.mle), prob = cont.p.hat)

# Comparing total plant establishmet by treatment ####

# Likelihood function under the null
lnLo <- function(burn,control,burn.trials,con.trials){
	
	both.counts <- c(burn,control)
	Ntrials <- round((burn.trials+con.trials)/2)
	p.hat <- sum(c(burn,control))/(Ntrials*length(both.counts))
	llikevec <- dbinom(x=both.counts, size=Ntrials, prob=p.hat, log=TRUE)
	return(sum(llikevec))
}

# Likelihood function under the alternative
lnL1 <- function(burn,control,burn.trials,con.trials){
	
	phat.burn <- sum(burn)/(burn.trials*length(burn))
	phat.cont <- sum(control)/(con.trials*length(control))
	
	llike.burn <- sum(dbinom(x=burn,size=round(burn.trials), prob=phat.burn, log=TRUE))
	llike.cont <- sum(dbinom(x=cont.est,size=round(control.lamb.mle), prob=phat.cont, log=TRUE))
	
	return(sum(c(llike.burn,llike.cont)))
}

# Computing Gsq = -2log(Lo/L1)
lnLo.hat <- lnLo(burn=burn.est, control=cont.est, 
								 burn.trials = burn.lamb.mle, con.trials = control.lamb.mle)
lnL1.hat <- lnL1(burn=burn.est, control=cont.est, 
								 burn.trials = burn.lamb.mle, con.trials = control.lamb.mle)

Gsq <- -2*(lnLo.hat-lnL1.hat)
alpha <- 0.001
Gsq.crit <- qchisq(p=1-alpha, df=3-1)
pvalue <- 1-pchisq(q=Gsq, df=3-1)

# Parametric bootstrapping ####

boot <- 5000
both <- c(burn.est,cont.est)
phat.Ho <- sum(both)/(N*length(both)); # MLE under H0 = sum(xis)/sum(nis)
Gsq.vec <- rep(0,boot)

	Ntrials <- round((burn.lamb.mle+control.lamb.mle)/2)

for(i in 1:boot){
	# Simulate data like the one observed but under the Null hypothesis
	burn.boot <- rbinom(n=burn.reps, size=round(burn.lamb.mle), prob=phat.Ho)
	cont.boot <- rbinom(n=control.reps, size=round(control.lamb.mle), prob=phat.Ho)
	
	lnLo.boot <- lnLo(burn=burn.boot, control=cont.boot,
										 burn.trials = burn.lamb.mle, con.trials = control.lamb.mle)	
	lnL1.boot <- lnL1(burn=burn.boot, control=cont.boot,
										 burn.trials = burn.lamb.mle, con.trials = control.lamb.mle)		
	
	Gsq.vec[i] <- -2*(lnLo.boot-lnL1.boot)
	
}

boot.pval <- sum(Gsq.vec > Gsq)/2000 # Proportion of boostrap Gsq's that are bigger than the observed Gsq
print(boot.pval) # compare to Chisquare pvalue
print(pvalue)

hist(Gsq.vec)
abline(v=Gsq, col="red")
