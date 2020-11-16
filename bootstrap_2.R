# 11 November 2020—Bootstra
# Set up the workspace ####
library(tidyverse)

# Bring in the data
seeds <- read.csv("seed_trap.csv", header = T)
# Subet the data
seeds_long <- pivot_longer(seeds, cols = 4:9,names_to = "SPECIES", 
													 values_to = "SEEDS")

# Subset the data by  treatment

burn <- seeds_long %>% 
filter(TREATMENT == "BURN", SPECIES == SPECIES)
control <- seeds_long %>% 
filter(TREATMENT == "CONTROL", SPECIES == SPECIES)

# Create the function for comparing seed arrival in burn and control traps ####
treatment <- function(burn, control ,plot.it=TRUE){

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
		"1"}else{outcome<- 
			"0"}
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
# Write a loop for running the function 1000 times and computing % success ####

# Choose the difference in mean seed arrival
pct.dec  <- 0.05 
pct.dec  <- 0.15 
pct.dec  <- 0.30 
pct.dec  <- 0.50 
pct.dec  <- 0.90 

# Create vectors for lambda values
lam.ctl.vec  <- c(5,10,50,250,500)
lam.burn.vec <- (1-pct.dec)*lam.ctl 

# Select one to run with the for loop
lam.ctl <- lam.ctl[3]
lam.burn <- lam.burn[3]

# Create a vector for the number of seed traps
samp.sizes <- seq(1:50, by = 5)

# Create empty objects for the for loop to fill
vec <- vector() # empty vector for ith outcomes
output.mat <- matrix(0, nrow = 50, ncol = 20) # empty matrix for the final output


for(i in 1:10){
	
	n.samples <- samp.sizes[i]
	
		for(j in 1:20){
				
				for(z in 1:1000){
					# generate samples from poisson distribution
					burn.numbers <- rpois(n = n.samples, lam.burn)
					control.numbers <- rpois(n = n.samples, lam.ctl)
			
					# turn it into a dataframe the function recognizes
					burn <- matrix(burn.numbers, ncol = 1)
					colnames(burn) <- c("SEEDS")
					burn <- as.data.frame(burn)
	
					control <- matrix(control.numbers, ncol = 1)
					colnames(control) <- c("SEEDS")
					control <- as.data.frame(control)
			
					# run the function
					ith.test <- treatment(burn=burn,control=control ,plot.it=FALSE)
					
					# calculate the proportion of successes
					vec[z] <- as.integer(ith.test$outcome)
					
					# put the proportion of successes in each cell 
					output.mat[i,j] <- sum(vec)/length(vec)
					}
		}
}


