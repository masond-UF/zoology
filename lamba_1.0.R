# 23 October 2020—Simulating lambda—David Mason ####
library(tidyverse)
seeds <- read.csv("seed_trap.csv", header = T)
seeds_long <- pivot_longer(seeds, cols = 4:9,names_to = "SPECIES", 
													 values_to = "SEEDS")
species <- c("RUBUS", "PRSE2", "PIPA2", "SMILAX", "PHAM4", "PIEC2")

# Subset the data by species & treatment

burn <- seeds_long %>% 
filter(TREATMENT == "BURN", SPECIES == SPECIES)
control <- seeds_long %>% 
filter(TREATMENT == "CONTROL", SPECIES == SPECIES)

treatment <- function(burn,control ,plot.it=TRUE){

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

treatment(burn = burn, control = control, plot.it = TRUE)
trial.test <- treatment(SPECIES="RUBUS", burn=burn, control=control, plot.it=TRUE)

# Subset the data by species & site
site <- function(SPECIES, seeds_long){
	# Computing Likelihood for Observed Data ####
	llh_poisson <- function(lambda, y){
  # log(likelihood) by summing 
  llh <- sum(dpois(y, lambda, log=TRUE))
  return(llh)
	}

	lambdas <- seq(0,600, by=1)

	# compute log-likelihood for all lambda values
	one_ll <- sapply(lambdas,function(x){llh_poisson(x,site_1$SEEDS)})
	max.one.ll <- max(one_ll)
	one.mle <- lambdas[one_ll==max.one.ll]
	one.rel.l <- exp(one_ll - max.one.ll)
	one.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	one.cut.m.rell <- one.ci.cuttoff - one.rel.l
	one.ci.indices <- which(one.cut.m.rell<0, arr.ind=TRUE)
	one.len.ciind <- length(one.ci.indices)
	one.left.CL <- lambdas[one.ci.indices[1]]
	one.right.CL <-lambdas[one.ci.indices[one.len.ciind]]
	
	four_ll <- sapply(lambdas,function(x){llh_poisson(x,site_4$SEEDS)})
	max.four.ll <- max(four_ll)
	four.mle <- lambdas[four_ll==max.four.ll]
	four.rel.l <- exp(four_ll - max.four.ll)
	four.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	four.cut.m.rell <- four.ci.cuttoff - four.rel.l
	four.ci.indices <- which(four.cut.m.rell<0, arr.ind=TRUE)
	four.len.ciind <- length(four.ci.indices)
	four.left.CL <- lambdas[four.ci.indices[1]]
	four.right.CL <-lambdas[four.ci.indices[four.len.ciind]]
	
	six_ll <- sapply(lambdas,function(x){llh_poisson(x,site_6$SEEDS)})
	max.six.ll <- max(six_ll)
	six.mle <- lambdas[six_ll==max.six.ll]
	six.rel.l <- exp(six_ll - max.six.ll)
	six.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	six.cut.m.rell <- six.ci.cuttoff - six.rel.l
	six.ci.indices <- which(six.cut.m.rell<0, arr.ind=TRUE)
	six.len.ciind <- length(six.ci.indices)
	six.left.CL <- lambdas[six.ci.indices[1]]
	six.right.CL <-lambdas[six.ci.indices[six.len.ciind]]
	
	seven_ll <- sapply(lambdas,function(x){llh_poisson(x,site_7$SEEDS)})
	max.seven.ll <- max(seven_ll)
	seven.mle <- lambdas[seven_ll==max.seven.ll]
	seven.rel.l <- exp(seven_ll - max.seven.ll)
	seven.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	seven.cut.m.rell <- seven.ci.cuttoff - seven.rel.l
	seven.ci.indices <- which(seven.cut.m.rell<0, arr.ind=TRUE)
	seven.len.ciind <- length(seven.ci.indices)
	seven.left.CL <- lambdas[seven.ci.indices[1]]
	seven.right.CL <-lambdas[seven.ci.indices[seven.len.ciind]]
	
	eight_ll <- sapply(lambdas,function(x){llh_poisson(x,site_8$SEEDS)})
	max.eight.ll <- max(eight_ll)
	eight.mle <- lambdas[eight_ll==max.eight.ll]
	eight.rel.l <- exp(eight_ll - max.eight.ll)
	eight.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	eight.cut.m.rell <- eight.ci.cuttoff - eight.rel.l
	eight.ci.indices <- which(eight.cut.m.rell<0, arr.ind=TRUE)
	eight.len.ciind <- length(eight.ci.indices)
	eight.left.CL <- lambdas[eight.ci.indices[1]]
	eight.right.CL <-lambdas[eight.ci.indices[eight.len.ciind]]
	
	sixteen_ll <- sapply(lambdas,function(x){llh_poisson(x,site_16$SEEDS)})
	max.sixteen.ll <- max(sixteen_ll)
	sixteen.mle <- lambdas[sixteen_ll==max.sixteen.ll]
	sixteen.rel.l <- exp(sixteen_ll - max.sixteen.ll)
	sixteen.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	sixteen.cut.m.rell <- sixteen.ci.cuttoff - sixteen.rel.l
	sixteen.ci.indices <- which(sixteen.cut.m.rell<0, arr.ind=TRUE)
	sixteen.len.ciind <- length(sixteen.ci.indices)
	sixteen.left.CL <- lambdas[sixteen.ci.indices[1]]
	sixteen.right.CL <-lambdas[sixteen.ci.indices[sixteen.len.ciind]]
	
	seventeen_ll <- sapply(lambdas,function(x){llh_poisson(x,site_17$SEEDS)})
	max.seventeen.ll <- max(seventeen_ll)
	seventeen.mle <- lambdas[seventeen_ll==max.seventeen.ll]
	seventeen.rel.l <- exp(seventeen_ll - max.seventeen.ll)
	seventeen.ci.cuttoff <- exp(-qchisq(p=0.95, df=1)/2)
	seventeen.cut.m.rell <- seventeen.ci.cuttoff - seventeen.rel.l
	seventeen.ci.indices <- which(seventeen.cut.m.rell<0, arr.ind=TRUE)
	seventeen.len.ciind <- length(seventeen.ci.indices)
	seventeen.left.CL <- lambdas[seventeen.ci.indices[1]]
	seventeen.right.CL <-lambdas[seventeen.ci.indices[seventeen.len.ciind]]

	# save the lambdas and log-likelihoods in a data frame
	one_df <- data.frame(ll=one_ll, lambda=lambdas)
	one_df$SITE <- c("1")
	four_df <- data.frame(ll=four_ll, lambda=lambdas)
	four_df$SITE <- c("4")
	six_df <- data.frame(ll=six_ll, lambda=lambdas)
	six_df$SITE <- c("6")
	seven_df <- data.frame(ll=seven_ll, lambda=lambdas)
	seven_df$SITE <- c("7")
	eight_df <- data.frame(ll=eight_ll, lambda=lambdas)
	eight_df$SITE <- c("8")
	sixteen_df <- data.frame(ll=sixteen_ll, lambda=lambdas)
	sixteen_df$SITE <- c("16")
	seventeen_df <- data.frame(ll=seventeen_ll, lambda=lambdas)
	seventeen_df$SITE <- c("17")

site <- rbind(one_df, four_df, six_df,
										seven_df, eight_df, sixteen_df, seventeen_df)
return(site)
}
