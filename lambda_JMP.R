# 9 October 2020—Simulating lambda—David Mason ####
library(tidyverse)
library(ggpubr)
seeds <- read.csv("seed_trap.csv", header = T)
seeds_long <- pivot_longer(seeds, cols = 4:9,names_to = "SPECIES", values_to = "SEEDS")
species <- c("RUBUS", "PRSE2", "PIPA2", "SMILAX", "PHAM4", "PIEC2")

# Subset the data by species & treatment ####
treatment <- function(SPECIES, seeds_long){
	burn <- seeds_long %>% 
	filter(TREATMENT == "BURN", SPECIES == SPECIES)
	control <- seeds_long %>% 
	filter(TREATMENT == "CONTROL", SPECIES == SPECIES)

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
	
	# save the lambdas and log-likelihoods in a data frame
	burn_df <- data.frame(ll=burn_ll, lambda=lambdas)
	burn_df$TREATMENT <- c("BURN")
	control_df <- data.frame(ll=control_ll, lambda=lambdas)
	control_df$TREATMENT <- c("CONTROL")

	treatment <- rbind(burn_df, control_df)
	return(treatment)
	
}
treatment_figure <- function(DATAFRAME){

d <- DATAFRAME	
	
burn_filt <- d %>% 
	filter(TREATMENT == "BURN")
control_filt <- d %>% 
	filter(TREATMENT == "CONTROL")

# Maximum Likelihood Estimate from Observed Data ####

lambdas <- seq(0,100, by=1)

burn_plot <- burn_filt %>% 
  ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(burn_filt$ll)], color="red",size=2) +
	scale_x_continuous(breaks = seq(0, 100, by = 10))
control_plot <- control_filt %>% 
  ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(control_filt$ll)], color="red",size=2) +
	scale_x_continuous(breaks = seq(0, 100, by = 10))

combo_plot <- ggarrange(burn_plot,control_plot)
}
# Subset the data by species & site ####
site <- function(SPECIES, seeds_long){
	site_1 <- seeds_long %>% 
		filter(SITE == "1", SPECIES == SPECIES)
	site_4 <- seeds_long %>% 
		filter(SITE == "4", SPECIES == SPECIES)
	site_6 <- seeds_long %>% 
		filter(SITE == "6", SPECIES == SPECIES)
	site_7 <- seeds_long %>% 
		filter(SITE == "7", SPECIES == SPECIES)
	site_8 <- seeds_long %>% 
		filter(SITE == "8", SPECIES == SPECIES)
	site_16 <- seeds_long %>% 
		filter(SITE == "16", SPECIES == SPECIES)
	site_17 <- seeds_long %>% 
		filter(SITE == "17", SPECIES == SPECIES)
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
site_figure <- function(DATAFRAME){
# Maximum Likelihood Estimate from Observed Data ####
lambdas <- seq(0,100, by=1)
d <- DATAFRAME

one_filt <- d %>% 
	filter(SITE == "1")
four_filt <- d %>% 
	filter(SITE == "4")
six_filt <- d %>% 
	filter(SITE == "6")
seven_filt <- d %>% 
	filter(SITE == "7")
eight_filt <- d %>% 
	filter(SITE == "8")
sixteen_filt <- d %>% 
	filter(SITE == "16")
seventeen_filt <- d %>% 
	filter(SITE == "17")

one_plot <- one_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(one_filt$ll)], color="red",size=2)

four_plot <- four_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(four_filt$ll)], color="red",size=2)

six_plot <- six_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(six_filt$ll)], color="red",size=2)

seven_plot <- seven_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(seven_filt$ll)], color="red",size=2)

eight_plot <- eight_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(eight_filt$ll)], color="red",size=2)

sixteen_plot <- sixteen_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(sixteen_filt$ll)], color="red",size=2)

seventeen_plot <- seventeen_filt %>% 
	ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(seventeen_filt$ll)], color="red",size=2)

combo_plot <- ggarrange(one_plot,four_plot, six_plot, seven_plot,
												eight_plot, sixteen_plot, seventeen_plot)
} 
