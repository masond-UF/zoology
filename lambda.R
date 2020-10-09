# 9 October 2020—Simulating lambda—David Mason ####
library(tidyverse)
library(ggpubr)
seeds <- read.csv("seed_trap.csv", header = T)
seeds_long <- pivot_longer(seeds, cols = 4:9, 
													 names_to = "SPECIES", values_to = "SEEDS")
# Subset the data by species & treatment ####
treatment <- function(SPECIES){
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

lambdas <- seq(0,100, by=1)

# compute log-likelihood for all lambda values
burn_ll <- sapply(lambdas,function(x){llh_poisson(x,burn$SEEDS)})
control_ll <- sapply(lambdas,function(x){llh_poisson(x,control$SEEDS)})
	
# save the lambdas and log-likelihoods in a data frame
burn_df <- data.frame(ll=burn_ll, lambda=lambdas)
burn_df$TREATMENT <- c("BURN")
control_df <- data.frame(ll=control_ll, lambda=lambdas)
control_df$TREATMENT <- c("CONTROL")

treatment <- rbind(burn_df, control_df)
return(treatment)
}
treatment_figure <- function(DATAFRAME){

burn_filt <- DATAFRAME %>% 
	filter(TREATMENT == "BURN")
control_filt <- DATAFRAME %>% 
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
site <- function(SPECIES){
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
four_ll <- sapply(lambdas,function(x){llh_poisson(x,site_4$SEEDS)})
six_ll <- sapply(lambdas,function(x){llh_poisson(x,site_6$SEEDS)})
seven_ll <- sapply(lambdas,function(x){llh_poisson(x,site_7$SEEDS)})
eight_ll <- sapply(lambdas,function(x){llh_poisson(x,site_8$SEEDS)})
sixteen_ll <- sapply(lambdas,function(x){llh_poisson(x,site_16$SEEDS)})
seventeen_ll <- sapply(lambdas,function(x){llh_poisson(x,site_17$SEEDS)})

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
rubus_site <- site("RUBUS")
# Maximum Likelihood Estimate from Observed Data ####
lambdas <- seq(0,100, by=1)

one_filt <- rubus_site %>% 
	filter(SITE == "1")
four_filt <- rubus_site %>% 
	filter(SITE == "4")
six_filt <- rubus_site %>% 
	filter(SITE == "6")
seven_filt <- rubus_site %>% 
	filter(SITE == "7")
eight_filt <- rubus_site %>% 
	filter(SITE == "8")
sixteen_filt <- rubus_site %>% 
	filter(SITE == "16")
seventeen_filt <- rubus_site %>% 
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
