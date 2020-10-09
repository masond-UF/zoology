# 9 October 2020—Simulating lambda—David Mason ####
library(tidyverse)
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

lambdas <- seq(0,600, by=1)

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
rubus_treatment <- treatment("RUBUS")
# Maximum Likelihood Estimate from Observed Data ####
rubus_treatment %>% 
  ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(ll)], color="red",size=2) +
	scale_x_continuous(breaks = seq(0, 600, by = 50)) +
	facet_wrap(~TREATMENT)
# Subset the data by species & site ####
rubus_site_1 <- seeds_long %>% 
	filter(SITE == "1", SPECIES == "RUBUS")
rubus_site_4 <- seeds_long %>% 
	filter(SITE == "4", SPECIES == "RUBUS")
rubus_site_6 <- seeds_long %>% 
	filter(SITE == "6", SPECIES == "RUBUS")
rubus_site_7 <- seeds_long %>% 
	filter(SITE == "7", SPECIES == "RUBUS")
rubus_site_8 <- seeds_long %>% 
	filter(SITE == "8", SPECIES == "RUBUS")
rubus_site_16 <- seeds_long %>% 
	filter(SITE == "16", SPECIES == "RUBUS")
rubus_site_17 <- seeds_long %>% 
	filter(SITE == "17", SPECIES == "RUBUS")
# Computing Likelihood for Observed Data ####
llh_poisson <- function(lambda, y){
  # log(likelihood) by summing 
  llh <- sum(dpois(y, lambda, log=TRUE))
  return(llh)
}

lambdas <- seq(0,600, by=1)

# compute log-likelihood for all lambda values
rubus_1_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_1$SEEDS)})
rubus_4_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_4$SEEDS)})
rubus_6_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_6$SEEDS)})
rubus_7_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_7$SEEDS)})
rubus_8_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_8$SEEDS)})
rubus_16_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_16$SEEDS)})
rubus_17_ll <- sapply(lambdas,function(x){llh_poisson(x,rubus_site_17$SEEDS)})

# save the lambdas and log-likelihoods in a data frame
rubus_1_df <- data.frame(ll=rubus_1_ll, lambda=lambdas)
rubus_1_df$SITE <- c("1")
rubus_4_df <- data.frame(ll=rubus_4_ll, lambda=lambdas)
rubus_4_df$SITE <- c("4")
rubus_6_df <- data.frame(ll=rubus_6_ll, lambda=lambdas)
rubus_6_df$SITE <- c("6")
rubus_7_df <- data.frame(ll=rubus_7_ll, lambda=lambdas)
rubus_7_df$SITE <- c("7")
rubus_8_df <- data.frame(ll=rubus_8_ll, lambda=lambdas)
rubus_8_df$SITE <- c("8")
rubus_16_df <- data.frame(ll=rubus_16_ll, lambda=lambdas)
rubus_16_df$SITE <- c("16")
rubus_17_df <- data.frame(ll=rubus_17_ll, lambda=lambdas)
rubus_17_df$SITE <- c("17")

rubus_site <- rbind(rubus_1_df, rubus_4_df, rubus_6_df,
										rubus_7_df, rubus_8_df, rubus_16_df, rubus_17_df)

# Maximum Likelihood Estimate from Observed Data ####
rubus_site %>% 
  ggplot(aes(x=lambda,y=ll))+
  geom_point(size=4,color="dodgerblue")+
  xlab("Lambda") +
  ylab("Log Likelihood")+
  theme_bw(base_size = 16) +
  geom_vline(xintercept = lambdas[which.max(ll)], color="red",size=2) +
	scale_x_continuous(breaks = seq(0, 600, by = 50)) +
	facet_wrap(~SITE)
