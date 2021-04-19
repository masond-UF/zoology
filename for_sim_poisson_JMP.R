# 23 October 2020—Simulating data—David Mason ####
# Subset the data ####
burn_sub <- as.data.frame(seeds_long %>% 
filter(TREATMENT == "BURN"))
control_sub <- as.data.frame(seeds_long %>% 
filter(TREATMENT == "CONTROL"))

site_1 <- rubus %>% 
		filter(SITE == "1")
site_4 <- rubus %>% 
		filter(SITE == "4")
site_6 <- rubus %>% 
		filter(SITE == "6")
site_7 <- rubus %>% 
		filter(SITE == "7")
site_8 <- rubus %>% 
		filter(SITE == "8")
site_16 <- rubus %>%
		filter(SITE == "16")
site_17 <- rubus %>% 
		filter(SITE == "17")

# Create the empty matrices ####
pois.test.results <- matrix(,nrow = 21, ncol = 2)
colnames(pois.test.results) <- c('samples', 'outcome')
pois.test.results <- as.data.frame(pois.test.results)

sample.test.results <- matrix(,nrow = 21, ncol = 2)
colnames(sample.test.results) <- c('samples', 'outcome')
sample.test.results <- as.data.frame(sample.test.results)

# Using samples drawn from a poisson distribution V1 ####
pct.dec  <- 0.05 
pct.dec  <- 0.15 
pct.dec  <- 0.30 
pct.dec  <- 0.50 
pct.dec  <- 0.90 

lam.ctl  <- c(2,4,6,8,10,30,50,100, 500)
lam.burn <- (1-pct.dec)*lam.ctl 
nlams <- length(lam.ctl)

list.out <-list()

# treatment
for(i in 1:nlams){
	
	# generate samples from poisson distribution
	burn.numbers <- rpois(n = 21, lam.burn[i])
	control.numbers <- rpois(n = 21, lam.ctl[i])
	# turn it into a dataframe the function reocgnizes
	burn <- matrix(burn.numbers, ncol = 1)
	colnames(burn) <- c("SEEDS")
	burn <- as.data.frame(burn)
	
	control <- matrix(control.numbers, ncol = 1)
	colnames(control) <- c("SEEDS")
	control <- as.data.frame(control)
	
	ith.test <- treatment(burn=burn,control=control ,plot.it=FALSE)
	list.out[[i]] <- ith.test
	
}
names(list.out) <- paste0(rep("lambda ctl = ",nlams), lam.ctl)

# site
pct.dec.a  <- 0.05 
pct.dec.b  <- 0.15 
pct.dec.c  <- 0.30 

lam.site.1  <- c(2,4,6,8,10,30,50,100)
lam.site.4 <- (1-pct.dec.a)*lam.ctl 
lam.site.6 <- (1-pct.dec.b)*lam.ctl 
lam.site.7 <- (1-pct.dec.c)*lam.ctl 
lam.site.8 <- (1-pct.dec.a)*lam.ctl 
lam.site.16 <- (1-pct.dec.b)*lam.ctl 
lam.site.17 <- (1-pct.dec.c)*lam.ctl 


nlams <- length(lam.site.1)

list.out <-list()

for(i in 1:nlams){
	
	# generate samples from poisson distribution
	site.1.numbers <- rpois(n = 6, lam.site.1[i])
	site.4.numbers <- rpois(n = 6, lam.site.4[i])
	site.6.numbers <- rpois(n = 6, lam.site.6[i])
	site.7.numbers <- rpois(n = 6, lam.site.7[i])
	site.8.numbers <- rpois(n = 6, lam.site.8[i])
	site.16.numbers <- rpois(n = 6, lam.site.16[i])
	site.17.numbers <- rpois(n = 6, lam.site.17[i])


	# turn it into a dataframe the function reocgnizes
	site_1 <- matrix(site.1.numbers, ncol = 1)
	colnames(site_1) <- c("SEEDS")
	site_1 <- as.data.frame(site_1)
	
	site_4 <- matrix(site.4.numbers, ncol = 1)
	colnames(site_4) <- c("SEEDS")
	site_4 <- as.data.frame(site_4)
	
	site_6 <- matrix(site.6.numbers, ncol = 1)
	colnames(site_6) <- c("SEEDS")
	site_6 <- as.data.frame(site_6)
	
	site_7 <- matrix(site.7.numbers, ncol = 1)
	colnames(site_7) <- c("SEEDS")
	site_7 <- as.data.frame(site_7)
	
	site_8 <- matrix(site.8.numbers, ncol = 1)
	colnames(site_8) <- c("SEEDS")
	site_8 <- as.data.frame(site_8)
	
	site_16 <- matrix(site.16.numbers, ncol = 1)
	colnames(site_16) <- c("SEEDS")
	site_16 <- as.data.frame(site_16)
	
	site_17 <- matrix(site.17.numbers, ncol = 1)
	colnames(site_17) <- c("SEEDS")
	site_17 <- as.data.frame(site_17)
	
	ith.test <- site(site_1=site_1,site_4=site_4, site_6=site_6, 
									 site_7=site_7,site_8=site_8, site_16=site_16,
								 	 site_17=site_17,plot.it=FALSE)
	list.out[[i]] <- ith.test
	
}
names(list.out) <- paste0(rep("lambda site 1 = ",nlams), lam.site.1)

# Using samples from the data ####
n <- 1
for(i in 1:21){
	
	# Build the counter
	n <- n + 1
	
	# Pull n random rows from the data
	burn <- sample_n(burn_sub, n)
	control <- sample_n(control_sub, n)

	test <- treatment(burn=burn, control=control, plot.it=FALSE)
	sample.test.results$samples[i] <- i
	sample.test.results$outcome[i] <- test$outcome
}

