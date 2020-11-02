# 23 October 2020—Simulating data—David Mason ####
# Subset the data ####
burn_sub <- as.data.frame(seeds_long %>% 
filter(TREATMENT == "BURN"))
control_sub <- as.data.frame(seeds_long %>% 
filter(TREATMENT == "CONTROL"))

# Create the empty matrices ####
pois.test.results <- matrix(,nrow = 21, ncol = 2)
colnames(pois.test.results) <- c('samples', 'outcome')
pois.test.results <- as.data.frame(pois.test.results)

sample.test.results <- matrix(,nrow = 21, ncol = 2)
colnames(sample.test.results) <- c('samples', 'outcome')
sample.test.results <- as.data.frame(sample.test.results)

# Using samples drawn from a poisson distribution ####
pct.dec  <- 0.15 
lam.ctl  <- c(2,4,6,8,10,30,50,100)
lam.burn <- (1-pct.dec)*lam.ctl 
nlams <- length(lam.ctl)

list.out <-list()
	
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
	
	ith.test <- treatment(burn,control ,plot.it=TRUE)
	list.out[[i]] <- ith.test
	
}

names(list.out) <- paste0(rep("lambda ctl = ",nlams), lam.ctl)

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

