# Choose the difference in mean seed arrival
pct.dec  <- 0.50 

# Create vectors for lambda values
lam.ctl.vec  <- c(5,50,500)
lam.burn.vec <- (1-pct.dec)*lam.ctl.vec 

# Create a vector for the number of seed traps
samp.sizes <- seq(1:50, by = 10)

# Create empty objects for the for loop to fill
vec <- vector() # empty vector for ith outcomes
output.mat <- matrix(0, nrow = 30, ncol = 10) # empty matrix for the final output

for (p in 1:3){
	# iterate through the lambda values
	lam.ctl <- lam.ctl.vec[p]
	lam.burn <- lam.burn.vec[p]
	
	for(i in 1:5){
		# iterate through the sample sizes 
		n.samples <- samp.sizes[i]
	
			for(j in 1:10){
				# create the columns
				
				for(z in 1:500){
					# in each column
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
}
