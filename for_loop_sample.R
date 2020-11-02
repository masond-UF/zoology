for(i in 1:42){
	# Filter so its just the species indicated
	burn <- filter(burn, SPECIES == SPECIES)
	control <- filter(control, SPECIES == SPECIES)
	
	x[i] <- x + 1
	
	# select random number of rows equal to iterations
	burn <- sample_n(burn, n)
	control <- sample_n(control, n)
	
	test <- treatment(burn=burn, control=control, plot.it=FALSE)
	test$outcome == "Separate model is best"
}