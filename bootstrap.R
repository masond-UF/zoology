pct.dec  <- 0.50 
lam.ctl  <- c(2,4,6,8,10,30,50,100)
lam.burn <- (1-pct.dec)*lam.ctl 
nlams <- length(lam.ctl)

list.out <-list()


poisson.regsim <- function(lam.burn, lam.ctl){
	
	# generate samples from poisson distribution
	burn.numbers <- rpois(n = 21, lam.burn)
	control.numbers <- rpois(n = 21, lam.ctl)
	# turn it into a dataframe the function reocgnizes
	burn <- matrix(burn.numbers, ncol = 1)
	colnames(burn) <- c("SEEDS")
	burn <- as.data.frame(burn)
	
	control <- matrix(control.numbers, ncol = 1)
	colnames(control) <- c("SEEDS")
	control <- as.data.frame(control)
	
	test <- treatment(burn=burn,control=control ,plot.it=FALSE)
	return(test)
}
poisson.regsim(lam.burn = lam.burn, lam.ctl = lam.ctl)

B <- 2000

for(i in 1:B){
	ith.test <- poisson.regsim(lam.burn = lam.burn, lam.ctl = lam.ctl)
	list.out[[i]] <- ith.test
}

df <- data.frame(matrix(unlist(list.out), nrow=length(list.out), byrow=T))

