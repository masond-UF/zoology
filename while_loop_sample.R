n <- 21
test <- treatment(burn=burn, control=control, plot.it=FALSE)

while(test$outcome == "Separate model is best"){
	burn <- sample_n(burn, n)
	control <- sample_n(control, n)
	test <- treatment(burn=burn, control=control, plot.it=FALSE)
	print(n)
 	n  <-  n-1
}
