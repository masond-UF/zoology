# To do a for loop, we need 2 things:
# 1. how many times I am going to do a calculation
# 2. which calculation I am going to do


# Each time we create a for loop, r creates a counter variable
# that keeps track of the iteration number
# This counter is user-defined (in a smart way)!!!

ts.data <- c(200,180,179,166,157,145,135,120,112,105)

len <- length(ts.data)

for(i in 1:(len-1)){
	
	ni <- ts.data[i]
	nip1 <- ts.data[(i+1)]
	
	print(c("i","ni","nip1"))
	print(c(i, ni, nip1))
	
}


for(i in 2:len){
	
	nim1 <- ts.data[(i-1)]
	ni <- ts.data[i]
	if(i==2){print(c("i nim1 ni"))}
	print(c(i, nim1, ni))
	
}


for(i in len:2){
	
	ni <- ts.data[i]
	nim1 <- ts.data[(i-1)]
	if(i==len){print(c("i nim1 ni"))}
	print(c(i, nim1, ni))
	
}


multiple.tsdata <- rbind(ts.data,ts.data,ts.data)
row.names(multiple.tsdata) <- c("rep1", "rep2", "rep3" )
colnames(multiple.tsdata) <- paste(rep("Time", 10), 1:10 )

# Selecting and printing one entire column at a time


for(i in 2:len){
	
	ni.vec <- multiple.tsdata[,i]
	nim1.vec <- multiple.tsdata[,(i-1)]
	if(i==len){print(c("nim1.vec ni.vec"))}
	print(cbind(nim1.vec, ni.vec))
	
}


output.mat <- matrix(0,nrow=10,ncol=20)
means.mat <- matrix(0,nrow=10,ncol=20)
samps.mat <- matrix(0,nrow=10,ncol=20)

musNsamps <- matrix("NA",nrow=10,ncol=20)



mu.vec <- (1:10)*10
samp.sizes <- seq(from=100,to=2000, by=100)

for(i in 1:10){

	ith.mu <- mu.vec[i]
	
	for(j in 1:20){
		
		jth.samp.s <- samp.sizes[j]
		
		means.mat[i,j] <- ith.mu
		samps.mat[i,j] <- jth.samp.s
		musNsamps[i,j] <- paste0(ith.mu,",",jth.samp.s)
		
		
		output.mat[i,j] <- mean(rnorm(n=jth.samp.s, mean=ith.mu, sd=1))
		
		
	}
}











