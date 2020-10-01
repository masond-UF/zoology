# 30 September 2020—Simulating a data set
# David Mason
library(tidyverse)
seeds <- read.csv("seed_trap.csv", header = T)
seeds <- seeds[,2:6]

spp.probs <- c(0.43,   0.24,   0.58,   0.81,  0.13)
survival <- matrix(rep(spp.probs,nrow(seeds)), nrow=nrow(seeds), ncol=5, byrow=TRUE)

colnames(survival) <- c('RUBUS', 'PRSE2', 'PIPA2', 'SMILAX', 'PHAM4')


germ <- matrix(nrow = 42, ncol = 5)
seeds <- as.matrix(seeds)
for(i in 1:42){
  germ[i,] <-  rbinom(n = 5, size = seeds[i,], prob = survival[i,])
}


