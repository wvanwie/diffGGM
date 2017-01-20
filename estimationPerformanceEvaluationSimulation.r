#############################################################################################################################################
#
# READ ME:
#
# This R-script serves to investigate the performance of the 'average' differential conditional dependence parameter. 
# Hereto data are simulated in accordance with the proposed model.
# Data are drawn for models of different sizes: p=25, 50, 75, 100.
# In the simulation design parameters (n1, n2 and gamma) are set (and need to be modified manually to investigate other settings).
# From simulated data then the model parameters are estimated for various choices of lambda. 
# The above is repeated 100 times. 
# The estimation is then followed by aggregation of the results per model-size.
# Now modify the design parameters (n1, n2 and gamma). 
# For the use of the plot script, modify gamma to be equal to 0.0, 0.1, 0.2, 0.3, 0.4.
# Gather the aggregated output files and move to plotting. 
#
#############################################################################################################################################

# erase history
rm(list=ls())

# activate libraries
library(mvtnorm)
library(Matrix)
library(diffGGM)

# set parameters (modify)
n1 <- n2 <- 100
gamma <- 0.5

# set other parameters
nInit <- 100
ps <- seq(25, 100, 25)
lambda <- exp(seq(log(0.01), log(1), length.out=100)) 

# generate precision matrix
diags <- list(rep(1, max(ps)), rep(0.5, max(ps)-1), rep(0.25, max(ps)-2), rep(0.1, max(ps)-3))
Omega <- as.matrix(bandSparse(max(ps), k = -c(0:3), diag = c(diags), symm=TRUE))

for (j in 1:nInit){
	# draw gene expression data
	gHats1 <- list()
	gHats2 <- list()	
	gHats3 <- list()
	gHats4 <- list()	

	Sigma1 <- solve(Omega[1:ps[1], 1:ps[1]])
	Sigma2 <- solve((1-gamma) * Omega[1:ps[1], 1:ps[1]] + gamma * diag(ps[1]))
	Y1 <- rmvnorm(n1, sigma=Sigma1)
	Y2 <- rmvnorm(n2, sigma=Sigma2)
	S1a <- t(Y1) %*% Y1 / n1
	S2a <- t(Y2) %*% Y2 / n2

	Sigma1 <- solve(Omega[1:ps[2], 1:ps[2]])
	Sigma2 <- solve((1-gamma) * Omega[1:ps[2], 1:ps[2]] + gamma * diag(ps[2]))
	Y1 <- rmvnorm(n1, sigma=Sigma1)
	Y2 <- rmvnorm(n2, sigma=Sigma2)
	S1b <- t(Y1) %*% Y1 / n1
	S2b <- t(Y2) %*% Y2 / n2

	Sigma1 <- solve(Omega[1:ps[3], 1:ps[3]])
	Sigma2 <- solve((1-gamma) * Omega[1:ps[3], 1:ps[3]] + gamma * diag(ps[3]))
	Y1 <- rmvnorm(n1, sigma=Sigma1)
	Y2 <- rmvnorm(n2, sigma=Sigma2)
	S1c <- t(Y1) %*% Y1 / n1
	S2c <- t(Y2) %*% Y2 / n2

	Sigma1 <- solve(Omega[1:ps[4], 1:ps[4]])
	Sigma2 <- solve((1-gamma) * Omega[1:ps[4], 1:ps[4]] + gamma * diag(ps[4]))
	Y1 <- rmvnorm(n1, sigma=Sigma1)
	Y2 <- rmvnorm(n2, sigma=Sigma2)
	S1d <- t(Y1) %*% Y1 / n1
	S2d <- t(Y2) %*% Y2 / n2

	for (k in 1:length(lambda)){
		print(paste("iteration: ",j, "; lambda: ", lambda[k], sep=""))	
		gHats1[[k]] <- diffGGMml(S1a, S2a, n1, n2, lambda[k])
		gHats2[[k]] <- diffGGMml(S1b, S2b, n1, n2, lambda[k])
		gHats3[[k]] <- diffGGMml(S1c, S2c, n1, n2, lambda[k])
		gHats4[[k]] <- diffGGMml(S1d, S2d, n1, n2, lambda[k])
	}
	save(ps, n1, n2, gHats1, gHats2, gHats3, gHats4, j, nInit, lambda, Sigma1, Sigma2, Omega, gamma, file=paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", j, "_diffGGMml.Rdata", sep=""))
}


# aggregate results for p=25
it <- 1
load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", 1, "_diffGGMml.Rdata", sep=""))
pcor1 <- matrix(NA, nrow=nInit, ncol=length(gHats1))
pcor2 <- matrix(NA, nrow=nInit, ncol=length(gHats1))
pcor3 <- matrix(NA, nrow=nInit, ncol=length(gHats1))
pcorr <- matrix(NA, nrow=nInit, ncol=length(gHats1))
u1 <- matrix(NA, nrow=nInit, ncol=length(gHats1))
u2 <- matrix(NA, nrow=nInit, ncol=length(gHats1))
gammaHat <- matrix(NA, nrow=nInit, ncol=length(gHats1))
for (it in 1:nInit){
	load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", it, "_diffGGMml.Rdata", sep=""))
	print(it)
	for (k in 1:length(gHats1)){
		gammaHat[it,k] <- gHats1[[k]]$gamma
		p <- ncol(gHats1[[k]]$Ps)
		pcor1[it,k] <- sum(band(gHats1[[k]]$Ps, 1, 1)) / (p - 1)
		pcor2[it,k] <- sum(band(gHats1[[k]]$Ps, 2, 2)) / (p - 2)
		pcor3[it,k] <- sum(band(gHats1[[k]]$Ps, 3, 3)) / (p - 3)
		pcorr[it,k] <- sum(band(gHats1[[k]]$Ps, 4, p)) / (p*(p-1)/2 - (3*p - 6))
		u1[it,k] <- mean(gHats1[[k]]$U1)
		u2[it,k] <- mean(gHats1[[k]]$U2)
	}
}
save(lambda, gamma, gammaHat, pcor1, pcor2, pcor3, pcorr, u1, u2, p, file=paste("simuRes_n", n1, "_gamma", gamma, "_p_", p, "_diffGGMml.Rdata", sep=""))

# aggregate results for p=50
it <- 1
load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", 1, "_diffGGMml.Rdata", sep=""))
pcor1 <- matrix(NA, nrow=nInit, ncol=length(gHats2))
pcor2 <- matrix(NA, nrow=nInit, ncol=length(gHats2))
pcor3 <- matrix(NA, nrow=nInit, ncol=length(gHats2))
pcorr <- matrix(NA, nrow=nInit, ncol=length(gHats2))
u1 <- matrix(NA, nrow=nInit, ncol=length(gHats2))
u2 <- matrix(NA, nrow=nInit, ncol=length(gHats2))
gammaHat <- matrix(NA, nrow=nInit, ncol=length(gHats2))
for (it in 1:nInit){
	load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", it, "_diffGGMml.Rdata", sep=""))
	print(it)
	for (k in 1:length(gHats2)){
		gammaHat[it,k] <- gHats2[[k]]$gamma
		p <- ncol(gHats2[[k]]$Ps)
		pcor1[it,k] <- sum(band(gHats2[[k]]$Ps, 1, 1)) / (p - 1)
		pcor2[it,k] <- sum(band(gHats2[[k]]$Ps, 2, 2)) / (p - 2)
		pcor3[it,k] <- sum(band(gHats2[[k]]$Ps, 3, 3)) / (p - 3)
		pcorr[it,k] <- sum(band(gHats2[[k]]$Ps, 4, p)) / (p*(p-1)/2 - (3*p - 6))
		u1[it,k] <- mean(gHats2[[k]]$U1)
		u2[it,k] <- mean(gHats2[[k]]$U2)        	
	}
}
save(lambda, gamma, gammaHat, pcor1, pcor2, pcor3, pcorr, u1, u2, p, file=paste("simuRes_n", n1, "_gamma", gamma, "_p_", p, "_diffGGMml.Rdata", sep=""))

# aggregate results for p=75
it <- 1
load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", 1, "_diffGGMml.Rdata", sep=""))
pcor1 <- matrix(NA, nrow=nInit, ncol=length(gHats3))
pcor2 <- matrix(NA, nrow=nInit, ncol=length(gHats3))
pcor3 <- matrix(NA, nrow=nInit, ncol=length(gHats3))
pcorr <- matrix(NA, nrow=nInit, ncol=length(gHats3))
u1 <- matrix(NA, nrow=nInit, ncol=length(gHats3))
u2 <- matrix(NA, nrow=nInit, ncol=length(gHats3))
gammaHat <- matrix(NA, nrow=nInit, ncol=length(gHats3))
for (it in 1:nInit){
	load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", it, "_diffGGMml.Rdata", sep=""))
	print(it)
	for (k in 1:length(gHats3)){
		gammaHat[it,k] <- gHats3[[k]]$gamma
		p <- ncol(gHats3[[k]]$Ps)
		pcor1[it,k] <- sum(band(gHats3[[k]]$Ps, 1, 1)) / (p - 1)
		pcor2[it,k] <- sum(band(gHats3[[k]]$Ps, 2, 2)) / (p - 2)
		pcor3[it,k] <- sum(band(gHats3[[k]]$Ps, 3, 3)) / (p - 3)
		pcorr[it,k] <- sum(band(gHats3[[k]]$Ps, 4, p)) / (p*(p-1)/2 - (3*p - 6))
		u1[it,k] <- mean(gHats3[[k]]$U1)
		u2[it,k] <- mean(gHats3[[k]]$U2)
	}
}
save(lambda, gamma, gammaHat, pcor1, pcor2, pcor3, pcorr, u1, u2, p, file=paste("simuRes_n", n1, "_gamma", gamma, "_p_", p, "_diffGGMml.Rdata", sep=""))

# aggregate results for p=100
it <- 1
load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", 1, "_diffGGMml.Rdata", sep=""))
pcor1 <- matrix(NA, nrow=nInit, ncol=length(gHats4))
pcor2 <- matrix(NA, nrow=nInit, ncol=length(gHats4))
pcor3 <- matrix(NA, nrow=nInit, ncol=length(gHats4))
pcorr <- matrix(NA, nrow=nInit, ncol=length(gHats4))
u1 <- matrix(NA, nrow=nInit, ncol=length(gHats4))
u2 <- matrix(NA, nrow=nInit, ncol=length(gHats4))
gammaHat <- matrix(NA, nrow=nInit, ncol=length(gHats4))
for (it in 1:nInit){
	load(paste("simuEstimation_n", n1, "_gamma", gamma, "_it_", it, "_diffGGMml.Rdata", sep=""))
	print(it)
	for (k in 1:length(gHats4)){
		gammaHat[it,k] <- gHats4[[k]]$gamma
		p <- ncol(gHats4[[k]]$Ps)
		pcor1[it,k] <- sum(band(gHats4[[k]]$Ps, 1, 1)) / (p - 1)
		pcor2[it,k] <- sum(band(gHats4[[k]]$Ps, 2, 2)) / (p - 2)
		pcor3[it,k] <- sum(band(gHats4[[k]]$Ps, 3, 3)) / (p - 3)
		pcorr[it,k] <- sum(band(gHats4[[k]]$Ps, 4, p)) / (p*(p-1)/2 - (3*p - 6))
		u1[it,k] <- mean(gHats4[[k]]$U1)
		u2[it,k] <- mean(gHats4[[k]]$U2)		
	}
}
save(lambda, gamma, gammaHat, pcor1, pcor2, pcor3, pcorr, u1, u2, p, file=paste("simuRes_n", n1, "_gamma", gamma, "_p_", p, "_diffGGMml.Rdata", sep=""))


