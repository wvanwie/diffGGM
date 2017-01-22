# erase history
rm(list=ls())

# activate libraries
library(mvtnorm)
library(Matrix)
library(rags2ridges)
library(diffGGM)

# set number of iterations
nInit <- 100

# set differential conditional dependence multiplier
gammaHa <- 0.0

# set sample size and dimension
n1 <- n2 <- 50
p <- 25

# generate precision matrix
diags <- list(rep(1, p), rep(0.5, p-1), rep(0.25, p-2), rep(0.1, p-3))
Omega <- as.matrix(bandSparse(p, k = -c(0:3), diag = c(diags), symm=TRUE))

# set significance levels
alphas <- c(1:50)/100

# set number of observation in null distribution
nNulldist <- 1000

# declare variables for results storing
calphas <- list()
nullDist <- list()
optLambdas <- list()
power <- list()
gammaObsDist <- list()

for (j in 1:nInit){

	print(paste("start iteration : ", j, sep=""), quote=FALSE)

	print("-> draw null data", quote=FALSE)
	Sigma1 <- solve(Omega)
	Y1 <- rmvnorm(n1, sigma=Sigma1)
	Y2 <- rmvnorm(n2, sigma=Sigma1)
	stage <- c(rep(0, n1), rep(1, n2))
	Yboth <- rbind(Y1, Y2)

	print("-> determine penalty parameter through LOOCV, using null data, repetitively", quote=FALSE)
	optLambda <- numeric()
	for (u in 1:10){
		Yshuffled <- Yboth[sample(1:nrow(Yboth), nrow(Yboth)),]
		optLambda <- c(optLambda, optPenalty.diffGGMnull(Yshuffled[which(stage==0),], Yshuffled[which(stage==1),], 10^(-10), 1, lambdaInit=0.1))
	}
	optLambdas[[j]] <- mean(optLambda)

	print("-> generate the null distribution through permutation for the given lambda", quote=FALSE)
	nullSlh <- numeric()
	for (u in 1:nNulldist){
		Yshuffled <- Yboth[sample(1:nrow(Yboth), nrow(Yboth)),]
		nullSlh <- c(nullSlh, diffGGMml(covML(Yshuffled[which(stage==0),]), covML(Yshuffled[which(stage==1),]), n1, n2, lambda=mean(optLambda), Uequal=TRUE)$gamma)
	}
	nullDist[[j]] <- nullSlh

	print("-> determine critical values", quote=FALSE)
	calphasSlh <- numeric(length=length(alphas))
	for (a in 1:length(alphas)){
		calphasSlh[a] <- sort(nullSlh)[ceiling(nNulldist * (1-alphas[a]))]
	}
	calphas[[j]] <- calphasSlh

	print("-> set Ha through specification of gamma, sample from Ha, estimate gamma", quote=FALSE)
	gammaObs <- numeric(length=nNulldist)
	for (u in 1:nNulldist){
		Sigma1 <- solve(Omega)
		Sigma2 <- solve((1-gammaHa) * Omega + gammaHa * diag(p))
		Y1 <- rmvnorm(n1, sigma=Sigma1)
		Y2 <- rmvnorm(n2, sigma=Sigma2)
		gammaObs[u] <- diffGGMml(covML(Y1), covML(Y2), n1, n2, lambda=mean(optLambda), Uequal=TRUE)$gamma
		
	}
	gammaObsDist[[j]] <- gammaObs

	print("-> assess test performance", quote=FALSE)
	powerSlh <- numeric(length=length(alphas))
	for (a in 1:length(alphas)){
		powerSlh[a] <- sum(gammaObs > calphas[[j]][a]) / nNulldist
	}
	power[[j]] <- powerSlh

	print("-> store results", quote=FALSE)
	save(nullDist, optLambdas, gammaObsDist, calphas, power, file=paste("power_n", n1,  "_p", p, "_gamma", gammaHa, ".Rdata", sep=""))
}


