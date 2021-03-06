##########################################################################################
#
# READ ME
#
# Analysis hedgehog signalling pathway data from the GEO prostate cancer study GSE3933.
# The script loads the data, included in the diffGGM package.
# The penalty parameter is chosen by LOOCV using `null data'.
# With the penalty parameter at hand, the observed test statistic is calculated.
# Followed by the construction of the null distribution by permutation. 
# Using the null distribution, the p-value is calculated.
# The test is performed for four different situation:
# - Conditional dependencies are assumed to be weaker in the cancer group.
# - Conditional dependencies are assumed to be weaker in the normal group.
# Both analyses are also done on Gaussianized data.
# 
##########################################################################################

# clear history 
rm(list=ls())

# load packages
library(Biobase)
library(diffGGM)
library(rags2ridges)
data(pcaHHs)

# extract data
Y <- t(exprs(pcaHHgse3933))
stage <- as.numeric(pData(pcaHHgse3933)[,1])

# determine penalty parameter through LOOCV, repetitively
# for both weakened conditional dependencies in the cancer stage and the reverse (i.e. weaker in the normal)
optLambdaCweak <- optLambdaCstrong <- numeric()
for (u in 1:10){
	shuffle <- sample(1:length(stage), length(stage))
	optLambdaCweak <- c(optLambdaCweak, optPenalty.diffGGMnull(Y[which(stage[shuffle]==0),], Y[which(stage[shuffle]==1),], 10^(-10), 1, lambdaInit=0.1))
	optLambdaCstrong <- c(optLambdaCstrong, optPenalty.diffGGMnull(Y[which(stage[shuffle]==1),], Y[which(stage[shuffle]==0),], 10^(-10), 1, lambdaInit=0.1))
}

# estimate gamma (= observed test statistic)
gammaObs <- rbind(c("gse3933", sum(stage==0), sum(stage==1), ncol(Y), "no", "no", "n.a.", median(optLambdaCweak), 
			diffGGMml(covML(Y[which(stage==0),]), covML(Y[which(stage==1),]), sum(stage==0), sum(stage==1), lambda=median(optLambdaCweak), Uequal=TRUE)$gamma, NA), 
		   c("gse3933", sum(stage==0), sum(stage==1), ncol(Y), "yes", "no", "n.a.", median(optLambdaCstrong), 
			diffGGMml(covML(Y[which(stage==1),]), covML(Y[which(stage==0),]), sum(stage==1), sum(stage==0), lambda=median(optLambdaCstrong), Uequal=TRUE)$gamma, NA))
colnames(gammaObs) <- c("dataset", "n1", "n2", "p", "c>n", "gauss", "outlier", "lambda", "gamma", "p-value")

# construct permutation null distrition
nullDist <- numeric()
for (u in 1:1000){
	shuffle <- sample(1:length(stage), length(stage))
	nullDist <- rbind(nullDist, c(diffGGMml(covML(Y[which(stage[shuffle]==0),]), covML(Y[which(stage[shuffle]==1),]), sum(stage==0), sum(stage==1), lambda=median(optLambdaCweak), Uequal=TRUE)$gamma, 
		                      diffGGMml(covML(Y[which(stage[shuffle]==1),]), covML(Y[which(stage[shuffle]==0),]), sum(stage==1), sum(stage==0), lambda=median(optLambdaCstrong), Uequal=TRUE)$gamma))
}

# calculate p-values
gammaObs[,10] <- c(sum(as.numeric(gammaObs[1,9]) <= c(nullDist[,1], as.numeric(gammaObs[1,9]))) / (nrow(nullDist) + 1), 
				sum(as.numeric(gammaObs[2,9]) <= c(nullDist[,2], as.numeric(gammaObs[2,9]))) / (nrow(nullDist) + 1))



#########################################
# same analysis with Gaussianized data
#########################################

# gaussianize
Y[which(stage==0),] <- apply(Y[which(stage==0),], 2, function(Z){ qnorm(length(Z) * (ecdf(Z)(Z)) / (length(Z) + 1)) })
Y[which(stage==1),] <- apply(Y[which(stage==1),], 2, function(Z){ qnorm(length(Z) * (ecdf(Z)(Z)) / (length(Z) + 1)) })

# determine penalty parameter through LOOCV, repetitively
# for both weakened conditional dependencies in the cancer stage and the reverse (i.e. weaker in the normal)
optLambdaCweak <- optLambdaCstrong <- numeric()
for (u in 1:10){
	shuffle <- sample(1:length(stage), length(stage))
	optLambdaCweak <- c(optLambdaCweak, optPenalty.diffGGMnull(Y[which(stage[shuffle]==0),], Y[which(stage[shuffle]==1),], 10^(-10), 1, lambdaInit=0.1))
	optLambdaCstrong <- c(optLambdaCstrong, optPenalty.diffGGMnull(Y[which(stage[shuffle]==1),], Y[which(stage[shuffle]==0),], 10^(-10), 1, lambdaInit=0.1))
}

# estimate gamma (= observed test statistic)
gammaObs <- rbind(gammaObs, c("gse3933", sum(stage==0), sum(stage==1), ncol(Y), "no", "yes", "n.a.", median(optLambdaCweak), 
			diffGGMml(covML(Y[which(stage==0),]), covML(Y[which(stage==1),]), sum(stage==0), sum(stage==1), lambda=median(optLambdaCweak), Uequal=TRUE)$gamma, NA), 
		   c("gse3933", sum(stage==0), sum(stage==1), ncol(Y), "yes", "yes", "n.a.", median(optLambdaCstrong), 
			diffGGMml(covML(Y[which(stage==1),]), covML(Y[which(stage==0),]), sum(stage==1), sum(stage==0), lambda=median(optLambdaCstrong), Uequal=TRUE)$gamma, NA))


# construct permutation null distrition
nullDist <- numeric()
for (u in 1:1000){
	shuffle <- sample(1:length(stage), length(stage))
	nullDist <- rbind(nullDist, c(diffGGMml(covML(Y[which(stage[shuffle]==0),]), covML(Y[which(stage[shuffle]==1),]), sum(stage==0), sum(stage==1), lambda=median(optLambdaCweak), Uequal=TRUE)$gamma, 
		                      diffGGMml(covML(Y[which(stage[shuffle]==1),]), covML(Y[which(stage[shuffle]==0),]), sum(stage==1), sum(stage==0), lambda=median(optLambdaCstrong), Uequal=TRUE)$gamma))
}

# calculate p-values
gammaObs[3:4,10] <- c(sum(as.numeric(gammaObs[3,9]) <= c(nullDist[,1], as.numeric(gammaObs[3,9]))) / (nrow(nullDist) + 1), 
				sum(as.numeric(gammaObs[4,9]) <= c(nullDist[,2], as.numeric(gammaObs[4,9]))) / (nrow(nullDist) + 1))


write.table(gammaObs, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", file="results_pcaHHgse3933.txt")

