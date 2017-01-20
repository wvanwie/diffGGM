#############################################################################################################################################
#
# READ ME:
#
# This R-script serves to plot the estimation performance of the 'average' differential conditional dependence parameter. 
# Empirical confidence intervals (determined by quantiles) of the regularization paths are plotted for gamma equal to 0.0, 0.1, 0.2, 0.3, 0.4, 0.5.
#
#############################################################################################################################################


# erase history
rm(list=ls())

# set sample size and dimensions
n <- 100
p <- 25

# calculate true condition number of precision matrix used
library(Matrix)
diags <- list(rep(1, p), rep(0.5, p-1), rep(0.25, p-2), rep(0.1, p-3))
Omega <- as.matrix(bandSparse(p, k = -c(0:3), diag = c(diags), symm=TRUE))
evs <- eigen(Omega, only.values=TRUE)$values
cnTrue <- max(evs) / min(evs)

# load results for gamma=0.5
load(paste("simuRes_n", n, "_gamma", 0.5, "_p_", p, "_diffGGMml.Rdata", sep=""))

# set quantiles of confidence for regularization paths
probs <- c(0.10, 0.5, 0.90)
quants <- apply(gammaHat, 2, quantile, probs=probs)

# assess for which choices of lambda yield a reasonable well-conditioned Omega
ids <- 1:100
lambda <- log(lambda[ids])
quants <- quants[,ids]
quantsCN <- apply(cns[,ids], 2, quantile, probs=0.5)
ids <- which(quantsCN < 100) 
lambdaCNtrue <- mean(max(lambda[which(quantsCN > cnTrue)]), min(lambda[which(quantsCN < cnTrue)]), na.rm=TRUE)
ids <- intersect(ids, which(lambda < -0.4))
if (length(ids) > 0){
	lambda <- lambda[ids]
	quants <- quants[,ids]
}

# plot curves for gamma=0.5
plot(quants[2,] ~ lambda, type="l", col="blue", lwd=2, ylab="gamma", xlab="log(lambda)",  main=paste("estimated gamma vs. lambda; n=", n, "; p=", p, sep=""), ylim=c(0, ceiling(10*max(quants))/10))
yCN <- (ceiling(10*max(quants))/10) * c(0:100)/100
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="blue", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="blue", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="blue", lwd=2)

# plot curves for gamma=0.4
load(paste("simuRes_n", n, "_gamma", 0.4, "_p_", p, "_diffGGMml.Rdata", sep=""))
quants <- apply(gammaHat, 2, quantile, probs=probs)
lambda <- log(lambda[ids])
quants <- quants[,ids]
lines(quants[2,] ~ lambda, lty=1, col="purple", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="purple", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="purple", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="purple", lwd=2)

# plot curves for gamma=0.3
load(paste("simuRes_n", n, "_gamma", 0.3, "_p_", p, "_diffGGMml.Rdata", sep=""))
quants <- apply(gammaHat, 2, quantile, probs=probs)
lambda <- log(lambda[ids])
quants <- quants[,ids]
lines(quants[2,] ~ lambda, lty=1, col="red", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="red", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="red", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="red", lwd=2)

# plot curves for gamma=0.2
load(paste("simuRes_n", n, "_gamma", 0.2, "_p_", p, "_diffGGMml.Rdata", sep=""))
quants <- apply(gammaHat, 2, quantile, probs=probs)
lambda <- log(lambda[ids])
quants <- quants[,ids]
lines(quants[2,] ~ lambda, lty=1, col="orange", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="orange", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="orange", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="orange", lwd=2)

# plot curves for gamma=0.1
load(paste("simuRes_n", n, "_gamma", 0.1, "_p_", p, "_diffGGMml.Rdata", sep=""))
quants <- apply(gammaHat, 2, quantile, probs=probs)
lambda <- log(lambda[ids])
quants <- quants[,ids]
lines(quants[2,] ~ lambda, lty=1, col="yellow", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="yellow", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="yellow", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="yellow", lwd=2)

# plot curves for gamma=0.0
load(paste("simuRes_n", n, "_gamma", 0.0, "_p_", p, "_diffGGMml.Rdata", sep=""))
quants <- apply(gammaHat, 2, quantile, probs=probs)
lambda <- log(lambda[ids])
quants <- quants[,ids]
lines(quants[2,] ~ lambda, lty=1, col="wheat1", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="wheat1", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="wheat1", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="wheat1", lwd=2)

# add condition number line
lines(yCN ~ rep(lambdaCNtrue, 101), lwd=0.8, col="black", lty=3)

# add legened
legend("topright", c(expression(paste(gamma, "=0.5", sep="")), expression(paste(gamma, "=0.4", sep="")), expression(paste(gamma, "=0.3", sep="")), expression(paste(gamma, "=0.2", sep="")), expression(paste(gamma, "=0.1", sep="")), expression(paste(gamma, "=0.0", sep=""))), lwd=2, col=c("blue", "purple", "red", "orange", "yellow", "wheat1"))


