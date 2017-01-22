#############################################################################################################################################
#
# READ ME:
#
# This R-script serves to plot the power curves of for testing  whether the 'average' differential conditional dependence parameter equals zero or not. 
# Empirical confidence intervals (determined by quantiles) of the power curves are plotted for H_a with gamma equal to 0.0, 0.1, 0.2, 0.3, 0.4, 0.5.
#
#############################################################################################################################################

# erase history
rm(list=ls())

# set number of genes
n <- 25
gammaAll <- seq(0.1, 0.5, 0.1)
probs <- c(0.10, 0.5, 0.90)
p <- 25

op <- par(pty="s")
load(paste("power_n", n, "_p" , p, "_gamma0.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
plot(quants[2,] ~ lambda, type="l", xlim=c(0, 0.5), ylim=c(0,1), xlab=expression(alpha), ylab=expression(paste("power = %(", hat(gamma) > c[alpha], ")", sep="")), lwd=2, col="wheat1", lty=1, main=expression(paste("correct ", H[a], "; power vs ", alpha, "; n=", 25, "; p=", 25, sep="")))
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="wheat1", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="wheat1", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="wheat1", lwd=2)


load(paste("power_n", n, "_p" , p, "_gamma0.5.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
lines(quants[2,] ~ lambda, lty=1, col="blue", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="blue", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="blue", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="blue", lwd=2)

load(paste("power_n", n, "_p" , p, "_gamma0.4.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
lines(quants[2,] ~ lambda, lty=1, col="purple", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="purple", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="purple", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="purple", lwd=2)

load(paste("power_n", n, "_p" , p, "_gamma0.3.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
lines(quants[2,] ~ lambda, lty=1, col="red", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="red", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="red", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="red", lwd=2)

load(paste("power_n", n, "_p" , p, "_gamma0.2.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
lines(quants[2,] ~ lambda, lty=1, col="orange", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="orange", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="orange", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="orange", lwd=2)

load(paste("power_n", n, "_p" , p, "_gamma0.1.Rdata", sep=""))
quants <- apply(matrix(unlist(power), ncol=length(power)), 1, quantile, probs=probs)
lambda <- c(1:50)/100
lines(quants[2,] ~ lambda, lty=1, col="yellow", lwd=2)
revOrder <- sort(1:length(lambda), decreasing=TRUE)
polygon(c(lambda, lambda[revOrder]), c(quants[1,], quants[3,revOrder]), col="yellow", density=50, border=NA)
lines(quants[1,] ~ lambda, lty=2, col="yellow", lwd=2)
lines(quants[3,] ~ lambda, lty=2, col="yellow", lwd=2)


legend("bottomright", c(expression(paste(gamma, "=0.5", sep="")), expression(paste(gamma, "=0.4", sep="")), expression(paste(gamma, "=0.3", sep="")), expression(paste(gamma, "=0.2", sep="")), expression(paste(gamma, "=0.1", sep="")), expression(paste(gamma, "=0.0", sep=""))), col=c("blue", "purple", "red", "orange", "yellow", "wheat1"), lwd=2, lty=c(1,2,3,4,5, 6))

par(op)




