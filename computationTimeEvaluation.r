# erase history
rm(list=ls())

# activate libraries
library(mvtnorm)
library(Matrix)
library(diffGGM)
library(microbenchmark)

# number of samples n1=n2=10, 50, 100, 500
n1 <- n2 <- 50

# penalty parameters
lambda <- 0.1

# set number of genes
p <- 75

# gamma parameter
gamma <- 0.5

# generate precision and covariance matrices
diags <- list(rep(1, p), rep(0.5, p-1), rep(0.25, p-2), rep(0.1, p-3))
Omega <- as.matrix(bandSparse(p, k = -c(0:3), diag = c(diags), symm=TRUE))
Sigma1 <- solve(Omega[1:p, 1:p])
Sigma2 <- solve((1-gamma) * Omega[1:p, 1:p] + gamma * diag(p))

# sample data
Y1 <- rmvnorm(n1, sigma=Sigma1)
Y2 <- rmvnorm(n2, sigma=Sigma2)
S1 <- t(Y1) %*% Y1 / n1
S2 <- t(Y2) %*% Y2 / n2

# perform timing (reported in milliseconds)
mean(microbenchmark(diffGGMml(S1, S2, n1, n2, lambda))$time * 10^(-9))



