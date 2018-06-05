##############################################################################
#### This Script defines the data generating function
#### Author : Adam Elder
##############################################################################
library(MASS)
## Function used to generate simulated data. The
## function generates ss observations from one of
## the three data generating mechanism (specified)
## by model, dim covariates, and rho correlation
## between each x for the first three figures.  

make_data <- function(ss, dim, rho, model = 1, b){
  x_cov <- matrix(rho, nrow = dim, ncol = dim)
  diag(x_cov) <- 1
  X <- mvrnorm(n = ss, mu = rep(0, dim), Sigma = x_cov)
  epsilon <- rnorm(ss)
  if (model == 1){
    data <- cbind(epsilon, X)
  }
  if (model == 2){
    data <- cbind(X[, 1]/4 + epsilon, X)
  }
  if (model == 3){
    beta <- c(rep(c(0.15, -0.1), each = 5), rep(0, dim - 10))
    data <- cbind(X %*% beta + epsilon, X)
  }
  if (model == 4){
    browser()
    theta_n <- c(b, rep(0, dim - 1))/sqrt(ss)
    data <- cbind(X %*% theta_n + epsilon, X)
  }
  return(data)
}

make_data_loc <- function(ss, mat, b){
  dim <- nrow(mat)
  eps <- matrix(rnorm(n = ss * dim), nrow = ss)
  X <- eps %*% mat
  epsilon <- rnorm(ss)
  theta_n <- c(b, rep(0, dim - 1))/sqrt(ss)
  data <- cbind(X %*% theta_n + epsilon, X)
  return(data)
}

cov_mat <- matrix(0.5, 10, 10)
diag(cov_mat) <- 1
cov_sqrt <- sqrtm(cov_mat)
for(i in 1:1000){
  max_betas[i] <- find_max_cor_beta(make_data_loc(1000, cov_sqrt, 0), find_sd = FALSE)
}

