##############################################################################
#### This Script defines the data generating function
#### Author : Adam Elder
##############################################################################

library(MASS)
make.data <- function(sample_size, dim, rho, model = 1){
  x_cov <- matrix(rho, nrow = dim, ncol = dim)
  diag(x_cov) <- 1
  X <- mvrnorm(n = sample_size, mu = rep(0, dim), Sigma = x_cov)
  epsilon <- rnorm(sample_size)
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
  return(data)
}
