################################################
#### This script is written to gives functions
#### which facilitate simulation studies based
#### on the ART test (specified in ART.R),
#### specifically, this script simulates
#### draws from the true sampling distribution
#### for local alternatives. 
################################################
library(MASS)

get_sample <- function(sims, dim, b){
  norm_var_mat <- matrix(0.5, nrow = dim, ncol = dim)
  diag(norm_var_mat) <- 1
  norm_data <- mvrnorm(n = sims, mu = rep(0, dim), Sigma = norm_var_mat)
  distr <- apply(norm_data, 1, FUN = one_draw, b = b, dim = dim)
  return(distr + b)
}

one_draw <- function(vec, b, dim){
  cov_b_0 <- c(1, rep(0.5, dim - 1))
  K <- which.max((vec + b * cov_b_0)^2)
  draw <- vec[K] + (0.5 + 0.5 * as.numeric(K == 1) - 1) * b
  return(draw)
}



