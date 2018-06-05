#################################################
#### This script contains the alternative testing
#### functions to be compared to the ART.
#################################################
#################################################
#### Author : Adam Elder
#################################################

test.once.lr <- function(data, alpha){
  ful.lm <- lm(data[, 1] ~ data[, -1])
  p.val <- anova(ful.lm)$`Pr(>F)`[1]
  return(as.numeric(p.val < alpha))
}

test.once.bc <- function(data, alpha){
  covs <- ncol(data) - 1
  obs <- nrow(data)
  p_vals <- rep(NA, covs)
  X <- as.matrix(data[, -1])
  Y <- data[, 1]
  var_y <- var(Y)
  cov_xy <- cov(X, Y)
  sum_x <- colSums(X)
  x_vars <- apply(data[, -1], 2, var)
  sum_x_2 <- colSums(X^2)
  se_betas <- sqrt(var_y/(sum_x_2 - sum_x^2/(obs)))
  betas <- cov_xy/x_vars
  t_stats <- abs(betas/se_betas)
  for(jj in 2:(covs + 1)){
   p_vals[jj - 1] <- summary(lm(data[, 1] ~ data[, jj]))$coefficients[2, 4] 
  }
  signif <- as.numeric(p_vals <= alpha/covs)
  return(as.numeric(sum(signif) > 0))
}

test.once.cpb <- function(data, boot_samples, alpha){
  covs <- ncol(data) - 1
  n_obs <- nrow(data)
  boot_distr <- rep(NA, boot_samples)
  theta_est <- find_max_cor_beta(data, find_sd = FALSE)
  sample_indx <- 1:n_obs
  for(boot_indx in 1:boot_samples){
    boot_est <- find_max_cor_beta(data[sample(sample_indx, replace = TRUE), ], find_sd = FALSE)
    boot_distr[boot_indx] <- sqrt(n_obs) * (boot_est - theta_est) 
  }
  quantile <- mean(as.numeric(sqrt(n_obs) * theta_est > boot_distr))
  return(as.numeric(quantile < (alpha/2) | quantile > (1 - alpha/2)))
}
