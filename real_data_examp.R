setwd("~/Dropbox/Spring 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_1.R")
source("sortlambda.R")
source("Testing.R")
library(MASS)


real_data <- read.csv("~/Desktop/real_data_examp.csv")
sub_data <- as.matrix(real_data[, -c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12)])

## This function estimates the limiting distribution of
## the sample statistic using a niave bootstrap.
## This procedure is known to be anti-conservative.

est_distr_cpb <- function(data, boot_samples){
  covs <- ncol(data) - 1
  n_obs <- nrow(data)
  boot_distr <- rep(NA, boot_samples)
  theta_est <- find_max_cor_beta(data, find_sd = FALSE)
  sample_indx <- 1:n_obs
  for(boot_indx in 1:boot_samples){
    boot_est <- find_max_cor_beta(data[sample(sample_indx, replace = TRUE), ], find_sd = FALSE)
    boot_distr[boot_indx] <- sqrt(n_obs) * (boot_est - theta_est) 
  }
  return(boot_distr)
}

test_stat <- sqrt(nrow(sub_data)) * find_max_cor_beta(sub_data)
ART_distr <- art_find_distr(sub_data, 1000, lambda = "db", 1000, alpha = FALSE)
naive_dist <- est_distr_cpb(sub_data, 3000)

## Creating Plots

par(mfrow = c(1, 1), mai = c(1.02, 0.82, 0.2, 0.42))
data_range <- c(min(c(ART_distr, naive_dist)),
                max(c(ART_distr, naive_dist)))
p1 <- hist(ART_distr, breaks = seq(data_range[1], data_range[2], length.out = 20))
p2 <- hist(naive_dist, breaks = seq(data_range[1], data_range[2], length.out = 20))

plot(p1, col = rgb(1, 0, 0, 1/4), xlim = data_range, freq = FALSE,
     border = FALSE, xlab = expression("Bootstrap Estimates of " * theta[0]), main = "")
plot(p2, col = rgb(0, 0, 1, 1/4), xlim = data_range, freq = FALSE, add = T, 
     border = FALSE, xlab = expression("Bootstrap Estimates of " * theta[0]), main = "")
abline(v = test_stat, lwd = 3)
legend(x = 450, y = 0.008, pch = c(15, 15, 15), 
       col = c("black", rgb(1, 0, 0, 1/4), rgb(0, 0, 1, 1/4)),
       legend = c("Test Statistic", "ART", "Naive Bootstrap"), bty = "n",
       cex = 1.2, pt.cex = 4)

## Finding p-values
find_p_bc <- function(data){
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
  min_indx <- which.min(p_vals) + 1
  print(summary(lm(data[, 1] ~ data[, min_indx])))
  return(min(p_vals) * covs)
}

find_p_bc(sub_data)
quantile(ART_distr, 0.025)
mean(test_stat > ART_distr)
mean(test_stat > naive_dist)
