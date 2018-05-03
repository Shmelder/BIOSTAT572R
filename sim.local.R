##################################################
#### Simulation Study of Local Alternatives
#### Used to create slides, not part of the
#### paper itself.
#### 
#### Author : Adam Elder
#################################################

library(MASS)

## Function used to generate simulated data.  The
## function generates ss observations from data 
## with covs covariates, and a local parameter of
## dimension less than covs.

gen.data <- function(ss, local, covs){
  x <- mvrnorm(n = ss, mu = rep(0, covs), Sigma = diag(covs) )
  true.beta <- c(local/sqrt(ss), rep(0, covs - length(local)))
  y <- rnorm(ss, sd = 1) + x %*% true.beta
  return(data.frame("y" = y, "x" = x))
}

gen_est <- function(ss, local, covs){
  obs_data <- gen.data(ss, local, covs)
  return(find_max_cor_param(obs_data))
}

get.est <- function(ss, local, covs){
  obs.data <- gen.data(ss, local, covs)
  reg <- lm(y ~., data = obs.data)
  return(which.max(coefficients(reg)[-1]))
}

find.distr <- function(sims, ss, local, covs, fixed = TRUE){
  sim_ests <- rep(NA, sims)
  if (fixed == TRUE){ lp = local * sqrt(ss) }else{lp = local}
  for(sim_indx in 1:sims){
    sim_ests[sim_indx] <- sqrt(ss) * (gen_est(ss = ss, local = lp, covs = covs) - max(abs(lp/sqrt(ss))))
  }
  return(sim_ests)
}

est_distr_boot <- function(data, b_s){
  obs <- nrow(data)
  beta_hat <- find_max_cor_param(data)
  boot.distr <- rep(NA, b_s)
  for(boot.indx in 1:b_s){
    boot_data <- data[sample(1:obs, replace = TRUE), ]
    boot_beta_hat <- find_max_cor_param(boot_data)
    boot.distr[boot.indx] <- sqrt(obs) * (boot_beta_hat - beta_hat)
  }
  return(boot.distr)
}

make_some_plots <- function(sample_size, sims, fixed = TRUE, pertrb){
  if(fixed == TRUE){lp = sqrt(sample_size) * pertrb}else{lp = pertrb}
  sample_distr <- find.distr(sims, sample_size, lp, 10, fixed = fixed)
  bootst_distr <- est_distr_boot(gen.data(sample_size, lp, 10), sims)
  
  boot_cdf <- ecdf(bootst_distr)
  samp_cdf <- ecdf(sample_distr)
  
  est.range <- range(c(sample_distr, bootst_distr))
  
  dat.2 <- data.frame(boot = boot_cdf(seq(from = est.range[1], to = est.range[2], length.out = 10000)),
                      samp = samp_cdf(seq(from = est.range[1], to = est.range[2], length.out = 10000)))
  
  dat <- data.frame("val" = c(sample_distr, bootst_distr), 
                    "Generation Method" = rep(c("Sample", "Bootstrap"), each = sims))
  
  A <- ggplot(dat.2, aes(boot, samp)) + geom_abline(slope = 1, intercept = 0, col = "red") + 
    geom_line(colour = "blue", size =1) + 
    labs(x = "\n Bootstrap Distribution Quantile", y = "Sampling Distribution Quantile \n") + 
    theme_minimal() + theme(text = element_text(size=20))
  
  B <- ggplot(dat, aes(x = val, fill = Generation.Method)) + 
    geom_density(alpha = .5, position="identity") + 
    scale_fill_discrete(name = "Generation Method : ") + 
    labs(x = "\n Beta", y = "\n Density") + 
    theme_minimal() + theme(legend.position="top", text = element_text(size=20))
  browser()
  cowplot::plot_grid(A, B)
}

null_sample_distr <- sample_distr <- find.distr(5000, 4000, 0, 10, fixed = TRUE)
make_some_plots_2 <- function(sample_size, sims, fixed = TRUE, pertrb){
  if(fixed == TRUE){lp = sqrt(sample_size) * pertrb}else{lp = pertrb}
  sample_distr <- find.distr(sims, sample_size, lp, 10, fixed = fixed)
  
  dat <- data.frame("val" = c(null_sample_distr, sample_distr), 
                    "Generation Method" = rep(c( "Local Sample", "Null Sample"), each = sims))
  dat <- subset(dat, val <= 5 & -7.5 <= val)
  A <- ggplot(dat, aes(x = val, fill = Generation.Method)) + 
    scale_fill_discrete(name = "Generation Method : ",
                        labels = c("Null Sampling Distribution", "Local Sampling Distribution")) + 
    geom_density(alpha = .5, position="identity") + 
    labs(x = "\n Beta", y = "\n Density") + xlim(c(-7.5 , 5)) + ylim(c(0.0, 0.45)) + 
    theme_minimal() + theme(legend.position="top", text = element_text(size=20))
  A
}

make_some_plots_2(4000, 5000, fixed = FALSE, 7)
make_some_plots_2(4000, 5000, fixed = FALSE, 3)
make_some_plots_2(4000, 5000, fixed = FALSE, 1)
make_some_plots_2(4000, 5000, fixed = FALSE, 0.5)
make_some_plots_2(4000, 5000, fixed = FALSE, 0.1)
