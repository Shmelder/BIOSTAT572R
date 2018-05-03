##############################################################################
#### This Script defines the ART proceedure
#### Author : Adam Elder
#### Requires get_max_beta.R
##############################################################################
library(MASS)
                                     
art_find_distr <- function(data, boot.s, lambda = "db", n_db){
  n.obs  <- nrow(data)
  n.covs <- ncol(data) - 1 
  sample_summary <- find_max_cor_beta(data, find_sd = TRUE)
  thnh <- sample_summary[1]
  t_stat <- sample_summary[2] 
  resids <- sample_summary[-c(1, 2)]
  boot.distr <- rep(NA, boot.s)
  for(boot1_indx in 1:boot.s){
    rs_indx <- sample(1:n.obs, replace =TRUE)
    boot_data <- data[rs_indx, ]
    boot_summary <- find_max_cor_beta(boot_data, find_sd = TRUE)
    if (lambda == "db"){
      lambda_n <- find_lambda(boot_data,
                              boot_summary[1], 
                              thnh, 
                              t_stat = t_stat, 
                              n_db = n_db, 
                              boot_resids = boot_summary[-c(1, 2)], 
                              n.obs = n.obs)
    }else{
      lambda_n <- lambda
    }
    if (lambda_n <= sample_summary[2] & lambda_n <= t_stat){
      boot.distr[boot1_indx] <- sqrt(n.obs) * (boot_summary[1] - thnh)
    }else{
      boot.distr[boot1_indx] <- calc_bbv(boot_data, data, 
                                         boot_resids = resids[rs_indx], 
                                         resids = resids)
    }
  }
  return(boot.distr)
}

art.test <- function(data, boot.s, lambda){
  n.obs <- nrow(data)
  test.stat <- sqrt(n.obs) * est_max_t(data, keep.resid = FALSE)$theta_n
  limit.distr <- art_find_distr_tunegiven(data, boot.s, lambda)
  z.crit <- quantile(limit.distr, probs = c(0.025, 0.975))
  print(test.stat)
  if (test.stat < z.crit[1] | test.stat > z.crit[2]){
    return(1)
  }else{
    return(0)
  }
}


