## This function is used to find the value lambda
## inside of the double bootstrap proceedure. 
## This function takes as arguments the 
## Bootstrapped sample, and the estimated test 
## statistic sqrt(n)(theta^\# - \theta)

find_lambda <- function(boot_data, boot_thn, thn, t_stat, n_db, boot_resids, n.obs){
  both_sample_df <- matrix(NA, nrow = n_db, ncol = 3)
  for(boot2_indx in 1:n_db){
    ns_indx <- sample(1:n.obs, replace = TRUE)
    nest_boot_data <- boot_data[ns_indx, ]
    nest_est <- find_max_cor_beta(nest_boot_data, find_sd = TRUE)
    nest_bbv <- calc_bbv(nest_boot_data, boot_data, 
                         boot_resids = boot_resids[ns_indx], resids = boot_resids)
    both_sample_df[boot2_indx, ] <- c(nest_est[2], 
                                      sqrt(n.obs) * (nest_est[1] - boot_thn),
                                      nest_bbv)
  }
  both_sample_df[, 1] <- max(both_sample_df[, 1], t_stat)
  found_lambda <- sort_lambda(both_sample_df, sqrt(n.obs) * (boot_thn - thn))
  return(found_lambda)
}




