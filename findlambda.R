## This function is used to find the value lambda
## inside of the double bootstrap proceedure. 
## This function takes as arguments the 
## Bootstrapped sample, and the estimated test 
## statistic sqrt(n)(theta^\# - \theta)

find_lambda <- function(boot_data, boot_thn, thn, t_stat, n_db, boot_resids, n.obs){
  both_sample_df <- matrix(NA, nrow = n_db, ncol = 3)
  ## The above data frame saves the values for both
  ## the centered percentile bootstrap and V_n^#
  ## and the associated T-statistic for n_db
  ## nested bootstrap draws. 
  for(boot2_indx in 1:n_db){
    ns_indx <- sample(1:n.obs, replace = TRUE)
    #Taking nested bootstrap draws, and recording the CBP and 
    # V_n^# values for each bootstrap draw
    nest_boot_data <- boot_data[ns_indx, ]
    nest_est <- find_max_cor_beta(nest_boot_data, find_sd = TRUE)
    nest_bbv <- calc_bbv(nest_boot_data, boot_data, 
                         boot_resids = boot_resids[ns_indx], resids = boot_resids)
    both_sample_df[boot2_indx, ] <- c(nest_est[2], 
                                      sqrt(n.obs) * (nest_est[1] - boot_thn),
                                      nest_bbv)
  }
  # Now, we pass this matrix to another function. The function (sort_lambda) 
  # finds the value of lambda_n for which 95 percent of all values fall under, 
  # or above the test statistic.  Increasing lambda will cause us to use 
  # V_n^* values more frequently, and thus cause our bootstraped distribution
  # to be closer to zero.  If no lambda_n value can be found to satisfy the
  # required condition, the lambda_n closest giving a rejection rate closest
  # to 0.05 is chosen.  
  both_sample_df[, 1] <- max(both_sample_df[, 1], t_stat)
  found_lambda <- sort_lambda(both_sample_df, sqrt(n.obs) * (boot_thn - thn))
  return(found_lambda)
}




