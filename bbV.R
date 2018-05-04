## This function is used to estimate the V^#_n
## when t_n < lambda. This function will use both
## the observed data and the bootstrapped data.

calc_bbv <- function(boot_data, data, boot_resids, resids){
  covs <- ncol(data) - 1
  ss <- length(resids)
  mean_cent <- rep(colMeans(boot_data), each = ss) 
  c_data      <- data - mean_cent
  c_data_boot <- boot_data - mean_cent
  z_nk <- sqrt(ss) * (crossprod(c_data_boot, boot_resids) - crossprod(c_data, resids))/ss
  var_k <- colSums(c_data_boot^2)/(ss - 1)
  max_k <- which.max(z_nk^2/var_k)
  max_z <- z_nk[max_k]/var_k[max_k]
  return(max_z)
}




