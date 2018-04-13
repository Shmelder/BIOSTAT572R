##############################################################################
#### This Script defines the ART proceedure
#### Author : Adam Elder
##############################################################################

library(MASS)

est_max_t <- function(data, keep.resid = FALSE){
  n.obs <- nrow(data)
  marginal.cor <- apply(data[, -1], 2, cor, y = data[, 1])
  max_cor_indx <- which.max(abs(marginal.cor))
  lm.max <- lm(data[, 1] ~ data[, max_cor_indx + 1])
  lm.coef <- as.numeric(lm.max$coefficients[2])
  lm.resids <- lm.max$residuals
  std.err.max <- sqrt(mean((lm.resids)^2) * n.obs/(n.obs - 1))
  base.t.stat <- lm.coef/std.err.max
  if(keep.resid == FALSE){
    return(list("tstat" = base.t.stat, "theta_n" = lm.coef))
  }else{
    return(list("tstat" = base.t.stat, "theta_n" = lm.coef, 
                "resids" = lm.resids))
    }
}
                                     
art_find_distr_tunegiven <- function(data, boot.s, lambda){
  n.obs  <- nrow(data)
  n.covs <- ncol(data) - 1 
  sample.summary <- est_max_t(data, keep.resid = TRUE)
  easy <- (abs(sample.summary$tstat) > lambda) # Tells us if the indicator for A^*_n is always on or not
  boot.distr <- rep(0, boot.s)
  if (easy == TRUE){
    for(boot.indx in 1:boot.s){
      boot.data <- data[sample(1:n.obs, replace = TRUE), ]
      boot.max.theta <- est_max_t(boot.data, keep.resid = FALSE)$theta_n
      boot.distr[boot.indx] <- sqrt(n.obs) * (boot.max.theta - sample.summary$theta_n)
    }
  }else{
    for(boot.indx in 1:boot.s){
      boot.sample <- sample(1:n.obs, replace = TRUE)
      boot.data <- data[boot.sample, ]
      boot.summary <- est_max_t(boot.data, keep.resid = FALSE)
      if (abs(boot.summary$tstat) > lambda){
        boot.distr[boot.indx] <- sqrt(n.obs) * (boot.summary$theta_n - sample.summary$theta_n)
      }else{
        # Here we are calculating Z^*_{n, k} to calculate V^*_n(0)
        gn.ep.xk <- pn.xk <- rep(0, n.covs)
        gn.ep <- sqrt(n.obs) * mean(sample.summary$resids[boot.sample])
        
        for(cov.indx in 1:n.covs){
          gn.ep.xk[cov.indx] <- sqrt(n.obs) * (mean(boot.data[, cov.indx] * sample.summary$resids[boot.sample]) - 
                                                 mean(   data[, cov.indx] * sample.summary$resids))
          pn.xk[cov.indx] <- mean(boot.data[, cov.indx])
        }
        znk <- gn.ep.xk - gn.ep * pn.xk
        boot.vars <- apply(boot.data[, -1], 2, var)
        vn.zero <- znk/boot.vars
        vn.zero.max <- which.max(abs(vn.zero))
        boot.distr[boot.indx] <- vn.zero[vn.zero.max]
      }
    }
  }
  hist(boot.distr)
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


