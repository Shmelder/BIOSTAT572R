## This function is used to find the maximally correlated 
## coefficient.  While this could be computed simply using
## built in functions, this function computes what we need
## slightly faster.   

find_max_cor_beta <- function(data, find_sd = FALSE){
  covs <- ncol(data) - 1
  ss <- nrow(data)
  s.data <- data - rep(colMeans(data), each = ss)
  s.data <- s.data / rep(sqrt(colSums(s.data^2)), each = ss)
  cors <- crossprod(s.data[, -1], s.data[, 1])
  max.k <- which.max(abs(cors)) + 1
  if(find_sd == TRUE){
    y.vec <- data[, 1]
    TX <- rbind(1, data[, max.k])
    X <- t(TX) 
    XTXinv <- solve(TX %*% X)
    betas <- XTXinv %*% TX %*% y.vec
    resids <- y.vec - as.numeric(X %*% betas)
    rse <- sum(resids^2)/(ss - 2)
    est <- c(betas[2], betas[2]/sqrt(XTXinv[2, 2] * rse), resids)
  }else{
    est <- cors[max.k - 1] / sd(data[, max.k])
  }
  return(est)
}
