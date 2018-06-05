## This function is used to find the value lambda
## inside of the double bootstrap proceedure. 
## This function takes as arguments the S 
## bootstrapped observations for both a CPB and 
## the estimate V_n for each sample. The funciton
## also takes the absolute value of the t-test 
## value. These arguments are given as a matrix
## with the first column as abs(t-test), second
## column is V_n values, and third column is CBP.

sort_lambda <- function(observations, test_stat){
  if (length(unique(observations[, 1])) == 1){
    return(observations[1, 1])
  }else{
    ord.obs <- observations[order(observations[, 1]), ]
    ord.vnb <- as.numeric(ord.obs[, 2] > test_stat) 
    ord.cpb <- as.numeric(ord.obs[, 3] > test_stat)
    n_obs <- length(ord.cpb)
    ord.dif <- (ord.vnb - ord.cpb)/n_obs
    num.abv <- mean(ord.cpb)
    step <- 1
    if(num.abv < 0.05){
      while(num.abv < 0.05 & step < n_obs){
        num.abv <- num.abv + ord.dif[step]
        step <- step + 1
      }
    }else{
      while(num.abv > 0.95 & step < n_obs){
        num.abv <- num.abv + ord.dif[step]
        step <- step + 1
      }
    }
    if(step == n_obs){
      ## If the precise quantile where we get 95% acceptance rate cannot
      ## be found, we find the lambda_n that gives us the value closest 
      ## to a 95% acceptence rate.  
      quants <- cumsum(c(mean(ord.cpb), ord.dif))
      quant.diffs <- c(abs(quants - 0.95), abs(quants - 0.05))
      browser()
      step = which.min(quant.diffs) %% n_obs 
    }
    return(ord.obs[step, 1]) 
  }
}
















