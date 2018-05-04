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
  return(ord.obs[step, 1])
}
















