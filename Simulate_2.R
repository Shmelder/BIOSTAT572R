################################################
#### This script is written to gives functions
#### which facilitate simulation studies based
#### on the ART test (specified in ART.R) and
#### other competing proceedures (specified in 
#### Testing.R). This script specifically is used
#### to recreate figures 4 and 5 
################################################
################################################

simulate_2 <- function(b, dim, sims, lim_sims, ss, lam, bs_sims){
  start_time <- Sys.time()
  limiting_distr <- get_sample(lim_sims, dim, b)
  test_res <- rep(NA, sims)
  for(sim_indx in 1:sims){
    sim_data <- make_data(ss = ss, dim = dim, rho = 0.5, model = 4, b = b)
    art_distr <- art_find_distr(data = sim_data, bs_sims, lambda = lam, bs_sims, alpha = FALSE)
    quants <- quantile(art_distr, c(0.025, 0.975))
    rej_rate <- (sum(limiting_distr < quants[1]) + sum(limiting_distr > quants[2]))/(lim_sims)
    test_res[sim_indx] <- rej_rate
    print(Sys.time() - start_time)
    cat(sim_indx, "Simulations have been completed \n")
  }
  print(Sys.time() - start_time)
  time.dif <- as.numeric(Sys.time() - start_time)/(60**2)
  test_res <- c(test_res, time.dif)
  names(test_res) <- c(1:sims, "Time")
  return(test_res)
}


