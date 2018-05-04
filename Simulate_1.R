################################################
#### This script is written to gives functions
#### which facilitate simulation studies based
#### on the ART test (specified in ART.R) and
#### other competing proceedures (specified in 
#### Testing.R). This script specifically is used
#### to recreate figures 1, 2, and 3. 
################################################
################################################

simulate <- function(model, n, dim, rho, sims, alpha){
  start_time <- Sys.time()
  test_res <- matrix(NA, nrow = sims, ncol = 4)
  for(sim_indx in 1:sims){
    sim_data <- make_data(ss = n, dim = dim, rho = rho, model = model)
    lr_test <-  test.once.lr(sim_data, alpha = alpha)
    bc_test <-  test.once.bc(sim_data, alpha = alpha)
    cpb_test <- test.once.cpb(sim_data, 1000, alpha = alpha)
    art_test <- art_find_distr(sim_data, 1000, lambda = "db", 1000, alpha = alpha)
    test_res[sim_indx, ] <- c(lr_test, bc_test, cpb_test, art_test)
    if(sim_indx%%100 == 0){
      time.dif <- as.numeric(Sys.time() - start_time)/(60**2)
      cat(sim_indx, " Simulations have been completed taking ", time.dif, "hours \n")
      }
  }
  rej.rate <- c(apply(sim_indx, 2, mean), time.dif)
  names(rej.rate) <- c("LR", "BC", "CPB", "ART", "Time")
  return(rej.rate)
}


