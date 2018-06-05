################################################
#### This script is written to gives to find run
#### simulations after being called by a script
#### command.  
################################################
################################################
#### Author : Adam Elder
################################################
################################################
setwd("~/Dropbox/Spring 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
library(MASS)
library(expm)
source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_1.R")
source("sortlambda.R")
source("Testing.R")

## This function takes in the simulation settings and
## runs a the ART test based on these settings.  sims
## gives the number of times the test will be run.

simulate_test <- function(b0, n, dim, sims, alpha, lamb){
  x_cov <- matrix(0.5, nrow = dim, ncol = dim)
  diag(x_cov) <- 1
  sqrt_mat <- sqrtm(x_cov)
  start_time <- Sys.time()
  test_res <- rep(NA, nrow = sims)
  for(sim_indx in 1:sims){
    sim_data <- make_data_loc(ss = n, mat = sqrt_mat, b = b0)
    test_res[sim_indx] <- art_find_distr(sim_data, 500, lambda = lamb, 500, alpha = alpha)
    if(sim_indx%%100 == 0){
      time.dif <- as.numeric(Sys.time() - start_time)/(60**2)
      cat(sim_indx, " Simulations have been completed taking ", time.dif, "hours \n")
    }
  }
  print(Sys.time() - start_time)
  print(as.numeric(Sys.time() - start_time))
  rej.rate <- mean(test_res)
  names(rej.rate) <- c("ART")
  return(rej.rate)
}

test_1 <- c(0, 10000, 1000, 10, 0.05)
test_2 <- c(3, 5000, 10, 10, 0.05)

## The following is used so the rscript can be run
## inside of the cluster, and used different simulation
## settings.  

syst.envr <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if(!is.na(syst.envr)){
  job.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}else{
  cat("No Job ID found, choosing the answer to life the universe and everything \n")
  job.id <- 42
}

if (job.id %% 2 == 0){
  job.settings <- test_1
  lamb.settings <- 4.3
}else{
  job.settings <- test_2
  lamb.settings <- "db"
}

cat("Job ID = ", job.id, "simulation Settings : ")
print(job.settings)

results <- simulate_test(b0 = job.settings[1],
                         n     = job.settings[2],
                         dim   = job.settings[3],
                         sims  = job.settings[4],
                         alpha = 0.05,
                         lamb = lamb.settings)

write.csv(results, paste0("simres4/test",job.id %% 2, "job_", floor(job.id/2) + 1, ".csv"))


