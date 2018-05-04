################################################
#### This script is written to gives to find run
#### simulations after being called by a script
#### command.  
################################################
################################################
#### Author : Adam Elder
################################################
################################################

source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_1.R")
source("sortlambda.R")
source("Testing.R")
library(MASS)

rho <- rep(c(0.8, 0.5, 0), each = 30)
sample_size <- rep(rep(c(100, 200), each = 15), times = 3)
dims <- rep(rep(c(10, 50, 100, 150, 200), each = 3), times = 6)
model <- rep(c(1, 2, 3), times = 30)

simulate_df <- data.frame("rho" = rho, "n" = sample_size,
                           "dims" = dims, "model" = model) 

syst.envr <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if(!is.na(syst.envr)){
  job.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}else{
  cat("No Job ID found, choosing the answer to life the universe and everything \n")
  job.id <- 42
}

job.settings <- simulate_df[job.id, ]

cat("Job ID = ", job.id, "simulation Settings : ")
print(job.settings)

results <- simulate(model = job.settings$model,
                    n     = job.settings$n,
                    dim   = job.settings$dims,
                    rho   = job.settings$rho,
                    sims  = 1,
                    alpha = 0.05)



