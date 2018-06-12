################################################
#### This script is written to gives to find run
#### simulations after being called by a script
#### command.  
################################################
################################################
#### Author : Adam Elder
################################################
################################################
#setwd("~/Dropbox/School_UW/2017-2018/Spring 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_1.R")
source("sortlambda.R")
source("Testing.R")
library(MASS)

## This is to create a matrix to give all possible combinations
## of simulation settings. 

rho <- rep(c(0.8, 0.5, 0), each = 30)
sample_size <- rep(rep(c(100, 200), each = 15), times = 3)
dims <- rep(rep(c(10, 50, 100, 150, 200), each = 3), times = 6)
model <- rep(c(1, 2, 3), times = 30)

simulate_df <- data.frame("rho" = rho, "n" = sample_size,
                           "dims" = dims, "model" = model)

## We want to run each setting five times.
simulate_df <- rbind(simulate_df, simulate_df, simulate_df,
                     simulate_df, simulate_df)

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

job.settings <- simulate_df[job.id, ]

cat("Job ID = ", job.id, "simulation Settings : ")
print(job.settings)

## Run a single simulation

results <- simulate(model = job.settings$model,
                    n     = job.settings$n,
                    dim   = job.settings$dims,
                    rho   = job.settings$rho,
                    sims  = 100,
                    bs    = 500,
                    alpha = 0.05)

write.csv(results, paste0("simres/n", job.settings$n, 
                          "mod", job.settings$model,
                          "dim", job.settings$dim,
                          "rho", job.settings$rho, 
                          "iter", floor(job.id/90) + 1, 
                          ".csv"))


