################################################
#### This script is written to gives to find run
#### simulations after being called by a script
#### command.  
################################################
################################################
#### Author : Adam Elder
################################################
################################################
#setwd("~/Dropbox/Spring 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_2.R")
source("find_limit_distr.R")
source("sortlambda.R")
source("Testing.R")
library(MASS)

## This is to create a matrix to give all possible combinations
## of simulation settings. 

local_param <- rep(seq(0, 5, 0.5), times = 6)
sample_size <- c(rep(5000, 22), rep(10000, 44))
dims <- c(rep(c(10, 50, 1000, 1000, 1000, 1000), each = 11))
sims <- c(rep(10, 22), rep(20, 44))/2
db <- rep(c("db", "db", "4.3", "6.1", "6.8", "8.6"), each = 11)

simulate_df <- data.frame("b_0" = local_param, "n" = sample_size,
                           "dims" = dims, "sims" = sims, "db" = as.character(db))

simulate_df$db <- db
simulate_df <- rbind(simulate_df, simulate_df)

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

cat("Job ID = ", job.id, "simulation Settings : \n")
print(job.settings)
lam_setting <- as.character(job.settings$db)
if (job.settings$db != "db"){
  lam_setting <- as.numeric(lam_setting)
}

## Run a single simulation

results <- simulate_2(b = job.settings$b_0,
                    dim   = job.settings$dims,
                    sims  = job.settings$sims,
                    lim_sims = 500000,
                    ss = job.settings$n,
                    lam = lam_setting,
                    bs_sims = 1000)

## Write results to a .csv file.

write.csv(results, paste0("simres2/", "sim_dim_", job.settings$dims,
                          "ss_", job.settings$n,
                          "b0_", job.settings$b_0, 
                          "iter_", floor((job.id - 1)/66) + 1, 
                          "db", job.settings$db,
                          ".csv"))


