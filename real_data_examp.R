setwd("~/Dropbox/Spring 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
source("ART.R")
source("bbV.R")
source("findlambda.R")
source("get_max_beta.R")
source("Sim_Study.R")
source("Simulate_1.R")
source("sortlambda.R")
source("Testing.R")
library(MASS)


real_data <- read.csv("~/Desktop/real_data_examp.csv")
sub_data <- real_data[, -c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
