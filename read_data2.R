#########################################################
#### Author : Adam Elder
#### This script contains the code used to read in 
#### simulation results for the fourth and fifth figure.
#########################################################
# setwd("/Users/adamelder/Dropbox/Spring\ 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
library(expm)
library(ggplot2)

read_one <- function(params, file_name){
  est_rej <- read.csv(file_name)[,2]
  est_rej <- est_rej[ -length(est_rej)]
  num_reps <- length(est_rej)
  param_mat <- data.frame("b0" = rep(params[1], num_reps),
                          "dims" = rep(params[2], num_reps),
                          "db" = rep(params[3], num_reps))
  rej_df <- as.data.frame(cbind(param_mat, est_rej))
  colnames(rej_df) <- c("b0", "dims", "db", "rejr")
  return(rej_df)
}

read_params <- function(ss, b, dim, db){
  file_name_start <- paste0("simres2/sim_dim_", dim, "ss_", ss, "b0_", b, "iter_")
  file_name_end   <- paste0("lambda", db)
  for(i in 1:2){
    iter_i <- paste0(file_name_start, i, file_name_end, ".csv")
    if(file.exists(iter_i)){
      res_i <- read_one(c(b, dim, db), iter_i)
      if (exists("results")){
        results <- rbind(results, res_i)
      }else{
        results <- res_i
      }
    }
  }
  return(results)
}

local_param <- rep(seq(0, 5, 0.5), times = 6)
sample_size <- c(rep(5000, 22), rep(10000, 44))
dims <- c(rep(c(10, 50, 1000, 1000, 1000, 1000), each = 11))
sims <- c(rep(10, 22), rep(20, 44))/2
db <- rep(c("db", "db", "4.3", "6.1", "6.8", "8.6"), each = 11)

simulate_df <- data.frame("b_0" = local_param, "n" = sample_size,
                          "dims" = dims, "sims" = sims, "db" = as.character(db))

simulate_df$db <- db
n_ests <- sum(simulate_df$sims) * 2
mast_indx <- cumsum(simulate_df$sims) * 2
master_mat <- as.data.frame(matrix(0, nrow = n_ests, ncol = 4))
settings <- simulate_df[1, ]
master_mat <- read_params(ss = settings$sims,
                          b = settings$b_0,
                          dim = settings$dims,
                          db = settings$db)

for(i in 2:66){
  settings <- simulate_df[i, ]
  i_mat <- read_params(ss = settings$sims,
                       b = settings$b_0,
                       dim = settings$dims,
                       db = settings$db)
  master_mat <- rbind(master_mat, i_mat)
}

master_mat$b0 <- as.numeric(as.character(master_mat$b0))

ggplot(subset, aes(x = b0, group = b0, y = rejr)) + geom_boxplot() +
  geom_point(data = points, colour = "red", alpha = 0.5)

sub_data <- subset(master_mat, dims == 10)



## This code is used to approximates the rejection 
## rate for a bonferroni based test to caompare 
## to ART


source("Sim_Study.R")
source("Testing.R")

gen_bonf_one <- function(b0, ss, sqrt_mat){
  gen_data <- make_data_loc(ss = ss, mat = sqrt_mat, b = b0)
  test <- test.once.bc(data = gen_data, alpha = 0.05)
  return(test)
}

gen_many_bonf <- function(b0, dim, ss, sims){
  x_cov <- matrix(0.5, nrow = dim, ncol = dim)
  diag(x_cov) <- 1
  sqrt_mat <- sqrtm(x_cov)
  rej_vec <- rep(NA, sims)
  for(rej_indx in 1:sims){
    if(rej_indx %% 50 == 0){cat(rej_indx/50)}
    rej_vec[rej_indx] <- gen_bonf_one(b0 = b0, ss = ss, sqrt_mat = sqrt_mat)
  }
  return(mean(rej_vec))
}

bzrs <- seq(0, 5, 0.5)
for(ii in 1:1){
  rej_rate_dim_10[ii] <- gen_many_bonf(bzrs[ii], 10, 5000, 1000)
}
rej_rate_dim_50 <- rep(0, 11)
for(ii in 1:11){
  rej_rate_dim_50[ii] <- gen_many_bonf(bzrs[ii], 50, 5000, 1000)
}

b0 <-  bzrs[job.id]
  
res_bo <- gen_many_bonf(b0, 1000, 10000, 1000)
res_and_bo <- c(b0, res_bo)
write.csv(res_and_bo, file = paste0("simres3/b0_", b0,".csv"))


## Here, we compare the Bonferroni and ART based methods

make_box_plot <- function(sub_data, bonf_data, labs){
  xlb <- labs[1]
  ylb <- labs[2]
  boxplot(rejr ~ b0, data = sub_data, main = labs[3], ylab = ylb, 
          xlab = xlb, pch = "+", cex.axis = 2, cex.lab = 2, cex.main = 2)
  points(1:11, bonf_data, col = "blue")
}

## These are the values for bonferoni based methods
## in the setting when p = 1000.  I read them from
## the clusters.  

bonf_p_1000 <- c(0.023, 0.033, 0.051,
                 0.102, 0.157, 0.250,
                 0.409, 0.555, 0.687,
                 0.832, 0.912)

## This code is used to recreate figure 4
layout(matrix(c(1, 2,3), ncol=3, byrow=TRUE), widths = c(4, 4, 3))

par(mai=rep(0.75, 4))
make_box_plot(sub_data = subset(master_mat, dims == 10), rej_rate_dim_10,
              labs = c(expression("b"[0]), "Rejection Rate", "p = 10"))
make_box_plot(sub_data = subset(master_mat, dims == 50), rej_rate_dim_50,
              labs = c(expression("b"[0]), "", "p = 50"))

par(mai=c(0,0,0,0))
plot.new()
legend("left",legend=c("Adaptive Resampling Test","Bonferoni Correction"), 
       col =c("black","blue"), title="Testing Method                      "
       , pch = c(0, 1) ,bty = "n", cex = 2)

## This code is used to recreate figure 5
layout(matrix(c(1, 2, 3, 4, 5, 3), ncol=3, byrow=TRUE), widths = c(4, 4, 3))

par(mai = c(0.25, 0.65, 0.5, 0.2))
make_box_plot(sub_data = subset(master_mat, db == 4.3), bonf_p_1000,
              labs = c(expression(""), "\n Rejection Rate", "a = 2"))
make_box_plot(sub_data = subset(master_mat, db == 6.1), bonf_p_1000,
              labs = c(expression(""), "", "a = 4"))

par(mai=c(0,0,0,0))
plot.new()
legend("left",legend=c("Adaptive Resampling Test","Bonferoni Correction"), 
       col =c("black","blue"), title="Testing Method                              "
       , pch = c(0, 1) ,bty = "n", cex = 2)

par(mai = c(0.5, 0.65, 0.5, 0.2))
make_box_plot(sub_data = subset(master_mat, db == 6.8), bonf_p_1000,
              labs = c(expression("b"[0]), "\n Rejection Rate", "a = 5"))
make_box_plot(sub_data = subset(master_mat, db == 8.6), bonf_p_1000,
              labs = c(expression("b"[0]), "", "a = 8"))

