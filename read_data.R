#########################################################
#### Author : Adam Elder
#### This script contains the code used to read in 
#### simulation results.
#########################################################
setwd("/Users/adamelder/Dropbox/Spring\ 2018/BIOSTAT572/SimulationStudy/BIOSTAT572R")
read_params <- function(n, mod, dim, rho){
  file_name <- paste0("simres/n", n, "mod", mod, "dim", dim, "rho", rho, "iter")
  for(i in 1:5){
    iter_i <- paste0(file_name, i, ".csv")
    if(file.exists(iter_i)){
      if (exists("results")){
        results <- rbind(results, t(read.csv(iter_i)[, 2]))
      }else{
        results <- t(read.csv(iter_i)[, 2]) 
      }
    }
  }
  if(exists("results")){
    ave_res <- apply(results, 2, mean)
  }else{
    ave_res <- rep(NA, 5)
  }
  names(ave_res) <- c("LRT", "Bonferroni", "CPB", "ART", "TIME")
  return(round(ave_res, 4))
}

rho <- rep(c(0.8, 0.5, 0), each = 30)
sample_size <- rep(rep(c(100, 200), each = 15), times = 3)
dims <- rep(rep(c(10, 50, 100, 150, 200), each = 3), times = 6)
model <- rep(c(1, 2, 3), times = 30)

simulate_df <- data.frame("rho" = rho, "n" = sample_size,
                          "dims" = dims, "model" = model)
simulate_df$LRT <- simulate_df$Bonferroni <- simulate_df$CPB <- simulate_df$ART <- NA

for( i in 1:nrow(simulate_df)){
  simulate_df[i, 5:8] <- read_params(n = simulate_df[i, 2],
                                     mod = simulate_df[i, 4],
                                     dim = simulate_df[i, 3],
                                     rho = simulate_df[i, 1])[c(4, 3, 2, 1)]
}

single_mod <- function(data, title, xlab, ylab, n_val, mod_val, rho_val, legty = "None", base.s){
  if (mod_val == 1){ylim <- 0.2}else{ylim <- 1}
  if (mod_val == 1){lsize <- 1}else{lsize <- 0}
  sub_data <- subset(data, rho == rho_val & n == n_val & model == mod_val, 
                     select = c("dims", "LRT", "Bonferroni", "CPB", "ART"))
  gg_data <- gather(sub_data, key = "test_type", value = "reject_rate", 2:5)
  if (legty == "A"){
    A <- ggplot(gg_data[complete.cases(gg_data), ], 
                aes(x = dims, y = reject_rate, 
                    group = test_type, shape = test_type, colour = test_type)) + 
      geom_point(size = 3) + geom_line() +
      coord_cartesian(ylim = c(0, ylim)) +
      labs(x = xlab, y = ylab, title = title) + 
      geom_abline(aes(intercept=0.05, slope=0, color = "red"), size = lsize,
                  show.legend = FALSE, linetype="dotted") + 
      scale_colour_discrete(name = "Test Type:   ",
                            breaks = c("LRT", "Bonferroni", "CPB", "ART"),
                            labels = c("Likelihood Ratio Test   ", "Bonferroni   ", "Naive Bootstrap   ", "ART ")) +
      scale_shape_discrete(name = "Test Type:   ",
                           breaks = c("LRT", "Bonferroni", "CPB", "ART"),
                           labels = c("Likelihood Ratio Test   ", "Bonferroni   ", "Naive Bootstrap   ", "ART "))+
      theme_minimal(base_size = base.s)
  }else{
    A <- ggplot(gg_data[complete.cases(gg_data), ], 
                aes(x = dims, y = reject_rate, 
                    group = test_type, shape = test_type, colour = test_type)) + 
      geom_point(size = 3) + geom_line() +
      coord_cartesian(ylim = c(0, ylim)) +
      labs(x = xlab, y = ylab, title = title) + 
      geom_abline(aes(intercept=0.05, slope=0, color = "red"), size = lsize,
                  show.legend = FALSE, linetype="dotted") + 
      theme_minimal(base_size = base.s) + theme(legend.position = "None")
  }
  
  return(A)
}

make_rho_plot <- function(data, rho_val, basesize){
  A <- single_mod(data, "Null Model", "", "Rejection Rate (n = 100)",
                  n_val = 100, mod_val = 1, rho_val = rho_val, base.s = basesize)
  B <- single_mod(data, "One Significant  Predictor", "Dimension", "",
                  n_val = 100, mod_val = 2, rho_val = rho_val, base.s = basesize)
  C <- single_mod(data, "Ten Significant Predictors", "", "",
                  n_val = 100, mod_val = 3, rho_val = rho_val, 
                  legty = "A", base.s = basesize)
  D <- single_mod(data, "", "", "Rejection Rate (n = 200)",
                  n_val = 200, mod_val = 1, rho_val = rho_val, base.s = basesize)
  E <- single_mod(data, "", "Dimension", "",
                  n_val = 200, mod_val = 2, rho_val = rho_val, base.s = basesize)
  EF <- single_mod(data, "", "", "",
                  n_val = 200, mod_val = 3, rho_val = rho_val, 
                  legty = "A", base.s = basesize)
  plot_1 <- plot_grid(A, B, C + theme(legend.position = "none"),
                      D, E, EF + theme(legend.position = "none"),
                      ncol = 3, rel_widths = c(1, 1, 1, 1, 1, 1))
  plot_leg <- get_legend(C + theme_minimal(base_size = basesize * 1.15) +theme(legend.position = "bottom", legend.key.size = unit(3, 'lines')))
  plot_grid(plot_1, plot_leg, ncol = 1, rel_heights = c(1, .1))
}



