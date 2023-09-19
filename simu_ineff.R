# Empirically verify the asymptomatic optimality of CORE-GCV
rm(list = ls())
gc()

## Import R Packages ##
list.of.packages <- c("mgcv", "MASS", "mvtnorm", "corpcor", "SparseM", "pracma", 
                      "Rfast", "Rcpp", "RcppEigen", "Matrix", "nabor", "lhs",
                      "ReIns", "ggplot2", "reshape2", "patchwork", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
library(mgcv)
library(MASS)
library(mvtnorm)
library(corpcor)
library(SparseM)
library(pracma)
library(Rfast)
library(Rcpp)
library(RcppEigen)
library(Matrix)
library(nabor)
library(lhs)
library(ReIns)
library(ggplot2)
library(reshape2)
library(patchwork)
library(gridExtra)

## Import the source file and utility file ##
setwd(getwd())
sourceCpp("source/funcs.cpp")
source("utility.R")

## Beginning ##
plot_ineff <- list()
n_list <- round(c(10^3, 10^3.5, 10^4, 10^4.5))  # full sample size n 
sub_list <- c(100, 200, 400, 800)  # subsample size r

for (func_type in 1:6) {
  print(paste0("function ", func_type))
  
  ## Initialization ##
  set.seed(1234)
  if (func_type <= 4) {
    k <- 40  # number of knots for univariate variables
  } else {
    k <- 30
  }
  k_ti <- 6  # number of knots for each marginal of 2D interactions
  nloop <- 100  # number of replicates
  SNR <- 5  # signal-to-noise ratio
  ineff <- matrix(0, nloop, length(n_list))
  
  ## Start Calculation ##
  for (i in 1:nloop) {
    for (j in 1:length(n_list)) {  
      n <- n_list[j]
      sub <- sub_list[j]
      set.seed(100 + 432 * j + 123 * i)
      
      ## Generate the covariate(s) and response ##
      if (func_type <= 4) {
        x1 <- runif(n, 0, 1)
        X <- cbind(x1)
      } else if (func_type == 5) {
        x1 <- runif(n, 0, 1); x2 <- runif(n, 0, 1)
        X <- cbind(x1, x2)
      } else if (func_type == 6) {
        x1 <- runif(n, 0, 1); x2 <- runif(n, 0, 1); x3 <- runif(n, 0, 1); x4 <- runif(n, 0, 1)
        X <- cbind(x1, x2, x3, x4)
      }
      y_raw <- yGenerate(X, func_type)  # response without noise
      y <- y_raw + rnorm(n, 0, sd(y_raw)/sqrt(SNR))  # response with noise
      
      ## Calculate design matrix and penalty matrix ##
      if (func_type <= 4) {
        kg1 <- kernel_generate(x1, k, id = 1:n, type_bs = "cr")
        XX <- kg1$X 
        D_list = list(D1 = kg1$D)
      } else if (func_type == 5) {
        kg1 <- kernel_generate_tp(x1, x2, k_ti, id = 1:n, type_bs = "te")
        XX <- cbind(kg1$X)
        D_list <- list(D1 = kg1$D)
      } else if (func_type == 6) {
        kg1 <- kernel_generate(x1, k, id = 1:n, type_bs = "cr")
        kg2 <- kernel_generate(x2, k, id = 1:n, type_bs = "cr")
        kg3 <- kernel_generate(x3, k, id = 1:n, type_bs = "cr")
        kg4 <- kernel_generate(x4, k, id = 1:n, type_bs = "cr")
        XX <- cbind(kg1$X, kg2$X, kg3$X, kg4$X)
        D_list <- list(D1 = kg1$D, D2 = kg2$D, D3 = kg3$D, D4 = kg4$D)
      } 
      dim <- length(D_list)
      p <- ncol(XX)
      
      ## Core-NAM ##
      # select the core-elements
      XXsub <- matrix(0, n, p)
      thres <- colnth(abs(XX), rep(sub, p), descending = TRUE, parallel = TRUE)
      for (kk in 1:p) {
        col <- abs(XX[, kk])
        id <- which(col > thres[kk])
        id <- c(id, which(col == thres[kk])[1:(sub - length(id))])
        XXsub[id, kk] <- XX[id, kk]
      }
      id_nz <- (rowSums(abs(XXsub)) > 0)
      XXrow <- XX[id_nz, ]
      XXsub <- XXsub[id_nz, ]
      ysub <- y[id_nz]
      y_raw_sub <- y_raw[id_nz]
      
      # estimation based on the smoothing parameter chose by CORE-GCV
      gcv <- optim(rep(1, dim), am.core.gcv, ysub = ysub, XXrow = XXrow,
                   XXsub = XXsub, D_list = D_list,
                   control = list(maxit = 200, warn.1d.NelderMead = FALSE))
      sp_core <- exp(gcv$par)
      fit <- am.core.fit(ysub, XXrow, XXsub, D_list, sp_core)
      b_core <- fit$b
      
      # estimation based on the smoothing parameter chose by minimizing the loss function
      loss <- optim(rep(1, dim), am.core.loss, y_raw_sub = y_raw_sub, ysub = ysub, 
                    XXrow = XXrow, XXsub = XXsub, D_list = D_list,
                    control = list(maxit = 200, warn.1d.NelderMead = FALSE))
      sp_core_opt <- exp(loss$par)
      fit <- am.core.fit(ysub, XXrow, XXsub, D_list, sp_core_opt)
      b_core_opt <- fit$b
      
      ineff[i, j] <- am.mse(y_raw_sub, XXrow, b_core)/am.mse(y_raw_sub, XXrow, b_core_opt)
    }
  }  
  
  ineff_mat <- data.frame(mean = apply(ineff, 2, mean),
                          sd = apply(ineff, 2, sd),
                          sub = log10(n_list))
  
  ## Plot the results ##
  if (func_type == 1) { func_name <- "F1"
  } else if (func_type == 2) { func_name <- "F2"
  } else if (func_type == 3) { func_name <- "F3"
  } else if (func_type == 4) { func_name <- "F4"
  } else if (func_type == 5) { func_name <- "F5"
  } else if (func_type == 6) { func_name <- "F6"
  }
  
  width <- (max(ineff_mat$sub) - min(ineff_mat$sub))/20
  
  p3 <- ggplot(ineff_mat, aes(x = sub, y = mean))
  p3 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_rect(colour = "black"),
                                axis.text = element_text(size = 18),
                                axis.title = element_text(size = 18),
                                legend.position = "",
                                legend.title = element_blank(),
                                legend.key.width = unit(3, "line"),
                                legend.key.height = unit(1, "line")) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = "#e41a1c"), 
                  width = width, show.legend = FALSE) +
    geom_line(aes(color = "#e41a1c"), alpha = 1, lwd = 0.9) +
    geom_point(aes(color = "#e41a1c"), alpha = 1, size = 2) +
    geom_hline(yintercept = 1, linetype = 'dashed', lwd = 0.7) + 
    scale_color_manual(values = c("#e41a1c")) +
    labs(title = func_name, x = expression(log[10](n)), y = "Inefficiency") +
    theme(plot.title = element_text(hjust = 0.5, size = 18))
  plot_ineff[[func_type]] <- p3
  print(p3)
}

pdf("simu_ineff.pdf", width = 13, height = 8)
grid.arrange(plot_ineff[[1]], plot_ineff[[2]], plot_ineff[[3]],
             plot_ineff[[4]], plot_ineff[[5]], plot_ineff[[6]],
             nrow = 2, ncol = 3)
dev.off()
