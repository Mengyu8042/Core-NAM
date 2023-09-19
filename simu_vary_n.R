# Compare FULL, UNIF, LowCon, and CORE w.r.t. MSE and PMSE versus increasing n
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
plot_mse <- list()
plot_pmse <- list()

n_list <- c(1e3, 1e4, 1e5, 1e6, 1e7)  # training set size n
n_test <- 1e+6  # testing set size
sub <- 200  # subsample size r

for (func_type in 1:6) {
  print(paste0("function ", func_type))
  
  ## Initialization ##
  set.seed(1234)
  use_bam <- TRUE
  if (func_type <= 4) {
    k <- 40  # number of knots for univariate variables
  } else {
    k <- 30
  }
  k_ti <- 6  # number of knots for each marginal of 2D interactions
  nloop <- 100  # number of replicates
  SNR <- 5  # signal-to-noise ratio
  nn_list <- round(c(1, 2^0.5, 2^1, 2^1.5, 2^2) * 5 * sub)
  num_method <- 4
  mse_meta <- array(0, dim = c(nloop, num_method, length(n_list)))
  mse_temp <- matrix(0, num_method, length(n_list))
  pmse_meta <- array(0, dim = c(nloop, num_method, length(n_list)))
  pmse_temp <- matrix(0, num_method, length(n_list))
  
  ## Start Calculation ##
  for (i in 1:nloop) {
    for (j in 1:length(n_list)) {  
      n <- n_list[j]
      set.seed(100 + 432 * j + 123 * i)
      
      ## Generate the covariate(s) and response ##
      if (func_type <= 4) {
        x1 <- runif(n, 0, 1)
        X <- cbind(x1)
        x1_t <- (1:n_test)/n_test
        X_t <- cbind(x1_t)
      } else if (func_type == 5) {
        x1 <- runif(n, 0, 1); x2 <- runif(n, 0, 1)
        X <- cbind(x1, x2)
        X_t <- grid_points(n_test, 2) 
        x1_t <- X_t[, 1]; x2_t <- X_t[, 2]
      } else if (func_type == 6) {
        x1 <- runif(n, 0, 1); x2 <- runif(n, 0, 1); x3 <- runif(n, 0, 1); x4 <- runif(n, 0, 1)
        X <- cbind(x1, x2, x3, x4)
        X_t <- grid_points(n_test, 4) 
        x1_t <- X_t[, 1]; x2_t <- X_t[, 2]; x3_t <- X_t[, 3]; x4_t <- X_t[, 4]
      }
      y_raw <- yGenerate(X, func_type)  # response without noise
      y <- y_raw + rnorm(n, 0, sd(y_raw)/sqrt(SNR))  # response with noise
      
      y_raw_t <- yGenerate(X_t, func_type)  # response without noise
      y_t <- y_raw_t + rnorm(n_test, 0, sd(y_raw_t)/sqrt(SNR))  # response with noise
      
      ## Calculate design matrix and penalty matrix ##
      if (func_type <= 4) {
        kg1 <- kernel_generate(x1, k, id = 1:n, type_bs = "cr")
        kg1_t <- kernel_generate(x1_t, k, id = 1:n_test, type_bs = "cr")
        XX <- kg1$X
        XX_t <- kg1_t$X
        D_list <- list(D1 = kg1$D)
      } else if (func_type == 5) {
        kg1 <- kernel_generate_tp(x1, x2, k_ti, id = 1:n, type_bs = "te")
        kg1_t <- kernel_generate_tp(x1_t, x2_t, k_ti, id = 1:n_test, type_bs = "te")
        XX <- cbind(kg1$X)
        XX_t <- cbind(kg1_t$X)
        D_list <- list(D1 = kg1$D)
      } else if (func_type == 6) {
        kg1 <- kernel_generate(x1, k, id = 1:n, type_bs = "cr")
        kg2 <- kernel_generate(x2, k, id = 1:n, type_bs = "cr")
        kg3 <- kernel_generate(x3, k, id = 1:n, type_bs = "cr")
        kg4 <- kernel_generate(x4, k, id = 1:n, type_bs = "cr")
        kg1_t <- kernel_generate(x1_t, k, id = 1:n_test, type_bs = "cr")
        kg2_t <- kernel_generate(x2_t, k, id = 1:n_test, type_bs = "cr")
        kg3_t <- kernel_generate(x3_t, k, id = 1:n_test, type_bs = "cr")
        kg4_t <- kernel_generate(x4_t, k, id = 1:n_test, type_bs = "cr")
        XX <- cbind(kg1$X, kg2$X, kg3$X, kg4$X)
        XX_t <- cbind(kg1_t$X, kg2_t$X, kg3_t$X, kg4_t$X)
        D_list <- list(D1 = kg1$D, D2 = kg2$D, D3 = kg3$D, D4 = kg4$D)
      } 
      dim <- length(D_list)
      p <- ncol(XX)
      
      ## FULL ##
      if (use_bam == FALSE) {
        gcv <- optim(rep(1, dim), am.gcv, y = y, XX = XX, D_list = D_list,
                     control = list(maxit = 100, warn.1d.NelderMead = FALSE))
        sp_full <- exp(gcv$par)  # best smoothing parameters
        fit <- am.fit(y, XX, D_list, sp_full)
        b <- fit$b  # estimated coefficients
        mse <- am.mse(y_raw, XX, b)  # MSE
        pmse <- am.mse(y_raw_t, XX_t, b)  # PMSE
        
      } else {
        if (func_type <= 4) {
          df <- as.data.frame(cbind(y, x1)) 
          colnames(df) <- c("y", "x1")
          df_test <- as.data.frame(cbind(y_t, x1_t)) 
          colnames(df_test) <- c("y", "x1")
          fit <- bam(y ~ s(x1, k = k, bs = "cr"), data = df, method = "GCV.Cp")
        } else if (func_type == 5) {
          df <- as.data.frame(cbind(y, x1, x2))
          colnames(df) <- c("y", "x1", "x2")
          df_test <- as.data.frame(cbind(y_t, x1_t, x2_t))
          colnames(df_test) <- c("y", "x1", "x2")
          fit <- bam(y ~ te(x1, x2, k = k_ti, bs = "ps"), 
                     data = df, method = "GCV.Cp")
        } else if (func_type == 6) {
          df <- as.data.frame(cbind(y, x1, x2, x3, x4))
          colnames(df) <- c("y", "x1", "x2", "x3", "x4")
          df_test <- as.data.frame(cbind(y_t, x1_t, x2_t, x3_t, x4_t))
          colnames(df_test) <- c("y", "x1", "x2", "x3", "x4")
          fit <- bam(y ~ s(x1, k = k, bs = "cr") + s(x2, k = k, bs = "cr") + 
                       s(x3, k = k, bs = "cr") + s(x4, k = k, bs = "cr"), 
                     data = df, method = "GCV.Cp")
        }
        mse <- sum((y_raw - fit$fitted.values)^2)/n
        pred <- predict(fit, df_test)
        pmse <- sum((y_raw_t - pred)^2)/n_test
      }

      ## Subsampling ##
      ## UNIF ##
      id_unif <- sort(sample(n, sub, replace = TRUE))
      XXsub <- XX[id_unif, ]
      ysub <- y[id_unif]
  
      gcv <- optim(rep(1, dim), am.gcv, y = ysub, XX = XXsub, D_list = D_list,
                   control = list(maxit = 100, warn.1d.NelderMead = FALSE))
      sp_unif <- exp(gcv$par)
      fit <- am.fit(ysub, XXsub, D_list, sp_unif)
      b_unif <- fit$b
      mse_unif <- am.mse(y_raw, XX, b_unif)  # MSE
      pmse_unif <- am.mse(y_raw_t, XX_t, b_unif)  # PMSE
      
      ## LowCon ##
      theta <- 1
      x_standardize <- function(x) {
        x <- (x - min(x))/(max(x) - min(x))
        return(x)
      }
      lhd_data <- apply(X, 2, x_standardize)
      lowcon_space <- apply(lhd_data, 2, quantile, probs = c(theta/100, 1 - theta/100))
      des = randomLHS(sub, dim(X)[2])
      for (dd in 1:dim(X)[2]) {
        des[, dd] = (des[, dd]) * (lowcon_space[2, dd] - lowcon_space[1, dd]) + lowcon_space[1, dd]
      }
      id_lowcon <- nabor::knn(lhd_data, des, k = 1, eps = 0.1)$nn.idx[, 1]
      XXsub <- XX[id_lowcon, ]
      ysub <- y[id_lowcon]
      
      gcv <- optim(rep(1, dim), am.gcv, y = ysub, XX = XXsub, D_list = D_list,
                   control = list(maxit = 100, warn.1d.NelderMead = FALSE))
      sp_lowcon <- exp(gcv$par)
      fit <- am.fit(ysub, XXsub, D_list, sp_lowcon)
      b_lowcon <- fit$b
      mse_lowcon <- am.mse(y_raw, XX, b_lowcon)  # MSE
      pmse_lowcon <- am.mse(y_raw_t, XX_t, b_lowcon)  # PMSE
      
      ## Core-NAM ##
      # select the core-elements
      n2 <- min(nn_list[j], n)
      id <- sort(sample(n, n2, replace = FALSE))
      XX2 <- XX[id, ]
      y2 <- y[id]
      XXsub <- matrix(0, n2, p)
      
      thres <- colnth(abs(XX2), rep(sub, p), descending = TRUE, parallel = TRUE)
      for (kk in 1:p) {
        col <- abs(XX2[, kk])
        id <- which(col > thres[kk])
        id <- c(id, which(col == thres[kk])[1:(sub - length(id))])
        XXsub[id, kk] <- XX2[id, kk]
      }
      id_nz <- (rowSums(abs(XXsub)) > 0)
      XXrow <- XX2[id_nz, ]
      XXsub <- XXsub[id_nz, ]
      ysub <- y2[id_nz]
      
      # choose the smoothing parameter via core-elements GCV
      gcv <- optim(rep(1, dim), am.core.gcv, ysub = ysub, XXrow = XXrow,
                   XXsub = XXsub, D_list = D_list,
                   control = list(maxit = 100, warn.1d.NelderMead = FALSE))
      sp_core <- exp(gcv$par)
      
      # core-elements estimation based on the chosen smoothing parameter
      fit <- am.core.fit(ysub, XXrow, XXsub, D_list, sp_core)
      b_core <- fit$b
      mse_core <- am.mse(y_raw, XX, b_core)  # MSE
      pmse_core <- am.mse(y_raw_t, XX_t, b_core)  # PMSE
      
      ######################################################
      mse_temp[, j] <- c(mse, mse_unif, mse_lowcon, mse_core)
      pmse_temp[, j] <- c(pmse, pmse_unif, pmse_lowcon, pmse_core)
    }
      
    mse_meta[i, , ] <- mse_temp
    pmse_meta[i, , ] <- pmse_temp
  }
  
  FULL.mse <- mse_meta[, 1, ]
  UNIF.mse <- mse_meta[, 2, ]
  LOWCON.mse <- mse_meta[, 3, ]
  CORE.mse <- mse_meta[, 4, ]
  
  FULL.pmse <- pmse_meta[, 1, ]
  UNIF.pmse <- pmse_meta[, 2, ]
  LOWCON.pmse <- pmse_meta[, 3, ]
  CORE.pmse <- pmse_meta[, 4, ]
  
  mse_mat <- data.frame(mean = c(apply(FULL.mse, 2, mean), apply(UNIF.mse, 2, mean), 
                                 apply(LOWCON.mse, 2, mean), apply(CORE.mse, 2, mean)),
                        sd = c(apply(FULL.mse, 2, sd), apply(UNIF.mse, 2, sd), 
                               apply(LOWCON.mse, 2, sd), apply(CORE.mse, 2, sd)),
                        Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                            each = length(log10(n_list)))),
                        sub = rep(log10(n_list), num_method))
  mse_mat$Method <- factor(mse_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))
  
  pmse_mat <- data.frame(mean = c(apply(FULL.pmse, 2, mean), apply(UNIF.pmse, 2, mean), 
                                  apply(LOWCON.pmse, 2, mean), apply(CORE.pmse, 2, mean)),
                        sd = c(apply(FULL.pmse, 2, sd), apply(UNIF.pmse, 2, sd), 
                               apply(LOWCON.pmse, 2, sd), apply(CORE.pmse, 2, sd)),
                        Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                            each = length(log10(n_list)))),
                        sub = rep(log10(n_list), num_method))
  pmse_mat$Method <- factor(pmse_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))
  
  ## Plot the results ##
  if (func_type == 1) { func_name <- "F1"
  } else if (func_type == 2) { func_name <- "F2"
  } else if (func_type == 3) { func_name <- "F3"
  } else if (func_type == 4) { func_name <- "F4"
  } else if (func_type == 5) { func_name <- "F5"
  } else if (func_type == 6) { func_name <- "F6"
  }
  
  width <- (max(log10(n_list)) - min(log10(n_list)))/10
  pd <- position_dodge(width/3)
  if (func_type == 6) {
    leg_pos <- c(0.75, 0.8)
  } else {
    leg_pos <- ""
  }
  
  # MSE
  p3 <- ggplot(mse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
  p3 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_rect(colour = "black"),
                                axis.text = element_text(size = 18),
                                axis.title = element_text(size = 18),
                                legend.position = leg_pos,
                                legend.title = element_blank(),
                                legend.background = element_rect(fill = "white", color = "black"),
                                legend.key.width = unit(3, "line"),
                                legend.key.height = unit(1, "line"),
                                legend.text = element_text(size = 13)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
    geom_line(aes(linetype = Method, size = Method), position = pd) +
    geom_point(position = pd, aes(shape = Method), size = 2) +
    scale_shape_manual(values = c(4, 7, 9, 1)) +
    scale_linetype_manual(values = c(4, 2, 6, 1)) +
    scale_size_manual(values = c(1, 1, 1, 1)) +
    scale_color_manual(values = c("#377eb8","#999999","#984ea3","#e41a1c")) +
    labs(title = func_name, x = expression(log[10](n)), y = "MSE") +
    theme(plot.title = element_text(hjust = 0.5, size = 18))
  plot_mse[[func_type]] <- p3
  print(p3)
  
  # PMSE
  p4 <- ggplot(pmse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
  p4 <- p4 + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_rect(colour = "black"),
                                axis.text = element_text(size = 18),
                                axis.title = element_text(size = 18),
                                legend.position = leg_pos,
                                legend.title = element_blank(),
                                legend.background = element_rect(fill = "white", color = "black"),
                                legend.key.width = unit(3, "line"),
                                legend.key.height = unit(1, "line"),
                                legend.text = element_text(size = 13)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
    geom_line(aes(linetype = Method, size = Method), position = pd) +
    geom_point(position = pd, aes(shape = Method), size = 2) +
    scale_shape_manual(values = c(4, 7, 9, 1)) +
    scale_linetype_manual(values = c(4, 2, 6, 1)) +
    scale_size_manual(values = c(1, 1, 1, 1)) +
    scale_color_manual(values = c("#377eb8","#999999","#984ea3","#e41a1c")) +
    labs(title = func_name, x = expression(log[10](n)), y = "PMSE") +
    theme(plot.title = element_text(hjust = 0.5, size = 18))
  plot_pmse[[func_type]] <- p4
  print(p4)
}

pdf("simu_n_mse.pdf", width = 13, height = 8)
grid.arrange(plot_mse[[1]], plot_mse[[2]], plot_mse[[3]], 
             plot_mse[[4]], plot_mse[[5]], plot_mse[[6]],
             nrow = 2, ncol = 3)
dev.off()

pdf("simu_n_pmse.pdf", width = 13, height = 8)
grid.arrange(plot_pmse[[1]], plot_pmse[[2]], plot_pmse[[3]], 
             plot_pmse[[4]], plot_pmse[[5]], plot_pmse[[6]],
             nrow = 2, ncol = 3)
dev.off()
