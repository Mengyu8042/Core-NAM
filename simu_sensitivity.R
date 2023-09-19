# Sensitivity analysis by varying the signal-to-noise ratio (SNR) and basis dimension (q)
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
func_type <- 1
plot_num <- 0

for (k in c(30, 40, 50)) {
  for (SNR in c(2, 5, 10)) {
    plot_num <- plot_num + 1
    print(plot_num)
    
    ## Initialization ##
    set.seed(1234)
    n <- 1e+4  # training set size 
    n_test <- 1e+6  # testing set size 
    nloop <- 100  # number of replicates
    num_method <- 4
    sub_meta <- (1:5) * 100  # subsample size r
    mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
    mse_temp <- matrix(0, num_method, length(sub_meta))
    pmse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
    pmse_temp <- matrix(0, num_method, length(sub_meta))
    
    ## Start Calculation ##
    for (i in 1:nloop) {
      set.seed(100 + 123 * i)
      
      ## Generate the covariate(s) and response ##
      x1 <- runif(n, 0, 1)
      X <- cbind(x1)
      x1_t <- (1:n_test) / n_test
      X_t <- cbind(x1_t)
      
      y_raw <- yGenerate(X, func_type)  # without noise
      y <- y_raw + rnorm(n, 0, sd(y_raw)/sqrt(SNR))  # observations with noise
      y_raw_t <- yGenerate(X_t, func_type)  # without noise
      y_t <- y_raw_t + rnorm(n_test, 0, sd(y_raw_t)/sqrt(SNR))  # observations with noise
      
      ## Calculate design matrix and penalty matrix ##
      kg1 <- kernel_generate(x1, k, id = 1:n, type_bs = "cr")
      kg1_t <- kernel_generate(x1_t, k, id = 1:n_test, type_bs = "cr")
      XX <- kg1$X 
      XX_t = kg1_t$X
      D_list = list(D1 = kg1$D)
      dim <- length(D_list)
      p <- ncol(XX)
      
      ## FULL ##
      gcv <- optim(rep(1, dim), am.gcv, y = y, XX = XX, D_list = D_list,
                   control = list(maxit = 500, warn.1d.NelderMead = FALSE))
      sp_full <- exp(gcv$par)  # best smoothing parameter
      fit <- am.fit(y, XX, D_list, sp_full)
      b <- fit$b  # estimated coefficients
      mse <- am.mse(y_raw, XX, b)  # MSE
      pmse <- am.mse(y_raw_t, XX_t, b)  # PMSE
      
      ## Subsampling ##
      for (j in 1:length(sub_meta)) {
        set.seed(100 + 432 * j + 123 * i)
        sub <- sub_meta[j]
        
        ## UNIF ##
        id_unif <- sort(sample(n, sub, replace = TRUE))
        XXsub <- XX[id_unif, ]
        ysub <- y[id_unif]
        
        gcv <- optim(rep(1, dim), am.gcv, y = ysub, XX = XXsub, D_list = D_list,
                     control = list(maxit = 500, warn.1d.NelderMead = FALSE))
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
                     control = list(maxit = 500, warn.1d.NelderMead = FALSE))
        sp_lowcon <- exp(gcv$par)
        fit <- am.fit(ysub, XXsub, D_list, sp_lowcon)
        b_lowcon <- fit$b
        mse_lowcon <- am.mse(y_raw, XX, b_lowcon)  # MSE
        pmse_lowcon <- am.mse(y_raw_t, XX_t, b_lowcon)  # PMSE
        
        ## Core-NAM ##
        # select the core-elements
        XXsub <- matrix(0, n, p)
        thres <- colnth(abs(XX), rep(sub, p), descending = TRUE, parallel = FALSE)
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
        
        # choose the smoothing parameter via core-elements GCV
        gcv <- optim(rep(1, dim), am.core.gcv, ysub = ysub, XXrow = XXrow,
                     XXsub = XXsub, D_list = D_list,
                     control = list(maxit = 500, warn.1d.NelderMead = FALSE))
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
    
    FULL.mse <- log10(mse_meta[, 1, ])
    UNIF.mse <- log10(mse_meta[, 2, ])
    LOWCON.mse <- log10(mse_meta[, 3, ])
    CORE.mse <- log10(mse_meta[, 4, ])
    
    FULL.pmse <- log10(pmse_meta[, 1, ])
    UNIF.pmse <- log10(pmse_meta[, 2, ])
    LOWCON.pmse <- log10(pmse_meta[, 3, ])
    CORE.pmse <- log10(pmse_meta[, 4, ])
    
    mse_mat <- data.frame(mean = c(apply(FULL.mse, 2, mean), apply(UNIF.mse, 2, mean), 
                                   apply(LOWCON.mse, 2, mean), apply(CORE.mse, 2, mean)),
                          sd = c(apply(FULL.mse, 2, sd), apply(UNIF.mse, 2, sd), 
                                 apply(LOWCON.mse, 2, sd), apply(CORE.mse, 2, sd)),
                          Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                              each = length(sub_meta))),
                          sub = rep(sub_meta, num_method))
    mse_mat$Method <- factor(mse_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))
    
    pmse_mat <- data.frame(mean = c(apply(FULL.pmse, 2, mean), apply(UNIF.pmse, 2, mean), 
                                    apply(LOWCON.pmse, 2, mean), apply(CORE.pmse, 2, mean)),
                           sd = c(apply(FULL.pmse, 2, sd), apply(UNIF.pmse, 2, sd), 
                                  apply(LOWCON.pmse, 2, sd), apply(CORE.pmse, 2, sd)),
                           Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                               each = length(sub_meta))),
                           sub = rep(sub_meta, num_method))
    pmse_mat$Method <- factor(pmse_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))
    
    ## Plot the results ##
    func_name <- "F1"
    width <- (max(sub_meta) - min(sub_meta))/10
    pd <- position_dodge(width/3)
    if (plot_num == 9) {
      leg_pos <- c(0.73, 0.81)
    } else {
      leg_pos <- ""
    }
    plot_title <- paste0("q = ", k, ", SNR = ", SNR)
    
    # MSE
    p3 <- ggplot(mse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
    p3 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_rect(colour = "black"),
                                  axis.text = element_text(size = 19),
                                  axis.title = element_text(size = 19),
                                  legend.position = leg_pos,
                                  legend.title = element_blank(),
                                  legend.background = element_rect(fill = "white", color = "black"),
                                  legend.key.width = unit(3, "line"),
                                  legend.key.height = unit(1, "line"),
                                  legend.text = element_text(size = 15)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
      geom_line(aes(linetype = Method, size = Method), position = pd) +
      geom_point(position = pd, aes(shape = Method), size = 2) +
      scale_shape_manual(values = c(4, 7, 9, 1)) +
      scale_linetype_manual(values = c(4, 2, 6, 1)) +
      scale_size_manual(values = c(1, 1, 1, 1)) +
      scale_color_manual(values = c("#377eb8","#999999","#984ea3","#e41a1c")) +
      labs(title = plot_title, x = "r", y = expression(log[10](MSE))) +
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    plot_mse[[plot_num]] <- p3
    print(p3)
    
    # PMSE
    p4 <- ggplot(pmse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
    p4 <- p4 + theme_bw() + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_rect(colour = "black"),
                                  axis.text = element_text(size = 19),
                                  axis.title = element_text(size = 19),
                                  legend.position = leg_pos,
                                  legend.title = element_blank(),
                                  legend.background = element_rect(fill = "white", color = "black"),
                                  legend.key.width = unit(3, "line"),
                                  legend.key.height = unit(1, "line"),
                                  legend.text = element_text(size = 15)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
      geom_line(aes(linetype = Method, size = Method), position = pd) +
      geom_point(position = pd, aes(shape = Method), size = 2) +
      scale_shape_manual(values = c(4, 7, 9, 1)) +
      scale_linetype_manual(values = c(4, 2, 6, 1)) +
      scale_size_manual(values = c(1, 1, 1, 1)) +
      scale_color_manual(values = c("#377eb8","#999999","#984ea3","#e41a1c")) +
      labs(title = plot_title, x = "r", y = expression(log[10](PMSE))) +
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    plot_pmse[[plot_num]] <- p4
    print(p4)
  }
}

pdf("simu_sens_mse.pdf", width = 13, height = 12)
grid.arrange(plot_mse[[1]], plot_mse[[2]], plot_mse[[3]], 
             plot_mse[[4]], plot_mse[[5]], plot_mse[[6]],
             plot_mse[[7]], plot_mse[[8]], plot_mse[[9]],
             nrow = 3, ncol = 3)
dev.off()

pdf("simu_sens_pmse.pdf", width = 13, height = 12)
grid.arrange(plot_pmse[[1]], plot_pmse[[2]], plot_pmse[[3]], 
             plot_pmse[[4]], plot_pmse[[5]], plot_pmse[[6]],
             plot_pmse[[7]], plot_pmse[[8]], plot_pmse[[9]],
             nrow = 3, ncol = 3)
dev.off()
