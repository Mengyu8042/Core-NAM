# Compare CORE based on smoothness-adaptive knots and uniform knots
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
sourceCpp("../source/funcs.cpp")
source("utility.R")

## Beginning ##
plot_mse <- list()
plot_pmse <- list()
func_type <- 4

## Initialization ##
set.seed(1234)
n <- 1e+4  # training set size 
n_test <- 1e+6  # testing set size 
k_range <- c(30, 40, 50)
nloop <- 10  # number of replicates
SNR <- 10  # signal-to-noise ratio
num_method <- 6
sub_meta <- (1:5) * 100  # subsample size r
mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
mse_temp <- matrix(0, num_method, length(sub_meta))
pmse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
pmse_temp <- matrix(0, num_method, length(sub_meta))

## Start Calculation ##
for (i in 1:nloop) {
  print(paste0("Replication ", i, "/", nloop))
  set.seed(100 + 123 * i)
  
  ## Generate the covariate(s) and response ##
  X <- matrix(runif(n), ncol = 1)
  X_t <- grid_points(n_test, 1) 
  
  y_raw <- yGenerate(X, func_type)  # response without noise
  y <- y_raw + rnorm(n, 0, sd(y_raw)/sqrt(SNR))  # response with noise
  
  y_raw_t <- yGenerate(X_t, func_type)  # response without noise
  y_t <- y_raw_t + rnorm(n_test, 0, sd(y_raw_t)/sqrt(SNR))  # response with noise
  
  # Adaptive knot selection
  adaptive_knots <- function(x, y, k, window_size = 201) {
    sy <- filter(y, rep(1/window_size, window_size), sides = 2)
    y_diff <- diff(sy) / diff(x)  # estimated roughness
    y_diff <- c(0, y_diff)
    y_diff[is.na(y_diff)] <- 0
    xk <- sort(sample(x, k - 2, prob = abs(y_diff), replace = FALSE))
    xk <- c(min(x), xk, max(x))
    xk <- unique(xk)
    return(xk)
  }
  
  # Uniform knot selection
  uniform_knots <- function(x, k) {
    xk <- sort(sample(x, k - 2, replace = FALSE))
    xk <- c(min(x), xk, max(x))
    xk <- unique(xk)
    return(xk)
  }
  
  # Generate univariate spline bases given knots
  generate_kernel_knots <- function(x, knots) {
    k <- length(knots)
    knots <- data.frame(x = knots)
    dat <- data.frame(x)
    sspec <- s(x, bs = "cr", k = k) # cubic regression spline
    smv <- smoothCon(sspec, data = dat, knots = knots)[[1]]
    S <- smv$S[[1]]
    D <- t(mroot(S, method = "svd"))
    X <- smv$X
    list(X = X, D = D)
  }
  
  for (k in k_range) {
    aknots = adaptive_knots(X[,1], y, k)
    uknots = uniform_knots(X[,1], k)
    
    ## Calculate design matrix and penalty matrix ##
    generate_kernels <- function(x, x_t, knots) {
      kg <- generate_kernel_knots(x, knots)
      kg_t <- generate_kernel_knots(x_t, knots)
      list(X = kg$X, X_t = kg_t$X, D = kg$D)
    }
    
    kernels <- generate_kernels(X[,1], X_t[,1], uknots)
    XX <- kernels$X
    XX_t <- kernels$X_t
    D_list <- list(D1 = kernels$D)
    dim <- length(D_list)
    
    ada_kernels <- generate_kernels(X[,1], X_t[,1], aknots)
    XXa <- ada_kernels$X
    XXa_t <- ada_kernels$X_t
    Da_list <- list(D1 = ada_kernels$D)
    
    ## Subsampling ##
    for (j in 1:length(sub_meta)) {
      set.seed(100 + 432 * j + 123 * i)
      sub <- sub_meta[j]
      
      ## Core-NAM ##
      # select the core-elements
      d <- ncol(XX)
      XXsub <- matrix(0, n, d)
      thres <- colnth(abs(XX), rep(sub, d), descending = TRUE, parallel = FALSE)
      for (kk in 1:d) {
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
      
      ## Core-NAM (adaptive) ##
      # select the core-elements
      d <- ncol(XXa)
      XXasub <- matrix(0, n, d)
      thres <- colnth(abs(XXa), rep(sub, d), descending = TRUE, parallel = FALSE)
      for (kk in 1:d) {
        col <- abs(XXa[, kk])
        id <- which(col > thres[kk])
        id <- c(id, which(col == thres[kk])[1:(sub - length(id))])
        XXasub[id, kk] <- XXa[id, kk]
      }
      id_nz <- (rowSums(abs(XXasub)) > 0)
      XXarow <- XXa[id_nz, ]
      XXasub <- XXasub[id_nz, ]
      ysub <- y[id_nz]
      
      # choose the smoothing parameter via core-elements GCV
      gcva <- optim(rep(1, dim), am.core.gcv, ysub = ysub, XXrow = XXarow,
                    XXsub = XXasub, D_list = Da_list,
                    control = list(maxit = 500, warn.1d.NelderMead = FALSE))
      spa_core <- exp(gcva$par)
      
      # core-elements estimation based on the chosen smoothing parameter
      fita <- am.core.fit(ysub, XXarow, XXasub, Da_list, spa_core)
      ba_core <- fita$b
      msea_core <- am.mse(y_raw, XXa, ba_core)  # MSE
      pmsea_core <- am.mse(y_raw_t, XXa_t, ba_core)  # PMSE      
      
      ######################################################
      mse_temp[((k/10-2)*2-1):((k/10-2)*2), j] <- c(mse_core, msea_core)
      pmse_temp[((k/10-2)*2-1):((k/10-2)*2), j] <- c(pmse_core, pmsea_core)
    }
  }

  mse_meta[i, , ] <- mse_temp
  pmse_meta[i, , ] <- pmse_temp
}

process_data <- function(data) {
  log_data <- log10(data)
  mean_sd_data <- do.call(rbind, lapply(1:ncol(log_data), function(ii) {
    data.frame(mean = apply(log_data[, ii, ], 2, mean, na.rm = T), 
               sd = apply(log_data[, ii, ], 2, sd, na.rm = T))
  }))
  return(mean_sd_data)
}

mse_processed <- process_data(mse_meta)
pmse_processed <- process_data(pmse_meta)

methods <- c("CORE-U (q=30)", "CORE-A (q=30)",
             "CORE-U (q=40)", "CORE-A (q=40)",
             "CORE-U (q=50)", "CORE-A (q=50)")

create_data_frame <- function(processed_data) {
  data.frame(mean = processed_data$mean,
             sd = processed_data$sd,
             Method = factor(rep(methods, each = length(sub_meta))),
             sub = rep(sub_meta, times = length(methods)))
}

mse_mat <- create_data_frame(mse_processed)
pmse_mat <- create_data_frame(pmse_processed)
mse_mat$Method <- factor(mse_mat$Method, levels = methods)
pmse_mat$Method <- factor(pmse_mat$Method, levels = methods)

## Plot the results ##
func_name <- paste0("F", func_type)
width <- (max(sub_meta) - min(sub_meta))/10
pd <- position_dodge(width/3)

# MSE
p3 <- ggplot(mse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
p3 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 19),
                              axis.title = element_text(size = 19),
                              legend.position = "",
                              legend.title = element_blank(),
                              legend.background = element_rect(fill = "white", color = "black"),
                              legend.key.width = unit(3, "line"),
                              legend.key.height = unit(1, "line"),
                              legend.text = element_text(size = 15)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
  geom_line(aes(linetype = Method, size = Method), position = pd) +
  geom_point(position = pd, aes(shape = Method), size = 2) +
  scale_shape_manual(values = c(4, 4, 9, 9, 1, 1)) +
  scale_linetype_manual(values = c(2, 1, 2, 1, 2, 1)) +
  scale_size_manual(values = c(1, 1, 1, 1, 1, 1)) +
  scale_color_manual(values = c("#9ecae1", "#3182bd", "#a1d99b", "#31a354", 
                                "#bdbdbd", "#636363")) +
  labs(title = func_name, x = "r", y = expression(log[10](MSE))) +
  theme(plot.title = element_text(hjust = 0.5, size = 19))
plot_mse[[func_type]] <- p3
print(p3)

# PMSE
p4 <- ggplot(pmse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
p4 <- p4 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 19),
                              axis.title = element_text(size = 19),
                              legend.position = "right",
                              legend.title = element_blank(),
                              legend.background = element_rect(fill = "white", color = "black"),
                              legend.key.width = unit(3, "line"),
                              legend.key.height = unit(1, "line"),
                              legend.text = element_text(size = 15,
                                                         margin = margin(t = 3, b = 3))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = width, position = pd, show.legend = FALSE) +
  geom_line(aes(linetype = Method, size = Method), position = pd) +
  geom_point(position = pd, aes(shape = Method), size = 2) +
  scale_shape_manual(values = c(4, 4, 9, 9, 1, 1)) +
  scale_linetype_manual(values = c(2, 1, 2, 1, 2, 1)) +
  scale_size_manual(values = c(1, 1, 1, 1, 1, 1)) +
  scale_color_manual(values = c("#9ecae1", "#3182bd", "#a1d99b", "#31a354", 
                                "#bdbdbd", "#636363")) +
  labs(title = func_name, x = "r", y = expression(log[10](PMSE))) +
  theme(plot.title = element_text(hjust = 0.5, size = 19))
plot_pmse[[func_type]] <- p4
print(p4)

pdf("simu_ada.pdf", width = 10.5, height = 4)
grid.arrange(plot_mse[[4]], plot_pmse[[4]],
             nrow = 1, ncol = 2, widths = c(4, 6.5))
dev.off()
