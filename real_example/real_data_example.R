# Real data example: total column ozone dataset
rm(list = ls())
gc()

## Import R Packages ##
list.of.packages <- c("mgcv", "MASS", "mvtnorm", "corpcor", "SparseM", "pracma", 
                      "Rfast", "Rcpp", "RcppEigen", "Matrix", "nabor", "lhs",
                      "ReIns", "ggplot2", "reshape2", "patchwork", "gridExtra",
                      "caret", "GGally", "ncdf4")
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
library(caret)
library(GGally)
library(ncdf4)

## Import the source file and utility file ##
setwd(getwd())
sourceCpp("../source/funcs.cpp")
source("utility.R")

## Import TCO data ## 
print("Load data...")
year_list <- 1978:2019
for (year in year_list) {
  file_name <- paste0("../data/NIWA-BS_CombinedTCO_V3.4.1_", year, "_Daily_Unpatched.nc")
  nc_data <- nc_open(file_name)
  # read the variable values
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude")
  doy <- 1:length(ncvar_get(nc_data, "time"))
  tco.array <- ncvar_get(nc_data, "tco")
  # view the fill values of missing data
  fillvalue <- ncatt_get(nc_data, "tco", "_FillValue")
  nc_close(nc_data)
  tco.array[tco.array == fillvalue$value] <- NA
  # convert tco.array to an n*p observation matrix
  dat_mat <- matrix(NA, prod(dim(tco.array)), 5)
  dat_mat[, 1] <- rep(year, prod(dim(tco.array)))
  dat_mat[, 2] <- rep(doy, each = dim(tco.array)[1] * dim(tco.array)[2])
  dat_mat[, 3] <- rep(rep(lat, each = dim(tco.array)[1]), times = dim(tco.array)[3])
  dat_mat[, 4] <- rep(rep(lon, times = dim(tco.array)[2]), times = dim(tco.array)[3])
  dat_mat[, 5] <- as.vector(tco.array)
  
  if (year == year_list[1]) {
    df0 <- dat_mat
  } else {
    df0 <- rbind(df0, dat_mat)
  }
}
rm(nc_data, tco.array, dat_mat, year, doy, lat, lon)
# remove NA
df0[is.na(df0) | is.nan(df0) | is.infinite(df0)] <- NA
df0 <- na.omit(df0)  
df0 <- as.data.frame(df0)
colnames(df0) <- c("year", "doy", "lat", "lon", "y")
df0$y <- log(df0$y + 1)
# standardization
df0[, 1:4] <- min_max_scale(df0[, 1:4])
df0[, 1:4] <- copula.trans(df0[, 1:4])
print("Data is prepared!")


## Beginning ##
ozone_kern_gen <- function(year, doy, lat, lon, id, k, k_ti) {
  # calculate design matrix and penalty matrix
  kg1 <- kernel_generate(year, k, id, type_bs = "cr")
  kg2 <- kernel_generate(doy, k, id, type_bs = "cc")
  kg3 <- kernel_generate_tp(lat, lon, k_ti, id, type_bs = "te")
  XX <- cbind(kg1$X, kg2$X, kg3$X) 
  D_list <- list(D1 = kg1$D, D2 = kg2$D, D3 = kg3$D1, D4 = kg3$D2)
  list(XX = XX, D_list = D_list)
}

## Initialization ##
set.seed(1234)
use_bam <- TRUE
n_test <- floor(0.01 * nrow(df0))  # testing set size 
n <- nrow(df0) - n_test  # training set size 
k <- 20  # number of knots for univariate variables
k_ti <- 6  # number of knots for each marginal of 2D interactions
nloop <- 10  # number of replicates
num_method <- 4 
sub_meta <- (1:5) * 200  # subsample size
mse_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
mse_temp <- matrix(0, num_method, length(sub_meta))
mae_meta <- array(0, dim = c(nloop, num_method, length(sub_meta)))
mae_temp <- matrix(0, num_method, length(sub_meta))

## Start Calculation ##
for (i in 1:nloop) {
  print(paste0("Replication ", i, "/", nloop))
  set.seed(100 + 123 * i)
  
  # training-test sets partition
  id <- sample(1:nrow(df0), n_test, replace = FALSE)
  df_test <- df0[id, ]
  df <- df0[-id, ]
  rm(id)

  res <- ozone_kern_gen(df_test$year, df_test$doy, df_test$lat, df_test$lon,
                        1:n_test, k, k_ti)
  XX_t <- res$XX; D_list <- res$D_list
  dim <- length(D_list)
  
  ## FULL ##
  print(" FULL")
  fit <- bam(y ~ s(year, k = k, bs = "cr") + s(doy, k = k, bs = "cc") + 
               te(lat, lon, k = k_ti, bs = "ps"), 
             data = df, method = "GCV.Cp")
  pred <- predict(fit, df_test)
  mse <- sum((df_test$y - pred)^2)/n_test
  mae <- sum(abs(df_test$y - pred))/n_test
  rm(fit, pred)
  
  ## Subsampling ##
  for (j in 1:length(sub_meta)) {
    set.seed(100 + 432 * j + 123 * i)
    sub <- sub_meta[j]
    print(paste0(" Subsample size r ", sub))
    
    ## UNIF ##
    print("  UNIF")
    id_unif <- sort(sample(n, sub, replace = TRUE))
    ysub <- df$y[id_unif]
    
    if (use_bam == TRUE) {
      df_unif <- df[id_unif, ]
      fit <- bam(y ~ s(year, k = k, bs = "cr") + s(doy, k = k, bs = "cc") + 
                   te(lat, lon, k = k_ti, bs = "ps"), 
                 data = df_unif, method = "GCV.Cp")
      pred <- predict(fit, df_test)
      mse_unif <- sum((df_test$y - pred)^2)/n_test
      mae_unif <- sum(abs(df_test$y - pred))/n_test
      rm(id_unif, ysub, df_unif, fit, pred)
    } else {
      df_unif <- df[id_unif, ]
      XXsub <- ozone_kern_gen(df_unif$year, df_unif$doy, df_unif$lat, 
                              df_unif$lon, 1:nrow(df_unif), k, k_ti)$XX
      gcv <- optim(rep(1, dim), am.gcv, y = ysub, XX = XXsub, D_list = D_list,
                   control = list(maxit = 200, warn.1d.NelderMead = FALSE))
      sp_unif <- exp(gcv$par)
      fit <- am.fit(ysub, XXsub, D_list, sp_unif)
      b_unif <- fit$b
      mse_unif <- am.mse(df_test$y, XX_t, b_unif)
      mae_unif <- am.mse(df_test$y, XX_t, b_unif, "MAE")
      rm(id_unif, XXsub, ysub, gcv, sp_unif, fit, b_unif)
    }
    
    ## LowCon ##
    print("  LowCon")
    theta <- 1
    x_standardize <- function(x) {
      x <- (x - min(x))/(max(x) - min(x))
      return(x)
    }
    lhd_data <- cbind(df$year, df$doy, df$lat, df$lon)
    lhd_data <- apply(lhd_data, 2, x_standardize)
    lowcon_space <- apply(lhd_data, 2, quantile, probs = c(theta/100, 1 - theta/100))
    des <- randomLHS(sub, ncol(lhd_data))
    for (dd in 1:ncol(lhd_data)) {
      des[, dd] <- (des[, dd]) * (lowcon_space[2, dd] - lowcon_space[1, dd]) + lowcon_space[1, dd]
    }
    id_lowcon <- nabor::knn(lhd_data, des, k = 1, eps = 0.1)$nn.idx[, 1]
    ysub <- df$y[id_lowcon]
    rm(lhd_data, lowcon_space, des)
    
    if (use_bam == TRUE) {
      df_lowcon <- df[id_lowcon, ]
      fit <- bam(y ~ s(year, k = k, bs = "cr") + s(doy, k = k, bs = "cc") + 
                   te(lat, lon, k = k_ti, bs = "ps"), 
                 data = df_lowcon, method = "GCV.Cp")
      pred <- predict(fit, df_test)
      mse_lowcon <- sum((df_test$y - pred)^2)/n_test 
      mae_lowcon <- sum(abs(df_test$y - pred))/n_test 
      rm(id_lowcon, ysub, df_lowcon, fit, pred)
    } else {
      df_lowcon <- df[id_lowcon, ]
      XXsub <- ozone_kern_gen(df_lowcon$year, df_lowcon$doy, df_lowcon$lat, 
                              df_lowcon$lon, 1:nrow(df_lowcon), k, k_ti)$XX
      gcv <- optim(rep(1, dim), am.gcv, y = ysub, XX = XXsub, D_list = D_list,
                   control = list(maxit = 200, warn.1d.NelderMead = FALSE))
      sp_lowcon <- exp(gcv$par)
      fit <- am.fit(ysub, XXsub, D_list, sp_lowcon)
      b_lowcon <- fit$b
      mse_lowcon <- am.mse(df_test$y, XX_t, b_lowcon)
      mae_lowcon <- am.mse(df_test$y, XX_t, b_lowcon, "MAE")
      rm(id_lowcon, XXsub, ysub, gcv, sp_lowcon, fit, b_lowcon)
    }
    
    ## Core-NAM ##
    print("  CORE")
    # select the core-elements
    n2 <- min(n, 20 * sub)
    id_reduc <- sort(sample(n, n2, replace = FALSE))
    df_reduc <- df[id_reduc, ]
    XXreduc <- ozone_kern_gen(df_reduc$year, df_reduc$doy, df_reduc$lat, 
                              df_reduc$lon, 1:nrow(df_reduc), k, k_ti)$XX
    p <- ncol(XXreduc)
    XXsub <- matrix(0, n2, p)
    thres <- colnth(abs(XXreduc), rep(sub, p), descending = TRUE, parallel = TRUE)
    for (kk in 1:p) {
      col <- abs(XXreduc[, kk])
      id <- which(col > thres[kk])
      id <- c(id, which(col == thres[kk])[1:(sub - length(id))])
      XXsub[id, kk] <- XXreduc[id, kk]
    }
    id_nz <- (rowSums(abs(XXsub)) > 0)
    XXrow <- XXreduc[id_nz, ]
    XXsub <- XXsub[id_nz, ]
    ysub <- df$y[id_reduc][id_nz]
    rm(id_reduc, XXreduc, thres, col, id, id_nz)
    
    # choose the smoothing parameter via core-elements GCV
    gcv <- optim(rep(1, dim), am.core.gcv, ysub = ysub, XXrow = XXrow,
                 XXsub = XXsub, D_list = D_list,
                 control = list(maxit = 200, warn.1d.NelderMead = FALSE))
    sp_core <- exp(gcv$par)
    
    # core-elements estimation based on the chosen smoothing parameter
    fit <- am.core.fit(ysub, XXrow, XXsub, D_list, sp_core)
    b_core <- fit$b
    mse_core <- am.mse(df_test$y, XX_t, b_core)  
    mae_core <- am.mse(df_test$y, XX_t, b_core, "MAE") 
    rm(XXsub, XXrow, ysub, gcv, sp_core, fit, b_core)
    
    ######################################################
    mse_temp[, j] <- c(mse, mse_unif, mse_lowcon, mse_core)
    mae_temp[, j] <- c(mae, mae_unif, mae_lowcon, mae_core)
    
    print("  mse_full, mse_unif, mse_lowcon, mse_core")
    print(mse_temp[, j])
    print("  mae_full, mae_unif, mae_lowcon, mae_core")
    print(mae_temp[, j])
  }
  
  rm(df, df_test, XX_t, D_list)
  mse_meta[i, , ] <- mse_temp
  mae_meta[i, , ] <- mae_temp
}

FULL.mse <- log10(mse_meta[, 1, ])
UNIF.mse <- log10(mse_meta[, 2, ])
LOWCON.mse <- log10(mse_meta[, 3, ])
CORE.mse <- log10(mse_meta[, 4, ])

FULL.mae <- log10(mae_meta[, 1, ])
UNIF.mae <- log10(mae_meta[, 2, ])
LOWCON.mae <- log10(mae_meta[, 3, ])
CORE.mae <- log10(mae_meta[, 4, ])

mse_mat <- data.frame(mean = c(colMeans(FULL.mse, na.rm = TRUE),
                               colMeans(UNIF.mse, na.rm = TRUE), 
                               colMeans(LOWCON.mse, na.rm = TRUE), 
                               colMeans(CORE.mse, na.rm = TRUE)),
                      sd = c(sqrt(colVars(FULL.mse, na.rm = TRUE)),
                             sqrt(colVars(UNIF.mse, na.rm = TRUE)), 
                             sqrt(colVars(LOWCON.mse, na.rm = TRUE)), 
                             sqrt(colVars(CORE.mse, na.rm = TRUE))),
                      Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                          each = length(sub_meta))),
                      sub = rep(sub_meta, num_method))
mse_mat$Method <- factor(mse_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))

mae_mat <- data.frame(mean = c(colMeans(FULL.mae, na.rm = TRUE),
                               colMeans(UNIF.mae, na.rm = TRUE), 
                               colMeans(LOWCON.mae, na.rm = TRUE), 
                               colMeans(CORE.mae, na.rm = TRUE)),
                      sd = c(sqrt(colVars(FULL.mae, na.rm = TRUE)), 
                             sqrt(colVars(UNIF.mae, na.rm = TRUE)), 
                             sqrt(colVars(LOWCON.mae, na.rm = TRUE)), 
                             sqrt(colVars(CORE.mae, na.rm = TRUE))),
                      Method = factor(rep(c("FULL", "UNIF", "LowCon", "CORE"),
                                          each = length(sub_meta))),
                      sub = rep(sub_meta, num_method))
mae_mat$Method <- factor(mae_mat$Method, levels = c("FULL", "UNIF", "LowCon", "CORE"))

## Plot the results ##
# MSE
width <- (max(sub_meta) - min(sub_meta))/10
pd <- position_dodge(width/3)
p3 <- ggplot(mse_mat, aes(x = sub, y = mean, group = Method, colour = Method))
p3 <- p3 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 18),
                              axis.title = element_text(size = 18),
                              legend.position = "",
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
  labs(x = "r", y = expression(log[10](PMSE))) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))
print(p3)

# MAE
width <- (max(sub_meta) - min(sub_meta))/10
pd <- position_dodge(width/3)
p4 <- ggplot(mae_mat, aes(x = sub, y = mean, group = Method, colour = Method))
p4 <- p4 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black"),
                              axis.text = element_text(size = 18),
                              axis.title = element_text(size = 18),
                              legend.position = c(0.72, 0.8),
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
  labs(x = "r", y = expression(log[10](PMAE))) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))
print(p4)

pdf("ozone_pred_mse_mae.pdf", width = 8, height = 3.5)
grid.arrange(p3, p4, nrow = 1, ncol = 2)
dev.off()
