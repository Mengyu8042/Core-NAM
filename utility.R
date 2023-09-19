######### Generate grid points #########
grid_points <- function(n, d) {
  ## generate n grid points in [0,1]^d
  n_d <- ceil(n^(1/d))
  x <- seq(0, 1, length.out = n_d)
  gp <- expand.grid(replicate(d, x, simplify = FALSE))
  gp <- gp[sample(n_d^d, n, replace = FALSE), ]
  gp <- as.matrix(gp)
  return(gp)
}

######### Standardize data to [0,1] ######### 
min_max_scale <- function(df) {
  dim <- ncol(df)
  for (j in 1:dim) {
    df[, j] <- (df[, j] - min(df[, j]))/(max(df[, j]) - min(df[, j]))
  }
  return(df)
}

######### Make a copula transformation #########
library(randtoolbox)
copula.trans <- function(dendat) {
  if (is.null(dim(dendat))) {
    n <- length(dendat)
    dendat <- matrix(dendat, n, 1)
  }
  n <- nrow(dendat)
  d <- ncol(dendat)
  copdat <- dendat
  
  for (ii in 1:d) {
    or <- rank(dendat[, ii], ties.method = "random")
    mones <- rep(0, n)
    for (i in 1:n) mones[or[i]] <- i  
    copdat[, ii] <- mones/(n + 1)
  }
  if (d == 1) {
    copdat <- copdat[, 1]
  }
  return(copdat)
}

######### Generate the response #########
yGenerate <- function(X, func_type = 1) {
  if (func_type == 1) {  # "blocks"
    x <- X[, 1]
    n <- length(x)
    t <- c(0.15, 0.3, 0.55, 0.8)
    h <- c(4., -4, 2, -5)
    y_raw <- rep(0, n)
    for (i in seq(1, length(h))) {
      y_raw <- y_raw + (h[i] * (1 + sign(x - t[i]))) / 2
    }  
  } else if (func_type == 2) {  # "bumps"
    x <- X[, 1]
    n <- length(x)
    t <- c(0.1, 0.16, 0.45, 0.65, 0.85)
    h <- c(4, 2, 3, 4, 5)
    w <- c(0.04, 0.02, 0.04, 0.06, 0.05)
    y_raw <- rep(0, n)
    for (j in 1:length(t)) {
      y_raw <- y_raw + h[j] / (1 + abs((x - t[j])/w[j]))^4
    }
  } else if (func_type == 3) {  # "heavi"
    x <- X[, 1]
    n <- length(x)
    y_raw <- 4 * sin(4 * pi * x) - sign(x -  0.3) - sign(0.72 - x)
  } else if (func_type == 4) {  # "doppler"
    x <- X[, 1]
    n <- length(x)
    eps <- 0.2
    y_raw <- sqrt(x * (1 - x)) * sin((2 * pi * (1 + eps))/ (x + eps))
  } else if (func_type == 5) {  # 2-d (with 1 interaction)
    x1 <- X[, 1]; x2 <- X[, 2]
    s1 <- 0.3; s2 <- 0.4
    y_raw <- 0.75/(s1 * s2 * pi) * exp(-(x1 - 0.2)^2/s1^2 - (x2 - 0.3)^2/s2^2) +
      0.45/(s1 * s2 * pi) * exp(-(x1 - 0.7)^2/s1^2 - (x2 - 0.8)^2/s2^2)
  } else if (func_type == 6) {  # 4-d (without interaction)
    x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]; x4 <- X[, 4]
    y1 <- exp(2 * x1)
    y2 <- sin((2 * pi * (1 + 0.2)) / (x2 + 0.2))
    y3 <- sin((2 * pi * (1 + 0.1)) / (x3 + 0.1))
    y4 <- 10^5 * x4^11 * (1 - x4)^6 + 10^3 * x4^3 * (1 - x4)^10
    y_raw <- y1 + y2 + y3 + y4
  }
  return(y_raw)
}

yGenerate2 <- function(X, func_type = 1) {
  if (func_type == 1) {
    x1 <- X[, 1]
    y_raw <- sin(3 * pi * x1)
  } else if (func_type == 2) {
    x1 <- X[, 1]
    y_raw <- exp(2 * x1^2)
  } else if (func_type == 3) {
    x1 <- X[, 1]
    y_raw <- 10^6 * x1^11 * (1 - x1)^6 + 10^4 * x1^3 * (1 - x1)^10
  } else if (func_type == 4) {  # 2-d (with 1 interaction)
    x1 <- X[, 1]; x2 <- X[, 2]
    y_raw <- 0.75/(0.3 * 0.4 * pi) * exp(-(x1 - 0.2)^2/0.3^2 - (x2 - 0.3)^2/0.4^2) +
      0.45/(0.3 * 0.4 * pi) * exp(-(x1 - 0.7)^2/0.3^2 - (x2 - 0.8)^2/0.4^2)
  } else if (func_type == 5) {  # 2-d (with 1 interaction)
    x1 <- X[, 1]; x2 <- X[, 2]
    y_raw <- 1.9 * (1.45 + exp(x1) * sin(13*(x1-0.6)^2)) * exp(-x2) * sin(7 * x2)
  } 
  return(y_raw)
}

######### Generate the design matrix X and square root penalty matrix D #########
kernel_generate <- function(x, k, id, type_bs = "cr") {
  ## generate univariate spline bases
  ## x: univariate variable 
  ## k: number of knots 
  ## id: index set 
  ## type_bs: type of splines ("cr" or "ps" or "linear" or "cc")
  xk <- data.frame(seq(min(x), max(x), length = k))
  dat <- data.frame(x)
  if (type_bs == "cr") {
    sspec <- s(x, bs = "cr", k = k)  # cubic regression spline
  } else if (type_bs == "ps") {
    sspec <- s(x, bs = "ps", m = c(2, 2), k = k)  # P-spline (order = 3, cubic spline)
  } else if (type_bs == "linear") {
    sspec <- s(x, bs = "ps", m = c(0, 2), k = k)  # P-spline (order = 1, linear basis)
  } else if (type_bs == "cc") {
    sspec <- s(x, bs = "cc", k = k)  # cyclic cubic regression spline
  }
  smv <- smoothCon(sspec, data = dat, knots = xk)[[1]]
  S <- smv$S[[1]]
  D <- t(mroot(S, method = "svd"))

  dat_pred <- data.frame(x[id])
  names(dat_pred) <- "x"
  X <- PredictMat(smv, dat_pred)
  
  list(X = X, D = D)
}

kernel_generate_tp <- function(x1, x2, k, id, type_bs = "ti") {
  ## generate tensor product interactions
  dim1 <- seq(min(x1), max(x1), length = k)
  dim2 <- seq(min(x2), max(x2), length = k)
  xk <- data.frame(expand.grid(dim1, dim2))
  dat <- data.frame(x1, x2)
  if (type_bs == "ti") {
    sspec <- ti(x1, x2, k = k, bs = "ps")
  } else if (type_bs == "te") {
    sspec <- te(x1, x2, k = k, bs = "ps")
  }
  smv <- smoothCon(sspec, data = dat, knots = xk)[[1]]
  S1 <- smv$S[[1]]
  S2 <- smv$S[[2]]
  D1 <- t(mroot(S1, method = "svd"))
  D2 <- t(mroot(S2, method = "svd"))

  dat_pred <- data.frame(x1[id], x2[id])
  names(dat_pred) <- c("x1", "x2")
  X <- PredictMat(smv, dat_pred)
  
  list(X = X, D = D1 + D2)
}

######### Fitting functions #########
am.init <- function(D_list, sp) {
  ## create the augmented model matrix...
  sp <- sqrt(sp)
  dim <- length(D_list)
  if (dim == 1) {
    D1 <- D_list[[1]]
    nD <- nrow(D1)
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1)))
  } else if (dim == 2) {
    D1 <- D_list[[1]]; D2 <- D_list[[2]]
    nD <- nrow(D1) + nrow(D2) 
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1, sp[2] * D2)))
  } else if (dim == 3) {
    D1 <- D_list[[1]]; D2 <- D_list[[2]]; D3 <- D_list[[3]]
    nD <- nrow(D1) + nrow(D2) + nrow(D3)
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1, sp[2] * D2, sp[3] * D3)))
  } else if (dim == 4) {
    D1 <- D_list[[1]]; D2 <- D_list[[2]]; D3 <- D_list[[3]]; D4 <- D_list[[4]]
    nD <- nrow(D1) + nrow(D2) + nrow(D3) + nrow(D4)
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1, sp[2] * D2, sp[3] * D3, sp[4] * D4)))
  } else if (dim == 5) {
    D1 <- D_list[[1]]; D2 <- D_list[[2]]; D3 <- D_list[[3]]; D4 <- D_list[[4]]; D5 <- D_list[[5]]
    nD <- nrow(D1) + nrow(D2) + nrow(D3) + nrow(D4) + nrow(D5)
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1, sp[2] * D2, sp[3] * D3,
                                            sp[4] * D4, sp[5] * D5)))
  } else if (dim == 6) {
    D1 <- D_list$D1; D2 <- D_list$D2; D3 <- D_list$D3; D4 <- D_list$D4; D5 <- D_list$D5; D6 <- D_list$D6
    nD <- nrow(D1) + nrow(D2) + nrow(D3) + nrow(D4) + nrow(D5) + nrow(D6)
    DD <- cbind(rep(0, nD), as.matrix(bdiag(sp[1] * D1, sp[2] * D2, sp[3] * D3,
                                            sp[4] * D4, sp[5] * D5, sp[6] * D6)))
  }
  return(DD)
}

am.fit <- function(y, XX, D_list, sp) {
  DD <- am.init(D_list, sp)
  Y <- c(y, rep(0, nrow(DD)))
  Z <- rbind(cbind(rep(1, nrow(XX)), XX), DD)
  
  ## fit model...
  b <- coef(lm(Y ~ Z - 1, na.action = na.exclude))
  b[is.na(b)] <- 0
  list(b = b)
}

am.mse <- function(y_raw, XX, b, metric = "MSE") {
  ## compute MSE or MAE...
  n <- length(y_raw)
  rsd <- y_raw - b[1] - as.vector(XX %*% b[2:(ncol(XX) + 1)])  ## residuals
  if (metric == "MSE") {
    return(sum(rsd^2)/n)  # MSE
  } else if (metric == "MAE") {
    return(sum(abs(rsd))/n)  # MAE
  }
}

am.gcv_score <- function(y, XX, D_list, sp, b) {
  ## compute GCV score...
  DD <- am.init(D_list, sp)
  Z <- rbind(cbind(rep(1, nrow(XX)), XX), DD)
  n <- length(y)
  
  ## compute trace of hat matrix...
  N <- nrow(Z)
  z <- base::qr(Z)
  hat <- rowSums(qr.qy(z, diag(1, nrow = N, ncol = z$rank))^2)
  trA <- sum(hat[1:n])
  
  y_fit <- Z %*% b
  rsd <- y - y_fit[1:n]  # residuals
  rss <- sum(rsd^2)  # residual SS
  gcv <- rss * n/(n - trA)^2  # GCV score
  list(gcv = gcv)
}

am.gcv <- function(lsp, y, XX, D_list) {
  ## function suitable for GCV optimization by optim
  ## return GCV-score from above functions
  b <- am.fit(y, XX, D_list, exp(lsp))$b
  am.gcv_score(y, XX, D_list, exp(lsp), b)$gcv
}

######### Fitting functions for Core-NAM #########
am.core.fit <- function(ysub, XXrow, XXsub, D_list, sp) {
  DD <- am.init(D_list, sp)
  Y <- c(ysub, rep(0, nrow(DD)))
  Z <- rbind(cbind(rep(1, nrow(XXrow)), XXrow), DD)
  Zsub <- rbind(cbind(rep(1, nrow(XXsub)), XXsub), DD)
  Zsub1t <- Matrix(t(Zsub), sparse = TRUE)
  
  ## fit model...
  b <- eigenMultSolveSp(Z, Zsub1t, Y)
  list(b = b)
}

am.core.gcv_score <- function(ysub, XXrow, XXsub, D_list, sp, b) {
  ## compute GCV score...
  n <- length(ysub)
  DD <- am.init(D_list, sp)
  Z <- rbind(cbind(rep(1, nrow(XXrow)), XXrow), DD)
  Zsub <- rbind(cbind(rep(1, nrow(XXsub)), XXsub), DD)
  Zsub1t <- Matrix(t(Zsub), sparse = TRUE)
  
  ## compute trace of hat matrix...
  trA <- eigenTraceHatSp(Z, Zsub1t, n)

  y_fit <- Z %*% b
  rsd <- ysub - y_fit[1:n]  # residuals
  rss <- sum(rsd^2)  # residual SS
  gcv <- rss * n/(n - trA)^2  # GCV score
  list(gcv = gcv)
}

am.core.gcv <- function(lsp, ysub, XXrow, XXsub, D_list) {
  ## function suitable for GCV optimization by optim
  ## return GCV-score from above functions
  b <- am.core.fit(ysub, XXrow, XXsub, D_list, exp(lsp))$b
  am.core.gcv_score(ysub, XXrow, XXsub, D_list, exp(lsp), b)$gcv
}

am.core.mse <- function(y_raw_sub, ysub, XXrow, XXsub, D_list, sp) {
  ## return the value of loss function
  DD <- am.init(D_list, sp)
  Y <- c(ysub, rep(0, nrow(DD)))
  Z <- rbind(cbind(rep(1, nrow(XXrow)), XXrow), DD)
  Zsub <- rbind(cbind(rep(1, nrow(XXsub)), XXsub), DD)
  Zsub1t <- Matrix(t(Zsub), sparse = TRUE)
  
  ## fit model..
  b <- eigenMultSolveSp(Z, Zsub1t, Y)
  
  ## compute mse...
  rsd <- y_raw_sub - b[1] - as.vector(XXrow %*% b[2:(ncol(XX)+1)])
  mse <- sum(rsd^2) / length(y_raw_sub)
  list(b = b, mse = mse)
}

am.core.loss <- function(lsp, y_raw_sub, ysub, XXrow, XXsub, D_list) {
  ## function suitable for optimization by optim
  am.core.mse(y_raw_sub, ysub, XXrow, XXsub, D_list, exp(lsp))$mse
}
