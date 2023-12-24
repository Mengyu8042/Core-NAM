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
  } else if (func_type == 7) {  # 12-d (with interaction)
    f1 <- function(x) x
    f2 <- function(x) (2*x - 1)^2
    f3 <- function(x) sin(2*pi*x) / (2 - sin(2*pi*x))
    f4 <- function(x) 0.1*sin(2*pi*x) + 0.2*cos(2*pi*x) + 0.3*sin(2*pi*x)^2 +
      0.4*cos(2*pi*x)^3 + 0.5*sin(2*pi*x)^3
    y_raw <- f1(X[,1]) + f2(X[,2]) + f3(X[,3]) + f4(X[,4]) +
      1.5*f1(X[,5]) + 1.5*f2(X[,6]) + 1.5*f3(X[,7]) + 1.5*f4(X[,8]) +
      2*f1(X[,9]) + 2*f2(X[,10]) + 2*f3(X[,11]) + 2*f4(X[,12])
  } else if (func_type == 8) {
    y_raw <- 0
    for (j in 1:18) {
      y_raw = y_raw + 10^4 * X[, j]^11 * (1 - X[, j])^6
    }
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
  
  list(X = X, D1 = D1, D2 = D2)
}

######### Fitting functions #########
am.init <- function(D_list, sp) {
  ## create the augmented penalty matrix
  sp <- sqrt(sp)
  nD <- sum(sapply(D_list, nrow)) # Total number of rows for all D matrices
  DD <- matrix(rep(0, nD), nrow = nD)
  # Create a list of scaled D matrices
  scaled_D_list <- lapply(seq_along(D_list), function(i) sp[i] * D_list[[i]])
  # Bind the column of zeros to the block-diagonal matrix
  DD <- cbind(DD, as.matrix(do.call(bdiag, scaled_D_list)))
  return(DD)
}

am.init2 <- function(D_list, sp) {
  ## create the augmented penalty matrix for models with interactions
  sp <- sqrt(sp)
  DD <- 0
  for (i in 1:length(D_list)) {
    DD <- DD + sp[i] * D_list[[i]]
  }
  nD <- nrow(DD) 
  DD <- cbind(rep(0, nD), DD) 
  return(DD)
}

am.fit <- function(y, XX, D_list, sp) {
  if (ncol(XX) == ncol(D_list[[1]])) {
    DD <- am.init2(D_list, sp)
  } else {
    DD <- am.init(D_list, sp)
  }
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
  if (ncol(XX) == ncol(D_list[[1]])) {
    DD <- am.init2(D_list, sp)
  } else {
    DD <- am.init(D_list, sp)
  }
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
  if (ncol(XXrow) == ncol(D_list[[1]])) {
    DD <- am.init2(D_list, sp)
  } else {
    DD <- am.init(D_list, sp)
  }
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
  if (ncol(XXrow) == ncol(D_list[[1]])) {
    DD <- am.init2(D_list, sp)
  } else {
    DD <- am.init(D_list, sp)
  }
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
  if (ncol(XXrow) == ncol(D_list[[1]])) {
    DD <- am.init2(D_list, sp)
  } else {
    DD <- am.init(D_list, sp)
  }
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
