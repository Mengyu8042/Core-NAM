######### Standardize data to [0,1] ######### 
min_max_scale <- function(df) {
  dim <- ncol(df)
  for (j in 1:dim) {
    df[, j] <- (df[, j] - min(df[, j]))/(max(df[, j]) - min(df[, j]))
  }
  return(df)
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

######### Copula transformation #########
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

######### Fitting functions #########
am.init <- function(D_list, sp) {
  ## create the augmented penalty matrix
  sp <- sqrt(sp)
  nD <- nrow(D_list[[1]]) + nrow(D_list[[2]]) + nrow(D_list[[3]])
  DD <- matrix(rep(0, nD), nrow = nD)
  # Create a list of scaled D matrices
  scaled_D_list <- list(D1 = sp[1] * D_list[[1]],
                        D2 = sp[2] * D_list[[2]],
                        D3 = sp[3] * D_list[[3]] + sp[4] * D_list[[4]])
  # Bind the column of zeros to the block-diagonal matrix
  DD <- cbind(DD, as.matrix(do.call(bdiag, scaled_D_list)))
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
