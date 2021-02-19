initialize_mcmc <- function(Y, tt, nt, B) {
  nsub <- nrow(Y) / nt
  nreg <- ncol(Y)
  L <- list()
  omega <- numeric(nreg)
  Yt <- matrix(0, nrow = nt, ncol = nreg * nsub)
  Ys <- matrix(0, nrow = nt, ncol = nreg * nsub)
  Yss <- matrix(0, nrow = nt * nsub, ncol = nreg)
  for (i in 1:nsub) {
    Yt[, ((i - 1) * nreg + 1):(i * nreg)] <- Y[((i - 1) * nt + 1):(i * nt),]
  }
  nr <- nsub * nreg
  Ly <- lapply(1:nr, function(i) Yt[,i])
  Lt <- lapply(1:nr, function(i) tt)
  fpca1 <- FPCA(Ly, Lt, optns = list(dataType = 'Dense', FVEthreshold = .99))
  for (i in 1:nr) {
    Ys[, i] <- fpca1$mu + fpca1$phi[,1:fpca1$selectK] %*% fpca1$xiEst[i,1:fpca1$selectK]
  }
  for (r in 1:nreg) {
    r_ind <- seq(from = r, length.out = nsub, by = nreg)
    for (i in 1:nsub) {
      Yss[, r] <- c(Ys[, r_ind])
    }
    omega[r] <- 1 / var(Y[,r] - Yss[,r], na.rm = TRUE)
  }
  
  svd1 <- svd(t(Ys))
  psi <- matrix(0, nt, fpca1$selectK)
  phi <- matrix(0, nreg * fpca1$selectK, nreg)
  xi <- matrix(0, nrow = nreg * nsub, ncol = fpca1$selectK)
  for (i in 1:nsub) {
    for (r in 1:nreg) {
      xi[(i - 1) * nreg + r, ] <- svd1$u[(i - 1) * nreg + r,1:fpca1$selectK] %*% diag(svd1$d[1:fpca1$selectK])
    }
  }
  lambda <- matrix(0, nrow = ncol(B), fpca1$selectK)
  eta <- matrix(0, nrow = nsub * nreg, ncol = fpca1$selectK)
  for (l in 1:fpca1$selectK) {
    psi[,l] <- svd1$v[,l]
    lambda[,l] <- solve(t(B) %*% B, t(B) %*% psi[,l])
    rxi <- t(matrix(xi[,l], nrow = nreg))
    svd2 <- svd(rxi)
    for (r in 1:nreg) {
      phi[((l - 1) * nreg + 1):(l * nreg),] <- svd2$v
      for (i in 1:nsub) {
        eta[(i - 1) * nreg + r, l] <- svd2$u[i, r] * svd2$d[r]
      }
    }
  }
  
  list(eta = eta, phi = phi, lambda = lambda, omega = omega)
}