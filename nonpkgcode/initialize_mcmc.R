library(fdapace)
initialize_mcmc <- function(Y, tt, nt, B, X, ldim = NULL, cumfve = NULL) {
  nsub <- nrow(Y) / nt
  nreg <- ncol(Y)
  L <- list()
  omega <- numeric(nreg)
  nr <- nsub * nreg
  Yt <- matrix(0, nrow = nt, ncol = nr)
  Ys <- matrix(0, nrow = nt, ncol = nr)
  Ysm <- matrix(0, nrow = nt, ncol = nr)
  Yss <- matrix(0, nrow = nt * nsub, ncol = nreg)
  for (i in 1:nsub) {
    Yt[, ((i - 1) * nreg + 1):(i * nreg)] <- Y[((i - 1) * nt + 1):(i * nt),]
  }
  Ly <- lapply(1:nr, function(i) Yt[,i])
  Lt <- lapply(1:nr, function(i) tt)
  fpca1 <- list()
  if (!is.null(cumfve)) {
    fpca1 <- FPCA(Ly, Lt, optns = list(dataType = 'Dense', FVEthreshold = cumfve))
  } else if (!is.null(ldim)) {
    fpca1 <- FPCA(Ly, Lt, optns = list(dataType = 'Dense',  methodSelectK = ldim))
  } else {
    fpca1 <- FPCA(Ly, Lt, optns = list(dataType = 'Dense',  FVEthreshold = .99))
  }
  for (i in 1:nr) {
    Ys[, i] <- fpca1$mu + fpca1$phi[,1:fpca1$selectK] %*% fpca1$xiEst[i,1:fpca1$selectK]
    Ysm[,i] <- Yt[, i] - fpca1$mu
  }
  pcomp1 <- prcomp(t(Ys))
  # pcomp2 <- prcomp(t(Yt))
  for (r in 1:nreg) {
    r_ind <- seq(from = r, length.out = nsub, by = nreg)
    for (i in 1:nsub) {
      Yss[, r] <- c(Ys[, r_ind])
    }
    omega[r] <- 1 / var(Y[,r] - Yss[,r], na.rm = TRUE)
  }
  mu1 <- solve(t(pcomp1$rotation[,1:fpca1$selectK]) %*% pcomp1$rotation[,1:fpca1$selectK], t(pcomp1$rotation[,1:fpca1$selectK]) %*% pcomp1$center)
  psi <- matrix(0, nt, fpca1$selectK)
  phi1 <- array(0, dim = c(nreg, nreg, fpca1$selectK))
  phi2 <- matrix(0, nreg * fpca1$selectK, nreg)
  lambda <- matrix(0, nrow = ncol(B), fpca1$selectK)
  eta <- matrix(0, nrow = nsub * nreg, ncol = fpca1$selectK)
  for (l in 1:fpca1$selectK) {
    psi[,l] <- pcomp1$rotation[, l]
    lambda[,l] <- solve(t(B) %*% B, t(B) %*% psi[,l])
    rxi <- t(matrix(pcomp1$x[,l], nrow = nreg))
    pcomp2 <- princomp(rxi)
    phi1[,,l] <- pcomp2$loadings[, 1:nreg]
    if (l > 1) {
      for (r in 1:nreg) {
        if (sum((phi1[, r, l] + phi1[, r, 1])^2) < sum((phi1[, r, l] - phi1[, r, 1])^2)) {
          phi1[, r, l] <- -phi1[, r, l]
        }
      }
    }
    for (r in 1:nreg) {
      phi2[((l - 1) * nreg + 1):(l * nreg),] <- phi1[,,l]
    }
    mu2 <- solve(pcomp2$loadings[, 1:nreg], pcomp2$center)
    for (i in 1:nsub) {
      eta[((i - 1) * nreg + 1):(i * nreg), l] <- pcomp2$scores[i, ] + mu2 + mu1[l]
    } 
  }
  phi_mat <- matrix(0, nrow = nreg * nreg, ncol = fpca1$selectK)
  for (r in 1:nreg) {
    for (rp in 1:nreg) {
      phi_mat[(r - 1) * nreg + rp, ] <- phi1[rp,r,]
    }
  }
  beta <- matrix(0, nreg * ncol(X), fpca1$selectK)
  preceta <- matrix(0, nreg, fpca1$selectK)
  d <- ncol(X)
  for (l in 1:fpca1$selectK) {
    for (r in 1:nreg) {
      seqr <- seq(from = r, to = nreg * nsub, by = nreg)
      etar <- eta[seqr, l]
      beta[((r - 1) * d + 1):(r * d), l] <- solve(t(X) %*% X, t(X) %*% etar)
      preceta[r, l] <- 1 / var(etar)
    }
  }
  rho <- mean(cor(phi_mat)[upper.tri(cor(phi_mat))])
  alpha <- mean(diag(cov(phi_mat)))
  C_rho_alpha <- mean(diag(cov(phi_mat)))
  delta_eta <- matrix(1, nreg, ldim)
  delta_eta[1, 1] <- preceta[1,1]
  for (l in 2:ldim) {
    delta_eta[1, l] <- preceta[1, l] / cumprod(delta_eta[1,])[l - 1]
  }
  for (r in 2:nreg) {
    delta_eta[r, 1] <- preceta[r, 1] / cumprod(delta_eta[,1])[r - 1]
  }
  for (r in 2:nreg) {
    for (l in 2:ldim) {
      delta_eta[r, l] <- preceta[r, l] / (cumprod(delta_eta[,l])[r - 1] * cumprod(delta_eta[1,])[l - 1])
    }
  }
  list(eta = eta, phi = phi2, phi_cube = phi1, lambda = lambda,
       omega = omega, rho = rho, alpha = alpha,
       phi_mat = phi_mat, beta = beta,
       preceta = preceta, C_rho_alpha = C_rho_alpha, delta_eta = delta_eta)
}
