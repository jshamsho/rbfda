library(refund)
initialize_mcmc_partial <- function(Y, tt, nt, B, X, pve = .99, ldim = NULL) {
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
  if (is.null(ldim)) {
    fpca1 <- refund::fpca.face(t(Yt), center = TRUE, pve = pve, knots = floor(length(tt) / 5))
  } else {
    fpca1 <- refund::fpca.face(t(Yt), center = TRUE, floor(length(tt) / 5), npc = ldim)
  }
  for (i in 1:nr) {
    Ys[, i] <- fpca1$mu + fpca1$efunctions[,1:fpca1$npc] %*% fpca1$scores[i,1:fpca1$npc]
    Ysm[,i] <- Yt[, i] - fpca1$mu
  }
  pcomp1 <- prcomp(t(Ys))
  for (r in 1:nreg) {
    r_ind <- seq(from = r, length.out = nsub, by = nreg)
    for (i in 1:nsub) {
      Yss[, r] <- c(Ys[, r_ind])
    }
    omega[r] <- 1 / var(Y[,r] - Yss[,r], na.rm = TRUE)
  }
  mu1 <- solve(t(pcomp1$rotation[,1:fpca1$npc]) %*% pcomp1$rotation[,1:fpca1$npc], t(pcomp1$rotation[,1:fpca1$npc]) %*% pcomp1$center)
  psi <- matrix(0, nt, fpca1$npc)
  phi1 <- array(0, dim = c(nreg, nreg, fpca1$npc))
  phi2 <- matrix(0, nreg * fpca1$npc, nreg)
  lambda <- matrix(0, nrow = ncol(B), fpca1$npc)
  eta <- matrix(0, nrow = nsub * nreg, ncol = fpca1$npc)
  for (l in 1:fpca1$npc) {
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
  
  phi_mat <- matrix(0, nrow = nreg * nreg, ncol = fpca1$npc)
  for (r in 1:nreg) {
    for (rp in 1:nreg) {
      phi_mat[(r - 1) * nreg + rp, ] <- phi1[rp,r,]
    }
  }
  beta <- matrix(0, nreg * ncol(X), fpca1$npc)
  preceta <- matrix(0, nreg, fpca1$npc)
  d <- ncol(X)
  
  for (l in 1:fpca1$npc) {
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
  delta_eta <- matrix(1, nreg, fpca1$npc)
  delta_eta[1, 1] <- preceta[1,1]
  for (l in 2:fpca1$npc) {
    delta_eta[1, l] <- preceta[1, l] / cumprod(delta_eta[1,])[l - 1]
  }
  for (r in 2:nreg) {
    delta_eta[r, 1] <- preceta[r, 1] / cumprod(delta_eta[,1])[r - 1]
  }
  for (r in 2:nreg) {
    for (l in 2:fpca1$npc) {
      delta_eta[r, l] <- preceta[r, l] / (cumprod(delta_eta[,l])[r - 1] * cumprod(delta_eta[1,])[l - 1])
    }
  }
  
  list(eta = eta, phi = phi2, phi_cube = phi1, lambda = lambda,
       omega = omega, rho = rho, alpha = alpha,
       phi_mat = phi_mat, beta = beta,
       preceta = preceta, C_rho_alpha = C_rho_alpha, delta_eta = delta_eta,
       npc = fpca1$npc)
}

initialize_mcmc_weak <- function(Y, tt, nt, B, X, pve = .99, ldim = NULL) {
  nsub <- nrow(Y) / nt
  Yt <- matrix(0, nsub, nt * nreg)
  for (i in 1:nsub) {
    for (j in 1:nreg) {
      idx1 <- (i - 1) * nt + 1
      idx2 <- i * nt
      idx3 <- (j - 1) * nt + 1
      idx4 <- j * nt
      Yt[i, idx3:idx4] <- Y[idx1:idx2, j]
    }
  }
  mean.Y <- apply(Yt, 2, mean, na.rm = TRUE)
  mean.Ytmat <- matrix(rep(mean.Y, nsub), nrow = nsub, ncol = nt * nreg, byrow = TRUE)    
  mean.Yt.c <- Yt - mean.Ytmat
  Yt.reg.c <- matrix(mean.Yt.c, nrow = nsub * nt, ncol = nreg, byrow = TRUE)
  Yt.func.c <- matrix(t(mean.Yt.c), nrow = nsub * nreg, ncol = nt, byrow = TRUE)
  
}
initialize_mcmc_weak <- function(Y, tt, nt, B, X, pve = .99, ldim = NULL) {
  
  nsub <- nrow(Y) / nt
  nreg <- ncol(Y)
  L <- list()
  omega <- numeric(nreg)
  nr <- nsub * nreg
  Yt <- array(0, dim = c(nt, nreg, nsub))
  for (i in 1:nsub) {
    idx1 <- (i - 1) * nt + 1
    idx2 <- i * nt
    Yt[,,i] <- Y[idx1:idx2,]
  }
  xbar <- apply(Yt, c(1, 2), mean)
  margf <- matrix(0, length(tt), length(tt))
  margr <- matrix(0, nreg, nreg)
  for (i in 1:nsub) {
    margf <- margf + 1 / nsub * (Yt[,,i] - xbar) %*% t(Yt[,,i] - xbar)
  }
  for (i in 1:nsub) {
    margr <- margr + 1 / nsub * t(Yt[,,i] - xbar) %*% (Yt[,,i] - xbar)
  }
  eigenf <- eigen(margf)
  eigenr <- eigen(margr)
  phi <- eigenr$vectors
  lambda <- matrix(0, nrow(B), )
  Yt <- matrix(0, nrow = nt, ncol = nr)
  Ys <- matrix(0, nrow = nt, ncol = nr)
  Ysm <- matrix(0, nrow = nt, ncol = nr)
  Yss <- matrix(0, nrow = nt * nsub, ncol = nreg)
  for (i in 1:nsub) {
    Yt[, ((i - 1) * nreg + 1):(i * nreg)] <- Y[((i - 1) * nt + 1):(i * nt),]
  }
  fpca1 <- refund::fpca.face(t(Yt), center = TRUE, pve = pve, npc = ldim,
                             knots = floor(length(tt) / 5))
  for (i in 1:nr) {
    Ys[, i] <- fpca1$mu + fpca1$efunctions[,1:fpca1$npc] %*%
      fpca1$scores[i,1:fpca1$npc]
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
  mu1 <- solve(t(pcomp1$rotation[,1:fpca1$npc]) %*% pcomp1$rotation[,1:fpca1$npc], t(pcomp1$rotation[,1:fpca1$npc]) %*% pcomp1$center)
  psi <- matrix(0, nt, fpca1$npc)
  phi1 <- array(0, dim = c(nreg, nreg, fpca1$npc))
  phi2 <- matrix(0, nreg * fpca1$npc, nreg)
  lambda <- matrix(0, nrow = ncol(B), fpca1$npc)
  eta <- matrix(0, nrow = nsub * nreg, ncol = fpca1$npc)
  for (l in 1:fpca1$npc) {
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
  
  phi_mat <- matrix(0, nrow = nreg * nreg, ncol = fpca1$npc)
  for (r in 1:nreg) {
    for (rp in 1:nreg) {
      phi_mat[(r - 1) * nreg + rp, ] <- phi1[rp,r,]
    }
  }
  beta <- matrix(0, nreg * ncol(X), fpca1$npc)
  preceta <- matrix(0, nreg, fpca1$npc)
  d <- ncol(X)
  
  for (l in 1:fpca1$npc) {
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
  delta_eta <- matrix(1, nreg, fpca1$npc)
  delta_eta[1, 1] <- preceta[1,1]
  for (l in 2:fpca1$npc) {
    delta_eta[1, l] <- preceta[1, l] / cumprod(delta_eta[1,])[l - 1]
  }
  for (r in 2:nreg) {
    delta_eta[r, 1] <- preceta[r, 1] / cumprod(delta_eta[,1])[r - 1]
  }
  for (r in 2:nreg) {
    for (l in 2:fpca1$npc) {
      delta_eta[r, l] <- preceta[r, l] / (cumprod(delta_eta[,l])[r - 1] * cumprod(delta_eta[1,])[l - 1])
    }
  }
  
  list(eta = eta, phi = phi2, phi_cube = phi1, lambda = lambda,
       omega = omega, rho = rho, alpha = alpha,
       phi_mat = phi_mat, beta = beta,
       preceta = preceta, C_rho_alpha = C_rho_alpha, delta_eta = delta_eta,
       npc = fpca1$npc)
}
