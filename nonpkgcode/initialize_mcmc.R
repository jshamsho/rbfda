library(refund)
library(magrittr)
new_partial_class <- function(Y, tt, B = NULL, X = NULL, pve = NULL, ldim = NULL) {
  nsub <- nrow(Y) / nt
  nreg <- ncol(Y)
  nr <- nsub * nreg
  omega <- numeric(nreg)
  npc <- NULL
  mu1 <- NULL
  mu2 <- NULL
  psi <- NULL
  phi_cube <- NULL
  phi_mat <- NULL
  lambda <- NULL
  eta <- NULL
  prec_eta <- NULL
  beta <- NULL
  delta_eta <- NULL
  Y.trans <- matrix(0, nrow = nt, ncol = nr)
  Y.trans.smoothed <- matrix(0, nrow = nt, ncol = nr)
  Y.trans.smoothed.centered <- matrix(0, nrow = nt, ncol = nr)
  if (is.null(B)) {
    B <- mgcv::smoothCon(s(tt, k = ceiling(length(tt) * 0.2), bs = "ps", m = 2),
                         data.frame(tt), absorb.cons = FALSE)[[1]]$X
  }
  if (is.null(X)) {
    X <- matrix(1, nrow = nsub, ncol = 1)
  }
  if (!is.null(ldim)) {
    npc <- ldim
  }
  if (is.null(pve)) {
    pve <- .99
  }
  members <- list(Y = Y, Y.trans = Y.trans, Y.trans.smoothed = Y.trans.smoothed,
                  Y.trans.smoothed.centered = Y.trans.smoothed.centered,
                  tt = tt, B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, mu1 = mu1, mu2 = mu2, psi = psi,
                  phi_cube = phi_cube, phi_mat = phi_mat, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta, pve = pve)
  return(structure(members, class = "partial_class"))
}

run_fpca <- function(partial_class) {
  for (v in 1:length(partial_class)) assign(names(partial_class)[v],
                                            partial_class[[v]])
  for (i in 1:nsub) {
    Y.trans[, ((i - 1) * nreg + 1):(i * nreg)] <- Y[((i - 1) * nt + 1):(i * nt),]
  }
  if (is.null(npc)) {
    fpca1 <- refund::fpca.face(t(Y.trans), center = TRUE,
                               pve = pve, knots = floor(length(tt) * 0.2))
  } else {
    fpca1 <- refund::fpca.face(t(Y.trans), center = TRUE, floor(length(tt) * 0.2),
                               npc = npc)
  }
  npc <- fpca1$npc
  
  for (i in 1:(nsub * nreg)) {
    Y.trans.smoothed[, i] <- fpca1$mu + fpca1$efunctions[,1:npc] %*%
      fpca1$scores[i, 1:npc]
    Y.trans.smoothed.centered[,i] <- Y.trans[, i] - fpca1$mu
  }
  members <- list(Y = Y, Y.trans = Y.trans, Y.trans.smoothed = Y.trans.smoothed,
                  Y.trans.smoothed.centered = Y.trans.smoothed.centered,
                  tt = tt, B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, mu1 = mu1, mu2 = mu2, psi = psi,
                  phi_cube = phi_cube, phi_mat = phi_mat, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta)
  return(structure(members, class = "partial_class"))
}

estimate_residual_error <- function(partial_class) {
  for (v in 1:length(partial_class)) assign(names(partial_class)[v],
                                            partial_class[[v]])
  Y.smoothed <- Y
  for (r in 1:nreg) {
    r_ind <- seq(from = r, length.out = nsub, by = nreg)
    for (i in 1:nsub) {
      Y.smoothed[, r] <- c(Y.trans.smoothed[, r_ind])
    }
    omega[r] <- 1 / var(Y[,r] - Y.smoothed[,r], na.rm = TRUE)
  }
  members <- list(Y = Y, Y.trans = Y.trans, Y.trans.smoothed = Y.trans.smoothed,
                  Y.trans.smoothed.centered = Y.trans.smoothed.centered,
                  tt = tt, B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, mu1 = mu1, mu2 = mu2, psi = psi,
                  phi_cube = phi_cube, phi_mat = phi_mat, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta)
  return(structure(members, class = "partial_class"))
}

set_size_param <- function(partial_class) {
  for (v in 1:length(partial_class)) assign(names(partial_class)[v],
                                            partial_class[[v]])
  psi <- matrix(0, nt, npc)
  phi_cube <- array(0, dim = c(nreg, nreg, npc))
  phi_mat <- matrix(0, nreg * npc, nreg)
  lambda <- matrix(0, nrow = ncol(B), npc)
  eta <- matrix(0, nrow = nsub * nreg, ncol = npc)
  beta <- matrix(0, nreg * ncol(X), npc)
  prec_eta <- matrix(0, nreg, npc)
  members <- list(Y = Y, Y.trans = Y.trans, Y.trans.smoothed = Y.trans.smoothed,
                  Y.trans.smoothed.centered = Y.trans.smoothed.centered,
                  tt = tt, B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, mu1 = mu1, mu2 = mu2, psi = psi,
                  phi_cube = phi_cube, phi_mat = phi_mat, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta)
  return(structure(members, class = "partial_class"))
}

set_param <- function(partial_class) {
  for (v in 1:length(partial_class)) assign(names(partial_class)[v],
                                            partial_class[[v]])
  pcomp1 <- prcomp(t(Y.trans.smoothed))
  mu1 <- solve(t(pcomp1$rotation[,1:npc]) %*% 
                 pcomp1$rotation[,1:npc],
               t(pcomp1$rotation[,1:npc]) %*% pcomp1$center)
  for (l in 1:npc) {
    psi[,l] <- pcomp1$rotation[, l]
    lambda[,l] <- solve(t(B) %*% B, t(B) %*% psi[,l])
    rxi <- t(matrix(pcomp1$x[,l], nrow = nreg))
    pcomp2 <- princomp(rxi)
    phi_cube[,,l] <- pcomp2$loadings[, 1:nreg]
    # if (l > 1) {
    #   for (r in 1:nreg) {
    #     if (sum((phi1[, r, l] + phi1[, r, 1])^2) < sum((phi1[, r, l] - phi1[, r, 1])^2)) {
    #       phi_cube[, r, l] <- -phi1[, r, l]
    #     }
    #   }
    # }
    for (r in 1:nreg) {
      phi_mat[((l - 1) * nreg + 1):(l * nreg),] <- phi_cube[,,l]
    }
    mu2 <- solve(pcomp2$loadings[, 1:nreg], pcomp2$center)
    for (i in 1:nsub) {
      eta[((i - 1) * nreg + 1):(i * nreg), l] <- pcomp2$scores[i, ] + mu2 + mu1[l]
    } 
  }
  
  # phi_mat <- matrix(0, nrow = nreg * nreg, ncol = npc)
  # for (r in 1:nreg) {
  #   for (rp in 1:nreg) {
  #     phi_mat[(r - 1) * nreg + rp, ] <- phi1[rp,r,]
  #   }
  # }
  
  d <- ncol(X)
  
  for (l in 1:npc) {
    for (r in 1:nreg) {
      seqr <- seq(from = r, to = nreg * nsub, by = nreg)
      beta[((r - 1) * d + 1):(r * d), l] <- 
        solve(t(X) %*% X, t(X) %*% eta[seqr, l])
      prec_eta[r, l] <- 1 / var(eta[seqr, l])
    }
  }
  
  delta_eta <- matrix(1, nreg, npc)
  delta_eta[1, 1] <- prec_eta[1,1]
  for (l in 2:npc) {
    delta_eta[1, l] <- prec_eta[1, l] / cumprod(delta_eta[1,])[l - 1]
  }
  if (nreg > 1) {
    for (r in 2:nreg) {
      delta_eta[r, 1] <- prec_eta[r, 1] / cumprod(delta_eta[,1])[r - 1]
    }
    for (r in 2:nreg) {
      for (l in 2:npc) {
        delta_eta[r, l] <- prec_eta[r, l] / (cumprod(delta_eta[,l])[r - 1] * 
                                               cumprod(delta_eta[1,])[l - 1])
      }
    }
  }
  
  
  members <- list(Y = Y, Y.trans = Y.trans, Y.trans.smoothed = Y.trans.smoothed,
                  Y.trans.smoothed.centered = Y.trans.smoothed.centered,
                  tt = tt, B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, mu1 = mu1, mu2 = mu2, psi = psi,
                  phi_cube = phi_cube, phi_mat = phi_mat, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta)
  return(structure(members, class = "partial_class"))
}

initialize_mcmc_partial <- function(Y, tt, B = NULL,
                                    X = NULL, pve = NULL, ldim = NULL) {
  partial_mcmc_init <- new_partial_class(Y = Y, tt = tt, B = B, X = X,
                                         pve = pve, ldim = ldim) %>%
    run_fpca %>%
    estimate_residual_error %>%
    set_size_param %>% 
    set_param
  return(partial_mcmc_init)
}

initialize_mcmc_weak <- function(Y, tt, B, X, pve = .99, npc = NULL) {
  nt <- length(tt)
  nsub <- nrow(Y) / nt
  Y.smoothed <- Y
  omega <- numeric(nreg)
  Y.trans <- matrix(0, nsub * nreg, ncol = nt)
  for (i in 1:nsub) {
    for (j in 1:nreg) {
      idx1 <- (i - 1) * nt + 1
      idx2 <- i * nt
      Y.trans[(i - 1) * nreg + j, ] <- Y[idx1:idx2, j]
    }
  }
 if (is.null(npc)) {
    fpca1 <- refund::fpca.face(Y.trans, center = TRUE,
                               pve = pve, knots = floor(length(tt) * 0.2))
  } else {
    fpca1 <- refund::fpca.face(Y.trans, center = TRUE, floor(length(tt) * 0.2),
                               npc = npc)
  }
  npc <- fpca1$npc
  psi <- fpca1$efunctions
  lambda <- matrix(0, ncol(B), npc)
  for (i in 1:npc) {
    lambda[, i] <- solve(t(B) %*% B) %*% t(B) %*% fpca1$efunctions[, i]
  }
  phi <- eigen(cov(Y))$vectors
  eta <- matrix(0, nsub * nreg, npc)
  prec_eta <- matrix(0, nreg, npc)
  beta <- matrix(0, nreg * ncol(X), npc)
  eigdesign <- matrix(0, nreg * nt, nreg * npc)
  i <- 0
  i <- 0
  for (j in 1:nreg) {
    for (l in 1:npc) {
      i <- i + 1
      eigdesign[, i] <- outer(psi[, l], phi[, j])
    }
  }
  for (i in 1:nsub) {
    idx1 <- (i - 1) * nreg + 1
    idx2 <- i * nreg
    idx3 <- (i - 1) * nt + 1
    idx4 <- i * nt
    eta[idx1:idx2, ] <- matrix(lm(c(Y[idx3:idx4,]) ~ eigdesign - 1)$coef,
                               nreg, byrow = TRUE)
  }
  d <- ncol(X)
  
  for (l in 1:npc) {
    for (r in 1:nreg) {
      seqr <- seq(from = r, to = nreg * nsub, by = nreg)
      beta[((r - 1) * d + 1):(r * d), l] <- 
        solve(t(X) %*% X, t(X) %*% eta[seqr, l])
      prec_eta[r, l] <- 1 / var(eta[seqr, l])
    }
  }
  for (i in 1:nsub) {
    for (l in 1:npc) {
      idx1 <- (i - 1) * nt + 1
      idx2 <- i * nt
      idx3 <- (i - 1) * nreg + 1
      idx4 <- i * nreg
      Y.smoothed[idx1:idx2, ] <- psi %*% t(eta[idx3:idx4, ]) %*% t(phi)
    }
  }
  for (r in 1:nreg) {
    omega[r] <- 1 / var(Y[,r] - Y.smoothed[,r], na.rm = TRUE)
  }
  
  delta_eta <- matrix(1, nreg, npc)
  delta_eta[1, 1] <- prec_eta[1,1]
  for (l in 2:npc) {
    delta_eta[1, l] <- prec_eta[1, l] / cumprod(delta_eta[1,])[l - 1]
  }
  if (nreg > 1) {
    for (r in 2:nreg) {
      delta_eta[r, 1] <- prec_eta[r, 1] / cumprod(delta_eta[,1])[r - 1]
    }
    for (r in 2:nreg) {
      for (l in 2:npc) {
        delta_eta[r, l] <- prec_eta[r, l] / (cumprod(delta_eta[,l])[r - 1] * 
                                               cumprod(delta_eta[1,])[l - 1])
      }
    }
  }
  
  members <- list(Y = Y, Y.trans = Y.trans, Y.smoothed = Y.smoothed, tt = tt,
                  B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, psi = psi, phi = phi, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta, pve = pve)
  return(members)
}

new_weak_class <- function(Y, tt, B = NULL, X = NULL, pve = NULL, ldim = NULL) {
  nsub <- nrow(Y) / nt
  nreg <- ncol(Y)
  nr <- nsub * nreg
  omega <- numeric(nreg)
  npc <- NULL
  psi <- NULL
  phi <- NULL
  lambda <- NULL
  eta <- NULL
  prec_eta <- NULL
  beta <- NULL
  delta_eta <- NULL
  Y.trans <- matrix(0, nrow = nt, ncol = nr)
  Y.smoothed <- matrix(0, nrow = nt * nsub, ncol = nreg)
  if (is.null(B)) {
    B <- mgcv::smoothCon(s(tt, k = ceiling(length(tt) * 0.2), bs = "ps", m = 2),
                         data.frame(tt), absorb.cons = FALSE)[[1]]$X
  }
  if (is.null(X)) {
    X <- matrix(1, nrow = nsub, ncol = 1)
  }
  if (!is.null(ldim)) {
    npc <- ldim
  }
  if (is.null(pve)) {
    pve <- .99
  }
  members <- list(Y = Y, Y.trans = Y.trans, Y.smoothed = Y.smoothed, tt = tt,
                  B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, psi = psi, phi = phi, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta, pve = pve)
  return(members)
}

run_fpca_weak <- function(weak_class) {
  for (v in 1:length(partial_class)) assign(names(partial_class)[v],
                                            partial_class[[v]])
  for (i in 1:nsub) {
    Y.trans[, ((i - 1) * nreg + 1):(i * nreg)] <- Y[((i - 1) * nt + 1):(i * nt),]
  }
  if (is.null(npc)) {
    fpca1 <- refund::fpca.face(t(Y.trans), center = TRUE,
                               pve = pve, knots = floor(length(tt) * 0.2))
  } else {
    fpca1 <- refund::fpca.face(t(Y.trans), center = TRUE, floor(length(tt) * 0.2),
                               npc = npc)
  }
  npc <- fpca1$npc
  members <- list(Y = Y, Y.trans = Y.trans, Y.smoothed = Y.smoothed, tt = tt,
                  B = B, X = X, nsub = nsub, nreg = nreg, npc = npc,
                  omega = omega, psi = psi, phi = phi, lambda = lambda,
                  eta = eta, prec_eta = prec_eta, beta = beta,
                  delta_eta = delta_eta, pve = pve)
  return(members)
}

