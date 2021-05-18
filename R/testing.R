rn <- function(n, p) {
  r <- sqrt(-log(1 - p / n))
  return(r)
}

get_theta <- function(eta, phi, nreg) {
  nsub <- nrow(eta) / nreg
  ldim <- ncol(eta)
  eta <- reshape_nreg(eta, nsub, nreg)
  p <- nreg * ldim
  mu <- -rn(nsub, p)^2 * (p - nsub + .5)
  for (i in 1:nreg) {
    mu <- mu + rn(nsub, nreg)^2 * (nreg - nsub + .5)
  }
  sigsq <- 2 * rn(nsub, p)^2
  for (i in 1:nreg) {
    sigsq <- sigsq - 2 * rn(nsub, nreg)^2
  }
  theta <- matrix(0, nrow = nsub, ncol = nreg * ldim)
  for (i in 1:nsub) {
    for (l in 1:ldim) {
      idx1 <- (l - 1) * nreg + 1
      idx2 <- l * nreg
      theta[i, idx1:idx2] <- phi[,,l] %*% eta[i, idx1:idx2]
    }
  }
  return(theta)
}

chisqd <- function(x, nreg, ldim) {
  covx <- cov(x)
  n <- nrow(x)
  p <- ncol(x)
  loglambda <- log(det(covx))
  for (j in 1:ldim) {
    idx1 <- (j - 1) * nreg + 1
    idx2 <- j * nreg
    loglambda <- loglambda - log(det(covx[idx1:idx2, idx1:idx2]))
  }
  loglambda <- (n + 1) / 2 * loglambda
  rho <- 1 - (2 * (p^3 - ldim * nreg^3) + 9 * (p^2  - ldim * nreg^2)) /
    (6 * (nsub + 1) * (p^2 - ldim * nreg^2))
  f <- .5 * (p^2 - ldim * nreg^2)
  return(pchisq(-2 * rho * loglambda, df = f, lower.tail = FALSE))
}

normald <- function(x, nreg, ldim) {
  covx <- cov(x)
  p <- nreg * ldim
  n <- nrow(x)
  mu <- -rn(n, p)^2 * (p - n + .5) + ldim * rn(n, nreg)^2 * (nreg - n + .5)
  sigsq <- 2 * rn(n, p)^2 - 2 * ldim * rn(n, nreg)^2
  loglambda <- log(det(covx))
  for (j in 1:ldim) {
    idx1 <- (j - 1) * nreg + 1
    idx2 <- j * nreg
    loglambda <- loglambda - log(det(covx[idx1:idx2, idx1:idx2]))
  }
  q <- (loglambda - mu) / sqrt(sigsq)
  return(pnorm(q, lower.tail = TRUE))
}

#' @export
get_pvals_partial <- function(result, func = normald) {
  iters <- result$control$iterations
  burnin <- result$control$burnin + 1
  pvals <- sapply(burnin:iters, function(i) {theta <- 
    get_theta(result$samples$eta[,,i], result$samples$phi[[i]], nreg);
  func(theta, nreg, ldim)}
  )
  return(pvals)
}

get_pmin <- function(pvals) {
  J <- length(pvals)
  pvals_sorted <- sort(pvals)
  return(1 - min(sapply(1:J, function(i)
    J * (1 - pvals_sorted[i]) / (J - i + 1))))
}
