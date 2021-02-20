library(mgcv)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
nt <- 100
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 9
nsub <- 500
ndf <- 25
d <- 2
ldim <- 4
sqexp <- function(t, tp) {
  exp(-(t - tp)^2 / .2)
}
kern <- matrix(0, nrow = nt, ncol = nt)
for (i in 1:nt) {
  for (j in 1:nt) {
    kern[i, j] <- 5 * sqexp(tt[i], tt[j])
  }
}
eigenfuncs <- eigen(kern)$vectors
eta <- matrix(0, nsub * nreg, ldim)
sigma <- ldim:1
for (l in 1:ldim) {
  for (i in 1:nsub) {
    for (r in 1:nreg) {
      print((ldim - l + 1) * 1 / r)
      eta[(i - 1) * nreg + r, l] <- rnorm(1, sd = (ldim - l + 1) * 1 / r)
    }
  }
}
phi <- array(dim = c(nreg, nreg, ldim))
orthmat <- pracma::randortho(nreg)
for (r in 1:nreg) {
  if (orthmat[r,r] < 0) orthmat[,r] <- -orthmat[,r]
}
for (l in 1:ldim) {
  # phi[,,l] <- diag(nreg)
  phi[,,l] <- orthmat
  # phi[,,l] <- rWishart(1, 10, diag(nreg))
  # phi[,,l] <- eigen(phi[,,l])$vectors
  # phi[,,l] <- pracma::randortho(nreg)
  # phi[,,l] <- matrix(rnorm(nreg * nreg), nreg, nreg)
  # phi[,,l] <- Re(eigen(phi[,,l])$vectors)
}
scale <- .2
rho <- .9
C_rho <- scale * rho * rep(1, ldim) %*% t(rep(1, ldim)) + scale * (1 - rho) * diag(ldim)
phi_mat <- matrix(0, nrow = nreg^2, ldim)
log_phi <- 0
for (r in 1:nreg) {
  for (rp in 1:nreg) {
    x <- phi[r, rp, ]
    phi_mat[(r - 1) * nreg + rp, ] <- x
    log_phi <- log_phi + mvtnorm::dmvnorm(x, rep(0, ldim), sigma = C_rho, log = TRUE)
  }
}
sum(mvtnorm::dmvnorm(phi_mat, rep(0, ldim), sigma = C_rho, log = TRUE))
cov(phi_mat)
log_phi
Y <- matrix(0, nsub * nt, ncol = nreg)
for (r in 1:nreg) {
  for (i in 1:nsub) {
    for (l in 1:ldim) {
      Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + 
        c(eigenfuncs[,l] %*% t(eta[((i - 1) * nreg + 1):(i * nreg), l]) %*% t(phi[,,l])) + c(rnorm(nt * nreg, sd = .1))
      
    }
  }
}
X <- cbind(rep(1, nsub), matrix(rnorm(nsub * (d - 1)), nrow = nsub, ncol = d - 1))
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
matplot(tt, B, xlab = "time", ylab = "spline", type = "l")
init_mcmc <- initialize_mcmc(Y, tt, nt, B, X, ldim = 4)
microbenchmark::microbenchmark(result <- run_mcmc(Y, X, B, tt, penalty, l, 1000, 100, 1, init_mcmc), times = 1)
plot(B %*% result$samples$lambda[,1,1000], type = "l")
for (i in 501:1000) {
  lines(B %*% result$samples$lambda[,1,i])
}

plot(Y[1:nt,4])
lines(result$samples$fit[1:nt,4], col = "blue")
