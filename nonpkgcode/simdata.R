library(mgcv)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
nt <- 100
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 9
nsub <- 50
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
      eta[(i - 1) * nreg + r, l] <- rnorm(1, sd = sigma[l])
    }
  }
}
phi <- array(dim = c(nreg, nreg, ldim))
for (l in 1:ldim) {
  phi[,,l] <- diag(nreg)
}
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
init_mcmc <- initialize_mcmc(Y, tt, nt, B)
microbenchmark::microbenchmark(result <- run_mcmc(Y, X, B, tt, penalty, l, 1000, 100, 1, init_mcmc), times = 1)
plot(B %*% result$samples$lambda[,1,1000], type = "l")
for (i in 501:1000) {
  lines(B %*% result$samples$lambda[,1,i])
}

plot(Y[1:nt,4])
lines(result$samples$fit[1:nt,4], col = "blue")
