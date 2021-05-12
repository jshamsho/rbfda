library(mgcv)
library(rbfda)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/simulate_data.R")
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/testing.R")

nsub <- 100
nt <- 60
nreg <- 7
ldim <- 4
d <- 2
ndf <- 15
tt <- seq(from = 0, to = 1, length.out = nt)
sim_data <- sim_weak(nt, nsub, nreg, ldim)
sim_data <- sim_partial(nt, nsub, nreg, ldim)
sim_data <- sim_partial_cs(nt, nsub, nreg, ldim, rho1 = .8)
sim_data <- sim_non_partial(nt, nsub, nreg, ldim, rho1 = .8, rho2 = .2)
X <- cbind(rep(1, nsub), matrix(rnorm(nsub * (d - 1)), nrow = nsub, ncol = d - 1))
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
init_mcmc <- initialize_mcmc_weak(sim_data$Y, tt, B, X, ldim = ldim)
plot(init_mcmc$psi[,1])
result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
                   penalty = penalty, ldim = ldim, iter = 5000, burnin = 1000,
                   thin = 1, init_ = init_mcmc, covstruct = "weak")
cumprod(result$samples$delta_eta1[,5000]) %*% t(cumprod(result$samples$delta_eta2[,5000]))
cumprod(init_mcmc$delta_eta1) %*% t(cumprod(init_mcmc$delta_eta2))
init_mcmc$prec_eta
pvals <- get_pvals_weak(result)
hist(pvals)
abline(v = median(pvals))
plot(pvals, type = "l")
median(pvals)
get_pmin(pvals)


b <- numeric(ndf)
Q <- matrix(0, ndf, ndf)
for (i in 1:nsub) {
  b <- b + t(B) %*% sim_data$Y[((i - 1) * nt + 1):(i * nt),] %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1] -
    t(B) %*% B %*% init_mcmc$lambda[,2:ldim] %*% t(init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 2:ldim]) %*% t(init_mcmc$phi) %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1]
  Q <- Q + as.numeric(init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1] %*% t(init_mcmc$phi) %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1]) * t(B) %*% B
}
lambda1 <- MASS::mvrnorm(n = 1, mu = solve(Q) %*% b, Sigma = solve(Q))
plot(B %*% lambda1, type = "l")
r <- 3
seqr <- ((r - 1) * d):(r * d)
plot.new()
for (i in 201:5000) {
  tmpsum <- numeric(nt)
  for (l in 1:init_mcmc$npc) {
    tmpsum <- tmpsum + B %*% result$samples$lambda[,l, i] %*% 
      t(result$samples$phi[[i]][r,,l]) %*% 
      (result$samples$beta[seqr, l, i]
  }
  lines(tmpsum, col = "blue")
}

x_vec <- c(1, 1)
meanf <- numeric(nt)
for (l in 1:ldim) {
  meanf <- meanf + B %*% result$samples$lambda[, l, 1] %*%
    result$samples$phi[[1]][r,,l] %*% 
    t(matrix(result$samples$beta[,l,1], d)) %*% x_vec
}
plot(tt, meanf, type = "l", ylim = c(-2,2))
for (i in 500:5000) {
  meanf <- numeric(nt)
  for (l in 1:ldim) {
    meanf <- meanf + B %*% result$samples$lambda[, l, i] %*%
      result$samples$phi[[i]][r,,l] %*% 
      t(matrix(result$samples$beta[,l,i], d)) %*% x_vec
  }
  lines(tt, meanf)
}

efunc <- 1
plot(B %*% result$samples$lambda[,efunc,100], type = "l")
evec <- numeric(5000)
for (i in 1000:5000) {
  lines(B %*% result$samples$lambda[,efunc,i])
  evec[i] <- (B %*% result$samples$lambda[,efunc,i])[30]
}
lines(-sim_data$psi[,efunc], col = "red")
lines(init_mcmc$psi[,efunc], col = "green")
lines(B %*% apply(result$samples$lambda[,efunc,],1,mean), col = "blue")
sum((init_mcmc$psi[,efunc] - sim_data$psi[,efunc])^2)
sum((B %*% apply(result$samples$lambda[,efunc,], 1, median) - sim_data$psi[,efunc])^2)

r <- 5
i <- 3
plot(sim_data$Y[((i - 1) * nt + 1):(i * nt),r])
seqr <- ((i - 1) * nreg + 1):(i * nreg)
for (i in 201:5000) {
  est <- B %*% result$samples$lambda[,, i] %*% 
      t(result$samples$eta[seqr, , i]) %*%
      result$samples$phi[r,,i]
  lines(est, col = "blue")
}

delta_eta_cumprod <- array(0, dim = c(nreg, ldim, 5000))
for (i in 1:5000) {
  delta_eta_cumprod[,,i] <- cumprod(result$samples$delta_eta1[,i]) %*% 
    t(cumprod(result$samples$delta_eta2[,i]))
}
r <- 1
l <- 1
plot(1 / delta_eta_cumprod[r,l,201:5000])
abline(h = 1 / init_mcmc$prec_eta[r,l])
abline(h = ((ldim - l + 1) * 1 / r)^2, col = "red")
1 / init_mcmc$prec_eta[r, l]
quantile(1 / delta_eta_cumprod[r,l,1:5000], c(.025, .975))
((ldim - l + 1) * 1 / r)^2








get_cumprod_coef <- function(input) {
  delta <- numeric(length(input))
  delta[1] <- input[1]
  if (length(input) > 1) {
    for (i in 2:length(input)) {
      delta[i] <- input[i] / input[i - 1]
    }
  }
  return(delta)
}

mynmf <- nmf(init_mcmc$prec_eta, rank = 3)
delta_eta1 <- basis(mynmf)
delta_eta2 <- coef(mynmf)
a_basis <- sqrt(delta_eta1[1,1] * delta_eta2[1,1]) / delta_eta1[1,1]
a_coef <- sqrt(delta_eta1[1,1] * delta_eta2[1,1]) / delta_eta2[1,1]
delta_eta1 <- a_basis * delta_eta1
delta_eta2 <- a_coef * delta_eta2
for (i in 1:3) {
  delta_eta1[,i] <- get_cumprod_coef(delta_eta1[,i])
  delta_eta2[i,] <- get_cumprod_coef(delta_eta2[i,])
}
xi_eta <- matrix(1, nreg * nsub, ldim)
delta_eta_diff(delta_eta1, delta_eta2, init_mcmc$eta, init_mcmc$beta, xi_eta, nsub,
               nreg, ldim, ndf, X)
