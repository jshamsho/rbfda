library(mgcv)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
nt <- 60
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 5
nsub <- 100
ndf <- 15
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
  # phi[,,l] <- orthmat
  # phi[,,l] <- rWishart(1, 10, diag(nreg))
  phi[,,l] <- pracma::randortho(nreg)
  # phi[,,l] <- orthmat
  # phi[,,l] <- matrix(rnorm(nreg * nreg), nreg, nreg)
  # phi[,,l] <- Re(eigen(phi[,,l])$vectors)
}
# phi[,,3] <- pracma::randortho(nreg)
scale <- init_mcmc$alpha
rho <- .001
C_rho <- scale * rho * rep(1, ldim) %*% t(rep(1, ldim)) + scale * (1 - rho) * diag(ldim)
phi_mat <- matrix(0, nrow = nreg^2, ldim)
log_phi <- 0
for (r in 1:nreg) {
  for (rp in 1:nreg) {
    # x <- phi[r, rp, ]
    x <- phi[rp, r,]
    phi_mat[(r - 1) * nreg + rp, ] <- x
    log_phi <- log_phi + mvtnorm::dmvnorm(x, rep(0, ldim), sigma = C_rho, log = TRUE)
  }
}
sum(mvtnorm::dmvnorm(phi_mat, rep(0, ldim), sigma = C_rho, log = TRUE))
cov(phi_mat)
log_phi
Y <- matrix(0, nsub * nt, ncol = nreg)
for (i in 1:nsub) {
  for (l in 1:ldim) {
    Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + 
      c(eigenfuncs[,l] %*% t(eta[((i - 1) * nreg + 1):(i * nreg), l]) %*% t(phi[,,l])) 
  }
  Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + c(rnorm(nt * nreg, sd = .1))
}
X <- cbind(rep(1, nsub), matrix(rnorm(nsub * (d - 1)), nrow = nsub, ncol = d - 1))
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
matplot(tt, B, xlab = "time", ylab = "spline", type = "l")
init_mcmc <- initialize_mcmc(Y, tt, nt, B, X, ldim = ldim)
microbenchmark::microbenchmark(result <- run_mcmc(Y, X, B, tt, penalty, ldim, 10000, 100, 1, init_mcmc), times = 1)

scale <- init_mcmc$alpha
rho <- .1
C_rho <- scale * rho * rep(1, ldim) %*% t(rep(1, ldim)) + scale * (1 - rho) * diag(ldim)
phi_mat <- matrix(0, nrow = nreg^2, ldim)
log_phi <- 0
for (r in 1:nreg) {
  for (rp in 1:nreg) {
    # x <- phi[r, rp, ]
    x <- result$samples$phi[[1000]][rp,r,]
    phi_mat[(r - 1) * nreg + rp, ] <- x
    log_phi <- log_phi + mvtnorm::dmvnorm(x, rep(0, ldim), sigma = C_rho, log = TRUE)
  }
}
sum(mvtnorm::dmvnorm(phi_mat, rep(0, ldim), sigma = C_rho, log = TRUE))
pairs(phi_mat)

efunc <- 1
plot(B %*% result$samples$lambda[,efunc,5001], type = "l")
for (i in 5001:10000) {
  lines(B %*% result$samples$lambda[,efunc,i])
}
lines(eigenfuncs[,efunc], col = "red")

r <- 2
i <- 4
plot(Y[((i - 1) * nt + 1):(i * nt),r])
seqr <- ((i - 1) * nreg + 1):(i * nreg)
for (i in 1:1000) {
  tmpsum <- numeric(nt)
  for (l in 1:ldim) {
    tmpsum <- tmpsum + B %*% result$samples$lambda[,l, i] %*% t(result$samples$phi[[i]][r,,l]) %*% result$samples$eta[seqr, l, i]
  }
  lines(tmpsum, col = "blue")
}
hist(result$samples$rho)

iter <- 1000
l <- 1
r <- 1
seqr <- seq(from = r, to = nsub * nreg, by = nreg)
xb <- X %*% result$samples$beta[((r - 1) * d + 1):(r * d), l, iter]
# sum(dgamma(result$samples$xi_eta[,,500], 1, rate = 1, log = TRUE))
var(result$samples$eta[seqr,l,iter] - xb)
result$samples$nu[iter]
delta_eta_cumprod <- array(0, dim = c(nreg, ldim, 10000))
for (i in 1:10000) {
  initd <- cumprod(result$samples$delta_eta[1,,i])
  delta_eta_cumprod[,1,i] <- cumprod(result$samples$delta_eta[,1,i])
  for (l in 2:ldim) {
    delta_eta_cumprod[,l,i] <- cumprod(result$samples$delta_eta[,l,i]) * initd[l - 1]
    
  }
}
plot(1 / delta_eta_cumprod[1,1,])
abline(h = 1 / init_mcmc$preceta[1,1])
quantile(1 / delta_eta_cumprod[1,1,1:10000], c(.025, .975))

1 / delta_eta_cumprod[,,1000]

hist(1 / delta_eta_cumprod[5,2,])
abline(h = 1 / init_mcmc$preceta[5,2])


abline(v = 1 / var(result$samples$eta[seqr,1,1000]))


iter <- 2
l <- 1
r <- 1
delta_eta_cumprod <- matrix(0, nreg, ldim)
this_delta <- result$samples$delta_eta[,,iter]
this_delta[r, l] <- 1
initd <- cumprod(this_delta[1,])
delta_eta_cumprod[,1] <- cumprod(this_delta[,1])
for (l in 2:ldim) {
  delta_eta_cumprod[,l] <- cumprod(this_delta[,l]) * initd[l - 1]
}
# for (rp in 1:r) {
#   for (lp in 1:l) {
#     delta_eta_cumprod[rp, lp] <- 0
#   }
# }
tmpsum <- 0
for (r in 1:nreg) {
  for (l in 1:ldim) {
    seqr <- seq(from = r, to = nsub * nreg, by = nreg)
    add_this <- t(result$samples$eta[seqr,l,iter - 1]) %*% 
      (result$samples$eta[seqr,l,iter - 1])
    multiply_by <- delta_eta_cumprod[r, l]
    tmpsum <- tmpsum + add_this * multiply_by
    print(result$samples$eta[seqr,l,iter][1:4])
    print(paste0("r = ", r, "  l = ", l, "  add_this = ", add_this, "  multiply_by = ", multiply_by))
  }
}
# microbenchmark::microbenchmark(result <- run_mcmc(Y, X, B, tt, penalty, ldim, 1, 100, 1, init_mcmc), times = 1)
tmpsum
drate = 1 + .5 * tmpsum
dshape = 1 + .5 * nsub * nreg * ldim
rgamma(20, shape = dshape, rate = drate)

microbenchmark::microbenchmark(result <- run_mcmc(Y, X, B, tt, penalty, ldim, 1000, 100, 1, init_mcmc), times = 1)

### Checking posterior measure
iter <- 1
delta_eta_cumprod <- matrix(0, nreg, ldim)
initd <- cumprod(result$samples$delta_eta[1,,iter])
delta_eta_cumprod[,1] <- cumprod(result$samples$delta_eta[,1,iter])
for (l in 2:ldim) {
  delta_eta_cumprod[,l] <- cumprod(result$samples$delta_eta[,l,iter]) * initd[l - 1]
}

tmpsum <- 0
for (r in 1:1) {
  for (l in 1:1) {
    seqr <- seq(from = r, to = nsub * nreg, by = nreg)
    tmpsum <- tmpsum + sum(dnorm(result$samples$eta[seqr,l,iter],
          mean = 0, sd = 1 / sqrt(delta_eta_cumprod[r, l]), log = TRUE))
    tmpsum <- tmpsum + dgamma(result$samples$delta_eta[r, l, iter], shape = 2, rate = 1, log = TRUE)
  }
}
tmpsum
1 / delta_eta_cumprod
plot(1 / result$samples$delta_eta[1,1,])
abline(h = 1 / init_mcmc$delta_eta[1,1], col = "red")
quantile(1 / result$samples$delta_eta[1,1,], c(.025, .975))




plot(1 / result$samples$delta_eta[1,1,])
max_iter <- which.min(result$samples$delta_eta[1,1,])
iter <- max_iter
sum(dnorm(result$samples$eta[seqr,l,iter],
          mean = 0, sd = 1 / sqrt(result$samples$delta_eta[1,1,iter]), log = TRUE))


z1 <- rnorm(1000, sd = 1/sqrt(100))
z2 <- rnorm(1000, mean = .1, sd = 1 / sqrt(100))
1 / 1000 * sum(z2 > z1)
mean(pnorm(z2 - z1, sd = sqrt(1 / 100 + 1 / 100)))
pnorm(.1, mean = 0, sd = sqrt(2 / 100))
myt <- numeric(10000)
stde <- numeric(100)
for (i in 1:10000) {
  z1 <- rnorm(100, sd = 1)
  z2 <- rnorm(100, mean = .5, sd = 1)
  # myt[i] <- t.test(z1, z2, alternative = "greater")$p.value
  # stde[i] <- t.test(z1, z2, alternative = "greater")$stderr
  myt[i] <- pnorm(mean(z2) - mean(z1), mean = 0, sd = sqrt(2/100))
}
median(myt)
hist(1 - myt)
mean(stde)
