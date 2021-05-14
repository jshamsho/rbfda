library(mgcv)
library(rbfda)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/simulate_data.R")
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/testing.R")

nsub <- 100
nt <- 60
nreg <- 5
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
init_mcmc <- initialize_mcmc_partial(sim_data$Y, tt, B, X, ldim = ldim)
plot(init_mcmc$psi[,1])
result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
                   penalty = penalty, ldim = ldim, iter = 5000, burnin = 1000,
                   thin = 1, init_ = init_mcmc, covstruct = "partial")
result$samples$phi[[1]][,,4] %>% tcrossprod()
pvals <- get_pvals_partial(result)
hist(pvals)
abline(v = median(pvals))
plot(pvals, type = "l")
median(pvals)
get_pmin(pvals)


r <- 1
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
for (i in 1:5000) {
  lines(B %*% result$samples$lambda[,efunc,i])
  evec[i] <- (B %*% result$samples$lambda[,efunc,i])[30]
}
lines(sim_data$psi[,efunc], col = "red")
lines(init_mcmc$psi[,efunc], col = "green")
lines(B %*% apply(result$samples$lambda[,efunc,],1,mean), col = "blue")
sum((init_mcmc$psi[,efunc] - sim_data$psi[,efunc])^2)
sum((B %*% apply(result$samples$lambda[,efunc,], 1, median) - sim_data$psi[,efunc])^2)
r <- 2
i <- 2
plot(sim_data$Y[((i - 1) * nt + 1):(i * nt),r])
seqr <- ((i - 1) * nreg + 1):(i * nreg)
for (i in 201:5000) {
  tmpsum <- numeric(nt)
  for (l in 1:init_mcmc$npc) {
    tmpsum <- tmpsum + B %*% result$samples$lambda[,l, i] %*% 
      t(result$samples$phi[[i]][r,,l]) %*% 
      result$samples$eta[seqr, l, i]
  }
  lines(tmpsum, col = "blue")
}

delta_eta_cumprod <- array(0, dim = c(nreg, init_mcmc$npc, 5000))
for (i in 1:5000) {
  initd <- cumprod(result$samples$delta_eta[1,,i])
  delta_eta_cumprod[,1,i] <- cumprod(result$samples$delta_eta[,1,i])
  for (l in 2:init_mcmc$npc) {
    delta_eta_cumprod[,l,i] <- cumprod(result$samples$delta_eta[,l,i]) * initd[l - 1]
  }
}
r <- 1
l <- 1
plot(1 / delta_eta_cumprod[r,l,201:5000])
abline(h = 1 / init_mcmc$prec_eta[r,l], col = "red")
abline(h = ((ldim - l + 1) * 1 / r)^2, col = "blue")
1 / init_mcmc$prec_eta[r, l]
quantile(1 / delta_eta_cumprod[r,l,1:1000], c(.025, .975))
((ldim - l + 1) * 1 / r)^2


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
iter <- 2
tmpsum <- 0
for (r in 1:nreg) {
  for (l in 1:init_mcmc$npc) {
    seqr <- seq(from = r, to = nsub * nreg, by = nreg)
    add_this <- t(result$samples$eta[seqr,l,iter - 1]) %*% 
      (result$samples$eta[seqr,l,iter - 1])
    multiply_by <- delta_eta_cumprod[r, l, iter - 1]
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


sample_tau <- function(y, mu) {
  n <- length(y)
  y_centered <- y - mu
  return(1 / rgamma(1, shape = 2 + .5 * n, rate = 2 + .5 * t(y_centered) %*% y_centered))
}

sample_mu <- function(y, tau) {
  n <- length(y)
  Q <- (n / tau + 0)
  b <- sum(y) / tau
  mu <- rnorm(1, b / Q, 1 / sqrt(Q))
}
sdseq <- seq(from = .2, to = 5, by = .5)
nseq <- nseq <- rep(30, 30)
s1 <- numeric(length(nseq))
s2 <- numeric(length(nseq))
for (n in 1:length(nseq)) {
  y <- rnorm(nseq[n], mean = .5, sd = 2)
  iter <- 15000
  tau_c <- numeric(iter)
  mu_c <- numeric(iter)
  test1 <- numeric(iter)
  test2 <- numeric(iter)
  tau <- 1
  mu <- 1
  for (i in 1:iter) {
    # tau <- sample_tau(y, mu)
    tau <- var(y)
    # mu <- sample_mu(y, tau)
    mu <- mean(y)
    tau_c[i] <- tau
    mu_c[i] <- mu
  }
  
  for (i in 1:iter) {
    y_centered <- y - mu_c[i]
    fakey_h0 <- rnorm(nseq[n], mean = 0, sqrt(tau_c[i]))
    fakey_h0_centered <- fakey_h0
    test1[i] <- t(y_centered) %*% y_centered / sqrt(tau_c[i])
    test2[i] <- t(fakey_h0_centered) %*% fakey_h0_centered / sqrt(tau_c[i])
    # test3[i] <- t(y) %*% y
  }
  
  
  s1[n] <- sum(test1 > test2) / iter
  s2[n] <- t.test(y)$p.value
}
plot(s1, s2)
abline(a = 0, b = 1)

pdim <- 4
rho <- .5
C11 <- rho * rep(1, pdim) %*% t(rep(1, pdim)) + (1 - rho) * diag(pdim)

cov1 <- array(0, dim = c(nreg * init_mcmc$npc, nreg * init_mcmc$npc, 1000))
for (i in 1:1000) {
  A1 <- reshape_nreg(result$samples$eta[,,i], nsub, nreg)
  cov1[,,i] <- cov(A1)
}
hist(cov1[15,22,])

x <- seq(from = 0, to = 5, length.out = 1000)
f1 <- .001 * (x - 5)^4
f1[1:100] <- 0
f1[1000] <- 1e-6
f2 <- .95 * f1
f2[1:120] <- 0
f3 <- .9 * f2
f3[1:140] <- 0

plot(x, f1, type = "l")
lines(x, f2, col = "red")
lines(x, f3, col = "green")

A1 <- diag(f2 / f1)
v <- which(is.nan(diag(A1)))
A1[v, v] <- 0
A2 <- diag(f3 / f2)
v <- which(is.nan(diag(A2)))
A2[v, v] <- 0

plot(x, f1, type = "l")
lines(x, A1 %*% f1, col = "red")
lines(x, A2 %*% f2, col = "green")
lines(x, A2 %*% A1 %*% f1, col = "blue")


x <- seq(from = 0, to = 5, length.out = 1000)
f1 <- .001 * (x - 5)^4
f1[1:100] <- 0
f2 <- .00004 * (x - 5)^6
f2[1:100] <- 0
plot(x, f1, type = "l")
lines(x, f2, col = "red")

crit <- which.max(f1)
mp <- f2[crit:1000] / f1[crit:1000]
mp[length(mp)] <- 0
mp_ext <- c(rep(0, crit - 1), mp)
plot(x, f1, type = "l")
lines(x, diag(mp_ext) %*% f1, col = "red")

