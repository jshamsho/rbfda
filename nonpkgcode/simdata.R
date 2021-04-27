library(mgcv)
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/initialize_mcmc.R")
source("/Users/johnshamshoian/Documents/R_projects/rbfda/nonpkgcode/simulate_data.R")

nt <- 60
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 4
nsub <- 300
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
      eta[(i - 1) * nreg + r, l] <- rnorm(1, sd = (ldim - l + 1) * 1 / r)
    }
  }
}
phi <- array(dim = c(nreg, nreg, ldim))
for (l in 1:ldim) {
  # phi[,,l] <- diag(nreg)
  # phi[,,l] <- rWishart(1, 10, diag(nreg))
  phi[,,l] <- pracma::randortho(nreg)
  # phi[,,l] <- orthmat
  # phi[,,l] <- matrix(rnorm(nreg * nreg), nreg, nreg)
  # phi[,,l] <- Re(eigen(phi[,,l])$vectors)
}
# sum(mvtnorm::dmvnorm(phi_mat, rep(0, ldim), sigma = C_rho, log = TRUE))
# cov(phi_mat)
# cor(phi_mat)
# log_phi
# c_inv <- solve(cor(phi_mat))
# rate <- 0
# for (i in 1:(nreg * nreg)) {
#   rate <- rate + phi_mat[i,] %*% c_inv %*% phi_mat[i,]
# }
# 1 / rgamma(10, shape = .5 + .5 * nreg * nreg * 4, rate = .5 + .5 * rate)
theta1 <- matrix(0, nrow = nsub, ncol = nreg * ldim)
etalong <- reshape_nreg(eta, nsub, nreg)
for (i in 1:nsub) {
  for (j in 1:nreg) {
    idx1 <- (j - 1) * ldim + 1
    idx2 <- j * ldim
    theta1[i, idx1:idx2] <- phi[,,j] %*% etalong[i, idx1:idx2]
  }
}
Y <- matrix(0, nsub * nt, ncol = nreg)
for (i in 1:nsub) {
  for (l in 1:ldim) {
    Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + 
      eigenfuncs[,l] %*% t(eta[((i - 1) * nreg + 1):(i * nreg), l]) %*% t(phi[,,l])
  }
  Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + c(rnorm(nt * nreg, sd = .1))
}

rho <- .9
cov_rho <- rho * rep(1, nreg * ldim) %*% t(rep(1, nreg * ldim)) + 
  (1 - rho) * diag(nreg * ldim)
theta1 <- MASS::mvrnorm(nsub, mu = rep(0, nrow(cov_rho)), Sigma = cov_rho)
Y <- matrix(0, nsub * nt, ncol = nreg)
for (i in 1:nsub) {
  for (l in 1:ldim) {
    idx1 <- (l - 1) * nreg + 1
    idx2 <- l * nreg
    Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + 
      eigenfuncs[, l] %*% t(theta1[i, idx1:idx2])
  }
  Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + c(rnorm(nt * nreg, sd = .1))
}

rho1 <- .5
rho2 <- .2
cov_rho <- matrix(0, nreg * ldim, nreg * ldim)
for (i in 1:ldim) {
  idx1 <- (i - 1) * nreg + 1
  idx2 <- i * nreg
  cov_rho[idx1:idx2, idx1:idx2] <- 3 * (ldim - i + 1) * (rho1 * rep(1, nreg) %*% 
    t(rep(1, nreg)) + 
    (1 - rho1) * diag(nreg))
}
matrows <- 1:(ldim * nreg)
for (i in 1:(ldim - 1)) {
  matrows <- matrows[-(1:nreg)]
  matcols <- ((i - 1) * nreg + 1):(i * nreg)
  for (j in matrows) {
    for (jp in matcols) {
      cov_rho[j, jp] <- rho2 * sqrt(cov_rho[j, j]) * sqrt(cov_rho[jp, jp])
      cov_rho[jp, j] <- cov_rho[j, jp]
    }
  }
  # for (j in matrows) {
    # cov_rho[j, i] <- rho2 * sqrt(cov_rho[])
  # }
}
cov_rho
theta1 <- MASS::mvrnorm(nsub, mu = rep(0, nrow(cov_rho)), Sigma = cov_rho)
Y <- matrix(0, nsub * nt, ncol = nreg)
for (i in 1:nsub) {
  for (l in 1:ldim) {
    idx1 <- (l - 1) * nreg + 1
    idx2 <- l * nreg
    Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + 
      eigenfuncs[, l] %*% t(theta1[i, idx1:idx2])
  }
  Y[((i - 1) * nt + 1):((i) * nt), ] <- Y[((i - 1) * nt + 1):((i) * nt), ] + c(rnorm(nt * nreg, sd = .1))
}
# 
# eigenfunc <- matrix(0, 60, 1000)
# var11 <- numeric(1000)
# for (b in 1:1000) {
#   print(b)
#   Yboot <- matrix(0, nt * nsub, nreg)
#   for (i in 1:nsub) {
#     this_sample <- sample(1:nsub, replace = TRUE)
#     Yboot[((i -1) * nt + 1):(i * nt), ] <- Y[((this_sample[i] -1) * nt + 1):(this_sample[i] * nt), ]
#   }
#   this_init <- initialize_mcmc(Yboot, tt, nt, B, X, ldim = 4)
#   eigenfunc[,b] <- B %*% this_init$lambda[,1]
#   var11[b] <- 1 / this_init$preceta[2,1]
# }
# for (b in 1:1000) {
#   if (sum((eigenfunc[,b] + eigenfunc[,1])^2) < sum((eigenfunc[,b] - eigenfunc[,1])^2)) {
#     eigenfunc[,b] <- -eigenfunc[,b]
#   }
# }
# matlines(-eigenfunc)
# Y <- sim_partial(nt, nsub, nreg, ldim = 5)
nsub <- 200
nt <- 60
nreg <- 4
ldim <- 4
sim_data <- sim_partial(nt, nsub, nreg, ldim)
sim_data <- sim_partial_cs(nt, nsub, nreg, ldim, rho1 = .8)
# sim_data <- sim_non_partial(nt, nsub, nreg, ldim, rho1 = .8, rho2 = .2)
X <- cbind(rep(1, nsub), matrix(rnorm(nsub * (d - 1)), nrow = nsub, ncol = d - 1))
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
init_mcmc <- initialize_mcmc_partial(sim_data$Y, tt, B, X, ldim = ldim)
init_mcmc <- initialize_mcmc_weak(Y, tt, B, X)
plot(init_mcmc$psi[,1])
result <- run_mcmc(sim_data$Y, X, B, tt, penalty, init_mcmc$npc, 5000, 1000, 1, init_mcmc)

gamma1 <- function(a, b) {
  gamma((b + 1) / 2) * gamma(a / 2) / (sqrt(pi) * gamma((a + b) / 2))
}
get_pval <- function(eta, nreg, iter) {
  p <- ncol(eta) * nreg
  nsub <- nrow(eta) / nreg
  eta1 <- reshape_nreg(eta[,,iter], nsub, nreg)
  cormat <- cor(eta1)
  dp <- sum(abs(cormat[upper.tri(cormat)]))
  # return(dp)
  mu <- p * (p - 1) / 2 * gamma1(nsub - 1, 1)
  tausq <- p * (p - 1) / 2 * (1 / (nsub - 1) - gamma1(nsub - 1, 1)^2)
  teststat <- (dp - mu) / sqrt(tausq)
  return(1 - pnorm(teststat))
}

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

chisqversion <- function(x, nreg, ldim) {
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
  print(paste0("loglambda = ", loglambda))
  rho <- 1 - (2 * (p^3 - ldim * nreg^3) + 9 * (p^2  - ldim * nreg^2)) /
    (6 * (nsub + 1) * (p^2 - ldim * nreg^2))
  f <- .5 * (p^2 - ldim * nreg^2)
  print(paste0("df = ", f))
  print(paste0("teststat = ", -2 * rho * loglambda))
  return(pchisq(-2 * rho * loglambda, df = f, lower.tail = FALSE))
}

normalversion <- function(x, nreg, ldim) {
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
iter <- 2000
theta2 <- get_theta(result$samples$eta[,,iter], result$samples$phi[[iter]], nreg)
theta2 <- get_theta(eta, phi, nreg)

chisqversion(theta2, nreg, ldim)
normalversion(theta2, nreg, ldim)
pvals <- sapply(501:5000, function(iter) {theta2 <- get_theta(result$samples$eta[,,iter], result$samples$phi[[iter]], nreg);
                normalversion(theta2, nreg, ldim)}
)
hist(pvals)
plot(pvals, type = "l")
scale <- init_mcmc$alpha
rho <- .1
C_rho <- scale * rho * rep(1, init_mcmc$npc) %*% t(rep(1, init_mcmc$npc)) + scale * (1 - rho) * diag(init_mcmc$npc)
phi_mat <- matrix(0, nrow = nreg^ 2, init_mcmc$npc)
log_phi <- 0
for (r in 1:nreg) {
  for (rp in 1:nreg) {
    # x <- phi[r, rp, ]
    x <- result$samples$phi[[1000]][rp,r,]
    phi_mat[(r - 1) * nreg + rp, ] <- x
    log_phi <- log_phi + mvtnorm::dmvnorm(x, rep(0, init_mcmc$npc), sigma = C_rho, log = TRUE)
  }
}
sum(mvtnorm::dmvnorm(phi_mat, rep(0, init_mcmc$npc), sigma = C_rho, log = TRUE))
pairs(phi_mat)
hist(result$samples$rho)
eta_mat <- matrix(0, nrow = nsub, ncol = ldim * nreg)
counter <- 1
for (l in 1:ldim) {
  for (r in 1:nreg) {
    col_ind <- (l - 1) * nreg + r
    seqr <- seq(from = r, to = nsub * nreg, by = nreg)
    eta_mat[, col_ind] <- result$samples$eta[seqr,l,1000]
    counter <- counter + 1
  }
}

cov(eta_mat)
efunc <- 1
plot(B %*% result$samples$lambda[,efunc,100], type = "l")
evec <- numeric(500)
for (i in 1:5000) {
  lines(B %*% result$samples$lambda[,efunc,i])
  evec[i] <- (B %*% result$samples$lambda[,efunc,i])[30]
}
lines(-sim_data$eigenfuncs[,efunc], col = "red")
lines(B %*% init_mcmc$lambda[,efunc], col = "green")
lines(B %*% apply(result$samples$lambda[,efunc,],1,mean), col = "blue")
sum((B %*% init_mcmc$lambda[,efunc] + sim_data$eigenfuncs[,efunc])^2)
sum((B %*% apply(result$samples$lambda[,efunc,], 1, median) + sim_data$eigenfuncs[,efunc])^2)
r <- 3
i <- 2
plot(sim_data$Y[((i - 1) * nt + 1):(i * nt),r])
seqr <- ((i - 1) * nreg + 1):(i * nreg)
for (i in 201:5000) {
  tmpsum <- numeric(nt)
  for (l in 1:init_mcmc$npc) {
    tmpsum <- tmpsum + B %*% result$samples$lambda[,l, i] %*% t(result$samples$phi[[i]][r,,l]) %*% result$samples$eta[seqr, l, i]
  }
  lines(tmpsum, col = "blue")
}
# hist(result$samples$rho)

iter <- 1000 
l <- 4
r <- 3
seqr <- seq(from = r, to = nsub * nreg, by = nreg)
xb <- X %*% result$samples$beta[((r - 1) * d + 1):(r * d), l, iter]
# sum(dgamma(result$samples$xi_eta[,,500], 1, rate = 1, log = TRUE))
var(result$samples$eta[seqr,l,iter] - xb)
result$samples$nu[iter]
delta_eta_cumprod <- array(0, dim = c(nreg, init_mcmc$npc, 1000))
for (i in 1:1000) {
  initd <- cumprod(result$samples$delta_eta[1,,i])
  delta_eta_cumprod[,1,i] <- cumprod(result$samples$delta_eta[,1,i])
  for (l in 2:init_mcmc$npc) {
    delta_eta_cumprod[,l,i] <- cumprod(result$samples$delta_eta[,l,i]) * initd[l - 1]
  }
}
r <- 2
l <- 1
plot(1 / delta_eta_cumprod[r,l,201:1000])
abline(h = 1 / init_mcmc$preceta[r,l])
abline(h = ((ldim - l + 1) * 1 / r)^2, col = "red")
1 / init_mcmc$preceta[r, l]
quantile(1 / delta_eta_cumprod[r,l,1:1000], c(.025, .975))
((ldim - l + 1) * 1 / r)^2

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

