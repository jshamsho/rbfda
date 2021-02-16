library(mgcv)
nt <- 100
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 9
nsub <- 50
ndf <- 25
d <- 2
sqexp <- function(t, tp) {
  exp(-(t - tp)^2 / .2)
}
kern <- matrix(0, nrow = nt, ncol = nt)
for (i in 1:nt) {
  for (j in 1:nt) {
    kern[i, j] <- 5 * sqexp(tt[i], tt[j])
  }
}
Y <- matrix(rnorm(nsub * nt * nreg, sd = .05), nrow = nsub * nt, ncol = nreg)
Y1 <- mvrnorm(nsub * nreg, mu = rep(0, nt), Sigma = kern + .05 * diag(nt)) %>% t() %>% matrix(nrow = nsub * nt, ncol = nreg)
# Y[1,1] <- NA
# Y[1,1:nreg] <- NA
# Y[2,1:nreg] <- NA
X <- matrix(rnorm(nsub * d), nrow = nsub, ncol = d)
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "tp", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
matplot(tt, B, xlab = "time", ylab = "spline", type = "l")
microbenchmark::microbenchmark(result <- run_mcmc(Y1, X, B, tt, penalty, 4, 1000, 100, 1), times = 1)
pcomp <- princomp(Y1)
matplot(pcomp$loadings[,1:3], type = "l")
pracma::trapz(tt, pcomp$loadings[,1] * pcomp$loadings[,3])
pracma::trapz(tt, pcomp$loadings[,1] * pcomp$loadings[,1])

l2norms <- numeric(3)
u <- matrix(0, nrow = nt, ncol = 3)
proj <- matrix(0, nrow = nt, ncol = 3)
orthf <- matrix(0, nrow = nt, ncol = 3)
u[,1] <- pcomp$loadings[,1]
orthf[,1] <- u[,1] / sqrt(pracma::trapz(tt, u[,1] * u[,1]))
proj[,1] <- pracma::trapz(tt, pcomp$loadings[,2] * pcomp$loadings[,1]) /
  pracma::trapz(tt, u[,1] * u[,1]) * pcomp$loadings[,1]
u[,2] <- pcomp$loadings[,2] - proj[,1]
orthf[,2] <- u[,2] / sqrt(pracma::trapz(tt, u[,2] * u[,2]))
proj[,2] <- pracma::trapz(tt, u[,1] * pcomp$loadings[,3]) / pracma::trapz(tt, u[,1] * u[,1]) * u[,1]
proj[,3] <- pracma::trapz(tt, u[,2] * pcomp$loadings[,3]) / pracma::trapz(tt, u[,2] * u[,2]) * u[,2]

u[,3] <- pcomp$loadings[,3] - proj[,2] - proj[,3]
orthf[,3] <- u[,3] / sqrt(pracma::trapz(tt, u[,3] * u[,3]))
pracma::trapz(tt, orthf[,2] * orthf[,3])
matlines(orthf, type = "-o")

myfunc <- function(p, tau) {
  lambda_j <- abs(rt(p, 1))
  k_vec <- 1 / (1 + tau^2 * p *lambda_j^2)
  m_eff <- sum(1 - k_vec)
  m_eff
}
hist(sapply(1:10000, function(i) myfunc(20, .1)))

p <- 6
l <- 3
rho <- .8
corrmat <- rho * outer(rep(1, l), rep(1, l)) + diag(1-rho, l)
M <- mvrnorm(10, mu = rep(0, l), Sigma = corrmat)
M
plot(M[,1], M[,2])
kronecker(corrmat, diag(p))


sigma <- 1:5
kronecker(diag(sigma), diag(10)) # correct

kronecker(diag(10), diag(sigma))
diag(1:5)
