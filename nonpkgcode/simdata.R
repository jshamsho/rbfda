library(mgcv)
nt <- 100
tt <- seq(from = 0, to = 1, length.out = nt)
nreg <- 5
nsub <- 50
ndf <- 25
sqexp <- function(t, tp) {
  exp(-(t - tp)^2)
}
kern <- matrix(0, nrow = nt, ncol = nt)
for (i in 1:nt) {
  for (j in 1:nt) {
    kern[i, j] <- sqexp(tt[i], tt[j])
  }
}
Y <- matrix(0, nrow = nsub * nt, ncol = nreg)
Y[1,1] <- NA
Y[1,1:nreg] <- NA
Y[2,1:nreg] <- NA
X <- matrix(1, nrow = nsub)
basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "tp", m = 2), data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
matplot(tt, B, xlab = "time", ylab = "spline", type = "l")
run_mcmc(Y, X, B, tt, penalty, 4, 1000, 100, 1)
