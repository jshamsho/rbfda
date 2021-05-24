args = commandArgs(trailingOnly=TRUE)
set.seed(args[1])
print(paste("Running", args[1]))
library(mgcv)
library(rrbfda)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
foreach(myseed = 1:5, .packages = c("rrbfda", "mgcv")) %dopar% {
  print(paste0("Working on seed ", myseed))
  set.seed(myseed)
  nsub <- 200
  nt <- 60
  nreg <- 6
  ldim <- 4
  ndf <- 15
  iterations <- 1000
  thin <- 1
  burnin <- 0
  tt <- seq(from = 0, to = 1, length.out = ntt)
  sim_data <- rrbfda::sim_weak(ntt, nsub, nreg, ldim = ldim)
  X <- cbind(rep(1, nsub))
  basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2),
                              data.frame(tt), absorb.cons = FALSE)
  B <- basisobj[[1]]$X
  penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
  init_mcmc <- rrbfda::initialize_mcmc_weak(sim_data$Y, tt, B, X)
  ldim_est <- ncol(init_mcmc$lambda)
  microbenchmark::microbenchmark(result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
                                                    penalty = penalty, ldim = ldim_est, iter = iterations, burnin = burnin,
                                                    thin = thin, init_ = init_mcmc, covstruct = "weak"), times = 1)
  
  NULL
}
parallel::stopCluster(cl)
for (myseed in 1:300) {
  print(paste0("Working on seed ", myseed))
  set.seed(myseed)
  nsub <- 200
  nt <- 60
  nreg <- 6
  ldim <- 4
  ndf <- 15
  iterations <- 1000
  thin <- 1
  burnin <- 0
  tt <- seq(from = 0, to = 1, length.out = nt)
  sim_data <- sim_weak(nt, nsub, nreg, ldim = ldim)
  # sim_data <- sim_non_weak(nt, nsub, nreg, ldim)
  # sim_data <- sim_partial(nt, nsub, nreg, ldim = ldim)
  # sim_data <- sim_partial_cs(nt, nsub, nreg, ldim, rho1 = .8)
  # sim_data <- sim_non_partial(nt, nsub, nreg, ldim, rho1 = .8, rho2 = .6)
  X <- cbind(rep(1, nsub))
  basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2),
                              data.frame(tt), absorb.cons = FALSE)
  B <- basisobj[[1]]$X
  penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
  init_mcmc <- initialize_mcmc_weak(sim_data$Y, tt, B, X)
  ldim_est <- ncol(init_mcmc$lambda)
  microbenchmark::microbenchmark(result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
                                                    penalty = penalty, ldim = ldim_est, iter = iterations, burnin = burnin,
                                                    thin = thin, init_ = init_mcmc, covstruct = "weak"), times = 1)
  
  eigenfunc_summary <- array(0, dim = c(nt, ldim, 3))
  for (l in 1:ldim) {
    eigenfunc_summary[, l, ] <- get_posteigenfunc(result, l)
  }
  eigenvec_summary <- array(0, dim = c(nreg, nreg, 3))
  for (r in 1:nreg) {
    eigenvec_summary[, r, ] <- get_posteigenvec(result, r)
  }
  eigenval_summary <- array(0, dim = c(nreg, ldim_est, 3))
  for (r in 1:nreg) {
    for (l in 1:ldim_est) {
      eigenval_summary[r, l, 1:3] <- quantile(1 / result$samples$sigmasqeta[r, l, (burnin + 1):iterations],
                                              c(.5, .025, .975))
    }
  }
  mean_summary <- array(0, dim = c(nt, nreg, 3))
  for (r in 1:nreg) {
    mean_summary[,r,] <- get_posteriormean(result, r, 1)
  }
  noise_summary <- 1 / t(apply(result$samples$omega[, (burnin + 1):iterations],
                               1, quantile, c(.5, .975, .025)))
  colnames(noise_summary) <- c("50%", "2.5%", "97.5%")
  nu_summary <- quantile(result$samples$nu[(burnin + 1):iterations], quantile(c(.025, .5, .975)))
  
  pvals <- get_pvals_weak(result)
  
  simstats <- list(sim_data <- sim_data,
                   eigenfunc_summary = eigenfunc_summary,
                   eigenvec_summary = eigenvec_summary,
                   eigenval_summary = eigenval_summary,
                   mean_summary = mean_summary,
                   noise_summary = noise_summary,
                   pvals = pvals)
  
  outfile <- paste0("/u/home/j/jshamsh1/Documents/rrbfda/output/n200_simweak_fitweak
/n200_simweak_fitweakforeach", myseed, ".RData")
  save(simstats, file = outfile)
}

# 
library (foreach)

fn<-function(i)
{
  set.seed(i)
  y <- rnorm(1)
  return(y)
}

x<-foreach(i=1:10) %do% fn(i)
print(x)
# num <- 2
# eigenvec_summary <- get_posteigenvec(result, num)
# plot(eigenvec_summary[,1], type = "o", ylim = c(-.9, .9))
# points(-sim_data$phi[, num], type = "o", col = "blue")
# points(eigenvec_summary[,2], type = "o")
# points(eigenvec_summary[,3], type = "o")
# # 
# var_mats <- array(0, dim = c(nreg, ldim_est, iterations))
# for (i in 1:iterations) {
#   var_mats[,,i] <- 1 / result$samples$sigmasqeta[,,i]
# }
# r <- 1
# l <- 1
# plot(var_mats[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r, l])
# abline(h = sim_data$sigma_mat[r, l], col = "red")
# 
# num <- 3
# plot(eigenfunc_summary[,num,1], ylim = c(-.5, .5), type = "l")
# lines(sim_data$psi[,num], col = "red")
# lines(eigenfunc_summary[,num,2])
# lines(eigenfunc_summary[,num,3])
# sim_data$sigma_mat[r, l]
# quantile(var_mats[r,l,burnin:iterations], c(.025, .975))
# 
# par(mfrow = c(2,2))
# r <- 1; l <- 1
# plot(var_mats[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r, l])
# abline(h = sim_data$sigma_mat[r, l], col = "red")
# r <- 2; l <- 1
# plot(var_mats[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r, l])
# abline(h = sim_data$sigma_mat[r, l], col = "red")
# r <- 1; l <- 2
# plot(var_mats[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r, l])
# abline(h = sim_data$sigma_mat[r, l], col = "red")
# r <- 2; l <- 2
# plot(var_mats[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r, l])
# abline(h = sim_data$sigma_mat[r, l], col = "red")
# 
# eigenval_summary
# sim_data$sigma_mat
# pvals <- get_pvals_weak(result)
# hist(pvals)
# 
# cumprod(result$samples$delta_eta1[,5000]) %*% t(cumprod(result$samples$delta_eta2[,5000]))
# cumprod(init_mcmc$delta_eta1) %*% t(cumprod(init_mcmc$delta_eta2))
# init_mcmc$prec_eta
# pvals <- get_pvals_weak(result)
# hist(pvals)
# abline(v = median(pvals))
# plot(pvals, type = "l")
# median(pvals)
# get_pmin(pvals)
# 
# 
# b <- numeric(ndf)
# Q <- matrix(0, ndf, ndf)
# for (i in 1:nsub) {
#   b <- b + t(B) %*% sim_data$Y[((i - 1) * nt + 1):(i * nt),] %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1] -
#     t(B) %*% B %*% init_mcmc$lambda[,2:ldim] %*% t(init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 2:ldim]) %*% t(init_mcmc$phi) %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1]
#   Q <- Q + as.numeric(init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1] %*% t(init_mcmc$phi) %*% diag(init_mcmc$omega) %*% init_mcmc$phi %*% init_mcmc$eta[((i - 1) * nreg + 1):(i * nreg), 1]) * t(B) %*% B
# }
# lambda1 <- MASS::mvrnorm(n = 1, mu = solve(Q) %*% b, Sigma = solve(Q))
# plot(B %*% lambda1, type = "l")
# r <- 3
# seqr <- ((r - 1) * d):(r * d)
# plot.new()
# for (i in 201:5000) {
#   tmpsum <- numeric(nt)
#   for (l in 1:init_mcmc$npc) {
#     tmpsum <- tmpsum + B %*% result$samples$lambda[,l, i] %*% 
#       t(result$samples$phi[[i]][r,,l]) %*% 
#       (result$samples$beta[seqr, l, i]
#   }
#   lines(tmpsum, col = "blue")
# }
# 
# x_vec <- c(1, 1)
# meanf <- numeric(nt)
# for (l in 1:ldim) {
#   meanf <- meanf + B %*% result$samples$lambda[, l, 1] %*%
#     result$samples$phi[[1]][r,,l] %*% 
#     t(matrix(result$samples$beta[,l,1], d)) %*% x_vec
# }
# plot(tt, meanf, type = "l", ylim = c(-2,2))
# for (i in 500:5000) {
#   meanf <- numeric(nt)
#   for (l in 1:ldim) {
#     meanf <- meanf + B %*% result$samples$lambda[, l, i] %*%
#       result$samples$phi[[i]][r,,l] %*% 
#       t(matrix(result$samples$beta[,l,i], d)) %*% x_vec
#   }
#   lines(tt, meanf)
# }
# 
# efunc <- 1
# plot(B %*% result$samples$lambda[,efunc,100], type = "l", ylim = c(-.5,.5))
# evec <- numeric(5000)
# for (i in 500:5000) {
  # lines(B %*% result$samples$lambda[,efunc,i])
  # evec[i] <- (B %*% result$samples$lambda[,efunc,i])[30]
}
# lines(sim_data$psi[,efunc], col = "red")
# lines(init_mcmc$psi[,efunc], col = "green")
# lines(B %*% apply(result$samples$lambda[,efunc,],1,mean), col = "blue")
# sum((init_mcmc$psi[,efunc] - sim_data$psi[,efunc])^2)
# sum((B %*% apply(result$samples$lambda[,efunc,], 1, mean) - sim_data$psi[,efunc])^2)
# 
# r <- 4
# i <- 50
# plot(sim_data$Y[((i - 1) * nt + 1):(i * nt),r])
# seqr <- ((i - 1) * nreg + 1):(i * nreg)
# for (i in burnin:iterations) {
#   est <- B %*% result$samples$lambda[,, i] %*% 
#       t(result$samples$eta[seqr, , i]) %*%
#       result$samples$phi[r,,i]
#   lines(est, col = "green")
# }
# 
# num_iter <- 5000
# delta_eta_cumprod <- array(0, dim = c(nreg, ldim, num_iter))
# for (i in 1:num_iter) {
#   q <- c(result$samples$delta_eta1[,,i], t(result$samples$delta_eta2[,,i]))
#   
#   delta_eta_cumprod[,,i] <- full_var(q)
# }
# r <- 1
# l <- 1
# plot(1 / delta_eta_cumprod[r,l,1:num_iter], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r,l])
# abline(h = ((ldim - l + 1) * 1 / r)^2, col = "red")
# 1 / init_mcmc$prec_eta[r, l]
# quantile(1 / delta_eta_cumprod[r,l,1000:num_iter], c(.025, .975))
# ((ldim - l + 1) * 1 / r)^2
# 
# 
# num_iter <- 5000
# l <- 2
# r <- 1
# plot(1 / result$samples$sigmasqeta[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r,l])
# abline(h = ((ldim - l + 1) * 1 / r)^2, col = "red")
# 1 / init_mcmc$prec_eta[r, l]
# quantile(1 / result$samples$sigmasqeta[r,l,1:num_iter], c(.025, .975))
# ((ldim - l + 1) * 1 / r)^2
# shape_param <- 1 + .5 * result$samples$eta[seqr,1,5000] %*% result$samples$eta[seqr,1,5000]
# rate_param <- 1 + .5 * nsub
# rgamma(1, shape = shape_param, rate = rate_param)
# seqr <- seq(from = r, by = nreg, length.out = nsub)
# var(result$samples$eta[seqr,l,1])
# 
# 
# 
# 
# get_cumprod_coef <- function(input) {
#   delta <- numeric(length(input))
#   delta[1] <- input[1]
#   if (length(input) > 1) {
#     for (i in 2:length(input)) {
#       delta[i] <- input[i] / input[i - 1]
#     }
#   }
#   return(delta)
# }
# 
# mynmf <- nmf(init_mcmc$prec_eta, rank = 4, method = "lee")
# 
# delta_eta1 <- basis(mynmf)
# delta_eta2 <- coef(mynmf)
# init_mcmc <- initialize_mcmc_weak(sim_data$Y, tt, B, X, ldim = ldim)
# get_delta_eta_density(init_mcmc$delta_eta1, init_mcmc$delta_eta2,
#                       init_mcmc$eta,
#                       init_mcmc$beta,
#                       init_mcmc$xi_eta, X)
# init_mcmc$delta_eta1
# init_mcmc$delta_eta2
# 
# a_basis <- sqrt(delta_eta1[1,1] * delta_eta2[1,1]) / delta_eta1[1,1]
# a_coef <- sqrt(delta_eta1[1,1] * delta_eta2[1,1]) / delta_eta2[1,1]
# delta_eta1 <- a_basis * delta_eta1
# delta_eta2 <- a_coef * delta_eta2
# for (i in 1:4) {
#   delta_eta1[,i] <- get_cumprod_coef(delta_eta1[,i])
#   delta_eta2[i,] <- get_cumprod_coef(delta_eta2[i,])
# }
# xi_eta <- matrix(1, nreg * nsub, ldim)
# 
# cdim <- 4
# delta_eta1_cumprod <- matrix(0, nreg, cdim)
# delta_eta2_cumprod <- matrix(0, cdim, ldim)
# for (c in 1:cdim) {
#   delta_eta1_cumprod[, c] <- cumprod(init_mcmc$delta_eta1[,c])
#   delta_eta2_cumprod[c, ] <- cumprod(init_mcmc$delta_eta2[c,])
# }
# full_var <- delta_eta1_cumprod %*% delta_eta2_cumprod
# print(full_var)
# 
# get_density <- function(q) {
#   bigdim <- length(q)
#   density_vec <- numeric(bigdim)
#   cdim <- ncol(init_mcmc$delta_eta1)
#   delta_eta1 <- matrix(q[1:(nreg * cdim)], nreg, cdim)
#   delta_eta2 <- matrix(q[(nreg * cdim + 1):bigdim], cdim, ldim, byrow = TRUE)
#   delta_eta1_cumprod <- matrix(0, nreg, cdim)
#   delta_eta2_cumprod <- matrix(0, cdim, ldim)
#   for (c in 1:cdim) {
#     delta_eta1_cumprod[, c] <- cumprod(delta_eta1[,c])
#     delta_eta2_cumprod[c, ] <- cumprod(delta_eta2[c,])
#   }
#   full_var <- delta_eta1_cumprod %*% delta_eta2_cumprod
#   density <- 0
#   for (rp in 1:nreg) {
#     rseq <- seq(from = rp, by = nreg, length.out = nsub)
#     for (lp in 1:ldim) {
#       etavec <- init_mcmc$eta[rseq, lp]
#       betavec <- init_mcmc$beta[((rp - 1) * d + 1):(rp * d), lp]
#       etamean <- X %*% betavec
#       mysd <- full_var[rp, lp]^(-.5)
#       density <- density + sum(dnorm(etavec, etamean, sd = mysd, log = TRUE))
#       # density <- density + .5 * nsub * log(full_var[rp, lp]) -.5 *
#       # full_var[rp, lp] *
#       # t(etavec - etamean) %*% (etavec - etamean)
#     }
#   }
#   
#   return(-density)
# }
# q <- rep(1, min(ldim, nreg) * (nreg + ldim))
# get_density(q)
# z <- rnorm(35)
# .5 * z%*%z
# get_density_real <- function() {
#   density <- 0
#   for (r in 1:nreg) {
#     rseq <- seq(from = r, by = nreg, length.out = nsub)
#     for (l in 1:ldim) {
#       etavec <- init_mcmc$eta[rseq, l]
#       betavec <- init_mcmc$beta[((r - 1) * d + 1):(r * d), l]
#       etamean <- X %*% betavec
#       mysd <- init_mcmc$prec_eta[r, l]^(-.5)
#       density <- density + sum(dnorm(etavec, etamean, sd = mysd, log = TRUE))
#       
#     }
#   }
#   return(-density)
# }
# 
# num_iters <- 200
# samples <- matrix(0, length(q), num_iters)
# densities <- numeric(num_iters)
# # q <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# q <- rep(1, min(ldim, nreg) * (nreg + ldim))
# q_o <- q
# p_o <- rnorm(length(q), sd = sqrt(.1))
# for (iter in 1:num_iters) {
#   print(iter)
#   # p <- rnorm(length(q), sd = sqrt(.1))
#   p <- rnorm(1, sd = sqrt(.1))
#   # get_density_real()
#   step_size <- .01
#   for (i in 1:1) {
#     p <- p - step_size / 2 * numDeriv::grad(get_density, q, method = "Richardson")[1]
#     q[1] <- q[1] + step_size * p * .1
#     p <- p - step_size / 2 * numDeriv::grad(get_density, q, method = "Richardson")[1]
#   }
#   densities[iter] <- get_density(q)
#   if (runif(1) < exp(.1 * .5 * p %*% p - get_density(q) - .1 * .5 * p_o %*% p_o + get_density(q_o))) {
#     q_o <- q
#     p_o <- p
#   }
#   samples[, iter] <- q_o
# }
# plot(densities)
# full_var(samples[, 200])
# 
# 
# num_iters <- 25000
# samples <- matrix(0, length(q), num_iters)
# densities <- numeric(num_iters)
# MH_ratios <- numeric(num_iters)
# q <- rep(1, min(ldim, nreg) * (nreg + ldim))
# q0 <- q
# cdim <- min(ldim, nreg)
# bigdim <- length(q)
# for (iter in 1:num_iters) {
#   delta_eta1 <- matrix(q0[1:(nreg * cdim)], nreg, cdim)
#   delta_eta2 <- matrix(q0[(nreg * cdim + 1):bigdim], cdim, ldim, byrow = TRUE)
#   prop <- get_delta_eta_proposal(delta_eta1, delta_eta2)
#   q <- c(prop$delta_eta1, t(prop$delta_eta2))
#   MH_ratio <- exp(-get_density(q) + get_density(q0))
#   if (runif(1) < MH_ratio) {
#     q0 <- q
#   }
#   MH_ratios[iter] <- MH_ratio
#   densities[iter] <- get_density(q0)
#   samples[,iter] <- q0
#   print(densities[iter])
#   print(get_density(q))
# }
# plot(densities)
# full_var(samples[, 5000])
# print(-get_density(q))
# print(-get_density(q_o))
# dnorm(p, sd = sqrt(.1), log = TRUE) %>% sum()
# dnorm(p_o, sd = sqrt(.1), log = TRUE) %>% sum()
# q <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# init_mcmc$xi_eta <- matrix(1, nsub * nreg, ldim)
# get_delta_eta1_grad(init_mcmc$delta_eta1, init_mcmc$delta_eta2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, 0, 0)
# numDeriv::grad(get_density, q)
# get_density(q + c(-.01, rep(0, length(q) - 1)))
# result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
#                    penalty = penalty, ldim = ldim, iter = 1, burnin = 1000,
#                    thin = 1, init_ = init_mcmc, covstruct = "weak")
# q <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# get_density(q)
# prop <- get_delta_eta_proposal(init_mcmc$delta_eta1, init_mcmc$delta_eta2)
# q1 <- c(prop$delta_eta1_proposal, t(prop$delta_eta2_proposal))
# get_density(q1)
# full_var(q)
# full_var(q1)
# full_var <- function(q) {
#   r <- 1
#   l <- 1
#   c <- 1
#   bigdim <- length(q)
#   density_vec <- numeric(bigdim)
#   cdim <- max(ncol(init_mcmc$delta_eta1), nrow(init_mcmc$delta_eta2))
#   delta_eta1 <- matrix(q[1:(nreg * cdim)], nreg, cdim)
#   delta_eta2 <- matrix(q[(nreg * cdim + 1):bigdim], cdim, ldim, byrow = TRUE)
#   delta_eta1_cumprod <- matrix(0, nreg, cdim)
#   delta_eta2_cumprod <- matrix(0, cdim, ldim)
#   for (c in 1:cdim) {
#     delta_eta1_cumprod[, c] <- cumprod(delta_eta1[,c])
#     delta_eta2_cumprod[c, ] <- cumprod(delta_eta2[c,])
#   }
#   full_var <- delta_eta1_cumprod %*% delta_eta2_cumprod
#   return(full_var)
# }
# 
# full_prec <- function(delta_eta1, delta_eta2) {
#   nreg <- nrow(delta_eta1)
#   ldim <- ncol(delta_eta2)
#   cdim <- ncol(delta_eta1)
#   delta_eta1_cumprod <- matrix(0, nreg, cdim)
#   delta_eta2_cumprod <- matrix(0, cdim, ldim)
#   for (c in 1:cdim) {
#     delta_eta1_cumprod[, c] <- cumprod(delta_eta1[,c])
#     delta_eta2_cumprod[c, ] <- cumprod(delta_eta2[c,])
#   }
#   return(delta_eta1_cumprod %*% delta_eta2_cumprod)
# }
# deriv_full_var <- function(q) {
#   r <- 1
#   l <- 1
#   c <- 1
#   bigdim <- length(q)
#   density_vec <- numeric(bigdim)
#   cdim <- ncol(init_mcmc$delta_eta1)
#   delta_eta1 <- matrix(q[1:(nreg * cdim)], nreg, cdim)
#   delta_eta2 <- matrix(q[(nreg * cdim + 1):bigdim], cdim, ldim, byrow = TRUE)
#   delta_eta1[2,1] <- 1
#   delta_eta1_cumprod <- matrix(0, nreg, cdim)
#   delta_eta2_cumprod <- matrix(0, cdim, ldim)
#   for (c in 1:cdim) {
#     delta_eta1_cumprod[, c] <- cumprod(delta_eta1[,c])
#     delta_eta2_cumprod[c, ] <- cumprod(delta_eta2[c,])
#   }
#   small_var <- delta_eta1_cumprod[,1] %*% t(delta_eta2_cumprod[1,])
#   return(small_var)
# }
# full_var(q)
# numDeriv::grad(full_var, q)[2]
# deriv_full_var(q)
# result <- run_mcmc(response = sim_data$Y, design = X, basis = B, time = tt,
#                    penalty = penalty, ldim = ldim, iter = 1, burnin = 1000,
#                    thin = 1, init_ = init_mcmc, covstruct = "weak")
# cdim <- min(ldim, nreg)
# num_iter <- 1000
# num_steps <- 1
# samples <- matrix(0, cdim * (nreg + ldim), num_iter)
# step_size <- .05
# # d1 <- init_mcmc$delta_eta1
# # d2 <- init_mcmc$delta_eta2
# d1_container <- array(0, dim = c(nreg, cdim, num_iter))
# d2_container <- array(0, dim = c(nreg, cdim, num_iter))
# d1 <- matrix(1, nreg, cdim)
# d2 <- matrix(1, cdim, ldim)
# d1_prop <- d1
# d2_prop <- d2
# p0 <- list()
# p0[[1]] <- matrix(rnorm(nreg * cdim, sd = .1), nreg, cdim)
# p0[[2]] <- matrix(rnorm(cdim * ldim, sd = .1), cdim, ldim)
# densities <- numeric(num_iter)
# for (iter in 1:num_iter) {
#   for (c in 1:cdim) {
#     for (r in 1:nreg) {
#       q <- d1[r, c]
#       p <- rnorm(1, sd = .1)
#       d1_prop <- d1
#       for (i in 1:num_steps) {
#         mygrad <- get_delta_eta1_grad(d1, d2,
#                                       init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, r - 1, c - 1)
#         p <- p - step_size / 2 * mygrad
#         q <- q + step_size * p * .1
#         d1_prop[r, c] <- q
#         mygrad <- get_delta_eta1_grad(d1_prop, d2,
#                                       init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, r - 1, c - 1)
#         p <- p - step_size / 2 * mygrad
#       }
#       p <- -p
#       density_original <- get_delta_eta_density(
#         d1, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       density_new <- get_delta_eta_density(
#         d1_prop, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       H_old <- .1 * .5 * p0[[1]][r, c]^2 + density_original
#       H_new <- .1 * .5 * p^2 + density_new
#       if (runif(1) < exp(H_old - H_new)) {
#         p0[[1]][r, c] <- p
#         d1[r, c] <- d1_prop[r, c]
#       }
#       d1_container[r, c, iter] <- d1[r, c]
#     }
#   }
#   
#   for (c in 1:cdim) {
#     for (l in 1:ldim) {
#       q <- d2[c, l]
#       p <- rnorm(1, sd = .1)
#       d2_prop <- d2
#       for (i in 1:num_steps) {
#         mygrad <- get_delta_eta1_grad(d1, d2,
#                                       init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, c - 1, l - 1)
#         p <- p - step_size / 2 * mygrad
#         q <- q + step_size * p * .1
#         d2_prop[c, l] <- q
#         mygrad <- get_delta_eta1_grad(d1, d2_prop,
#                                       init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, c - 1, l - 1)
#         p <- p - step_size / 2 * mygrad
#       }
#       p <- -p
#       density_original <- get_delta_eta_density(
#         d1, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       density_new <- get_delta_eta_density(
#         d1, d2_prop, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       H_old <- .1 * .5 * p0[[2]][c, l]^2 + density_original
#       H_new <- .1 * .5 * p^2 + density_new
#       if (runif(1) < exp(H_old - H_new)) {
#         p0[[2]][c, l] <- p
#         d2[c, l] <- d2_prop[c, l]
#       }
#       d2_container[c, l, iter] <- d2[c, l]
#     }
#   }
#   
#   densities[iter] <- get_delta_eta_density(d1, d2,
#                                            init_mcmc$eta,
#                                            init_mcmc$beta,
#                                            init_mcmc$xi_eta, X)
#   print(densities[iter])
# }
# q0 <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# full_var(q0)
# prec_mats <- array(0, dim = c(nreg, ldim, num_iter))
# for (i in 1:num_iter) {
#   q1 <- c(d1_container[,,i], t(d2_container[,,i]))
#   prec_mats[,,i] <- full_var(q1)
# }
# init_mcmc$prec_eta
# plot(prec_mats[1,1,], type = "l")
# get_delta_eta_density(init_mcmc$delta_eta1, init_mcmc$delta_eta2,
#                       init_mcmc$eta,
#                       init_mcmc$beta,
#                       init_mcmc$xi_eta, X)
# 
# num_iter <- 5000
# samples <- matrix(0, cdim * (nreg + ldim), num_iter)
# d1_container <- array(0, dim = c(nreg, cdim, num_iter))
# d2_container <- array(0, dim = c(cdim, ldim, num_iter))
# d1 <- matrix(1, nreg, cdim)
# d2 <- matrix(1, cdim, ldim)
# d1_prop <- d1
# d2_prop <- d2
# densities <- numeric(num_iter)
# for (iter in 1:num_iter) {
#   for (c in 1:cdim) {
#     for (r in 1:nreg) {
#       q <- -1
#       while(q <= 0) {
#         q <- d1[r, c] + rnorm(1, .025)
#       }
#       d1_prop <- d1
#       d1_prop[r, c] <- q
#       density_original <- -get_delta_eta_density(
#         d1, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       density_new <- -get_delta_eta_density(
#         d1_prop, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       prior_original <- dgamma(d1[r, c], shape = 2, rate = 1, log = TRUE)
#       prior_new <- dgamma(d1_prop[r, c], shape = 2, rate = 1, log = TRUE)
#       p1 <- pnorm(d1[r, c])
#       p2 <- pnorm(d1_prop[r, c])
#       H_old <- density_original + prior_original - p1
#       H_new <- density_new + prior_new - p2
#       if (runif(1) < exp(H_new - H_old)) {
#         d1[r, c] <- d1_prop[r, c]
#       }
#       d1_container[r, c, iter] <- d1[r, c]
#     }
#   }
#   
#   for (c in 1:cdim) {
#     for (l in 1:ldim) {
#       q <- -1
#       while(q <= 0) {
#         q <- d2[c, l] + rnorm(1, .025)
#       }
#       d2_prop <- d2
#       d2_prop[c, l] <- q
#       density_original <- -get_delta_eta_density(
#         d1, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       density_new <- -get_delta_eta_density(
#         d1, d2_prop, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
#       prior_original <- dgamma(d2[c, l], shape = 2, rate = 1, log = TRUE)
#       prior_new <- dgamma(d1_prop[c, l], shape = 2, rate = 1, log = TRUE)
#       p1 <- pnorm(d2[c, l])
#       p2 <- pnorm(d2_prop[c, l])
#       H_old <- density_original + prior_original - p1
#       H_new <- density_new + prior_new - p2
#       if (runif(1) < exp(H_new - H_old)) {
#         d2[c, l] <- d2_prop[c, l]
#       }
#       d2_container[c, l, iter] <- d2[c, l]
#     }
#   }
#   densities[iter] <- get_delta_eta_density(d1, d2,
#                                            init_mcmc$eta,
#                                            init_mcmc$beta,
#                                            init_mcmc$xi_eta, X)
#   print(densities[iter])
# }
# 
# q0 <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# full_var(q0)
# prec_mats <- array(0, dim = c(nreg, ldim, num_iter))
# for (i in 1:num_iter) {
#   q1 <- c(d1_container[,,i], t(d2_container[,,i]))
#   prec_mats[,,i] <- full_var(q1)
# }
# init_mcmc$prec_eta
# plot(1 / prec_mats[1,1,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[1,1])
# prec_mats[,,5000]
# init_mcmc$xi_eta <- matrix(1, nreg * nsub, ldim)
# get_delta_eta_density(init_mcmc$delta_eta1, init_mcmc$delta_eta2,
#                       init_mcmc$eta,
#                       init_mcmc$beta,
#                       init_mcmc$xi_eta, X)
# 
# get_delta_eta_density(init_mcmc$delta_eta1, init_mcmc$delta_eta2,
#                       init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X)
# get_delta_eta_density(init_mcmc$delta_eta1, init_mcmc$delta_eta2,
#                       result$samples$eta[,,1000], result$samples$beta[,,1000], result$samples$xi_eta[,,1000], X)
# get_delta_eta_density(result$samples$delta_eta1[,,1000], result$samples$delta_eta2[,,1000],
#                       result$samples$eta[,,1000], result$samples$beta[,,1000], result$samples$xi_eta[,,1000], X)
# q0 <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# 
# q1 <- c(result$samples$delta_eta1[,,1000], t(result$samples$delta_eta2[,,1000]))
# full_var(q0)
# full_var(q1)
# init_mcmc$prec_eta
# 
# r <- 3
# l <- 2
# seqr <- seq(from = 1, by = nreg, length.out = nsub)
# sd0 <- 1 / sqrt(full_var(q0)[r,l])
# sd1 <- 1/sqrt(result$samples$xi_eta[seqr, l,1000] * full_var(q1)[r,l])
# sum(dnorm(result$samples$eta[seqr, 1, 1000], sd0, log = TRUE))
# sum(dnorm(result$samples$eta[seqr, 1, 1000], mean(sd1), log = TRUE))
# 
# mynmf <- nmf(1 / init_mcmc$prec_eta, rank = min(nreg, ldim))
# delta_eta1 <- basis(mynmf)
# delta_eta2 <- coef(mynmf)
# for (i in 1:min(nreg, ldim)) {
#   init_mcmc$delta_eta1[,i] <- get_cumprod_coef(delta_eta1[,i])
#   init_mcmc$delta_eta2[i,] <- get_cumprod_coef(delta_eta2[i,])
# }
# 
# 
# num_iter <- 10
# thin <- 1
# samples <- matrix(0, cdim * (nreg + ldim), num_iter)
# d1_container <- array(0, dim = c(nreg, cdim, num_iter))
# d2_container <- array(0, dim = c(cdim, ldim, num_iter))
# d1 <- matrix(1, nreg, cdim)
# d2 <- matrix(1, cdim, ldim)
# d1 <- init_mcmc$delta_eta1
# d2 <- init_mcmc$delta_eta2
# d1_prop <- matrix(1, nreg, cdim)
# d2_prop <- matrix(1, cdim, ldim)
# init_mcmc$xi_eta <- matrix(1, nsub * nreg, ldim)
# 
# for (i in 1:num_iter) {
#   print(i)
#   for (th in 1:thin) {
#     for (c in 1:4) {
#       
#       d1[,c] <- identify_delta(d1[,c])
#       d2[c,] <- identify_delta(d2[c,])
#       
#       sd1 <- .025 * d1[,c]
#       sd2 <- .025 * d2[c,]
#       d1_prop[,c] <- rtruncnorm(nreg, a = rep(0, nreg), mean = d1[,c], sd = sd1)
#       d2_prop[c,] <- rtruncnorm(ldim, a = rep(0, ldim), mean = d2[c,], sd = sd2)
#       
#       
#       
#       new_d <- compute_delta_eta_density_c(d1_prop, d2_prop, init_mcmc$eta, init_mcmc$beta,
#                                            init_mcmc$xi_eta, X, c - 1)
#       # print("")
#       old_d <- compute_delta_eta_density_c(d1, d2, init_mcmc$eta, init_mcmc$beta, init_mcmc$xi_eta, X, c - 1)
#       print(new_d)
#       print(old_d)
#       # mhr <- new_d - old_d - sum(log(c(d1_prop[,c], d2_prop[c,]))) + sum(log(c(d1[,c], d2[c,])))
#       p1 <- sum(pnorm(d1_prop[,c], mean = 0, sd = sd1, log.p = TRUE))
#       p2 <- sum(pnorm(d2_prop[c,], mean = 0, sd = sd2, log.p = TRUE))
#       p3 <- sum(pnorm(d1[,c], mean = 0, sd = sd1, log.p = TRUE))
#       p4 <- sum(pnorm(d2[c,], mean = 0, sd = sd2, log.p = TRUE))
#       print(p1)
#       print(p3)
#       print(p2)
#       print(p4)
#       mhr <- new_d - old_d + p3 + p4 - p1 - p1
#       print(mhr)
#       if (runif(1) < exp(mhr)) {
#         d1[,c] <- d1_prop[,c]
#         d2[c,] <- d2_prop[c,]
#         print("accepted")
#       }
#       d1_container[,,i] <- d1
#       d2_container[,,i] <- d2
#     }
#   }
# }
# 
# full_var_container <- array(0, dim = c(nreg, ldim, num_iter))
# for (i in 1:num_iter) {
#   q <- c(d1_container[,,i], t(d2_container[,,i]))
#   full_var_container[,,i] <- full_var(q)
# }
# r <- 3
# l <- 2
# 1 / full_var_container[,,num_iter]
# full_var_container[,,num_iter]
# plot(1 / full_var_container[r,l,], type = "l")
# abline(h = 1 / init_mcmc$prec_eta[r,l], col = "red")
# identify_delta <- function(delta) {
#   delta_cumprod <- cumprod(delta)
#   delta <- get_cumprod_coef(sort(delta_cumprod))
#   return(delta)
# }
# 
# q0 <- c(init_mcmc$delta_eta1, t(init_mcmc$delta_eta2))
# q1 <- c(d1_container[,,5000], t(d2_container[,,5000]))
# 1 / full_var(q0)
# 1 / full_var(q1)
# 
# 
# d1 <- init_mcmc$delta_eta1
# d2 <- init_mcmc$delta_eta2
# q0 <- c(d1, t(d2))
# full_var(q0)
# # d1 %*% t(d2)
# a_coef <- sqrt(d1[1,1] * d2[1,1])
# d1[,1] <- d1[,1] * a_coef^(1/7) / d1[1,1]
# d2[1,] <- d2[1,] * a_coef^(1/4) / d2[1,1]
# q0 <- c(d1, t(d2))
# full_var(q0)
# # d1 %*% t(d2)
# 
# q0 <- c(d1, t(d2))
# q1 <- c(d1_con, t(d2_identify))
# full_var(q0)
# full_var(q1)
# 
# 
# library("future")
# plan(multicore)
# 
# myfun <- function() {
#   future(fun2())
#   
#   return(1+1)
# }
# 
# ## Use multicore futures
# plan(multicore, workers = availableCores() - 1)
# # devtools::install_github("hathawayj/buildings")
# library(buildings) # remember that the 'permits' data object is created when the library is loaded.
# a <- 4
# ff <- function(x){
#   for (i in 1:1000){
#     i
#   }
#   
#   ggplot() + geom_point(x = permits[x, "value"])
# }
# list_object <- as.list(1:7500)
# tic()
# temp1 <- map(list_object, ff)
# toc()
# 
# 
# tic()
# temp1 <- future_map(list_object, ff)
# toc()
# 
# z_func <- function(seed) {
#   set.seed(seed)
#   print(rnorm(1))
# }
# z_func(2)
# library(future.apply)
# plan(multicore)
# future_lapply(1:4, z_func, future.seed = TRUE)
