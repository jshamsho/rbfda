#' @export
get_posteigenfunc <- function(result, num = 1, conf = .95) {
  ldim <- result$data$ldim
  tt <- result$data$time
  iterations <- result$control$iterations
  burnin <- result$control$burnin
  eigenmat <- matrix(0, length(tt), iterations - burnin)
  B <- result$data$basis
  for (iter in (burnin + 1):iterations) {
    index <- iter - burnin
    eigenmat[, index] <- B %*% result$samples$lambda[, num, iter]
  }
  mean_eigen <- apply(eigenmat, 1, mean)
  sd_eigen <- apply(eigenmat, 1, sd)
  mult <- -qnorm((1 - conf) / 2)
  summary_eigen <- matrix(0, length(tt), 3)
  summary_eigen[, 1] <- mean_eigen
  summary_eigen[, 2] <- mean_eigen - mult * sd_eigen
  summary_eigen[, 3] <- mean_eigen + mult * sd_eigen
  return(summary_eigen)
}

#' @export
get_posteigenvec <- function(result, num = 1, conf = .95) {
  nreg <- result$data$nreg
  iterations <- result$control$iterations
  burnin <- result$control$burnin
  eigenmat <- matrix(0, nreg, iterations - burnin)
  for (iter in (burnin + 1):iterations) {
    index <- iter - burnin
    eigenmat[, index] <- result$samples$phi[, num, iter]
  }
  mean_eigen <- apply(eigenmat, 1, mean)
  sd_eigen <- apply(eigenmat, 1, sd)
  mult <- -qnorm((1 - conf) / 2)
  summary_eigen <- matrix(0, nreg, 3)
  summary_eigen[, 1] <- mean_eigen
  summary_eigen[, 2] <- mean_eigen - mult * sd_eigen
  summary_eigen[, 3] <- mean_eigen + mult * sd_eigen
  return(summary_eigen)
}

#' @export
get_posteriormean <- function(result, region, x, conf = .95) {
  iterations <- result$control$iterations
  burnin <- result$control$burnin
  tt <- result$data$time
  nt <- length(tt)
  ldim <- result$data$ldim
  meanmat <- matrix(0, nt, iterations - burnin)
  B <- result$data$basis
  d <- length(x)
  for (iter in (burnin + 1):iterations) {
    index <- iter - burnin
    for (l in 1:ldim) {
      meanmat[, index] <- meanmat[, index] + 
        B %*% result$samples$lambda[, l, iter] %*%
        result$samples$phi[region,,iter] %*% 
        t(matrix(result$samples$beta[,l,iter], d)) %*% x
    }
  }
  mean_mean <- apply(meanmat, 1, mean)
  sd_mean <- apply(meanmat, 1, sd)
  mult <- -qnorm((1 - conf) / 2)
  summary_mean <- matrix(0, nt, 3)
  summary_mean[, 1] <- mean_mean
  summary_mean[, 2] <- mean_mean - mult * sd_mean
  summary_mean[, 3] <- mean_mean + mult * sd_mean
  return(summary_mean)
}
