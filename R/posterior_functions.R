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
get_posteigenvec <- function(result, conf = .95) {
  nreg <- result$data$nreg
  iterations <- result$control$iterations
  burnin <- result$control$burnin
  ldim <- result$data$ldim
  summary_eigen <- NULL
  if (result$control$covstruct == "weak") {
    summary_eigen <- array(0, dim = c(nreg, nreg, 3))
    for (r in 1:nreg) {
      eigenmat <- matrix(0, nreg, iterations - burnin)
      for (iter in (burnin + 1):iterations) {
        index <- iter - burnin
        eigenmat[, index] <- result$samples$phi[, r, iter]
      }
      mean_eigen <- apply(eigenmat, 1, mean)
      sd_eigen <- apply(eigenmat, 1, sd)
      mult <- -qnorm((1 - conf) / 2)
      summary_eigen[, r, 1] <- mean_eigen
      summary_eigen[, r, 2] <- mean_eigen - mult * sd_eigen
      summary_eigen[, r, 3] <- mean_eigen + mult * sd_eigen
    }
  } else if (result$control$covstruct == "partial") {
    summary_eigen <- array(0, dim = c(nreg, nreg, ldim, 3)) 
    for (l in 1:ldim) {
      for (r in 1:nreg) {
        eigenmat <- matrix(0, nreg, iterations - burnin)
        for (iter in (burnin + 1):iterations) {
          index <- iter - burnin
          eigenmat[, index] <- result$samples$phi[[iter]][, r, l]
        }
        mean_eigen <- apply(eigenmat, 1, mean)
        sd_eigen <- apply(eigenmat, 1, sd)
        mult <- -qnorm((1 - conf) / 2)
        summary_eigen[, r, l, 1] <- mean_eigen
        summary_eigen[, r, l, 2] <- mean_eigen - mult * sd_eigen
        summary_eigen[, r, l, 3] <- mean_eigen + mult * sd_eigen
      }
    }
  }
  return(summary_eigen)
}

#' @export
get_postmean <- function(result, x, conf = .95) {
  iterations <- result$control$iterations
  burnin <- result$control$burnin
  tt <- result$data$time
  nt <- length(tt)
  ldim <- result$data$ldim
  nreg <- result$data$nreg
  B <- result$data$basis
  d <- length(x)
  summary_mean <- array(0, dim = c(nt, nreg, 3))
  for (j in 1:nreg) {
    meanmat <- matrix(0, nt, iterations - burnin)
    for (iter in (burnin + 1):iterations) {
      index <- iter - burnin
      for (l in 1:ldim) {
        if (result$control$covstruct == "weak") {
          meanmat[, index] <- meanmat[, index] + 
            B %*% result$samples$lambda[, l, iter] %*%
            result$samples$phi[j,,iter] %*% 
            t(matrix(result$samples$beta[,l,iter], d)) %*% x
        } else if (result$control$covstruct == "partial") {
          meanmat[, index] <- meanmat[, index] + 
            B %*% result$samples$lambda[, l, iter] %*%
            result$samples$phi[[iter]][j,,l] %*% 
            t(matrix(result$samples$beta[,l,iter], d)) %*% x
        }
      }
    }
    mean_mean <- apply(meanmat, 1, mean)
    sd_mean <- apply(meanmat, 1, sd)
    mult <- -qnorm((1 - conf) / 2)
    summary_mean[, j, 1] <- mean_mean
    summary_mean[, j, 2] <- mean_mean - mult * sd_mean
    summary_mean[, j, 3] <- mean_mean + mult * sd_mean
  }
  return(summary_mean)
}
