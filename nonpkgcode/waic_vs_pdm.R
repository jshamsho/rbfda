library(tidyverse)
runthese <- paste0("n50_simweak_waic", 1:300, ".RData") %in% 
  dir(paste0("/Users/johnshamshoian/Documents/",
             "R_projects/rrbfda/output/",
             "n50_simweak_waic")) %>%
  which()

nsim <- length(runthese)
nreg <- 6
ldim <- 4
nt <- 60
nsub <- 200
M <- 200
iterations <- 7500
q_M <- seq(from = 0.005, to = .995, length.out = M)
order_pdm <- floor(quantile(1:iterations, q_M))
tt <- seq(from = 0, to = 1, length.out = 60)
true_eigenfunc <- fda::getbasismatrix(tt, fda::create.fourier.basis(nbasis = ldim))
mean_accuracy <- matrix(0, nreg, nsim)
mean_coverage <- matrix(0, nreg, nsim)
eigenfunc_accuracy <- matrix(0, ldim, nsim)
eigenfunc_coverage <- matrix(0, ldim, nsim)
eigenvec_accuracy <- matrix(0, nreg, nsim)
eigenvec_coverage <- matrix(0, nreg, nsim)
eigenval_accuracy <- numeric(nsim)
eigenval_coverage <- numeric(nsim)
mean_zscore <- numeric(nsim)
median_zscore <- numeric(nsim)
mean_pval <- numeric(nsim)
median_pval <- numeric(nsim)
decision_pdm <- numeric(nsim)
decision_waic <- numeric(nsim)
decision_abs_waic <- numeric(nsim)
yjpmin <- numeric(nsim)
i <- 1
for (num in runthese) {
  
  loadthis <- paste0("/Users/johnshamshoian/Documents/", 
                     "R_projects/rrbfda/output/",
                     "n50_simweak_waic/n50_simweak_waic", 
                     num, ".RData")
  load(loadthis)
  pdm_partial <- qnorm(simstats$pvals_partial)
  pdm_weak <- qnorm(simstats$pvals_weak)
  se_pdm <- sqrt(var(pdm_partial) + var(pdm_weak))
  if (abs(mean(pdm_partial) - mean(pdm_weak)) > 2 * se_pdm) {
    if (mean(pdm_partial) > mean(pdm_weak)) {
      decision_pdm[i] <- 1
    } else {
      decision_pdm[i] <- 2
    }
  } else {
    decision_pdm[i] <- 3
  }
  if (abs(simstats$waic_diff$diff) > 2 * simstats$waic_diff$diff_se) {
    if (simstats$waic_diff$diff > 0) {
      decision_waic[i] <- 1
    } else {
      decision_waic[i] <- 2
    }
  } else {
    decision_waic[i] <- 3
  }
  if (simstats$waic_diff$diff > 0) {
    decision_abs_waic[i] <- 1
  }
  i <- i + 1
}

round(table(decision_waic) / nsim * 100, 2)
round(table(decision_pdm) / nsim * 100, 2)
round(sum(decision_abs_waic) / nsim * 100, 2)