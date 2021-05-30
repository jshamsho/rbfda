library(tidyverse)
runthese <- paste0("n200_simnonpartial_fitpartial", 1:200, ".RData") %in% 
  dir(paste0("/Users/johnshamshoian/Documents/",
             "R_projects/rrbfda/output/",
             "n200_simnonpartial_fitpartial")) %>%
  which()

nsim <- length(runthese)
nreg <- 6
ldim <- 4
nt <- 60
nsub <- 200
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
mean_pval <- numeric(nsim)
median_pval <- numeric(nsim)
i <- 1
for (num in runthese) {
  
  loadthis <- paste0("/Users/johnshamshoian/Documents/", 
                     "R_projects/rrbfda/output/",
                     "n200_simnonpartial_fitpartial/n200_simnonpartial_fitpartial", 
                     num, ".RData")
  load(loadthis)
  true_eigenvec <- simstats$sim_data$phi
  for (r in 1:nreg) {
    mean_accuracy[r, i] <- pracma::trapz(tt, simstats$mean_summary[,r,1]^2)
    mean_coverage[r, i] <- mean(simstats$mean_summary[,r,2] < 0 & simstats$mean_summary[,r,3] > 0)
  }
  for (l in 1:ldim) {
    simstats$eigenfunc_summary[,l,] <- simstats$eigenfunc_summary[,l,] * sqrt(nt)
    prop_accuracy1 <- 
      pracma::trapz(tt, (simstats$eigenfunc_summary[,l,1] - true_eigenfunc[,l])^2)
    prop_accuracy2 <- 
      pracma::trapz(tt, (simstats$eigenfunc_summary[,l,1] + true_eigenfunc[,l])^2)
    prop_coverage1 <- mean((simstats$eigenfunc_summary[,l,2] <= true_eigenfunc[,l]) &
                             (simstats$eigenfunc_summary[,l,3] >= true_eigenfunc[,l]))
    prop_coverage2 <- mean((simstats$eigenfunc_summary[,l,2] <= -true_eigenfunc[,l]) &
                             (simstats$eigenfunc_summary[,l,3] >= -true_eigenfunc[,l]))
    eigenfunc_accuracy[l, i] <- min(prop_accuracy1, prop_accuracy2)
    eigenfunc_coverage[l, i] <- max(prop_coverage1, prop_coverage2)
  }
  
# 
  for (r in 1:nreg) {
    # prop_accuracy1 <- sum((simstats$eigenvec_summary[,r,1,1] - true_eigenvec[,r,1])^2)
    # prop_accuracy2 <- sum((simstats$eigenvec_summary[,r,1,1] + true_eigenvec[,r,1])^2)
    # prop_coverage1 <- mean((simstats$eigenvec_summary[,r,2] <= true_eigenvec[,r,1]) &
                             # (simstats$eigenvec_summary[,r,3] >= true_eigenvec[,r]))
    # prop_coverage2 <- mean((simstats$eigenvec_summary[,r,2] <= -true_eigenvec[,r]) &
                             # (simstats$eigenvec_summary[,r,3] >= -true_eigenvec[,r]))
    # eigenvec_accuracy[r, i] <- min(prop_accuracy1, prop_accuracy2)
    # eigenvec_coverage[r, i] <- max(prop_coverage1, prop_coverage2)
  }
#   
  # eigenval_accuracy[i] <- sum((simstats$eigenval_summary[,,1] - simstats$sim_data$sigma_mat)^2)
  # eigenval_coverage[i] <- mean((simstats$eigenval_summary[,,2] <= simstats$sim_data$sigma_mat) &
  #                                  (simstats$eigenval_summary[,,3] >= simstats$sim_data$sigma))
  # 
  mean_pval[i] <- mean(simstats$pvals)
  median_pval[i] <- median(simstats$pvals)
  i <- i + 1
}
round(sum(mean_pval < .05) / nsim * 100, 2)
round(sum(mean_pval < .10) / nsim * 100, 2)
round(mean(apply(mean_accuracy, 1, median)) * 1e3, digits = 2)
# apply(eigenfunc_accuracy, 1, quantile, c(.025, .5, .975))
round(apply(eigenfunc_accuracy, 1, median) * 1e3, digits = 2)
round(apply(eigenvec_accuracy, 1, median) * 1e3, digits = 2)
