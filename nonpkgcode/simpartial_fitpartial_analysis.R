# any(dir() == paste0("n50_simweak_fitweak", 1:2, ".RData"))
setwd("/Users/johnshamshoian/Documents/R_projects/rrbfda/output/n50_simpartial_fitpartial")
runthese <- which(sapply(1:300, function(i) any(dir() == paste0("n50_simpartial_fitpartial", i, ".RData"))))
nreg <- 6
ldim <- 4
nt <- 60
nsub <- 50
tt <- seq(from = 0, to = 1, length.out = 60)
nsim <- length(runthese)
true_eigenfunc <- fda::getbasismatrix(tt, fda::create.fourier.basis(nbasis = ldim))
mean_accuracy <- matrix(0, nreg, nsim)
mean_coverage <- matrix(0, nreg, nsim)
eigenfunc_accuracy <- matrix(0, ldim, nsim)
eigenfunc_coverage <- matrix(0, ldim, nsim)
eigenval_accuracy <- numeric(nsim)
eigenval_coverage <- numeric(nsim)
mean_pval <- numeric(nsim)
median_pval <- numeric(nsim)
for (num in runthese) {
  
  loadthis <- paste0("/Users/johnshamshoian/Documents/", 
                     "R_projects/rrbfda/output/",
                     "n50_simpartial_fitpartial/n50_simpartial_fitpartial", 
                     num, ".RData")
  load(loadthis)
  # true_eigenvec <- simstats$sim_data$phi
  for (r in 1:nreg) {
    mean_accuracy[r, num] <- pracma::trapz(tt, simstats$mean_summary[,r,1]^2)
    mean_coverage[r, num] <- mean(simstats$mean_summary[,r,2] < 0 & simstats$mean_summary[,r,3] > 0)
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
    eigenfunc_accuracy[l, num] <- min(prop_accuracy1, prop_accuracy2)
    eigenfunc_coverage[l, num] <- max(prop_coverage1, prop_coverage2)
  }
  
  
  eigenval_accuracy[num] <- sum((simstats$eigenval_summary[,,1] - simstats$sim_data$sigma_mat)^2)
  eigenval_coverage[num] <- mean((simstats$eigenval_summary[,,2] <= simstats$sim_data$sigma_mat) &
                                   (simstats$eigenval_summary[,,3] >= simstats$sim_data$sigma))

  mean_pval[num] <- mean(simstats$pvals)
  median_pval[num] <- median(simstats$pvals)
}

plot(mean_pval)
# plot(median_pval)
mean(mean_pval < .05)
mean(median_pval < .05)
hist(mean_pval, main = "Partial separability goondess of fit under null", xlab = "p-value")
