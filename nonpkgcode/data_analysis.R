library(tidyverse)
library(rrbfda)
library(mgcv)
load("/Users/johnshamshoian/Documents/R_projects/bfcr/ASD/pa.dat.RData")
demographics <- read.csv("/Users/johnshamshoian/Documents/R_projects/bfcr/ASD/demographic_data.csv")
source("/Users/johnshamshoian/Documents/R_projects/rrbfda/nonpkgcode/data_analysis_helper.R")
tt <- seq(from = 6, to = 14, by = .25)
chan_id <- c('Fp1', 'Fp2','F9','F7','F3','Fz','F4','F8','F10','T9','T7',
             'C3','Cz','C4','T8','T10','P9','P7','P3','Pz','P4','P8','P10','O1','O2')
left_temporal <- 10:11
right_temporal <- 15:16
frontal <- 1:9
central <- 12:14
occiptal_parietal <- 17:25
pa.dat.wide <- pa.dat %>% 
  mutate(logy = log(y)) %>%
  select(-y) %>% 
  spread(reg, logy)

colnames(pa.dat.wide) <- c("ID", "group", "func", "Age", chan_id)
pa.dat.wide$left_temporal <- rowMeans(pa.dat.wide[, chan_id[left_temporal]])
pa.dat.wide$right_temporal <- rowMeans(pa.dat.wide[, chan_id[right_temporal]])
pa.dat.wide$frontal <- rowMeans(pa.dat.wide[, chan_id[frontal]])
pa.dat.wide$central <- rowMeans(pa.dat.wide[, chan_id[central]])
pa.dat.wide$occiptal_parietal <- rowMeans(pa.dat.wide[, chan_id[occiptal_parietal]])
pa.dat.wide <- pa.dat.wide %>% 
  select(ID, group, left_temporal, right_temporal, frontal, central, occiptal_parietal)
pa.dat.asd <- pa.dat.wide %>% 
  filter(group == 2)
pa.dat.asd.y <- pa.dat.asd %>%
  select(-ID, -group)
pa.dat.td <- pa.dat.wide %>% 
  filter(group == 1)
pa.dat.td.y <- pa.dat.td %>%
  select(-ID, -group)

nsub.asd <- pa.dat.asd$ID %>% unique() %>% length()
indices.asd <- demographics$ParticipantNumber %in% pa.dat.asd$ID %>% which()
nsub.td <- pa.dat.td$ID %>% unique() %>% length()
indices.td <- demographics$ParticipantNumber %in% pa.dat.td$ID %>% which()

nreg <- 5
tt <- seq(from = 6, to = 14, by = .25)
ndf <- 12
X.asd <- cbind(rep(1, nsub.asd), demographics[indices.asd,])
X.asd <- cbind(X.asd$`rep(1, nsub.asd)`, X.asd$Age..Months.)
X.asd[,2] <- (X.asd[,2] - mean(X.asd[,2])) / sd(X.asd[,2])
X.td <- cbind(rep(1, nsub.td), demographics[indices.td,])
X.td <- cbind(X.td$`rep(1, nsub.td)`, X.td$Age..Months.)
X.td[,2] <- (X.td[,2] - mean(X.td[,2])) / sd(X.td[,2])

basisobj <- mgcv::smoothCon(s(tt, k = ndf, bs = "ps", m = 2),
                            data.frame(tt), absorb.cons = FALSE)
B <- basisobj[[1]]$X
penalty <- basisobj[[1]]$S[[1]] * basisobj[[1]]$S.scale
init_mcmc.asd <- rrbfda::initialize_mcmc_weak(as.matrix(pa.dat.asd.y), tt, B, X.asd, pve = .95)
ldim.asd <- init_mcmc.asd$lambda %>% ncol()
init_mcmc.td <- rrbfda::initialize_mcmc_weak(as.matrix(pa.dat.td.y), tt, B, X.td, pve = .975)
ldim.td <- init_mcmc.td$lambda %>% ncol()
iter <- 5000
thin <- 1
burnin <- 1000
result.asd <- run_mcmc(response = as.matrix(pa.dat.asd.y), design = X.asd, basis = B, time = tt, 
                   penalty = penalty, ldim = ldim.asd, iter = iter, burnin = burnin,
                   thin = thin, init_ = init_mcmc.asd, covstruct = "weak")
result.td <- run_mcmc(response = as.matrix(pa.dat.td.y), design = X.td, basis = B, time = tt, 
                      penalty = penalty, ldim = ldim.td, iter = iter, burnin = burnin,
                      thin = thin, init_ = init_mcmc.td, covstruct = "weak")
predictions_tibble <- tibble(response = numeric(4 * nt),
                             mean = numeric(4 * nt),
                             lower = numeric(4 * nt),
                             upper = numeric(4 * nt),
                             freq = rep(tt, 4),
                             subj = rep(1:4, each = nt))

sub <- 12
reg <- 2
predictions <- get_postprediction(result.asd, sub, reg)
predictions_tibble[1:nt, 1] <- pa.dat.asd.y[((sub - 1) * nt + 1):(sub * nt), reg]
predictions_tibble[1:nt, 2] <- predictions[,1]
predictions_tibble[1:nt, 3] <- predictions[,2]
predictions_tibble[1:nt, 4] <- predictions[,3]
sub <- 24
reg <- 4
predictions <- get_postprediction(result.asd, sub, reg)
predictions_tibble[(nt + 1):(2 * nt), 1] <- pa.dat.asd.y[((sub - 1) * nt + 1):(sub * nt), reg]
predictions_tibble[(nt + 1):(2 * nt), 2] <- predictions[,1]
predictions_tibble[(nt + 1):(2 * nt), 3] <- predictions[,2]
predictions_tibble[(nt + 1):(2 * nt), 4] <- predictions[,3]
sub <- 15
reg <- 3
predictions <- get_postprediction(result.td, sub, reg)
predictions_tibble[(2 * nt + 1):(3 * nt), 1] <- pa.dat.td.y[((sub - 1) * nt + 1):(sub * nt), reg]
predictions_tibble[(2 * nt + 1):(3 * nt), 2] <- predictions[,1]
predictions_tibble[(2 * nt + 1):(3 * nt), 3] <- predictions[,2]
predictions_tibble[(2 * nt + 1):(3 * nt), 4] <- predictions[,3]
sub <- 30
reg <- 5
predictions <- get_postprediction(result.td, sub, reg)
predictions_tibble[(3 * nt + 1):(4 * nt), 1] <- pa.dat.td.y[((sub - 1) * nt + 1):(sub * nt), reg]
predictions_tibble[(3 * nt + 1):(4 * nt), 2] <- predictions[,1]
predictions_tibble[(3 * nt + 1):(4 * nt), 3] <- predictions[,2]
predictions_tibble[(3 * nt + 1):(4 * nt), 4] <- predictions[,3]
title.labs <- c("ASD subject 12, region 2", "ASD subject 24, region 4",
                "TD subject 15, region 3", "TD subject 30, region 5")
names(title.labs) <- c("1", "2", "3", "4")
p1 <- predictions_tibble %>%
  ggplot(aes(freq, response)) +
  geom_point(color = "grey70") +
  geom_line(aes(freq, lower)) +
  geom_line(aes(freq, upper)) +
  facet_wrap(~subj, labeller = labeller(subj = title.labs), nrow = 2) +
  labs(x = "Frequency (Hz)", y = "Negative log spectral density") +
  ylim(-4, -.5) + 
  theme_Publication(base_family = "Arial")
outfile <- "/Users/johnshamshoian/Documents/R_projects/rrbfda/figures/sample_paths.eps"
ggsave(outfile, p1, scale = .7, width = 8, height = 7, units = "in")

mean_tibble <- tibble(mean = numeric(10 * nt),
                      lower = numeric(10 * nt),
                      upper = numeric(10 * nt),
                      freq = rep(tt, 10),
                      region = rep(rep(1:5, each = nt), 2),
                      group = rep(1:2, each = 5 * nt),
                      baseline = rep(0, 10 * nt))
means.asd <- get_postmean(result.asd, c(0, 1))
means.td <- get_postmean(result.td, c(0, 1))
for (r in 1:5) {
  mean_tibble$mean[((r - 1) * nt + 1):(r * nt)] <- means.asd[,r,1]
  mean_tibble$lower[((r - 1) * nt + 1):(r * nt)] <- means.asd[,r,2]
  mean_tibble$upper[((r - 1) * nt + 1):(r * nt)] <- means.asd[,r,3]
}
add_this <- 5 * nt
for (r in 1:5) {
  mean_tibble$mean[add_this + ((r - 1) * nt + 1):(r * nt)] <- means.td[,r,1]
  mean_tibble$lower[add_this + ((r - 1) * nt + 1):(r * nt)] <- means.td[,r,2]
  mean_tibble$upper[add_this + ((r - 1) * nt + 1):(r * nt)] <- means.td[,r,3]
}
group.labs <- c("ASD", "TD")
names(group.labs) <- c(1,2)
region.labs <- c("Left temporal", "Right temporal", "Frontal", "Central", "Occiptal-Pareital")
names(region.labs) <- 1:5
p2 <- mean_tibble %>% ggplot(aes(freq, mean)) +
  geom_line() +
  geom_line(aes(freq, lower), linetype = 2) +
  geom_line(aes(freq, upper), linetype = 2) +
  geom_line(aes(freq, baseline), linetype = 3) + 
  facet_grid(group ~ region, labeller = labeller(region = region.labs, group = group.labs)) +
  labs(y = "Negative log spectral density", x = "Frequency") +
  theme_Publication(base_family = "Arial")
outfile <- "/Users/johnshamshoian/Documents/R_projects/rrbfda/figures/means.eps"
ggsave(outfile, p2, width = 10, height = 4, units = "in")

eigen_tibble <- tibble(mean = numeric(6 * nt),
                       freq = rep(tt, 6),
                       group = rep(1:2, each = 3 * nt),
                       number = rep(rep(1:3, each = nt),2),
                       baseline = rep(0, 6 * nt))
add_this <- 3 * nt
for (l in 1:3) {
  eigenf.asd <- get_posteigenfunc(result.asd, l)
  eigen_tibble$mean[((l - 1) * nt + 1):(l * nt)] <- eigenf.asd[,1]
  eigenf.td <- get_posteigenfunc(result.td, l)
  tmp_eigen <- eigenf.td[,1]
  if (sum((tmp_eigen + eigenf.asd[,1])^2) < 
      sum((tmp_eigen - eigenf.asd[,1])^2)) {
    tmp_eigen <- -tmp_eigen
  }
  eigen_tibble$mean[add_this + ((l - 1) * nt + 1):(l * nt)] <- tmp_eigen
}
number.labs <- c("Eigenfunction 1", "Eigenfunction 2", "Eigenfunction 3")
names(number.labs) <- c(1, 2, 3)
p3 <- eigen_tibble %>%
  ggplot(aes(freq, mean)) +
  geom_line() +
  geom_line(aes(freq, baseline), lty = 3) + 
  ylim(-.4,.4) +
  facet_grid(group ~ number, labeller = labeller(group = group.labs, number = number.labs)) + 
  labs(y = "Negative log spectral power", x = "Frequency") +
  theme_Publication(base_family = "Arial")
outfile <- "/Users/johnshamshoian/Documents/R_projects/rrbfda/figures/eigen.eps"
ggsave(outfile, p3, width = 6, height = 4, units = "in", scale = 1.25)

eigenval_summary <- function(result, num) {
  burnin <- result$control$burnin
  iter <- result$control$iterations
  quantile(sapply((burnin + 1):iter, function(i) {
                  sum(1 / result$samples$sigmasqeta[,num, i]) / 
                    sum(1 / result$samples$sigmasqeta[,,i])}),
           c(.025, .975)) * 100
}

eigenval_summary(result.asd, 1)
eigenval_summary(result.asd, 2)
eigenval_summary(result.asd, 3)
eigenval_summary(result.td, 1)
eigenval_summary(result.td, 2)
eigenval_summary(result.td, 3)

pvals.asd <- get_pvals_weak(result.asd)
pnorm(mean(qnorm(pvals.asd)))
pvals.td <- get_pvals_weak(result.td)
pnorm(mean(qnorm(pvals.td)))


quantile(sapply((burnin + 1):iter, function(i) {
  sum(1 / result.asd$samples$sigmasqeta[,1, i]) / 
    sum(1 / result.asd$samples$sigmasqeta[,,i])
}), c(.025, .975)) * 100
quantile(sapply((burnin + 1):iter, function(i) {
  sum(1 / result.asd$samples$sigmasqeta[,2, i]) / 
    sum(1 / result.asd$samples$sigmasqeta[,,i])
}), c(.025, .975)) * 100
pvals <- get_pvals_partial(result)
hist(pvals)
pnorm(mean(qnorm(pvals)))
M <- 200
iterations <- 4000
q_M <- seq(from = 0.005, to = .995, length.out = M)
order_pdm <- floor(quantile(1:iterations, q_M))
pvals_sorted <- sort(pvals, decreasing = TRUE)
q_pval <- pvals_sorted[order_pdm]
min(1, iterations * (q_pval) / (iterations - order_pdm + 1))

eta_reshaped <- reshape_nreg(result$samples$eta[,,4000], nsub.asd, 5)
cor(eta_reshaped)
image(cor(eta_reshaped))

eta_reshaped <- reshape_nreg(init_mcmc$eta, nsub.asd, 5)
cor(eta_reshaped)
image(cor(eta_reshaped))

