num_cores <- 7
batch_size <- 7
myseq <- seq(from = 1, to = 300, by = num_cores)
for (i in myseq) {
  sysCall <- paste("Rscript sim_script_partial.R", num_cores, batch_size, i)
  system(sysCall)
}
