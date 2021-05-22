runs <- 1:100
for (run in runs) {
  sysCall <- paste("nohup Rscript sim_script_weak.R", run)
  system(sysCall, wait = FALSE)
}

# system("sleep 2", wait = TRUE)

runs <- 801:900
for (run in runs) {
  sysCall <- paste("nohup Rscript sim_script_weak.R", run)
  system(sysCall, wait = FALSE)
}
# 
# runs <- 1:16
# 
# 
# for (run in runs) {
#   if (run %% 4) {
#     system("")
#   }
#   sysCall <- paste("nohup Rscript runthis.R", run)
#   system(sysCall, wait = FALSE)
# }