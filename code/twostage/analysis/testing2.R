R0_vax = 0.5
assignment_mechanism = 0
cutoff_react = 120
cutoff_antici = 120
to_hist_antici <- c()
to_hist_react <- c()
cnt <- 0
for (sim_num in 0:999) {
  react <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism, '_sim', sim_num, '_SEIR.csv'))
  antici <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism, '_sim', sim_num, '_SEIR_antici.csv'))
  to_hist_react <- c(to_hist_react, length(which(react$assignment == 'na' & react$time2inf_trt < cutoff_react & react$node[1] != 'na')))
  to_hist_antici <- c(to_hist_antici, length(which(antici$assignment == 'na' & antici$time2inf_trt < cutoff_antici & antici$node[1] != 'na')))
}
hist(to_hist_react, col='blue')
hist(to_hist_antici, col='red', add=TRUE)
mean(to_hist_react)
mean(to_hist_antici)

# for (sim_num in 0:999) {
#   react <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism, '_sim', sim_num, '_SEIR.csv'))
#   if (length(which(react$assignment == 'na' & react$time2inf_trt < cutoff_react)) == 0) {
#     break
#   }
# }
