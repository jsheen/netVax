




R0_vax = 0
assignment_mechanism = 0.1
cutoff_react = 120
cutoff_antici = 120
to_hist_antici <- c()
to_hist_react <- c()
cnt <- 0
for (sim_num in 0:1999) {
  react <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism, '_sim', sim_num, '_SEIR.csv'))
  antici <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism, '_sim', sim_num, '_SEIR_antici.csv'))
  to_hist_react <- c(to_hist_react, length(which(react$assignment == 'na' & react$time2inf_trt < cutoff_react & react$node[1] != 'na')))
  to_hist_antici <- c(to_hist_antici, length(which(antici$assignment == 'na' & antici$time2inf_trt < cutoff_antici & antici$node[1] != 'na')))
}
sampled_react <- c()
sampled_antici <- c()
for (j in 1:10) {
  for (i in 1:2000) {
    sampled_react <- c(sampled_react, length(which(sample(c(rep(1, to_hist_react[i]), rep(0, 1000 - to_hist_react[i])), size=100) == 1)))
    sampled_antici <- c(sampled_antici, length(which(sample(c(rep(1, to_hist_antici[i]), rep(0, 1000 - to_hist_antici[i])), size=100) == 1)))
  }
}
threshold = 2
length(which(sampled_react < threshold)) / 20000
length(which(sampled_antici < threshold)) / 20000
#hist(to_hist_react, col='blue')
#hist(to_hist_antici, col='red', add=TRUE)
mean(to_hist_react) / 1000
mean(to_hist_antici) / 1000




R0_vax = 0
assignment_mechanism_1 = 0
assignment_mechanism_2 = 0.1
cutoff = 120
to_hist_mech1 <- c()
to_hist_mech2 <- c()
cnt <- 0
for (sim_num in 0:999) {
  mech1 <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism_1, '_sim', sim_num, '_SEIR_antici.csv'))
  mech2 <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff0.8_assign', assignment_mechanism_2, '_sim', sim_num, '_SEIR_antici.csv'))
  to_hist_mech1 <- c(to_hist_mech1, length(which(mech1$assignment == 'na' & mech1$time2inf_trt < cutoff & mech1$node[1] != 'na')))
  to_hist_mech2 <- c(to_hist_mech2, length(which(mech2$assignment == 'na' & mech2$time2inf_trt < cutoff & mech2$node[1] != 'na')))
}
sampled_mech1 <- c()
sampled_mech2 <- c()
for (i in 1:1000) {
  sampled_mech1 <- c(sampled_mech1, length(which(sample(c(rep(1, to_hist_mech1[i]), rep(0, 1000 - to_hist_mech1[i])), size=100) == 1)))
  sampled_mech2 <- c(sampled_mech2, length(which(sample(c(rep(1, to_hist_mech2[i]), rep(0, 1000 - to_hist_mech2[i])), size=100) == 1)))
}
threshold = 2
length(which(sampled_mech1 < threshold))
length(which(sampled_mech2 < threshold))


