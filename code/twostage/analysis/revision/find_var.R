


# Delta 1
assig_frac <- 0.05
over <- F
vaxEff <- 0.8
R0_vax <- 0.9
assign_type <- 'trad'
fs <- c()
for (sim_num in 0:(N_sims - 1)) {
  if (over) {
    test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff', vaxEff, '_assign', assign_frac, '_assigntype', assign_type, '_sim', sim_num, '_SEIR_antici_revision.csv'))
  } else {
    test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_R0wt2_R0vax', R0_vax, '_eit0.005_vaxEff', vaxEff, '_assign', assign_frac, '_assigntype', assign_type, '_sim', sim_num, '_SEIR_Pois_antici_revision.csv'))
  }
  if (test_cluster$node[1] != 'na') {
    fs <- c(fs, length(which(test_cluster$time2inf_trt < 150)))
  }
}
var(fs/1000)
