
# What percentage of simulations have less than the threshold number of infections
threshold <- 3
less_threshold_react <- 0
less_threshold_antici <- 0
n_infect_react <- c()
n_infect_antici <- c()
R0_vax = 0.5
assignment_mechanism = 0.1
N_sims = 2000
for (sim_num in 0:(N_sims - 1)) {
  test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', assignment_mechanism, '_sim', sim_num, '_SEIR.csv'))
  if (length(which(test_cluster$time2inf_trt < 1000)) < threshold) {
    less_threshold_react <- less_threshold_react + 1
  }
  n_infect_react <- c(n_infect_react, length(which(test_cluster$time2inf_trt < 1000)))
  test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', assignment_mechanism, '_sim', sim_num, '_SEIR_antici.csv'))
  if (length(which(test_cluster$time2inf_trt < 1000)) < threshold) {
    less_threshold_antici <- less_threshold_antici + 1
  }
  n_infect_antici <- c(n_infect_antici, length(which(test_cluster$time2inf_trt < 1000)))
}
less_threshold_react / N_sims
less_threshold_antici / N_sims
mean(n_infect_react)
mean(n_infect_antici)
sd(n_infect_react)
sd(n_infect_antici)
length(which(n_infect_antici == 4))


# Write this test to try to see what's going on
num_groups_arm = 2
res_react <- c()
res_antici <- c()
for (i in 1:1000) {
  # React
  r_cont <- c()
  r_treat <- c()
  for (j in 1:num_groups_arm) {
    cont_react <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', 0, '_sim', sample(0:1999, 1), '_SEIR.csv'))
    r_cont <- c(r_cont, length(which(cont_react$time2inf < 1000)))
    treat_react <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', 0.1, '_sim', sample(0:1999, 1), '_SEIR.csv'))
    r_treat <- c(r_treat, length(which(treat_react$time2inf < 1000)))
  }
  # Antici
  a_cont <- c()
  a_treat <- c()
  for (j in 1:num_groups_arm) {
    cont_antici <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', 0, '_sim', sample(0:1999, 1), '_SEIR_antici.csv'))
    a_cont <- c(a_cont, length(which(cont_antici$time2inf < 1000)))
    treat_antici <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', 0.1, '_sim', sample(0:1999, 1), '_SEIR_antici.csv'))
    a_treat <- c(a_treat, length(which(treat_antici$time2inf < 1000)))
  }
  # Permutation test to get power
  n_perm = 1000
  outcomes <- to_analyze$status
  N_clus_c_perm <- as.integer(length(which(to_analyze$cond == 0)) / 100)
  N_clus_t_perm <- as.integer(length(which(to_analyze$cond == 1)) / 100)
  perms_hist <- c()
  for (perm_dex in 1:n_perm) {
    temp_to_analyze <- to_analyze
    clus_assigns <- sample(c(rep(0, N_clus_c_perm), rep(1, N_clus_t_perm)), size = N_clus_c_perm + N_clus_t_perm, replace=F)
    count_clus <- 1
    for (f_clus_num in unique(temp_to_analyze$clus_num)) {
      clus_dex <- which(temp_to_analyze$clus_num == f_clus_num)
      temp_to_analyze$cond[clus_dex] <- rep(clus_assigns[count_clus], length(clus_dex))
      count_clus <- count_clus + 1
    }
    perms_hist <- c(perms_hist, calc_est(temp_to_analyze))
  }
  pval <- 1 - (length(which(perms_hist < est)) / length(perms_hist))
  pval_res <- c(pval_res, pval)
  
}







