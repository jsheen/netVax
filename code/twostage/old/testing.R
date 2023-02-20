
set.seed(0)
N_sims = 1000 # Total number of cluster simulations in simulation bank
N_trials = 1000 # Number of trial simulations to conduct
cutoff = 90
R0_vax = 0.25

realized_assign = 0.25
N_cluster = 5
sums = c()
tot_count_burnout = 0
for (trial_num in 1:N_trials) {
  if (trial_num %% 200 == 0) {
    print(trial_num)
  }
  sum_clus = 0
  count_burnout = 0
  for (clus in 1:N_cluster) {
    cluster_df_indic <- NA
    cluster_df <- NA
    while (is.na(cluster_df_indic)) {
      test <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_mort0.85_eit0.005_vaxEff0.6_assign', realized_assign, '_sim', sample(0:(N_sims - 1), 1), '.csv'))
      if (test$node[1] != 'na') {
        cluster_df <- test
        cluster_df_indic = 1
      }
    }
    denom_to_add <- length(which(cluster_df$assignment == 'na'))
    num_to_add <- length(which(cluster_df$assignment == 'na' & cluster_df$time2inf_trt < cutoff))
    cluster_average <- num_to_add / denom_to_add
    if (cluster_average < 1 / 1000) {
      print(cluster_average)
      count_burnout = count_burnout + 1
    }
    sum_clus = sum_clus + cluster_average
  }
  sums = c(sums, sum_clus / N_cluster)
  tot_count_burnout = tot_count_burnout + count_burnout
}
mean(sums)
hist(sums)
tot_count_burnout / N_trials




