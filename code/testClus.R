# Parameters -------------------------------------------------------------------
right_censor <- 120
N_clusters <- c(1000)
overdispersions <- c(1)
R0_wts <- c(3)
R0_vaxs <- c(0)
morts <- c(0.85)
vaxEffs <- c(0.6)
param_sets <- list()
param_sets_dex <- 1
for (i in N_clusters) {
  for (j in overdispersions) {
    for (k in R0_wts) {
      for (l in R0_vaxs) {
        for (m in morts) {
          for (n in vaxEffs) {
            param_sets[[param_sets_dex]] <- c(i, j, k, l, m, n)
            param_sets_dex <- param_sets_dex + 1
          }
        }
      }
    }
  }
}
for (param_set in param_sets) {
  print(param_set)
  N_cluster <- param_set[1]
  overdispersion <- param_set[2]
  R0_wt <- param_set[3]
  R0_vax <- param_set[4]
  mort <- param_set[5]
  vax_eff <- param_set[6]
  eit <- 0.005
  count_bad <- 0
  for (sim_num in c(0:2999)) {
    filename <- paste0('~/netVax_sim_results/CRT_N', N_cluster, '_k', overdispersion, 
                       '_R0wt', R0_wt, '_R0vax', R0_vax, '_mort', mort, '_eit', eit,
                       '_vaxEff', vax_eff, '_sim', sim_num, '.csv')
    single_clus <- read.csv(filename)
    if (length(unique(single_clus$time2inf_trt[which(single_clus$assignment != 'na')])) == 3 |
        length(unique(single_clus$time2inf_con[which(single_clus$assignment != 'na')])) == 3) {
      count_bad <- count_bad + 1
    }
  }
  print(count_bad / 3000)
}

