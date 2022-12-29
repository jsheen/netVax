# Script used to simulate two-stage trials -------------------------------------
set.seed(0)
N_trials = 1000
N_groups = 10
N_sims = 1000
cutoff = 90
assignment_mechanisms = c(0, 0.1, 0.25, 0.5, 0.75)
R0_vax = 0.25
final_ls <- list()
final_ls_dex <- 1
for (trial_num in 1:N_trials) {
  clusters_to_use <- sample(seq(0, N_sims - 1, 1), N_groups)
  if (N_groups %% length(assignment_mechanisms) != 0) {
    stop('The number of groups should be divisible by the number of assignment mechanisms.')
  }
  clusters_to_use_dex <- 1
  trial_dfs <- list()
  for (realized_assign in rep(assignment_mechanisms, N_groups / length(assignment_mechanisms))) {
    cluster_for_trial <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_mort0.85_eit0.005_vaxEff0.6_assign', realized_assign, '_sim', clusters_to_use[clusters_to_use_dex], '.csv'))
    trial_dfs[[clusters_to_use_dex]] <- cluster_for_trial
    clusters_to_use_dex <- clusters_to_use_dex + 1
  }
  # First, calculate percent infected within each cluster of different assignment
  ave_infect_non_enrolled_num <- rep(0, length(assignment_mechanisms))
  ave_infect_non_enrolled_denom <- rep(0, length(assignment_mechanisms))
  for (assignment_mech_dex in 1:length(assignment_mechanisms)) {
    for (trial_df_dex in seq(assignment_mech_dex, N_groups, length(assignment_mechanisms))) {
      cluster_df <- trial_dfs[[trial_df_dex]]
      if (cluster_df$node[1] != 'na') {
        denom_to_add <- length(which(cluster_df$assignment == 'na'))
        num_to_add <- length(which(cluster_df$assignment == 'na' & cluster_df$time2inf_trt < cutoff))
        ave_infect_non_enrolled_denom[assignment_mech_dex] <- denom_to_add
        ave_infect_non_enrolled_num[assignment_mech_dex] <- num_to_add
      }
    }
  }
  ave_infect_non_enrolled <- ifelse(ave_infect_non_enrolled_denom == 0, NA, ave_infect_non_enrolled_num / ave_infect_non_enrolled_denom)
  temp_df <- data.frame(matrix(ave_infect_non_enrolled[1:length(ave_infect_non_enrolled) - 1] - ave_infect_non_enrolled[2:length(ave_infect_non_enrolled)], nrow=1, ncol=length(ave_infect_non_enrolled) - 1))
  final_ls[[final_ls_dex]] <- temp_df
  final_ls_dex <- final_ls_dex + 1
}
final_df <- do.call(rbind, final_ls)

for (col_dex in 1:ncol(final_df)) {
  print(paste0('Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Mean: ', mean(final_df[,col_dex], na.rm=T)))
  #print(paste0('Variance: ', var(final_df[,col_dex], na.rm=T)))
}






