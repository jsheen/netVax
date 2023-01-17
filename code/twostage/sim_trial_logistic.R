# Script used to simulate two-stage trials -------------------------------------
# Libraries --------------------------------------------------------------------
library("stats")

# Parameters -------------------------------------------------------------------
set.seed(0)
N_sims = 1000 # Total number of cluster simulations in simulation bank
N_sample = 100 # Number sampled from each cluster
N_trials = 1000 # Number of trial simulations to conduct
cutoff = 90
assignment_mechanisms = c(0, 0.1, 0.25)
N_assignment_mechanism_sets = 1
N_groups = length(assignment_mechanisms) * N_assignment_mechanism_sets
R0_vax = 0.25

# Trial simulations ------------------------------------------------------------
final_ls <- list()
final_est_eff_ls <- list()
final_pval_ls <- list()
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
  sampled_ls <- list()
  for (assignment_mech_dex in 1:length(assignment_mechanisms)) {
    assignment_mech_ls <- list()
    assignment_mech_ls_dex <- 1
    for (trial_df_dex in seq(assignment_mech_dex, N_groups, length(assignment_mechanisms))) {
      cluster_df <- trial_dfs[[trial_df_dex]]
      if (cluster_df$node[1] != 'na') {
        # The following is to get the true effect
        denom_to_add <- length(which(cluster_df$assignment == 'na'))
        num_to_add <- length(which(cluster_df$assignment == 'na' & cluster_df$time2inf_trt < cutoff))
        ave_infect_non_enrolled_denom[assignment_mech_dex] <- denom_to_add
        ave_infect_non_enrolled_num[assignment_mech_dex] <- num_to_add
        # The following is prepare the dataframe for the cox-ph analysis
        mod_cluster_df <- cluster_df
        mod_cluster_df <- mod_cluster_df[which(mod_cluster_df$assignment == 'na'),]
        mod_cluster_df$assignment_mech <- assignment_mech_dex
        mod_cluster_df$clus_num <- paste0(assignment_mech_dex, '_', trial_df_dex)
        mod_cluster_df$status <- ifelse(mod_cluster_df$time2inf_trt < cutoff, 2, 1) # Where 2 is infected, 1 is censored
        mod_cluster_df <- mod_cluster_df[sample(1:nrow(mod_cluster_df), N_sample),] # Sample some number of this, assuming samp_num << number that are 'na'
        assignment_mech_ls[[assignment_mech_ls_dex]] <- mod_cluster_df
        assignment_mech_ls_dex <- assignment_mech_ls_dex + 1
      }
    }
    assignment_mech_df <- do.call(rbind, assignment_mech_ls) 
    if (is.null(assignment_mech_df)) {
      sampled_ls[[assignment_mech_dex]] <- data.frame(matrix(NA, nrow=1, ncol=1))
    } else {
      sampled_ls[[assignment_mech_dex]] <- assignment_mech_df
    }
  }
  # First, get information from sampled clusters
  est_eff_res <- c()
  pval_res <- c()
  for (sampled_dex in 1:(length(sampled_ls) - 1)) {
    low_df <- sampled_ls[[sampled_dex]]
    high_df <- sampled_ls[[sampled_dex + 1]]
    if (ncol(low_df) > 1 & ncol(high_df) > 1) {
      low_df$cond <- 1
      high_df$cond <- 2
      to_analyze <- rbind(low_df, high_df)
      res.cox <- coxph(Surv(time2inf_trt, status) ~ cond, data = to_analyze)
      est_eff_res <- c(est_eff_res, summary(res.cox)[7]$coefficients[2])
      pval_res <- c(pval_res, summary(res.cox)[7]$coefficients[5])
    } else {
      est_eff_res <- c(est_eff_res, NA)
      pval_res <- c(pval_res, NA)
    }
  }
  final_est_eff_ls[[final_ls_dex]] <- data.frame(matrix(est_eff_res, nrow=1, ncol=(length(assignment_mechanisms) - 1)))
  final_pval_ls[[final_ls_dex]] <- data.frame(matrix(pval_res, nrow=1, ncol=(length(assignment_mechanisms) - 1)))
  
  # Second, get true effect (estimand)
  ave_infect_non_enrolled <- ifelse(ave_infect_non_enrolled_denom == 0, NA, ave_infect_non_enrolled_num / ave_infect_non_enrolled_denom)
  temp_df <- data.frame(matrix(ave_infect_non_enrolled[1:length(ave_infect_non_enrolled) - 1] - ave_infect_non_enrolled[2:length(ave_infect_non_enrolled)], nrow=1, ncol=length(ave_infect_non_enrolled) - 1))
  final_ls[[final_ls_dex]] <- temp_df
  final_ls_dex <- final_ls_dex + 1
}
final_df <- do.call(rbind, final_ls)
final_est_eff_df <- do.call(rbind, final_est_eff_ls)
final_pval_df <- do.call(rbind, final_pval_ls)

# Q1) What is the the median and variance of the estimated effects under this trial design?
for (col_dex in 1:ncol(final_est_eff_df)) {
  print(paste0('Effect estimate: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Median: ', median(final_est_eff_df[,col_dex], na.rm=T)))
  print(paste0('Variance: ', var(final_est_eff_df[,col_dex], na.rm=T)))
}

# Q2) What is the power to detect an effect under this trial design using the cox proportional hazards model?
for (col_dex in 1:ncol(final_pval_df)) {
  print(paste0('Power: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(length(which(final_pval_df[,col_dex] < 0.05)) / length(which(!is.na(final_pval_df[,col_dex]))))
}

# Q3) What is the true effect across groups (mean and variance)?
for (col_dex in 1:ncol(final_df)) {
  print(paste0('True effect: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Mean: ', mean(final_df[,col_dex], na.rm=T)))
  print(paste0('Variance: ', var(final_df[,col_dex], na.rm=T)))
}

