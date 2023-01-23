# Script used to simulate two-stage trials -------------------------------------
# Libraries --------------------------------------------------------------------
library("stats")
library("lme4")

# Parameters -------------------------------------------------------------------
set.seed(0)
N_sims = 1000 # Total number of cluster simulations in simulation bank
N_sample = 200 # Number sampled from each cluster
N_trials = 200 # Number of trial simulations to conduct
cutoff = 90
assignment_mechanisms = c(0, 0.1, 0.25)
N_assignment_mechanism_sets = 5
N_groups = length(assignment_mechanisms) * N_assignment_mechanism_sets
R0_vax = 0.25
if (N_groups %% length(assignment_mechanisms) != 0) {
  stop('The number of groups should be divisible by the number of assignment mechanisms.')
}

# Get simulations to use for each assignment mechanism -------------------------
to_use_ls <- list()
to_use_ls_dex <- 1
for (assignment_mechanism in assignment_mechanisms) {
  to_use <- c()
  for (sim_num in 0:(N_sims - 1)) {
    test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_mort0.85_eit0.005_vaxEff0.6_assign', assignment_mechanism, '_sim', sim_num, '.csv'))
    if (test_cluster$node[1] != 'na') {
      to_use <- c(to_use, sim_num)
    }
  }
  to_use_ls[[to_use_ls_dex]] <- to_use
  to_use_ls_dex <- to_use_ls_dex + 1
}

# Trial simulations ------------------------------------------------------------
final_ls <- list()
final_est_eff_ls <- list()
final_pval_ls <- list()
final_ls_dex <- 1
for (trial_num in 1:N_trials) {
  print(trial_num)
  clusters_to_use <- c()
  for (assignment_mechanism_dex in 1:length(assignment_mechanisms)) {
    clusters_to_use <- c(clusters_to_use, sample(to_use_ls[[assignment_mechanism_dex]], N_assignment_mechanism_sets))
  }
  clusters_to_use_mat <- matrix(clusters_to_use, nrow=length(assignment_mechanisms), ncol=N_assignment_mechanism_sets, byrow=TRUE)
  clusters_to_use_final <- c(clusters_to_use_mat)
  clusters_to_use_dex <- 1
  trial_dfs <- list()
  for (realized_assign in rep(assignment_mechanisms, N_groups / length(assignment_mechanisms))) {
    cluster_for_trial <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_mort0.85_eit0.005_vaxEff0.6_assign', realized_assign, '_sim', clusters_to_use_final[clusters_to_use_dex], '.csv'))
    trial_dfs[[clusters_to_use_dex]] <- cluster_for_trial
    clusters_to_use_dex <- clusters_to_use_dex + 1
  }
  # First, calculate percent infected within each cluster of different assignment
  sum_cluster_average <- rep(0, length(assignment_mechanisms))
  sampled_ls <- list()
  for (assignment_mech_dex in 1:length(assignment_mechanisms)) {
    assignment_mech_ls <- list()
    assignment_mech_ls_dex <- 1
    for (trial_df_dex in seq(assignment_mech_dex, N_groups, length(assignment_mechanisms))) {
      cluster_df <- trial_dfs[[trial_df_dex]]
      if (cluster_df$node[1] == 'na') {
        stop('Error in enrolled cluster')
      } else {
        # The following is to get the true effect
        denom_to_add <- length(which(cluster_df$assignment == 'na'))
        num_to_add <- length(which(cluster_df$assignment == 'na' & cluster_df$time2inf_trt < cutoff))
        cluster_average <- num_to_add / denom_to_add
        sum_cluster_average[assignment_mech_dex] <- sum_cluster_average[assignment_mech_dex] + cluster_average
        # The following is prepare the dataframe for the logistic analysis
        mod_cluster_df <- cluster_df
        mod_cluster_df <- mod_cluster_df[which(mod_cluster_df$assignment == 'na'),]
        mod_cluster_df$assignment_mech <- assignment_mech_dex
        mod_cluster_df$clus_num <- paste0(assignment_mech_dex, '_', trial_df_dex)
        mod_cluster_df$status <- ifelse(mod_cluster_df$time2inf_trt < cutoff, 1, 0) # Where 1 is infected, 0 is censored
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
      low_df$cond <- 0
      high_df$cond <- 1
      to_analyze <- rbind(low_df, high_df)
      # Perform t-test
      
      
      
      
      res.logistic <- glm(formula=status ~ factor(cond), family=binomial(link = "logit"), data=to_analyze)
      pval_res <- c(pval_res, summary(res.logistic)[12]$coefficients[8])
      if (summary(res.logistic)[12]$coefficients[8] < 0.05) {
        est_eff_res <- c(est_eff_res, summary(res.logistic)[12]$coefficients[2])
      } else {
        est_eff_res <- c(est_eff_res, NA)
      }
    } else {
      est_eff_res <- c(est_eff_res, NA)
      pval_res <- c(pval_res, NA)
    }
  }
  final_est_eff_ls[[final_ls_dex]] <- data.frame(matrix(est_eff_res, nrow=1, ncol=(length(assignment_mechanisms) - 1)))
  final_pval_ls[[final_ls_dex]] <- data.frame(matrix(pval_res, nrow=1, ncol=(length(assignment_mechanisms) - 1)))
  
  # Second, get true effect (estimand)
  ave_cluster_average <- sum_cluster_average / N_assignment_mechanism_sets
  temp_df <- data.frame(matrix(ave_cluster_average[1:length(ave_cluster_average) - 1] - ave_cluster_average[2:length(ave_cluster_average)], nrow=1, ncol=length(ave_cluster_average) - 1))
  final_ls[[final_ls_dex]] <- temp_df
  final_ls_dex <- final_ls_dex + 1
}
final_df <- do.call(rbind, final_ls)
final_est_eff_df <- do.call(rbind, final_est_eff_ls)
final_pval_df <- do.call(rbind, final_pval_ls)

# Q1) What is the the mean and variance of the estimated effects under this trial design?
for (col_dex in 1:ncol(final_est_eff_df)) {
  print(paste0('Effect estimate: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Mean: ', mean(final_est_eff_df[,col_dex], na.rm=T)))
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

