# Script used to simulate two-stage trials -------------------------------------
# Libraries --------------------------------------------------------------------
library("stats")
library("lme4")

# Parameters -------------------------------------------------------------------
set.seed(0)
N_sims = 2000 # Total number of cluster simulations in simulation bank
N_sample = 100 # Number sampled from each cluster
N_trials = 1000 # Number of trial simulations to conduct
cutoff = 120
num_bootstrap_sample = 100
assignment_mechanisms = c(0, 0.1)
N_assignment_mechanism_sets = 2
N_groups = length(assignment_mechanisms) * N_assignment_mechanism_sets
R0_vax = 1.1
if (N_groups %% length(assignment_mechanisms) != 0) {
  stop('The number of groups should be divisible by the number of assignment mechanisms.')
}
threshold_inclusion = 3

# Get simulations to use for each assignment mechanism -------------------------
to_use_ls <- list()
to_use_ls_dex <- 1
for (assignment_mechanism in assignment_mechanisms) {
  to_use <- c()
  for (sim_num in 0:(N_sims - 1)) {
    test_cluster <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', assignment_mechanism, '_sim', sim_num, '_SEIR.csv'))
    if (test_cluster$node[1] != 'na') {
      to_use <- c(to_use, sim_num)
    }
  }
  if (length(to_use) != N_sims) {
    stop('Should use all sims.')
  }
  to_use_ls[[to_use_ls_dex]] <- to_use
  to_use_ls_dex <- to_use_ls_dex + 1
}

# Function to run trial simulation ---------------------------------------------
run_trial <- function(trial_num) {
  clusters_to_use <- c()
  for (assignment_mechanism_dex in 1:length(assignment_mechanisms)) {
    clusters_to_use <- c(clusters_to_use, sample(to_use_ls[[assignment_mechanism_dex]], N_assignment_mechanism_sets))
  }
  clusters_to_use_mat <- matrix(clusters_to_use, nrow=length(assignment_mechanisms), ncol=N_assignment_mechanism_sets, byrow=TRUE)
  clusters_to_use_final <- c(clusters_to_use_mat)
  clusters_to_use_dex <- 1
  trial_dfs <- list()
  for (realized_assign in rep(assignment_mechanisms, N_groups / length(assignment_mechanisms))) {
    cluster_for_trial <- read.csv(paste0('~/netVax/code_output/twostage/sims/2stg_N1000_k1_R0wt3_R0vax', R0_vax, '_eit0.005_vaxEff1_assign', realized_assign, '_sim', clusters_to_use_final[clusters_to_use_dex], '_SEIR.csv'))
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
  bs_est_eff_res <- c()
  pval_res <- c()
  for (sampled_dex in 1:(length(sampled_ls) - 1)) {
    low_df <- sampled_ls[[sampled_dex]]
    high_df <- sampled_ls[[sampled_dex + 1]]
    if (ncol(low_df) > 1 & ncol(high_df) > 1) {
      low_df$cond <- 0
      high_df$cond <- 1
      to_analyze <- rbind(low_df, high_df)
      # Get rid of those with a single outcome (either all the unvaccinated were infected, or all not infected)
      to_delete <- c()
      clus_num_f_cnt <- 0
      clus_num_f_cont_cnt <- 0
      for (clus_num_f in unique(to_analyze$clus_num)) {
        if (length(unique(to_analyze[which(to_analyze$clus_num == clus_num_f),]$status)) == 1 | 
            length(which(to_analyze$clus_num == clus_num_f & to_analyze$status == 1)) < threshold_inclusion) {
          to_delete <- c(to_delete, which(to_analyze$clus_num == clus_num_f))
          clus_num_f_cnt <- clus_num_f_cnt + 1
          if (to_analyze[which(to_analyze$clus_num == clus_num_f),]$cond[1] == 0) {
            clus_num_f_cont_cnt <- clus_num_f_cont_cnt + 1
          }
        }
      }
      if (!is.null(to_delete)) {
        to_analyze <- to_analyze[-to_delete,]
      }
      if (length(which(to_analyze$assignment_mech == 1)) < 101 | length(which(to_analyze$assignment_mech == 2)) < 101) { # make sure contrast error does not come up. think it needs at least two clusters per arm for the contrast.
        pval_res <- c(pval_res, NA)
        est_eff_res <- c(est_eff_res, NA)
        bs_est_eff_res <- c(bs_est_eff_res, NA)
      } else {
        res.logistic <- glmer(formula=status ~ factor(cond) + (1|clus_num), family=binomial(link = "logit"), data=to_analyze)
        pval <- summary(res.logistic)[10]$coefficients[8]
        pval_res <- c(pval_res, pval)
        if (pval < 0.05) {
          # Use predict function
          to_predict <- data.frame(matrix(c(0, 1), nrow=2, ncol=1))
          colnames(to_predict) <- c('cond')
          rownames(to_predict) <- c('low', 'high')
          predicted <- predict(res.logistic, to_predict, type='response',re.form=NA)
          if (length(unique(predicted)) != 2) {
            stop(unique(predicted))
            #stop('Error in number of predicted values.')
          }
          est_eff_res <- c(est_eff_res, unname(predicted[1] - predicted[2]))
          # Do bootstrap estimate
          bs_ests <- c()
          for (iter in 1:num_bootstrap_sample) {
            # 1) construct bootstrap sample:
            passed <- FALSE
            while (!passed) {
              bs_samp_ls <- list()
              bs_samp_ls_dex <- 1
              control_to_analyze <- to_analyze[which(to_analyze$cond == 0),]
              N_clus_control <- as.integer(nrow(control_to_analyze) / 100)
              treated_to_analyze <- to_analyze[which(to_analyze$cond == 1),]
              N_clus_treated <- as.integer(nrow(treated_to_analyze) / 100)
              if (nrow(control_to_analyze) %% 100 != 0 | nrow(treated_to_analyze) %% 100 != 0) {
                stop('Error in number of observations in either treatment or control.')
              }
              for (index1 in 1:N_clus_control) {
                clus_samp <- control_to_analyze[which(control_to_analyze$clus_num == sample(unique(control_to_analyze$clus_num), 1), 1),]
                # First, get of the infecteds, randomly sample at least the number of threshold inclusion
                infecteds_clus_samp <- which(clus_samp$status == 1)
                if (length(infecteds_clus_samp) < threshold_inclusion) {
                  stop('Error in threshold inclusion for bootstrap.')
                }
                not_infecteds_clus_samp <- which(clus_samp$status == 0)
                chosen_infecteds_clus_samp <- sample(infecteds_clus_samp, threshold_inclusion, replace=TRUE) # A little iffy here... since first we sample among infecteds, then everyone else #FIX THIS PUT IN A STOP IF TOO LONG SEARCH
                rest_bs_indices <- sample(1:nrow(clus_samp), nrow(clus_samp) - threshold_inclusion, replace = TRUE) # Randomly choose the rest of the indices
                bs_samp_ls[[bs_samp_ls_dex]] <- clus_samp[c(chosen_infecteds_clus_samp, rest_bs_indices),]
                bs_samp_ls_dex <- bs_samp_ls_dex + 1
              }
              for (index2 in 1:N_clus_treated) {
                clus_samp <- treated_to_analyze[which(treated_to_analyze$clus_num == sample(unique(treated_to_analyze$clus_num), 1), 1),]
                # First, get of the infecteds, randomly sample at least the number of threshold inclusion
                infecteds_clus_samp <- which(clus_samp$status == 1)
                if (length(infecteds_clus_samp) < threshold_inclusion) {
                  stop('Error in threshold inclusion for bootstrap.')
                }
                not_infecteds_clus_samp <- which(clus_samp$status == 0)
                chosen_infecteds_clus_samp <- sample(infecteds_clus_samp, threshold_inclusion, replace=TRUE) # A little iffy here... since first we sample among infecteds, then everyone else #FIX THIS PUT IN A STOP IF TOO LONG SEARCH
                rest_bs_indices <- sample(1:nrow(clus_samp), nrow(clus_samp) - threshold_inclusion, replace = TRUE) # Randomly choose the rest of the indices
                bs_samp_ls[[bs_samp_ls_dex]] <- clus_samp[c(chosen_infecteds_clus_samp, rest_bs_indices),]
                bs_samp_ls_dex <- bs_samp_ls_dex + 1
              }
              bs_samp <- do.call(rbind, bs_samp_ls)
              bs_samp <- bs_samp[order(bs_samp$cond),]
              if (nrow(bs_samp) != nrow(to_analyze)) {
                stop('Error in creating bootstrap. (1)')
              }
              if (length(unique(bs_samp$cond)) == 2) {# & 
                  #length(which(bs_samp$status == 1 & bs_samp$cond == 0)) / length(which(bs_samp$cond == 0)) != length(which(bs_samp$status == 1 & bs_samp$cond == 1)) / length(which(bs_samp$cond == 1))) { 
                passed <- TRUE
              } else if (length(unique(bs_samp$cond)) > 2) {
                stop('Error in creating bootstrap. (2)')
              }
            }
            # 2) run logistic and predict
            res.logistic_bs <- glmer(formula=status ~ factor(cond) + (1|clus_num), family=binomial(link = "logit"), data=bs_samp)
            bs_pval <- summary(res.logistic_bs)[10]$coefficients[8]
            if (bs_pval < 0.05) {
              to_predict <- data.frame(matrix(c(0, 1), nrow=2, ncol=1))
              colnames(to_predict) <- c('cond')
              rownames(to_predict) <- c('low', 'high')
              predicted_bs <- predict(res.logistic_bs, to_predict, type='response', re.form=NA)
              if (length(unique(predicted_bs)) != 2) {
                stop(unique(predicted_bs))
                #stop('Error in number of predicted values.')
              }
              bs_ests <- c(bs_ests, unname(predicted_bs[1] - predicted_bs[2]))
            }
          }
          bs_est_eff_res <- c(bs_est_eff_res, unname(abs(quantile(bs_ests, probs=c(.025)) - quantile(bs_ests, probs=c(.975)))))
        } else {
          est_eff_res <- c(est_eff_res, NA)
          bs_est_eff_res <- c(bs_est_eff_res, NA)
        }
      }
    } else {
      est_eff_res <- c(est_eff_res, NA)
      pval_res <- c(pval_res, NA)
    }
  }
  write.csv(c(clus_num_f_cnt / N_groups, (clus_num_f_cont_cnt / N_groups) / 2), paste0('~/netVax/scratch/', trial_num, '_overreactionary_tes.csv'))
  # Code to get true ASE effect (temp_df)
  ave_cluster_average <- sum_cluster_average / N_assignment_mechanism_sets
  temp_df <- data.frame(matrix(ave_cluster_average[1:length(ave_cluster_average) - 1] - ave_cluster_average[2:length(ave_cluster_average)], nrow=1, ncol=length(ave_cluster_average) - 1))
  # Code to store the true effect (est_eff_res) the bootstraps (bs_est_eff_res) and the p_values (p_val_res) and true ASE effect (temp_df)
  to_return_ls = list(data.frame(matrix(est_eff_res, nrow=1, ncol=(length(assignment_mechanisms) - 1))),
                      data.frame(matrix(bs_est_eff_res, nrow=1, ncol=(length(assignment_mechanisms) - 1))),
                      data.frame(matrix(pval_res, nrow=1, ncol=(length(assignment_mechanisms) - 1))),
                      temp_df)
  return(to_return_ls)
}

# Run in parallel --------------------------------------------------------------
library(foreach)
library(doParallel)
library(pracma)
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
final <- foreach(i=1:N_trials) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  library(lme4)
  res = run_trial(i)
  res
}
stopCluster(cl)
save(final, file = paste0("~/netVax/code_output/twostage/rData/final", R0_vax, "_SEIR.RData"))

# Load results -----------------------------------------------------------------
load(paste0("~/netVax/code_output/twostage/rData/final", R0_vax, "_SEIR.RData"))
final_est_eff_res <- list()
final_bs_est_eff_res <- list()
final_pval_res <- list()
final_ASE <- list()
list_dex <- 1
for (i in 1:length(final)) {
  final_est_eff_res[[list_dex]] <- final[[list_dex]][[1]]
  final_bs_est_eff_res[[list_dex]] <- final[[list_dex]][[2]]
  final_pval_res[[list_dex]] <- final[[list_dex]][[3]]
  final_ASE[[list_dex]] <- final[[list_dex]][[4]]
  list_dex <- list_dex + 1
}
final_est_eff_res_df <- do.call(rbind, final_est_eff_res)
final_bs_est_eff_res_df <- do.call(rbind, final_bs_est_eff_res)
final_pval_res_df <- do.call(rbind, final_pval_res)
final_ASE_df <- do.call(rbind, final_ASE)

# Answer questions -------------------------------------------------------------
# Q1) What is the the mean and variance of the estimated effects under this trial design?
for (col_dex in 1:ncol(final_est_eff_res_df)) {
  print(paste0('Effect estimate: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Mean: ', mean(final_est_eff_res_df[,col_dex] * 100, na.rm=T)))
  print(paste0('SD: ', sd(final_est_eff_res_df[,col_dex] * 100, na.rm=T)))
  print(paste0('Mean of the width of the 95% bootstrap CI: ', mean(final_bs_est_eff_res_df[,col_dex] * 100, na.rm=T)))
}

# Q2) What is the power to detect an effect under this trial design using the logistic model?
for (col_dex in 1:ncol(final_pval_res_df)) {
  print(paste0('Power: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(length(which(final_pval_res_df[,col_dex] < 0.05)) / length(which(!is.na(final_pval_res_df[,col_dex]))))
}

# Q3) What is the true effect across groups (mean and variance)?
for (col_dex in 1:ncol(final_ASE_df)) {
  print(paste0('True effect: Reduced incidence of treatment assignment from:', assignment_mechanisms[col_dex],' to: ', assignment_mechanisms[col_dex + 1]))
  print(paste0('Mean: ', mean(final_ASE_df[,col_dex] * 100, na.rm=T)))
  print(paste0('SD: ', sd(final_ASE_df[,col_dex] * 100, na.rm=T)))
}

# Average spillover effect histogram
hist(final_ASE_df[,2], main="From 0% to 10% Vaccinated", xlab="True Average Spillover Effect")
abline(v=0, col='red', lwd=3, lty=2)

# Type I error histogram
hist(final_ASE_df[,1], main="Comparing 0% to 0% Vaccinated", xlab="True Average Spillover Effect")
abline(v=0, col='red', lwd=3, lty=2)

