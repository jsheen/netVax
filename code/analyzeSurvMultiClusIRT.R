# ------------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script to analyze survival data of transmissible vaccine trial
# ------------------------------------------------------------------------------

# Import libraries -------------------------------------------------------------
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)
set.seed(0)

# Parameters -------------------------------------------------------------------
right_censor <- 120
N_clusters <- c(1000)
overdispersions <- c(1)
R0_wts <- c(3)
 vaxs <- c('R0=0_treat=0.5', 'R0=1.1_treat=0.5', 'R0=1.1_treat=0.1')
morts <- c(0.85)
vaxEffs <- c(0.6)
multi_clusters <- c(2, 4, 10, 20)
param_sets <- list()
param_sets_dex <- 1
for (i in N_clusters) {
  for (j in overdispersions) {
    for (k in R0_wts) {
      for (l in vaxs) {
        for (m in morts) {
          for (n in vaxEffs) {
            for (o in multi_clusters) {
              param_sets[[param_sets_dex]] <- c(i, j, k, l, m, n, o)
              param_sets_dex <- param_sets_dex + 1
            }
          }
        }
      }
    }
  }
}

# Go through simulation results and calculate cox proportional hazard ----------
res <- list()
res_dex <- 1
for (param_set in param_sets) {
  print(param_set)
  N_cluster <- param_set[1]
  overdispersion <- param_set[2]
  R0_wt <- param_set[3]
  vax <- param_set[4]
  mort <- param_set[5]
  vax_eff <- param_set[6]
  multi_cluster <- param_set[7]
  eit <- 0.005
  R0_vax <- gsub('R0=', '', strsplit(vax, "_")[[1]][1])
  vax_treat <- gsub('treat=', '', strsplit(vax, "_")[[1]][2])
  runSim <- function(N_cluster, overdispersion, R0_wt, R0_vax, mort, eit, vax_treat, vax_eff, sim_num, multi_cluster) {
    if (multi_cluster == 1) {
      filename <- paste0('~/netVax_sim_results/N', N_cluster, '_k', overdispersion, 
                         '_R0wt', R0_wt, '_R0vax', R0_vax, '_mort', mort, '_eit', eit, '_vaxTreat', vax_treat, 
                         '_vaxEff', vax_eff, '_sim', sim_num, '.csv')
      single_clus <- read.csv(filename)
      single_clus$clusNum <- 1
      if (single_clus$node[1] != 'na') {
        single_clus <- single_clus[-which(single_clus$assignment == 'na'),]
        # (Infection) Take care of censoring and status
        # Get end of simulation time for infections
        end_sim_time_inf <- max(single_clus$time2inf)
        # If equal to end_sim_time, then was never infected (status == 1), else was infected (status == 2)
        single_clus$status_inf <- ifelse(single_clus$time2inf == end_sim_time_inf, 1, 2)
        # For all that were never infected (status == 1), change to right censor. These nodes have end of simulation
        # time as their time2hazard. If end of sim time is less than right censor, should be replaced because we censor
        # them at right censor. If end of sim time is greater than right censor, should be replaced because we also 
        # stop observation at right censor. If equal, nothing will happen with the following code block. Those with
        # status 2, i.e. infected at some point of the simulation, will not be touched.
        single_clus$time2inf <- ifelse(single_clus$status_inf == 1, right_censor, single_clus$time2inf)
        # If end_sim_time_inf > right_censor, then possible that some that were infected should be censored as they
        # have infection times greater than right censor. If end_sim_time_inf <= right_censor, then none of the
        # time2infs need to change.
        if (end_sim_time_inf > right_censor) {
          single_clus$time2inf <- ifelse(single_clus$status_inf == 2 & single_clus$time2inf > right_censor,
                                         right_censor, single_clus$time2inf)
        }
        # (Death) Take care of censoring and status
        # Get end of simulation time for deaths
        end_sim_time_death <- max(single_clus$time2death)
        # If equal to end_sim_time, then was never infected (status == 1), else was infected (status == 2)
        single_clus$status_death <- ifelse(single_clus$time2death == end_sim_time_death, 1, 2)
        # For all that were never infected (status == 1), change to right censor. These nodes have end of simulation
        # time as their time2hazard. If end of sim time is less than right censor, should be replaced because we censor
        # them at right censor. If end of sim time is greater than right censor, should be replaced because we also 
        # stop observation at right censor. If equal, nothing will happen with the following code block. Those with
        # status 2, i.e. infected at some point of the simulation, will not be touched.
        single_clus$time2death <- ifelse(single_clus$status_death == 1, right_censor, single_clus$time2death)
        # If end_sim_time_death > right_censor, then possible that some that were infected should be censored as they
        # have infection times greater than right censor. If end_sim_time_death <= right_censor, then none of the
        # time2deathss need to change.
        if (end_sim_time_death > right_censor) {
          single_clus$time2death <- ifelse(single_clus$status_death == 2 & single_clus$time2death > right_censor,
                                           right_censor, single_clus$time2death)
        }
      }
      sim_res <- single_clus
    } else {
      sim_res_ls <- list()
      used <- c()
      for (cluster_dex in 1:multi_cluster) {
        is_empty_df <- TRUE
        while (is_empty_df) {
          to_join <- sample(c(0:2999), size=1)
          filename <- paste0('~/netVax_sim_results/N', N_cluster, '_k', overdispersion, 
                             '_R0wt', R0_wt, '_R0vax', R0_vax, '_mort', mort, '_eit', eit, '_vaxTreat', vax_treat, 
                             '_vaxEff', vax_eff, '_sim', to_join, '.csv')
          single_clus <- read.csv(filename)
          if (single_clus$node[1] != 'na' & !(to_join %in% used)) {
            is_empty_df <- FALSE
            used <- c(used, to_join)
          }
        }
        single_clus$clusNum <- cluster_dex
        single_clus <- single_clus[-which(single_clus$assignment == 'na'),]
        # (Infection) Take care of censoring and status
        end_sim_time_inf <- max(single_clus$time2inf)
        single_clus$status_inf <- ifelse(single_clus$time2inf == end_sim_time_inf, 1, 2)
        single_clus$time2inf <- ifelse(single_clus$status_inf == 1, right_censor, single_clus$time2inf)
        if (end_sim_time_inf > right_censor) {
          single_clus$time2inf <- ifelse(single_clus$status_inf == 2 & single_clus$time2inf > right_censor,
                                         right_censor, single_clus$time2inf)
        }
        # (Death) Take care of censoring and status
        end_sim_time_death <- max(single_clus$time2death)
        single_clus$status_death <- ifelse(single_clus$time2death == end_sim_time_death, 1, 2)
        single_clus$time2death <- ifelse(single_clus$status_death == 1, right_censor, single_clus$time2death)
        if (end_sim_time_death > right_censor) {
          single_clus$time2death <- ifelse(single_clus$status_death == 2 & single_clus$time2death > right_censor,
                                         right_censor, single_clus$time2death)
        }
        sim_res_ls[[cluster_dex]] <- single_clus
      }
      sim_res <- do.call(rbind, sim_res_ls)
    }
    if (sim_res$node[1] != 'na') {
      sim_res_infect <- sim_res
      sim_res_death <- sim_res
      # Infection
      sim_res_infect$time2death <- NULL
      sim_res_infect$status_death <- NULL
      res_cox_inf <- coxph(Surv(time2inf, status_inf) ~ assignment + strata(clusNum), data = sim_res_infect)
      # plot_df <- with(sim_res_infect, data.frame(assignment = c('c', 't')))
      # fit <- survfit(res_cox_inf, newdata = plot_df)
      # ggsurvplot(fit, data=sim_res_infect, conf.int = TRUE, legend.labs=c("Control", "Treatment"),
      #            ggtheme = theme_minimal())
      p_val_inf <- summary(res_cox_inf)[7]$coefficients[5]
      effect_estimate_inf <- 'na'
      if (p_val_inf <= 0.05) {
        effect_estimate_inf <- 1 - exp(unname(res_cox_inf$coefficients))
      }
      # Death outcome
      sim_res_death$time2inf <- NULL
      sim_res_death$status_inf <- NULL
      res_cox_death <- coxph(Surv(time2death, status_death) ~ assignment + strata(clusNum), data = sim_res_death)
      # plot_df <- with(sim_res_death, data.frame(assignment = c('c', 't')))
      # fit <- survfit(res_cox_death, newdata = plot_df)
      # ggsurvplot(fit, data=sim_res_death, conf.int = TRUE, legend.labs=c("Control", "Treatment"),
      #            ggtheme = theme_minimal())
      p_val_death <- summary(res_cox_death)[7]$coefficients[5]
      effect_estimate_death <- 'na'
      if (p_val_death <= 0.05) {
        effect_estimate_death <- 1 - exp(unname(res_cox_death$coefficients))
      }
      return(c(effect_estimate_inf, effect_estimate_death, p_val_inf, p_val_death))
    } else {
      return(c('na', 'na', 'na', 'na'))
    }
  }
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload computer
  registerDoParallel(cl)
  finalMatrix <- foreach(i=seq(0, 2999, 1), .combine=rbind) %dopar% {
    library(survival)
    library(survminer)
    library(foreach)
    library(doParallel)
    tempMatrix = runSim(N_cluster=N_cluster, overdispersion=overdispersion, R0_wt=R0_wt, R0_vax=R0_vax, 
                        mort=mort, eit=eit, vax_treat=vax_treat, vax_eff=vax_eff, multi_cluster=multi_cluster, sim_num=i)
    tempMatrix = data.frame(matrix(tempMatrix, ncol=4, nrow=1))
    colnames(tempMatrix) <- c('eff_est_inf', 'eff_est_death', 'p_val_inf', 'p_val_death')
    tempMatrix
  }
  stopCluster(cl)
  finalMatrix_copy <- finalMatrix
  if (length(which(finalMatrix_copy$p_val_inf == 'na')) != length(which(finalMatrix_copy$p_val_death == 'na'))) {
    stop("Incorrect number of zero_I.")
  }
  zero_I <- length(which(finalMatrix_copy$p_val_inf == 'na'))
  finalMatrix_copy <- finalMatrix_copy[which(finalMatrix_copy$p_val_inf != 'na'),] # Get rid of all that did not have enough infections at intervention, and thus no p-value
  if (nrow(finalMatrix_copy) + zero_I != nrow(finalMatrix)) {
    stop("Incorrect number of zero_I after deletion.")
  }
  power_inf <- length(which(as.numeric(finalMatrix_copy$p_val_inf) <= 0.05)) / 3000
  power_death <- length(which(as.numeric(finalMatrix_copy$p_val_death) <= 0.05)) / 3000
  effect_estimates_inf <- paste(finalMatrix_copy$eff_est_inf[which(finalMatrix_copy$eff_est_inf != 'na')], collapse='_')
  effect_estimates_death <- paste(finalMatrix_copy$eff_est_death[which(finalMatrix_copy$eff_est_death != 'na')], collapse='_')
  new_row <- data.frame(matrix(ncol=14, nrow=1))
  colnames(new_row) <- c('N', 'k', 'R0_wt', 'R0_vax', 'mort', 'eit', 'vaxTreat', 'vaxEff', 'multi_cluster',
                         'effect_estimates_inf', 'effect_estimates_death', 'power_inf', 'power_death',
                         'zero_I')
  new_row$N[1] <- N_cluster
  new_row$k[1] <- overdispersion
  new_row$R0_wt[1] <- R0_wt
  new_row$R0_vax[1] <- R0_vax
  new_row$mort[1] <- mort
  new_row$eit[1] <- eit
  new_row$vaxTreat[1] <- vax_treat
  new_row$vaxEff[1] <- vax_eff
  new_row$multi_cluster[1] <- multi_cluster
  new_row$effect_estimates_inf[1] <- effect_estimates_inf
  new_row$effect_estimates_death[1] <- effect_estimates_death
  new_row$power_inf[1] <- power_inf
  new_row$power_death[1] <- power_death
  new_row$zero_I[1] <- zero_I
  res[[res_dex]] <- new_row
  res_dex <- res_dex + 1
}
res_df <- do.call(rbind, res)

# Plotting ---------------------------------------------------------------------
l <- list()
l_dex <- 1
l2 <- list()
l2_dex <- 1
for (row_dex in 1:nrow(res_df)) {
  # Infection
  hist_df <- data.frame(as.numeric(strsplit(res_df$effect_estimates_inf[row_dex], '_')[[1]]))
  colnames(hist_df) <- c('effect_estimate')
  quant <- as.numeric(quantile(hist_df$effect_estimate, c(0.025, 0.975)))
  hist_df <- data.frame(hist_df[which(hist_df$effect_estimate > quant[1] & hist_df$effect_estimate < quant[2]),])
  colnames(hist_df) <- c('effect_estimate')
  l[[l_dex]] <- ggplot(hist_df, aes(x=effect_estimate)) + 
    ggtitle(paste0('Inf: R0_vax=', res_df$R0_vax[row_dex], '; vax_treat=', res_df$vaxTreat[row_dex],
                   '; Power=', round(res_df$power_inf[row_dex], digits=2), '; Clus=', res_df$multi_cluster[row_dex])) + 
    geom_histogram(bins=30) + geom_vline(xintercept=0.6, col='red') + theme(plot.title = element_text(size = 5))
  l_dex <- l_dex + 1
  # Death
  hist_df <- data.frame(as.numeric(strsplit(res_df$effect_estimates_death[row_dex], '_')[[1]]))
  colnames(hist_df) <- c('effect_estimate')
  quant <- as.numeric(quantile(hist_df$effect_estimate, c(0.025, 0.975)))
  hist_df <- data.frame(hist_df[which(hist_df$effect_estimate > quant[1] & hist_df$effect_estimate < quant[2]),])
  colnames(hist_df) <- c('effect_estimate')
  l2[[l2_dex]] <- ggplot(hist_df, aes(x=effect_estimate)) + 
    ggtitle(paste0('Death: R0_vax=', res_df$R0_vax[row_dex], '; vax_treat=', res_df$vaxTreat[row_dex], 
                   '; Power=', round(res_df$power_death[row_dex], digits=2),'; Clus=', res_df$multi_cluster[row_dex])) + 
    geom_histogram(bins=30) + geom_vline(xintercept=0.6, col='red') + theme(plot.title = element_text(size = 5))
  l2_dex <- l2_dex + 1
}
ggsave(filename="~/netVax/code_output/plots/N=1000_inf_multi.pdf", marrangeGrob(grobs = l, nrow=1, ncol=4),
       width=10, height=3, units='in', dpi=600)
ggsave(filename="~/netVax/code_output/plots/N=1000_death_multi.pdf", marrangeGrob(grobs = l2, nrow=1, ncol=4),
       width=10, height=3, units='in', dpi=600)


