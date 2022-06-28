# ------------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script to analyze survival data of transmissible vaccine trial
# ------------------------------------------------------------------------------

# Import libraries -------------------------------------------------------------
library('survival')
library('survminer')

# Parameters -------------------------------------------------------------------
right_censor <- 90
N_clusters <- c(1000)
overdispersions <- c(1)
R0_wts <- c(3)
vaxs <- c('R0=0_treat=0.5', 'R0=1.5_treat=0.5', 'R0=2_treat=0.5',
          'R0=0_treat=0.2', 'R0=1.5_treat=0.2', 'R0=2_treat=0.2')
morts <- c(0.85)
vaxEffs <- c(0.6)
param_sets <- list()
param_sets_dex <- 1
for (i in N_clusters) {
  for (j in overdispersions) {
    for (k in R0_wts) {
      for (l in vaxs) {
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

# Go through simulation results and calculate cox proportional hazard ----------
res <- list()
res_dex <- 1
for (param_set in param_sets) {
  N_cluster <- param_set[1]
  overdispersion <- param_set[2]
  R0_wt <- param_set[3]
  vax <- param_set[4]
  mort <- param_set[5]
  vax_eff <- param_set[6]
  eit <- 0.005
  R0_vax <- gsub('R0=', '', strsplit(vax, "_")[[1]][1])
  vax_treat <- gsub('treat=', '', strsplit(vax, "_")[[1]][2])
  zero_I <- 0
  effect_estimates_inf <- c()
  effect_estimates_death <- c()
  alpha_satis_inf <- 0
  alpha_satis_death <- 0
  for (sim_num in seq(0, 2999, 1)) {
    filename <- paste0('~/netVax/code_output/sim_results/N', N_cluster, '_k', overdispersion, 
                       '_R0wt', R0_wt, '_R0vax', R0_vax, '_mort', mort, '_eit', eit, '_vaxTreat', vax_treat, 
                       '_vaxEff', vax_eff, '_sim', sim_num, '.csv')
    sim_res <- read.csv(filename)
    if (sim_res$node[1] =='na') {
      sim_res <- sim_res[-which(sim_res$assignment == 'na'),]
      sim_res_infect <- sim_res
      sim_res_death <- sim_res
      # Infection
      sim_res_infect$time2death <- NULL
      sim_res_infect$time2inf <- ifelse(sim_res_infect$time2inf > right_censor, right_censor, sim_res_infect$time2inf)
      sim_res_infect$status <- ifelse(sim_res_infect$time2inf == right_censor, 1, 2) # 1=right censored, 2=infected
      res_cox_inf <- coxph(Surv(time2inf, status) ~ assignment, data = sim_res_infect)
      res_cox_inf
    
      # Death outcome
      sim_res_death$time2inf <- NULL
      sim_res_death$time2death <- ifelse(sim_res_death$time2death > right_censor, right_censor, sim_res_death$time2death)
      sim_res_death$status <- ifelse(sim_res_death$time2death == right_censor, 1, 2) # 1=right censored, 2=dead
      res_cox_inf <- coxph(Surv(time2death, status) ~ assignment, data = sim_res_death)
      
      
    } else {
      zero_I <- zero_I + 1
    }
  }
  new_row <- data.frame(matrix(ncol=12, nrow=1))
  colnames(new_row) <- c('N', 'k', 'R0_wt', 'R0_vax', 'mort', 'eit', 'vaxTreat', 'vaxEff', 
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
  new_row$effect_estimates_inf[1] <- paste0(effect_estimates_inf)
  new_row$effect_estimates_death[1] <- paste0(effect_estimates_death)
  new_row$power_inf[1] <- alpha_satis_inf / 3000
  new_row$power_death[1] <- alpha_satis_death / 3000
  new_row$zero_I[1] <- zero_I
  res[[res_dex]] <- new_row
  res_dex <- res_dex + 1
}

# Plotting ---------------------------------------------------------------------











