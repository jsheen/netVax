# ------------------------------------------------------------------------------
# @author: Justin Sheen
#
# @description: script to analyze survival data of transmissible vaccine trial
# ------------------------------------------------------------------------------

# Import libraries -------------------------------------------------------------
library('survival')
library('survminer')

# Parameters -------------------------------------------------------------------
right_censoring <- 90
N_clusters <- c(1000, 10000)
overdispersions <- c(1)
R0_wts <- c(3)
vaxs <- c('R0=0_treat=0.5', 'R0=1.5_treat=0.5', 'R0=2_treat=0.5',
          'R0=0_treat=0.2', 'R0=1.5_treat=0.2', 'R0=2_treat=0.2')
morts <- c(0.85)

# Go through simulation results and calculate cox proportional hazard ----------
for (sim_num in seq(0, 2999, 1)) {
  
}