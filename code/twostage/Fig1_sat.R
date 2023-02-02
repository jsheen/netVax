# Script to create figures showing saturation of both the indirect and overall
# effects of vaccination if the R0 of the vaccine is low vs. high --------------

# Libraries --------------------------------------------------------------------
library(ggplot2)

# Parameters -------------------------------------------------------------------
set.seed(0)
N_sims = 1000 # Total number of cluster simulations in simulation bank
N_trials = 1000 # Number of trial simulations to conduct
cutoff = 90
assignment_mechanisms = c(0, 0.1, 0.25)
N_assignment_mechanism_sets = 4
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
final_ls_baseline <- list()
final_ls_dex <- 1
for (trial_num in 1:N_trials) {
  if (trial_num %% 200 == 0) {
    print(trial_num)
  }
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
  for (assignment_mech_dex in 1:length(assignment_mechanisms)) {
    for (trial_df_dex in seq(assignment_mech_dex, N_groups, length(assignment_mechanisms))) {
      cluster_df <- trial_dfs[[trial_df_dex]]
      if (cluster_df$node[1] == 'na') {
        stop('Error in enrolled cluster')
      } else {
        # The following is to get the true indirect effect
        denom_to_add <- length(which(cluster_df$assignment == 'na'))
        num_to_add <- length(which(cluster_df$assignment == 'na' & cluster_df$time2inf_trt < cutoff))
        cluster_average <- num_to_add / denom_to_add
        sum_cluster_average[assignment_mech_dex] <- sum_cluster_average[assignment_mech_dex] + cluster_average
      }
    }
  }
  # First, get true ASE effect (estimand)
  ave_cluster_average <- sum_cluster_average / N_assignment_mechanism_sets
  temp_df <- data.frame(matrix(ave_cluster_average[1:length(ave_cluster_average) - 1] - ave_cluster_average[2:length(ave_cluster_average)], nrow=1, ncol=length(ave_cluster_average) - 1))
  final_ls[[final_ls_dex]] <- temp_df
  
  # Second, get true indirect effect (i.e. % not infected of unvaccinated)
  temp_df <- data.frame(matrix((1 - ave_cluster_average) * (1 - assignment_mechanisms), nrow=1, ncol=length(ave_cluster_average)))
  final_ls_baseline[[final_ls_dex]] <- temp_df

  # Add iterator
  final_ls_dex <- final_ls_dex + 1
}
final_df <- do.call(rbind, final_ls)
final_df_baseline <- do.call(rbind, final_ls_baseline)

# Plot 1: ASE ------------------------------------------------------------------
ASE_ls <- list()
ASE_ls_dex <- 1
for (col_dex in 1:ncol(final_df)) {
  new_row <- data.frame(matrix(c(paste0(assignment_mechanisms[col_dex] * 100, '% to ', assignment_mechanisms[col_dex + 1] * 100, '%'), mean(final_df[,col_dex] * 100, na.rm=T), sd(final_df[,col_dex] * 100, na.rm=T)), nrow=1, ncol=3))
  ASE_ls[[ASE_ls_dex]] <- new_row
  ASE_ls_dex <- ASE_ls_dex + 1
}
ASE_df <- do.call(rbind, ASE_ls)
colnames(ASE_df) <- c('perc_vax', 'mean', 'sd')
ASE_plot <- ggplot(ASE_df) +
  geom_bar(aes(x=perc_vax, y=as.numeric(mean)), stat="identity", fill="skyblue", alpha=0.7) +
  #geom_errorbar(aes(x=perc_vax, ymin=as.numeric(mean)-as.numeric(sd), ymax=as.numeric(mean)+as.numeric(sd)), width=.2, position=position_dodge(.9)) +
  ggtitle(expression("Average Spillover Effect (ASE) for Vax "*R[0]*"= 1.1")) + ylab('% Reduction in Incidence of Unvaccinated Pop. (ASE)') + xlab('Increase in Vaccination') +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x= element_text(size=14),
    axis.text.y= element_text(size=14)
  )
ASE_plot

# Plot 2: Indirect, direct, overall effects ------------------------------------
vaxEff_plot_ls <- list()
vaxEff_plot_ls_dex <- 1
for (col_dex in 1:ncol(final_df_baseline)) {
  new_row <- data.frame(matrix(c(paste0(assignment_mechanisms[col_dex] * 100, '%'),
                                 'a) vaccinated,\nuninfect',
                                 assignment_mechanisms[col_dex] * 100,
                                 paste0(assignment_mechanisms[col_dex] * 100, '%'),
                                 'b) unvaccinated,\nuninfect',
                                 mean(final_df_baseline[,col_dex], na.rm=T) * 100), nrow=2, ncol=3, byrow=T))
  vaxEff_plot_ls[[vaxEff_plot_ls_dex]] <- new_row
  vaxEff_plot_ls_dex <- vaxEff_plot_ls_dex + 1
}
vaxEff_plot_df <- do.call(rbind, vaxEff_plot_ls)
colnames(vaxEff_plot_df) <- c('perc_vax', 'cond', 'value')
vaxEff_plot <- ggplot(vaxEff_plot_df, aes(fill=cond, y=as.numeric(value), x=perc_vax)) +
  geom_bar(position="stack", stat="identity") +
  ggtitle(expression("Protective Effects for Vax "*R[0]*"= 1.1")) + ylab('% Not Infected') + xlab('% Vaccinated') +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x= element_text(size=14),
    axis.text.y= element_text(size=14),
    legend.position='none'
  )
vaxEff_plot




