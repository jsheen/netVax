# Script used to simulate two-stage trials -------------------------------------
set.seed(0)
N_trials = 1000
N_groups = 5
N_sims = 1000
assignment_mechanisms = c(0, 0.1, 0.25, 0.5, 0.75)
for (trial in 1:N_trials) {
  clusters_to_use <- sample(N_sims, N_groups)
  # Get the treatment assignment of this one
  
}

