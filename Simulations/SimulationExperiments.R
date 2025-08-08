###############################################################################
# Simulation study runner
# - Executes multiple configurations of the simulation pipeline
# - Collates detection-rate summaries and quick p-value checks
###############################################################################

library(tidyverse)
library(data.table)  # for rbindlist()

# Fix the RNG seed so runs are reproducible end-to-end
set.seed(42)

# -----------------------------------------------------------------------------
# 2-parent scenarios
# -----------------------------------------------------------------------------

# All parents up-regulated
result_2parentsUp <- run_simulation_pipeline(
  DEG_genes = c(1:2),          # true parent indices
  n_genes = 20,                # total genes simulated
  n_timepoints = 200,          # time points per experiment
  amplitude_range = c(0, 1000),# expression amplitude range
  n_repetitions = 5,           # batches
  experiments_per_rep = 100,   # experiments per batch
  regulation_type = "up"       # direction used to simulate target
)

# Mixed up/down regulation across parents
result_2parentsBoth <- run_simulation_pipeline(
  DEG_genes = c(1:2),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "both"
)

# All parents down-regulated
result_2parentsDown <- run_simulation_pipeline(
  DEG_genes = c(1:2),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "down"
)

# -----------------------------------------------------------------------------
# 3-parent scenarios
# -----------------------------------------------------------------------------

# All parents up-regulated
result_3parentsUp <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "up"
)

# Mixed up/down regulation across parents
result_3parentsBoth <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "both"
)

# All parents down-regulated
result_3parentsDown <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "down"
)

# -----------------------------------------------------------------------------
# Summaries & quick checks
# - calculate_detection_rates(): your metric aggregator
# - pval_minimal_top(): helper to inspect minimum/top p-values vs truth
# -----------------------------------------------------------------------------

############ 2 parents results ############

# 2 parents UP
aux <- calculate_detection_rates(
  df = result_2parentsUp[["summaries"]],
  n_true_parents = 2
)
aux
pval_minimal_top(aux, n_true = 2)

# 2 parents DOWN
aux <- calculate_detection_rates(
  df = result_2parentsDown[["summaries"]],
  n_true_parents = 2
)
aux
pval_minimal_top(aux, n_true = 2)

# 2 parents MIXED
aux <- calculate_detection_rates(
  df = result_2parentsBoth[["summaries"]],
  n_true_parents = 2
)
aux
pval_minimal_top(aux, n_true = 2)

############ 3 parents results ############

# 3 parents UP
aux <- calculate_detection_rates(
  df = result_3parentsUp[["summaries"]],
  n_true_parents = 3
)
aux
pval_minimal_top(aux, n_true = 3)

# 3 parents DOWN
aux <- calculate_detection_rates(
  df = result_3parentsDown[["summaries"]],
  n_true_parents = 3
)
aux
pval_minimal_top(aux, n_true = 3)

# 3 parents MIXED
aux <- calculate_detection_rates(
  df = result_3parentsBoth[["summaries"]],
  n_true_parents = 3
)
aux
pval_minimal_top(aux, n_true = 3)

