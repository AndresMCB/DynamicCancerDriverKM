

###############################################################################
library(tidyverse)
library(data.table)  # for rbindlist()

set.seed(42)

result_2parentsUp <- run_simulation_pipeline(
  DEG_genes = c(1:2),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "up"
)


result_2parentsBoth <- run_simulation_pipeline(
  DEG_genes = c(1:2),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "both"
)

result_2parentsDown <- run_simulation_pipeline(
  DEG_genes = c(1:2),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "down"
)


result_3parentsUp <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "up"
)


result_3parentsBoth <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "both"
)

result_3parentsDown <- run_simulation_pipeline(
  DEG_genes = c(1:3),
  n_genes = 20,
  n_timepoints = 200,
  amplitude_range = c(0, 1000),
  n_repetitions = 5,
  experiments_per_rep = 100,
  regulation_type = "down"
)

############2 parents results ############
# 2 parents UP
calculate_detection_rates(df =result_2parentsUp[["summaries"]]
                          ,n_true_parents = 2)
aux 
pval_minimal_top(aux, n_true = 2)
# 2 parents DOWN
calculate_detection_rates(df =result_2parentsDown[["summaries"]]
                          ,n_true_parents = 2)
aux 
pval_minimal_top(aux, n_true = 2)
# 2 parents MIXED
aux <- calculate_detection_rates(df =result_2parentsBoth[["summaries"]]
                          ,n_true_parents = 2)
aux 
pval_minimal_top(aux, n_true = 2)

############3 parents results ############
# 3 parents UP
aux <- calculate_detection_rates(df =result_3parentsUp[["summaries"]]
                          ,n_true_parents = 3)
aux 
pval_minimal_top(aux, n_true = 3)
# 3 parents DOWN
aux <- calculate_detection_rates(df =result_3parentsDown[["summaries"]]
                          ,n_true_parents = 3)
aux 
pval_minimal_top(aux, n_true = 3)
# 3 parents MIXED
aux <- calculate_detection_rates(df =result_3parentsBoth[["summaries"]]
                          ,n_true_parents = 3)
aux 
pval_minimal_top(aux, n_true = 3)



