run_simulation_pipeline <- function(
    DEG_genes = c(1:2),
    n_genes = 20,
    n_timepoints = 200,
    amplitude_range = c(0, 1000),
    n_repetitions = 5,
    experiments_per_rep = 100,
    regulation_type = "up"
) {
  
  parent_genes <- paste0("G", DEG_genes)
  total_experiments <- n_repetitions * experiments_per_rep
  
  res_list <- vector("list", total_experiments)
  summary_list <- vector("list", total_experiments)
  
  for (rep in seq_len(n_repetitions)) {
    for (i in seq_len(experiments_per_rep)) {
      experiment_index <- (rep - 1) * experiments_per_rep + i
      
      # Simulate gene expression
      df <- simulate_gene_expression(
        n_genes = n_genes,
        DEG_genes = DEG_genes,
        n_timepoints = n_timepoints,
        amplitude_range = amplitude_range,
        regulation_type = regulation_type
      )
      
      # Simulate target gene
      targetSim <- simulate_target_from_expression(
        parent_genes = parent_genes,
        expression_matrix = df,
        #seed = NULL
      )
      
      df_with_target <- targetSim$result_df
      
      # Run causal inference
      res <- run_dcdkm_test(df_with_target)
      res[["targetSim"]] <- targetSim$param
      res_list[[experiment_index]] <- res
      
      # Evaluate performance
      summary <- evaluate_inferred_drivers(res, true_parents = parent_genes)
      summary_list[[experiment_index]] <- summary
    }
  }
  
  all_summaries <- data.table::rbindlist(summary_list)
  
  return(list(
    results = res_list,
    summaries = all_summaries
  ))
}
