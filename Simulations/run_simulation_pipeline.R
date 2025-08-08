#' Run a simulated causal-inference pipeline over many experiments
#'
#' Repeats a simulation where (1) gene expression is generated, (2) a target
#' gene is simulated from specified parent genes, (3) a causal inference test
#' is run, and (4) performance is evaluated. Results and per-experiment
#' summaries are returned.
#'
#' @param DEG_genes Integer vector of indices for the true driver/parent genes (1-based).
#' @param n_genes Total number of genes to simulate per experiment.
#' @param n_timepoints Number of time points per simulation.
#' @param amplitude_range Numeric length-2 vector giving the min/max amplitude of signals.
#' @param n_repetitions Number of outer repetitions (batches).
#' @param experiments_per_rep Number of experiments to run in each repetition.
#' @param regulation_type Character string for direction of regulation (e.g., "up" or "down").
#'
#' @return A list with:
#' \describe{
#'   \item{results}{List of per-experiment outputs from \code{run_dcdkm_test()}, augmented
#'   with \code{$targetSim} parameters used to generate the target.}
#'   \item{summaries}{data.table row-bind of per-experiment evaluation summaries.}
#' }
#'
#' @examples
#' set.seed(1)
#' out <- run_simulation_pipeline(DEG_genes = 1:2, n_genes = 20, n_timepoints = 200)
#' str(out$summaries)
run_simulation_pipeline <- function(
    DEG_genes = c(1:2),          # indices (1-based) of true driver genes
    n_genes = 20,                # total genes to simulate
    n_timepoints = 200,          # time points per experiment
    amplitude_range = c(0, 1000),# amplitude range for simulated expression
    n_repetitions = 5,           # how many batches to run
    experiments_per_rep = 100,   # experiments per batch
    regulation_type = "up"       # direction of regulation for target simulation
) {
  # Convert driver indices to gene IDs that will exist in the simulated data
  parent_genes <- paste0("G", DEG_genes)

  # Total number of experiments across all repetitions
  total_experiments <- n_repetitions * experiments_per_rep

  # Preallocate containers for speed and clarity
  res_list <- vector("list", total_experiments)      # stores full inference results
  summary_list <- vector("list", total_experiments)  # stores performance summaries

  # Outer loop over repetitions (batches)
  for (rep in seq_len(n_repetitions)) {

    # Inner loop over experiments within a repetition
    for (i in seq_len(experiments_per_rep)) {

      # Linear index into the preallocated lists
      experiment_index <- (rep - 1) * experiments_per_rep + i

      # --- 1) Simulate multigene expression time series ---
      df <- simulate_gene_expression(
        n_genes        = n_genes,
        DEG_genes      = DEG_genes,
        n_timepoints   = n_timepoints,
        amplitude_range= amplitude_range,
        regulation_type= regulation_type
      )

      # --- 2) Simulate target gene as a (causal) function of parent genes ---
      # targetSim is expected to include:
      #   $result_df: expression with target appended
      #   $param: parameters used (e.g., parents, weights, noise)
      targetSim <- simulate_target_from_expression(
        parent_genes      = parent_genes,
        expression_matrix = df
        )

      # Combine original expression with the simulated target
      df_with_target <- targetSim$result_df

      # --- 3) Run causal inference on the combined data ---
      res <- run_dcdkm_test(df_with_target)

      # Attach target simulation parameters for traceability/repro
      res[["targetSim"]] <- targetSim$param

      # Store full result for this experiment
      res_list[[experiment_index]] <- res

      # --- 4) Evaluate how well the inferred drivers matched the truth ---
      summary <- evaluate_inferred_drivers(res, true_parents = parent_genes)

      # Store summary metrics (e.g., precision/recall/AUROC—depends on implementation)
      summary_list[[experiment_index]] <- summary
    }
  }

  # Row-bind all per-experiment summaries into a single data.table
  all_summaries <- data.table::rbindlist(summary_list)

  # Return both raw results and aggregated summaries
  return(list(
    results   = res_list,
    summaries = all_summaries
  ))
}
