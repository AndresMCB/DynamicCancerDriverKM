####################### simulate_gene_expression ###############################
#' Simulate multigene expression over pseudotime with DEGs
#'
#' Generates a pseudotime vector and builds expression for a set of genes,
#' where nominated DEGs follow sigmoid trajectories (up/down/mixed) and
#' non-DEGs follow simple background processes. Returns a data.frame sorted
#' by pseudotime with columns: Pseudotime, G1..Gn.
#'
#' @param n_genes Integer; total number of genes to simulate.
#' @param DEG_genes Integer vector; indices (1-based) of DEGs among genes.
#' @param n_timepoints Integer; number of pseudotime points.
#' @param amplitude_range Numeric length-2; min/max amplitude for DEGs.
#' @param regulation_type "up", "down", or "both" for DEG directions.
#' @return data.frame with pseudotime and simulated expression.
simulate_gene_expression <- function(
    n_genes = 20,
    DEG_genes = c(1),
    n_timepoints = 300,
    amplitude_range = c(10, 1000),
    regulation_type = "up"
){
  # Soft-depend on dplyr (only if needed). This is fine for scripts, but
  # avoid install-from-function in packages.
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)

  # Pseudotime ~ U(0,1), centred to mean 0 for nicer sigmoid midpoints
  pseudotime <- runif(n_timepoints, min = 0, max = 1)
  pseudotime <- pseudotime - mean(pseudotime)

  # Number of DEGs
  k <- length(DEG_genes)

  # Assign directions per DEG according to requested mode
  directions <- switch(
    regulation_type,
    "up"   = rep("up", k),
    "down" = rep("down", k),
    "both" = sample(c("up", "down"), k, replace = TRUE),
    stop("Invalid regulation_type")
  )

  # Sigmoid curve generator for DEGs with optional down-flip and noise
  sigmoid <- function(t, midpoint, slope, baseline = 0,
                      amplitude = 100, noise_prop = 0.05, direction = "up") {
    raw <- 1 / (1 + exp(-slope * (t - midpoint)))   # logistic
    expr <- if (direction == "up") raw else 1 - raw # flip for down-reg
    expr <- baseline + amplitude * expr             # scale + shift
    noise_sd <- amplitude * noise_prop              # noise ∝ amplitude
    pmax(rnorm(length(expr), mean = expr, sd = noise_sd), 0) # non-negative
  }

  # Simple background models for non-DEGs (varied but light-touch)
  random_expression <- function(t) {
    type <- sample(c("flat", "normal", "lognormal"), 1)

    if (type == "flat") {
      base <- runif(1, 50, 300)
      noise_sd <- runif(1, 2, 8)
      return(pmax(rep(base, length(t)) + rnorm(length(t), 0, noise_sd), 0))

    } else if (type == "normal") {
      mean_expr <- runif(1, 50, 300)
      sd_expr <- runif(1, 5, 20)
      return(pmax(rnorm(length(t), mean = mean_expr, sd = sd_expr), 0))

    } else {
      # Lightly varying log-normal baseline
      return(pmax(rlnorm(length(t), meanlog = log(80), sdlog = 0.3), 0))
    }
  }

  # Container: rows = timepoints, cols = genes
  gene_matrix <- matrix(NA, nrow = n_timepoints, ncol = n_genes)

  # DEG inflection points near centre to avoid edge artefacts
  midpoints <- runif(length(DEG_genes), -0.4, 0.4)
  names(midpoints) <- DEG_genes

  # Build each gene’s profile
  for (i in 1:n_genes) {
    if (i %in% DEG_genes) {
      # Parametrise DEG sigmoid
      idx <- which(DEG_genes == i)
      direction <- directions[idx]
      amp <- runif(1, amplitude_range[1], amplitude_range[2])
      midpoint <- midpoints[i]
      slope <- runif(1, 5, 50)               # steeper → sharper switch
      baseline <- runif(1, 10, 1000)

      gene_matrix[, i] <- sigmoid(
        pseudotime, midpoint, slope,
        amplitude = amp, baseline = baseline, direction = direction
      )
    } else {
      # Background expression for non-DEGs
      gene_matrix[, i] <- random_expression(pseudotime)
    }
  }

  # Label columns and return sorted by pseudotime (handy for plotting)
  colnames(gene_matrix) <- paste0("G", 1:n_genes)
  df <- data.frame(Pseudotime = pseudotime, gene_matrix)
  df <- df[order(df$Pseudotime), ]

  return(df)
}

################## simulate_target_from_expression #############################
#' Simulate a target (y) as a noisy linear mix of parent genes
#'
#' @param parent_genes Character vector; column names of parent genes in `expression_matrix`.
#' @param expression_matrix data.frame/matrix with gene columns.
#' @param noise_level Numeric in [0,1]; noise SD as a proportion of the signal SD.
#' @return list with $result_df (original + y), and $param (parents/weights).
simulate_target_from_expression <- function(parent_genes,
                                            expression_matrix,
                                            noise_level  = 0.05) {
  outcome <- list()

  # Basic input checks
  if (!all(parent_genes %in% colnames(expression_matrix))) {
    stop("Not all parent_genes are present in the expression matrix columns.")
  }
  if (!is.numeric(noise_level) || noise_level < 0 || noise_level > 1) {
    # (Note: error text said 'noise_proportion' in original; corrected here.)
    stop("noise_level must be a number between 0 and 1.")
  }

  # Pull parent expression as matrix
  parents_matrix <- as.matrix(expression_matrix[, parent_genes, drop = FALSE])

  # Random signed weights with magnitudes in [0.2, 1]
  magnitude <- runif(ncol(parents_matrix), min = 0.2, max = 1)
  weights <- magnitude * sample(c(-1, 1), ncol(parents_matrix), replace = TRUE)

  # Linear signal + Gaussian noise
  signal <- as.vector(parents_matrix %*% weights)
  signal_sd <- sd(signal)
  noise_sd <- noise_level * signal_sd
  noise <- rnorm(length(signal), mean = 0, sd = noise_sd)
  target <- signal + noise

  # Shift to be non-negative (convenient for “expression-like” targets)
  target <- target + abs(min(target))

  # Assemble output
  outcome$result_df <- as.data.frame(expression_matrix)
  outcome$result_df$y <- target
  outcome$param$GE <- as.data.frame(cbind(parents_matrix, target))
  outcome$param$weights <- weights

  return(outcome)
}

############################ evaluate_inferred_drivers ########################
#' Compare inferred drivers to ground-truth parents
#'
#' Expects a `result` from run_dcdkm_test() containing $geneScore and
#' $InferredDrivers. Summarises how many true parents appear in the top-N and
#' an “extended” set (ties included) and reports parent scores.
#'
#' @param result List; output from run_dcdkm_test().
#' @param true_parents Character vector; ground-truth parent gene IDs.
#' @return data.frame with detection counts, thresholds, and score strings.
evaluate_inferred_drivers <- function(result, true_parents) {

  # Structure checks (fail fast with helpful errors)
  if (!all(c("geneScore", "InferredDrivers") %in% names(result))) {
    stop("The result list must contain 'geneScore' and 'InferredDrivers'.")
  }

  # Extract
  gene_features <- result$geneScore$features
  inferred_df <- result$InferredDrivers

  # Exclude the target if present
  inferred_df <- subset(inferred_df, features != "y")

  # Partition truth into found/missing in the feature set
  found_parents <- true_parents[true_parents %in% gene_features]
  missing_parents <- setdiff(true_parents, gene_features)

  if (length(found_parents) == 0) {
    stop("None of the true parents were found in geneScore.")
  }
  if (length(missing_parents) > 0) {
    warning(paste("These true parents are missing from geneScore:",
                  paste(missing_parents, collapse = ", ")))
  }

  # Rank inferred drivers by decreasing score
  inferred_df <- inferred_df[order(-inferred_df$score), ]

  # Top-N = number of true parents
  n_parents <- length(true_parents)
  if (n_parents > nrow(inferred_df)) {
    stop("Number of true parents exceeds number of inferred drivers.")
  }

  # Threshold for top-N selection (handles ties at the boundary)
  top_cut_score <- inferred_df$score[n_parents]
  top_set <- inferred_df$features[inferred_df$score >= top_cut_score]

  # Extended set threshold: include everything down to the highest score
  # strictly below the top cut (effectively top set + one tie block below)
  lower_scores <- inferred_df$score[inferred_df$score < top_cut_score]
  extended_cut_score <- if (length(lower_scores) == 0) {
    min(inferred_df$score)
  } else {
    max(lower_scores)
  }
  extended_set <- inferred_df$features[inferred_df$score >= extended_cut_score]

  # Tally overlaps
  detected_top <- intersect(true_parents, top_set)
  detected_extended <- intersect(true_parents, extended_set)
  not_in_inferred <- setdiff(true_parents, inferred_df$features)

  # Collect per-parent scores (NA if absent)
  parent_scores <- sapply(true_parents, function(p) {
    score <- inferred_df$score[inferred_df$features == p]
    if (length(score) == 0) NA else score
  })

  # Human-readable strings
  score_string <- paste(paste(names(parent_scores),
                              round(parent_scores, 4), sep = "="),
                        collapse = "; ")
  top_hits_string <- paste(detected_top, collapse = ", ")
  extended_hits_string <- paste(detected_extended, collapse = ", ")

  # One-row summary
  summary_df <- data.frame(
    n_true_parents = length(true_parents),
    n_detected_in_top = length(detected_top),
    n_detected_in_extended = length(detected_extended),
    top_cut_score = top_cut_score,
    extended_cut_score = extended_cut_score,
    n_not_in_inferred = length(not_in_inferred),
    detected_top = top_hits_string,
    detected_extended = extended_hits_string,
    parent_scores = score_string,
    stringsAsFactors = FALSE
  )

  return(summary_df)
}

########################plot_gene_simulation##################################
#' Quick plot of selected simulated genes over pseudotime
#'
#' @param df data.frame from simulate_gene_expression().
#' @param genes Character vector of gene column names to plot.
#' @param title Title string.
#' @return ggplot object.
plot_gene_simulation <- function(df, genes,
                                 title = "Simulated Gene Expression") {
  # Check requested genes exist
  missing <- setdiff(genes, colnames(df))
  if (length(missing) > 0) {
    stop(paste("These genes are missing in the data frame:",
               paste(missing, collapse = ", ")))
  }

  # Long format for ggplot
  plot_data <- df[, c("Pseudotime", genes)]
  plot_data <- tidyr::pivot_longer(plot_data, -Pseudotime,
                                   names_to = "Gene", values_to = "Expression")

  # Line plot per gene across pseudotime
  ggplot(plot_data, aes(x = Pseudotime, y = Expression, colour = Gene)) +
    geom_line(size = 1) +
    labs(title = title, y = "Expression", x = "Pseudotime") +
    theme_minimal()
}

############################# run_dcdkm_test ############################
#' Run Dynamic Cancer Driver KM (DCDKM) on simulated data
#'
#' Bins pseudotime, enumerates/score models, and derives driver scores for
#' predictors relative to a target gene. Returns ranked models and driver tables.
#'
#' @param df data.frame; first column is pseudotime, others are features incl. target.
#' @param target Character; name of target column (default "y").
#' @param covariate Optional; feature used for time binning. Defaults to first feature.
#' @param nBins Integer; number of time bins.
#' @param parallel Logical; enable parallel scoring if supported.
#' @param verbose Logical; chatty progress.
#' @return list with topModels, geneScore, InferredDrivers, formulas, modelScore.
run_dcdkm_test <- function(df, target = "y", covariate = NULL,
                           nBins = 50, parallel = TRUE, verbose = TRUE) {
  # --- Dependency setup -------------------------------------------------------
  # Install on-the-fly if missing (OK for scripts; avoid inside packages/CI).
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
  if (!requireNamespace("AMCBGeneUtils", quietly = TRUE))
    devtools::install_github("AndresMCB/AMCBGeneUtils")
  if (!requireNamespace("DynamicCancerDriverKM", quietly = TRUE))
    devtools::install_github("AndresMCB/DynamicCancerDriverKM")

  library(DynamicCancerDriverKM)
  library(tidyverse)
  library(tictoc)

  # --- Features/covariate -----------------------------------------------------
  Features <- colnames(df)[-1]            # assume first col is pseudotime
  predictors <- setdiff(Features, target)  # everything except target

  if (is.null(covariate)) {
    covariate <- Features[1]               # default: first gene as covariate
  }

  # --- Time binning -----------------------------------------------------------
  binned <- DCDKM.BinTime(
    Mat1 = df[Features],
    covariate = covariate,
    Features = Features,
    nBins = nBins,
    pTime = df[, 1] # pseudotime vector
  )

  # --- Model scoring ----------------------------------------------------------
  results <- vector(mode = "list", length = 0)
  if (verbose) message(paste("Running DCDKM for target gene:", target))

  tic()
  features <- c(target, predictors)   # features vector starts with target
  k <- length(features)
  testModels <- DCDKM.GetModels(k)    # grouped model sets

  invariantScore <- list(score = NULL, formulas = NULL)

  # Score each group and accumulate
  for (j in seq_along(testModels)) {
    aux <- DCDKM.modelScoring(
      models = testModels[[j]],
      binned = binned,
      parallel = parallel,
      features = features,
      targetIndex = 1,   # target is first in 'features'
      num.folds = 2
    )
    invariantScore$score <- c(invariantScore$score, aux$modelScores)
    invariantScore$formulas <- c(invariantScore$formulas, aux$formulas)
  }

  # --- Rank/select top models and derive driver scores ------------------------
  index <- order(invariantScore$score)          # lower is better
  models <- unlist(testModels, recursive = FALSE)
  topModels <- models[index[1:k]]               # keep top k
  geneScore <- driverScore(topModels, features) # per-gene contribution
  InferredDrivers <- geneScore %>% dplyr::filter(score > 0)

  # --- Pack results -----------------------------------------------------------
  results$topModels <- topModels
  results$geneScore <- geneScore
  results$InferredDrivers <- InferredDrivers
  results$formulas <- invariantScore$formulas[index[1:k]]
  results$modelScore <- invariantScore$score[index[1:k]]
  toc()

  return(results)
}

#########################calculate_detection_rates##############################
#' Count detection frequencies across simulation runs
#'
#' For each possible detection count (n_true, n_true-1, ..., 0), tally how many
#' runs achieved exactly that count in the top and extended sets.
#'
#' @param df data.frame from evaluate_inferred_drivers() row-binds.
#' @param n_true_parents Integer; number of true parents (2 or 3 in your study).
#' @return data.frame with rows per detection level and counts per set.
calculate_detection_rates <- function(df, n_true_parents = 2) {
  # Schema check
  required_cols <- c("n_true_parents", "n_detected_in_top", "n_detected_in_extended")
  if (!all(required_cols %in% colnames(df))) {
    stop("Dataframe must contain columns: n_true_parents, n_detected_in_top, n_detected_in_extended")
  }

  # Tally exact counts for top and extended sets
  top_stats <- NULL
  ext_stats <- NULL

  for (i in n_true_parents:0) {
    # Exactly i detections in the top set
    aux <- sum(df$n_detected_in_top == i)
    top_stats <- c(top_stats, aux)

    # Exactly i detections in the extended set
    aux <- sum(df$n_detected_in_extended == i)
    ext_stats <- c(ext_stats, aux)
  }

  # Assemble result; values are counts (you can convert to proportions later)
  result <- data.frame(
    metric = paste0("detected ", n_true_parents:0),
    Top_Percentage = top_stats,
    Extended_Percentage = ext_stats
  )

  return(result)
}

########################pval_minimal_top#######################################
#' Binomial tail p-value for exact full recovery in the top set
#'
#' Tests whether the observed number of “perfect” recoveries (all true parents
#' in the minimal top set) exceeds what you'd expect from random ranking.
#'
#' @param res_df data.frame from calculate_detection_rates(); first row is “detected n”.
#' @param n_true Integer; number of true drivers (2 or 3).
#' @param G Integer; total candidate genes (default 20).
#' @param N Integer; total runs (default 500).
#' @return One-sided binomial upper-tail p-value, P[X >= k], X~Bin(N, p0).
pval_minimal_top <- function(res_df, n_true, G = 20, N = 500) {
  k  <- res_df$Top_Percentage[1]     # observed perfect recoveries
  p0 <- 1 / choose(G, n_true)        # null success prob: random top-n match
  pbinom(k - 1, size = N, prob = p0, # P[X >= k]
         lower.tail = FALSE)
}
