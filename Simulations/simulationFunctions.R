####################### simulate_gene_expression ###############################
# This function simulates gene expression data over pseudotime, including 
# a specified number of differentially expressed genes (DEGs) with either 
# up- or down-regulation using sigmoid curves, and background noise for 
# non-DEGs.

simulate_gene_expression <- function(
    n_genes = 20,                   # Total number of genes to simulate
    DEG_genes = c(1),               # Indices of DEGs among the total genes
    n_timepoints = 300,             # Number of pseudotime points to simulate
    amplitude_range = c(10, 1000),  # Range of expression amplitudes
    regulation_type = "up"          # Direction of DEGs: "up", "down", or "both"
){
  # Load 'dplyr' package if not already installed
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  
  # Generate pseudotime values from uniform distribution [0, 1], then centre
  pseudotime <- runif(n_timepoints, min = 0, max = 1)
  pseudotime <- pseudotime - mean(pseudotime)
  
  # Number of DEGs
  k <- length(DEG_genes)
  
  # Assign regulation directions for each DEG
  directions <- switch(
    regulation_type,
    "up" = rep("up", k),
    "down" = rep("down", k),
    "both" = sample(c("up", "down"), k, replace = TRUE),
    stop("Invalid regulation_type")  # Throw error if invalid option
  )
  
  # Define sigmoid function for modelling DEG expression over pseudotime
  sigmoid <- function(t, midpoint, slope, baseline = 0, 
                      amplitude = 100, noise_prop = 0.05, direction = "up") {
    # Standard sigmoid curve
    raw <- 1 / (1 + exp(-slope * (t - midpoint)))
    
    # Flip for downregulated genes
    expr <- if (direction == "up") raw else 1 - raw
    
    # Scale and shift with baseline and amplitude
    expr <- baseline + amplitude * expr
    
    # Add noise (Gaussian) proportional to amplitude
    noise_sd <- amplitude * noise_prop
    # Ensure non-negative values
    pmax(rnorm(length(expr), mean = expr, sd = noise_sd), 0) 
  }
  
  # Define random expression model for non-DEGs
  random_expression <- function(t) {
    # Random background type
    type <- sample(c("flat", "normal", "lognormal"), 1)  
    
    if (type == "flat") {
      # Constant baseline with small noise
      base <- runif(1, 50, 300)
      noise_sd <- runif(1, 2, 8)
      return(pmax(rep(base, length(t)) + rnorm(length(t), 0, noise_sd), 0))
      
    } else if (type == "normal") {
      # Normally distributed expression with lower SD
      mean_expr <- runif(1, 50, 300)
      sd_expr <- runif(1, 5, 20)
      return(pmax(rnorm(length(t), mean = mean_expr, sd = sd_expr), 0))
      
    } else {
      # Log-normal expression, mild variance
      return(pmax(rlnorm(length(t), meanlog = log(80), sdlog = 0.3), 0))
    }
  }
  
  # Initialise matrix to hold gene expression data
  gene_matrix <- matrix(NA, nrow = n_timepoints, ncol = n_genes)
  
  # Assign midpoints for sigmoid curves (DEG inflection points)
  midpoints <- runif(length(DEG_genes), -0.4, 0.4)
  names(midpoints) <- DEG_genes
  
  # Loop over all genes to simulate expression profiles
  for (i in 1:n_genes) {
    if (i %in% DEG_genes) {
      # For DEGs: use sigmoid function with randomised parameters
      idx <- which(DEG_genes == i)
      direction <- directions[idx]
      amp <- runif(1, amplitude_range[1], amplitude_range[2])
      midpoint <- midpoints[i]
      slope <- runif(1, 5, 50)
      baseline <- runif(1, 10, 1000)
      
      # Generate DEG expression
      gene_matrix[, i] <- sigmoid(pseudotime, midpoint, slope,
                                  amplitude = amp, baseline = baseline, 
                                  direction = direction)
    } else {
      # For non-DEGs, simulate background expression
      gene_matrix[, i] <- random_expression(pseudotime)
    }
  }
  
  # Convert expression matrix to data frame with pseudotim
  colnames(gene_matrix) <- paste0("G", 1:n_genes)
  df <- data.frame(Pseudotime = pseudotime, gene_matrix)
  
  # Sort the dataset by pseudotime for smoother visualisation
  df <- df[order(df$Pseudotime), ]
  
  return(df)
}

################## simulate_target_from_expression #############################
# This function simulates a target gene (or response variable) as a noisy 
# linear combination of selected parent genes from an existing expression matrix.

simulate_target_from_expression <- function(parent_genes, 
                                            expression_matrix, 
                                            noise_level  = 0.05
                                            ) {
  outcome <- list()
  
  # Validate input
  if (!all(parent_genes %in% colnames(expression_matrix))) {
    stop("Not all parent_genes are present in the expression matrix columns.")
  }
  if (!is.numeric(noise_level) || noise_level < 0 || noise_level > 1) {
    stop("noise_proportion must be a number between 0 and 1.")
  }

  
  # Extract parent gene expression data
  parents_matrix <- as.matrix(expression_matrix[, parent_genes, drop = FALSE])
  
  # Generate signal via weighted linear combination of parent genes
  magnitude <- runif(ncol(parents_matrix), min = 0.2, max = 1)
  weights <- magnitude * sample(c(-1, 1), ncol(parents_matrix), replace = TRUE)
  signal <- as.vector(parents_matrix %*% weights)
  
  # Add Gaussian noise proportional to the signal's SD
  signal_sd <- sd(signal)
  noise_sd <- noise_level*signal_sd
  noise <- rnorm(length(signal), mean = 0, sd = noise_sd)
  target <- signal + noise 
  
  # Ensure non-negative expression values
  target <- target + abs(min(target))
  
  # Store results
  outcome$result_df <- as.data.frame(expression_matrix)
  outcome$result_df$y <- target
  outcome$param$GE <- as.data.frame(cbind(parents_matrix,target))
  outcome$param$weights <- weights
  
  return(outcome)
}


############################ evaluate_inferred_drivers ########################
# This function evaluates how well inferred gene drivers match a 
# set of true parent genes.
# It returns a summary data frame showing how many true parents were recovered 
# in the top inferred features, and includes scores and ranks for comparison.

evaluate_inferred_drivers <- function(result, true_parents) {
  
  # Check that required fields exist in the result object
    if (!all(c("geneScore", "InferredDrivers") %in% names(result))) {
    stop("The result list must contain 'geneScore' and 'InferredDrivers'.")
  }
  
  # Extract gene features and inferred driver table
  gene_features <- result$geneScore$features
  inferred_df <- result$InferredDrivers
  
  # Remove 'y' from the inferred drivers if present
  inferred_df <- subset(inferred_df, features != "y")
  
  # Identify which true parents were found in the feature list
  found_parents <- true_parents[true_parents %in% gene_features]
  missing_parents <- setdiff(true_parents, gene_features)
  
  # Stop if no true parents were found in the feature list
  if (length(found_parents) == 0) {
    stop("None of the true parents were found in geneScore.")
  }
  
  # Issue a warning if any true parents were not found
  if (length(missing_parents) > 0) {
    warning(paste("These true parents are missing from geneScore:",
                  paste(missing_parents, collapse = ", ")))
  }
  
  # Sort inferred drivers by decreasing score (higher score = more likely driver)
  inferred_df <- inferred_df[order(-inferred_df$score), ]
  
  # Determine number of top candidates to compare based on number of true parents
  n_parents <- length(true_parents)
  if (n_parents > nrow(inferred_df)) {
    stop("Number of true parents exceeds number of inferred drivers.")
  }
  
  # Identify the score threshold for top-n selection
  top_cut_score <- inferred_df$score[n_parents]
  
  # Identify the score threshold for top-n selection
  top_set <- inferred_df$features[inferred_df$score >= top_cut_score]
  
  # Extended set: features with scores >= the max score *below* the cut-off
  lower_scores <- inferred_df$score[inferred_df$score < top_cut_score]
  extended_cut_score <- if (length(lower_scores) == 0) {
    min(inferred_df$score)
  } else {
    max(lower_scores)
  }
  
  extended_set <- inferred_df$features[inferred_df$score >= extended_cut_score]
  
  # Evaluate recovery: how many true parents are in the top/extended sets
  detected_top <- intersect(true_parents, top_set)
  detected_extended <- intersect(true_parents, extended_set)
  not_in_inferred <- setdiff(true_parents, inferred_df$features)
  
  # Extract scores for each true parent (or NA if missing)
  parent_scores <- sapply(true_parents, function(p) {
    score <- inferred_df$score[inferred_df$features == p]
    if (length(score) == 0) NA else score
  })
  
  # Format output strings for clarity
  score_string <- paste(paste(names(parent_scores), 
                              round(parent_scores, 4), sep = "="), 
                        collapse = "; ")
  top_hits_string <- paste(detected_top, collapse = ", ")
  extended_hits_string <- paste(detected_extended, collapse = ", ")
  
  #Create and return a summary data frame with evaluation metrics
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
plot_gene_simulation <- function(df, genes, 
                                 title = "Simulated Gene Expression") {
  # Check gene names exist
  missing <- setdiff(genes, colnames(df))
  if (length(missing) > 0) {
    stop(paste("These genes are missing in the data frame:", 
               paste(missing, collapse = ", ")))
  }
  
  # Prepare long-format data for ggplot
  plot_data <- df[, c("Pseudotime", genes)]
  plot_data <- tidyr::pivot_longer(plot_data, -Pseudotime, 
                                   names_to = "Gene", values_to = "Expression")
  
  # Plot
  ggplot(plot_data, aes(x = Pseudotime, y = Expression, colour = Gene)) +
    geom_line(size = 1) +
    labs(title = title, y = "Expression", x = "Pseudotime") +
    theme_minimal()
}

############################# run_dcdkm_test ############################
# This function runs the Dynamic Cancer Driver KM (DCDKM) test on a our simulated
# gene expression dataset over pseudotime. It identifies dynamic drivers  
# of a target gene using model-based scoring with time-binning.

run_dcdkm_test <- function(df, target = "y", covariate = NULL, 
                           nBins = 50, parallel = TRUE, verbose = TRUE) {
  # --- Dependency checks and installation ---
  # Ensure required packages are installed
  
  if (!requireNamespace("devtools", quietly = TRUE)) 
    install.packages("devtools")
  if (!requireNamespace("tidyverse", quietly = TRUE)) 
    install.packages("tidyverse")
  if (!requireNamespace("AMCBGeneUtils", quietly = TRUE)) 
    devtools::install_github("AndresMCB/AMCBGeneUtils")
  if (!requireNamespace("DynamicCancerDriverKM", quietly = TRUE)) 
    devtools::install_github("AndresMCB/DynamicCancerDriverKM")
  
  # Load required packages
  library(DynamicCancerDriverKM)
  library(tidyverse)
  library(tictoc)
  
  # --- Feature preparation ---
  Features <- colnames(df)[-1] # Exclude pseudotime (assumed first column)
  predictors <- setdiff(Features, target) # All features except target
  
  
  # Set covariate for time binning; default to first feature if not supplied
  if (is.null(covariate)) {
    covariate <- Features[1]  # use first gene as covariate if none given
  }
  
  # --- Bin pseudotime using DCDKM utility ---
  binned <- DCDKM.BinTime(
    Mat1 = df[Features],
    covariate = covariate,
    Features = Features,
    nBins = nBins,
    pTime = df[, 1] # pseudotime vector
  )
  
  # Initialise results container
  results <- vector(mode = "list", length = 0)
  
  if (verbose) {
    message(paste("Running DCDKM for target gene:", target))
  }
  # --- Main DCDKM model scoring ---
  tic() # Start timing
  features <- c(target, predictors)  # Target + all predictors
  k <- length(features)  # Total number of features
  testModels <- DCDKM.GetModels(k)  # Get all possible models
  
  # Prepare object to store model scores and formulas
  invariantScore <- list(score = NULL, formulas = NULL)
  
  # Evaluate each group of models returned by DCDKM.GetModels
  for (j in seq_along(testModels)) {
    aux <- DCDKM.modelScoring(
      models = testModels[[j]],
      binned = binned,
      parallel = parallel,
      features = features,
      targetIndex = 1,        # 'target' is the first in features vector
      num.folds = 2           # Cross-validation folds
    )
    # Store scores and formulas
    invariantScore$score <- c(invariantScore$score, aux$modelScores)
    invariantScore$formulas <- c(invariantScore$formulas, aux$formulas)
  }
  
  # --- Rank and extract top models ---
  index <- order(invariantScore$score)  # Rank by score (ascending)
  models <- unlist(testModels, recursive = FALSE)  # Flatten list of models
  topModels <- models[index[1:k]]  # Select top k models
  
  # Score each gene (predictor) by its contribution across top models
  geneScore <- driverScore(topModels, features)
  
  # Filter only predictors with non-zero score (likely drivers)
  InferredDrivers <- geneScore %>% filter(score > 0)
  
  # --- Store and return results ---
  results$topModels <- topModels
  results$geneScore <- geneScore
  results$InferredDrivers <- InferredDrivers
  results$formulas <- invariantScore$formulas[index[1:k]]
  results$modelScore <- invariantScore$score[index[1:k]]
  toc()  # End timing
  
  return(results)
}



#########################calculate_detection_rates##############################
# 
# This function calculates how often a certain number of true parent genes
# were detected in the 'top' and 'extended' inferred driver sets, based on
# summary data (e.g., from evaluate_inferred_drivers output).

calculate_detection_rates <- function(df, n_true_parents = 2) {
  # Ensure the input dataframe has required columns
  required_cols <- c("n_true_parents", "n_detected_in_top", "n_detected_in_extended")
  if (!all(required_cols %in% colnames(df))) {
    stop("Dataframe must contain columns: n_true_parents, n_detected_in_top, n_detected_in_extended")
  }
  
  # Initialise counters for each detection level (from n to 0)
  top_stats <- NULL
  ext_stats <- NULL
  
  for (i in n_true_parents:0) {
    # Count how many rows had exactly 'i' true parents detected in top set
    aux <- sum(df$n_detected_in_top == i)
    top_stats <- c(top_stats, aux)
    
    # Count how many rows had exactly 'i' true parents detected in extended set
    aux <- sum(df$n_detected_in_extended == i)
    ext_stats <- c(ext_stats, aux)
  }
  
  # Create result table: one row per detection count (n to 0)
  result <- data.frame(
    metric = paste0("detected ", n_true_parents:0),
    Top_Percentage = top_stats,
    Extended_Percentage = ext_stats
  )
  
  return(result)
}

########################pval_minimal_top############################

# Compute a one-sided binomial tail p-value for "exact full recovery (minimal top)"
# Assumptions:
# - res_df comes from your calculate_detection_rates(...), where:
#     * row 1 corresponds to "detected n" (i.e., exact full recovery),
#     * column Top_Percentage holds the count of successes across N runs.
# - G is the total number of candidate genes (default 20).
# - n_true is the number of true drivers (2 or 3).
# - N is the number of simulation runs (default 500).
# - Under the null (random ranking wrt the true drivers), the per-run success
#   probability is p0 = 1 / choose(G, n_true). We then test P(X >= k) with X~Bin(N, p0).

pval_minimal_top <- function(res_df, n_true, G = 20, N = 500) {
  k  <- res_df$Top_Percentage[1]      # observed successes (exact full recovery) in N runs
  p0 <- 1 / choose(G, n_true)         # null per-run success prob: 1 / C(G, n_true)
  pbinom(k - 1, size = N, prob = p0,  # one-sided binomial upper tail
         lower.tail = FALSE)          # returns P[X >= k]
}
