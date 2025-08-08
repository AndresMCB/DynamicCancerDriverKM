# ============================================================================ #
# save_simulation_results.R
#
# Purpose:
#   Save all simulation outputs from the current R session into the Simulations/
#   folder structure for sharing (e.g., on GitHub). This includes:
#     - Per-configuration summaries (CSV only)
#     - A combined master summary table
#     - Full result objects (.rds, xz-compressed)
#     - A manifest file of saved objects
#     - Session info for reproducibility
#
# Usage:
#   1. Run your experiment script to create result_* objects in memory.
#   2. Run: source("Simulations/save_simulation_results.R")
# ============================================================================ #

# ---- folder setup ------------------------------------------------------------
dirs <- list(
  root      = "Simulations",
  summaries = file.path("Simulations", "results", "summaries"),
  objects   = file.path("Simulations", "results", "objects")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

suppressPackageStartupMessages({
  library(data.table)   # fwrite, rbindlist
})

# ---- helper functions --------------------------------------------------------
write_summary <- function(df, path_no_ext) {
  # Save as CSV only (no extra dependencies)
  data.table::fwrite(df, paste0(path_no_ext, ".csv"))
}

add_cfg <- function(summaries, cfg) {
  data.table::as.data.table(summaries)[, config := cfg]
}

save_rds <- function(obj, path) {
  saveRDS(obj, path, compress = "xz")
}

# ---- check required result objects -------------------------------------------
required_objs <- c(
  "result_2parentsUp","result_2parentsDown","result_2parentsBoth",
  "result_3parentsUp","result_3parentsDown","result_3parentsBoth"
)
missing <- setdiff(required_objs, ls(envir = .GlobalEnv))
if (length(missing)) {
  stop("Missing simulation result objects in the environment: ",
       paste(missing, collapse = ", "))
}

# ---- save per-configuration summaries ----------------------------------------
write_summary(result_2parentsUp$summaries,   file.path(dirs$summaries, "summaries_2parents_up"))
write_summary(result_2parentsDown$summaries, file.path(dirs$summaries, "summaries_2parents_down"))
write_summary(result_2parentsBoth$summaries, file.path(dirs$summaries, "summaries_2parents_both"))
write_summary(result_3parentsUp$summaries,   file.path(dirs$summaries, "summaries_3parents_up"))
write_summary(result_3parentsDown$summaries, file.path(dirs$summaries, "summaries_3parents_down"))
write_summary(result_3parentsBoth$summaries, file.path(dirs$summaries, "summaries_3parents_both"))

# ---- save combined master summary --------------------------------------------
master <- data.table::rbindlist(list(
  add_cfg(result_2parentsUp$summaries,   "2p_up"),
  add_cfg(result_2parentsDown$summaries, "2p_down"),
  add_cfg(result_2parentsBoth$summaries, "2p_both"),
  add_cfg(result_3parentsUp$summaries,   "3p_up"),
  add_cfg(result_3parentsDown$summaries, "3p_down"),
  add_cfg(result_3parentsBoth$summaries, "3p_both")
), use.names = TRUE, fill = TRUE)

write_summary(master, file.path(dirs$summaries, "summaries_master"))

# ---- save full simulation objects --------------------------------------------
save_rds(result_2parentsUp,   file.path(dirs$objects, "result_2parentsUp.rds"))
save_rds(result_2parentsDown, file.path(dirs$objects, "result_2parentsDown.rds"))
save_rds(result_2parentsBoth, file.path(dirs$objects, "result_2parentsBoth.rds"))
save_rds(result_3parentsUp,   file.path(dirs$objects, "result_3parentsUp.rds"))
save_rds(result_3parentsDown, file.path(dirs$objects, "result_3parentsDown.rds"))
save_rds(result_3parentsBoth, file.path(dirs$objects, "result_3parentsBoth.rds"))

# ---- save manifest + session info --------------------------------------------
manifest <- data.table::data.table(
  file = c(
    "result_2parentsUp.rds","result_2parentsDown.rds","result_2parentsBoth.rds",
    "result_3parentsUp.rds","result_3parentsDown.rds","result_3parentsBoth.rds"
  ),
  config = c("2p_up","2p_down","2p_both","3p_up","3p_down","3p_both"),
  n_runs = c(
    nrow(result_2parentsUp$summaries),
    nrow(result_2parentsDown$summaries),
    nrow(result_2parentsBoth$summaries),
    nrow(result_3parentsUp$summaries),
    nrow(result_3parentsDown$summaries),
    nrow(result_3parentsBoth$summaries)
  )
)
data.table::fwrite(manifest, file.path(dirs$objects, "manifest.csv"))

writeLines(capture.output(sessionInfo()),
           file.path(dirs$root, "sessionInfo.txt"))

message("✓ Simulation results saved under `Simulations/results/`.")

