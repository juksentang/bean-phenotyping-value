# ============================================================================
# Layer 1 — Master Orchestrator
# ============================================================================
# Run the full Layer 1 simulation and save results.

source(file.path(getwd(), "layer1_genetic_sim/R/05_monte_carlo.R"))

run_layer1 <- function(n_reps = NULL, output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "layer1_genetic_sim/outputs")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load configs
  master_cfg <- load_config("master_config")
  gen_cfg    <- load_config("genetic_params")
  breed_cfg  <- load_config("breeding_params")

  # Override n_reps if provided (useful for smoke tests)
  if (!is.null(n_reps)) master_cfg$mc_reps <- n_reps

  # Run Monte Carlo
  mc <- run_monte_carlo(master_cfg, gen_cfg, breed_cfg,
                         n_reps = master_cfg$mc_reps)

  # Save results
  out_path <- file.path(output_dir, "mc_results.csv")
  fwrite(mc$results, out_path)
  cat("Results saved to:", out_path, "\n")

  # Save summary
  summary_path <- file.path(output_dir, "mc_summary.csv")
  fwrite(mc$summary, summary_path)
  cat("Summary saved to:", summary_path, "\n")

  # Try arrow format if available
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow_path <- file.path(output_dir, "mc_results.arrow")
    arrow::write_feather(mc$results, arrow_path)
    cat("Arrow results saved to:", arrow_path, "\n")
  }

  mc
}

# Run if executed directly
if (sys.nframe() == 0) {
  # Check for command line args: first arg = n_reps
  args <- commandArgs(trailingOnly = TRUE)
  n_reps <- if (length(args) > 0) as.integer(args[1]) else NULL
  run_layer1(n_reps = n_reps)
}
