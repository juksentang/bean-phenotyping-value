# ============================================================================
# Layer 2 — Master Orchestrator
# ============================================================================
# Runs the full economic evaluation on Layer 1 results.

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/04_sensitivity.R"))
source(file.path(getwd(), "layer2_economics/R/05_visualize.R"))

run_layer2 <- function(mc_results = NULL, output_dir = NULL,
                        run_mc_sens = TRUE, n_mc_draws = 1000) {
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  econ_cfg   <- load_config("economic_params")
  master_cfg <- load_config("master_config")

  # Load Layer 1 results if not provided
  if (is.null(mc_results)) {
    mc_path <- file.path(getwd(), "layer1_genetic_sim/outputs/mc_results.csv")
    if (!file.exists(mc_path)) stop("Layer 1 results not found at: ", mc_path)
    mc_results <- fread(mc_path)
    cat("Loaded", nrow(mc_results), "rows from Layer 1 results\n")
  }

  cat("=== Layer 2: Economic Evaluation ===\n\n")

  # --- Step 1: Base case economic evaluation ---
  cat("Step 1: Base case evaluation...\n")
  econ <- evaluate_all_reps(mc_results, econ_cfg, master_cfg)

  cat("\n=== Economic Summary ===\n")
  cat(sprintf("Mean NPV (tool vs traditional): $%s\n",
              formatC(econ$summary$mean_npv, format = "f", big.mark = ",", digits = 0)))
  cat(sprintf("Median NPV: $%s\n",
              formatC(econ$summary$median_npv, format = "f", big.mark = ",", digits = 0)))
  cat(sprintf("95%% CI: [$%s, $%s]\n",
              formatC(econ$summary$q025_npv, format = "f", big.mark = ",", digits = 0),
              formatC(econ$summary$q975_npv, format = "f", big.mark = ",", digits = 0)))
  cat(sprintf("P(NPV > 0): %.1f%%\n", econ$summary$prob_positive * 100))
  cat(sprintf("Mean IRR: %.1f%%\n", econ$summary$mean_irr * 100))
  cat(sprintf("Mean incremental genetic gain: %.2f kg/ha/yr\n",
              econ$summary$mean_incr_dg))

  # Save economic results
  fwrite(econ$results, file.path(output_dir, "economic_results.csv"))
  cat("\nEconomic results saved.\n")

  # --- Step 2: Tornado sensitivity analysis ---
  cat("\nStep 2: Tornado sensitivity analysis...\n")
  tornado <- run_tornado(mc_results, econ_cfg, master_cfg)
  cat("\nTornado Results:\n")
  print(tornado$tornado[, .(param, lo_val, hi_val, npv_lo, npv_hi, npv_range)])

  fwrite(tornado$tornado, file.path(output_dir, "tornado_results.csv"))
  plot_tornado(tornado, file.path(output_dir, "tornado_diagram.png"))

  # --- Step 3: Visualizations ---
  cat("\nStep 3: Generating plots...\n")
  plot_gain_trajectories(mc_results, file.path(output_dir, "gain_comparison.png"))
  plot_adoption_curves(econ_cfg = econ_cfg,
                        output_path = file.path(output_dir, "adoption_curves.png"))

  # --- Step 4: MC sensitivity (optional, slower) ---
  if (run_mc_sens) {
    cat("\nStep 4: Monte Carlo sensitivity (", n_mc_draws, "draws)...\n")
    mc_sens <- run_mc_sensitivity(mc_results, n_draws = n_mc_draws,
                                   econ_cfg = econ_cfg, master_cfg = master_cfg)
    cat(sprintf("MC Sensitivity: Mean NPV = $%s, P(NPV>0) = %.1f%%\n",
                formatC(mc_sens$mean_npv, format = "f", big.mark = ",", digits = 0),
                mc_sens$prob_positive * 100))

    plot_npv_distribution(mc_sens$npv_distribution,
                           file.path(output_dir, "npv_distribution.png"))

    saveRDS(mc_sens, file.path(output_dir, "mc_sensitivity.rds"))
  }

  cat("\n=== Layer 2 Complete ===\n")
  cat("Outputs saved to:", output_dir, "\n")

  list(economic = econ, tornado = tornado,
       mc_sensitivity = if (run_mc_sens) mc_sens else NULL)
}

# --- New: full Layer 2 sweep pipeline consuming cell_summaries.csv ---

run_layer2_sweep <- function(cell_summaries_path = NULL,
                              output_dir = NULL,
                              sobol_N = NULL) {
  source(file.path(getwd(), "layer2_economics/R/07_sweep_economics.R"))
  source(file.path(getwd(), "layer2_economics/R/08_break_even.R"))
  source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))

  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("╔══════════════════════════════════════════════╗\n")
  cat("║         Layer 2 Sweep Pipeline              ║\n")
  cat("╚══════════════════════════════════════════════╝\n\n")

  cat("━━━ Step 1: Lookup table ━━━\n")
  lookup <- build_lookup_table(cell_summaries_path = cell_summaries_path,
                                output_path = file.path(output_dir,
                                                          "lookup_table.parquet"))

  cat("\n━━━ Step 2: Break-even surface ━━━\n")
  be <- extract_break_even(lookup = lookup, output_dir = output_dir)

  cat("\n━━━ Step 3: Global Sobol (8D) ━━━\n")
  sa <- run_sobol(N = sobol_N, output_dir = output_dir)

  cat("\n━━━ Step 4: Figures ━━━\n")
  # Heatmap with break-even contour (stressed slice where contour is non-empty)
  stressed_err_cost <- be$stressed_marginals[slice == "err_vs_cost"]
  stressed_contour_ec <- be$stressed_contours[slice == "err_vs_cost"]
  if (nrow(stressed_err_cost) > 0) {
    plot_response_heatmap(
      stressed_err_cost,
      x_name = "x_val", y_name = "y_val",
      x_label = "Error reduction", y_label = "Cost reduction",
      title = "NPV response: err × cost (stressed economics)",
      contour_dt = stressed_contour_ec,
      output_path = file.path(fig_dir, "heatmap_err_cost_stressed.png")
    )
  }

  # Baseline heatmap — all positive NPV, no contour
  base_err_cost <- be$marginals[slice == "err_vs_cost"]
  if (nrow(base_err_cost) > 0) {
    plot_response_heatmap(
      base_err_cost,
      x_name = "x_val", y_name = "y_val",
      x_label = "Error reduction", y_label = "Cost reduction",
      title = "NPV response: err × cost (baseline economics)",
      output_path = file.path(fig_dir, "heatmap_err_cost_baseline.png")
    )
  }

  base_cyc_budget <- be$marginals[slice == "cyc_vs_budget"]
  if (nrow(base_cyc_budget) > 0) {
    plot_response_heatmap(
      base_cyc_budget,
      x_name = "x_val", y_name = "y_val",
      x_label = "Cycle compression (years)", y_label = "Total budget (USD)",
      title = "NPV response: cycle compression × budget (baseline)",
      output_path = file.path(fig_dir, "heatmap_cyc_budget_baseline.png")
    )
  }

  plot_sobol_bars(sa, output_path = file.path(fig_dir, "sobol_indices.png"))

  plot_profit_surface(be$profit_frac,
                       output_path = file.path(fig_dir, "profit_surface.png"))

  cat("\n╔══════════════════════════════════════════════╗\n")
  cat("║    Layer 2 Sweep Complete                    ║\n")
  cat("╚══════════════════════════════════════════════╝\n")
  cat("Outputs in:", output_dir, "\n")

  invisible(list(lookup = lookup, break_even = be, sobol = sa))
}


# Run if executed directly
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0 && args[1] == "--sweep") {
    run_layer2_sweep()
  } else {
    run_layer2(run_mc_sens = FALSE)
  }
}
