# ============================================================================
# Layer 2 — Conservative Analysis (cycle_compression = 0)
# ============================================================================
# The full 5D sweep shows that cycle_compression dominates NPV at 61% of
# total-order variance, which masks the contribution of error_reduction
# and cost_reduction. cycle_compression is also the LEAST verified of the
# five tool-defining parameters: tool accelerating phenotyping data capture
# by 70% (as reported for a smartphone tool) does not automatically shorten the full 10-yr
# variety-release pipeline, which involves national variety trials,
# farmer participatory steps, and biological season constraints.
#
# This module runs a conservative analysis that fixes cycle_compression
# at 0 (i.e., assumes tool gives zero release-time acceleration) and asks:
# under that worst case, is the tool still profitable? What conditions on
# error_reduction and cost_reduction are needed for break-even?
#
# Outputs:
#   layer2_economics/outputs/conservative_lookup.parquet  — 748 × 625 slice
#   layer2_economics/outputs/conservative_break_even.csv  — err×cost contour
#   layer2_economics/outputs/sobol_indices_cyc0.csv       — 7-dim Sobol
#   layer2_economics/outputs/figures/conservative_*.png

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/01_gain_to_value.R"))
source(file.path(getwd(), "layer2_economics/R/03_npv_irr.R"))
source(file.path(getwd(), "layer2_economics/R/05_visualize.R"))
source(file.path(getwd(), "layer2_economics/R/07_sweep_economics.R"))
source(file.path(getwd(), "layer2_economics/R/08_break_even.R"))
source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))


# --- Filter lookup to cycle_compression = 0 ---
conservative_lookup <- function(lookup = NULL) {
  if (is.null(lookup)) lookup <- load_lookup()
  sub <- lookup[cycle_compression == 0]
  cat(sprintf("Conservative subset: %d rows (%.1f%% of full lookup)\n",
              nrow(sub), 100 * nrow(sub) / nrow(lookup)))
  sub
}


# --- Marginal aggregator in per-ha/yr units ---
cons_marginal_eaa <- function(slice_dt, x, y,
                                stat = c("median", "mean", "min", "max")) {
  stat <- match.arg(stat)
  agg <- slice_dt[, .(eaa_per_ha_yr = get(stat)(eaa_per_ha_yr, na.rm = TRUE),
                       npv_median = median(npv, na.rm = TRUE)),
                   by = c(x, y)]
  setnames(agg, c(x, y, "eaa_per_ha_yr", "npv_median"))
  setorderv(agg, c(x, y))
  agg
}


# --- Main conservative analysis runner ---
run_conservative_analysis <- function(lookup = NULL,
                                        econ_cfg = NULL,
                                        N_sobol = NULL,
                                        output_dir = NULL) {
  if (is.null(lookup))   lookup   <- load_lookup()
  if (is.null(econ_cfg)) econ_cfg <- load_config("economic_params")
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  fig_dir <- file.path(output_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("╔══════════════════════════════════════════════╗\n")
  cat("║  Conservative analysis (cycle_compression=0) ║\n")
  cat("╚══════════════════════════════════════════════╝\n\n")

  cons <- conservative_lookup(lookup)

  # Write the subset for downstream inspection
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(cons, file.path(output_dir, "conservative_lookup.parquet"))
  } else {
    fwrite(cons, file.path(output_dir, "conservative_lookup.csv.gz"))
  }

  # --- Summary statistics in per-ha/yr terms ---
  cat("\n=== Conservative NPV & EAA summary ===\n")
  cat(sprintf("NPV range:     $%.1fM to $%.1fM  (median $%.1fM)\n",
              min(cons$npv)/1e6, max(cons$npv)/1e6, median(cons$npv)/1e6))
  cat(sprintf("EAA per ha/yr: $%.2f to $%.2f  (median $%.2f)\n",
              min(cons$eaa_per_ha_yr), max(cons$eaa_per_ha_yr),
              median(cons$eaa_per_ha_yr)))
  cat(sprintf("P(NPV > 0):    %.1f%%\n", 100 * mean(cons$npv > 0)))

  # --- Baseline slice: price=$0.50, r=8%, A_max=35% (nearest), area=5M ---
  base_price <- econ_cfg$bean_price_usd_per_kg
  base_r     <- econ_cfg$discount_rate
  base_A     <- econ_cfg$adoption_ceiling
  base_area  <- econ_cfg$total_area_ha
  nearest <- function(x, g) g[which.min(abs(g - x))]
  bp <- nearest(base_price, sort(unique(cons$bean_price)))
  br <- nearest(base_r,     sort(unique(cons$discount_rate)))
  bA <- nearest(base_A,     sort(unique(cons$adoption_ceiling)))
  bAr <- nearest(base_area, sort(unique(cons$total_area_ha)))
  cat(sprintf("\nBaseline econ slice: price=%.2f, r=%.2f, A_max=%.2f, area=%.0f\n",
              bp, br, bA, bAr))

  slice <- cons[bean_price == bp & discount_rate == br &
                 adoption_ceiling == bA & total_area_ha == bAr]
  cat(sprintf("  → %d cells (cyc=0 × single econ point)\n", nrow(slice)))
  cat(sprintf("  EAA/ha/yr range: $%.2f to $%.2f\n",
              min(slice$eaa_per_ha_yr), max(slice$eaa_per_ha_yr)))
  cat(sprintf("  P(NPV > 0): %.1f%%\n", 100 * mean(slice$npv > 0)))

  # --- 2D break-even on (err_reduction × cost_reduction) ---
  marg_ec_eaa <- cons_marginal_eaa(slice, "error_reduction", "cost_reduction")
  setnames(marg_ec_eaa, c("error_reduction", "cost_reduction"), c("x_val", "y_val"))

  # Extract NPV=0 contour from the EAA marginal (sign is preserved under EAA)
  x_levels <- sort(unique(marg_ec_eaa$x_val))
  y_levels <- sort(unique(marg_ec_eaa$y_val))
  z_mat <- matrix(NA_real_, nrow = length(x_levels), ncol = length(y_levels))
  for (i in seq_along(x_levels)) for (j in seq_along(y_levels)) {
    z_mat[i, j] <- marg_ec_eaa[x_val == x_levels[i] & y_val == y_levels[j], eaa_per_ha_yr][1]
  }
  cl <- grDevices::contourLines(x_levels, y_levels, z_mat, levels = 0)
  contour_dt <- if (length(cl) == 0) {
    data.table(x = numeric(), y = numeric(), segment = integer())
  } else {
    do.call(rbind, lapply(seq_along(cl), function(s) {
      data.table(x = cl[[s]]$x, y = cl[[s]]$y, segment = s)
    }))
  }
  fwrite(marg_ec_eaa, file.path(output_dir, "conservative_marginal_err_cost.csv"))
  fwrite(contour_dt,  file.path(output_dir, "conservative_break_even_contour.csv"))
  cat(sprintf("\nBreak-even contour: %d points\n", nrow(contour_dt)))

  # --- Figures ---
  # EAA heatmap with break-even contour
  plot_response_heatmap(
    marg_ec_eaa,
    x_name = "x_val", y_name = "y_val",
    fill_name = "eaa_per_ha_yr",
    title = "Conservative case: cycle_compression = 0",
    x_label = "Error reduction",
    y_label = "Cost reduction",
    contour_dt = contour_dt,
    output_path = file.path(fig_dir, "conservative_heatmap_err_cost.png"),
    fill_label = "USD/ha/yr (EAA)",
    fmt = "%.2f"
  )

  # NPV EAA histogram across full conservative subset
  p_hist <- ggplot(cons, aes(x = eaa_per_ha_yr)) +
    geom_histogram(bins = 60, fill = "steelblue", alpha = 0.75, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick", linewidth = 1) +
    labs(title = "Conservative case: EAA distribution across L2 × L1 grid",
         subtitle = "cycle_compression = 0 (no release-time speedup)",
         x = "Equivalent annual annuity (USD/ha/yr)",
         y = "Count",
         caption = sprintf("P(NPV > 0) = %.1f%%", 100 * mean(cons$npv > 0))) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "conservative_eaa_histogram.png"),
          p_hist, width = 9, height = 5.5, dpi = 150)
  cat("conservative_eaa_histogram.png saved\n")

  # --- Conditional Sobol: 7-dim (drop cycle_compression, fix at 0) ---
  cat("\n━━━ Conditional Sobol (7D, cycle_compression = 0) ━━━\n")
  sobol_cons <- run_sobol(
    N = N_sobol,
    output_dir = output_dir,
    fix_l1 = list(cycle_compression = 0),
    label = "cyc0"
  )

  plot_sobol_bars(
    sobol_cons,
    output_path = file.path(fig_dir, "sobol_indices_cyc0.png")
  )

  invisible(list(
    conservative = cons,
    slice        = slice,
    marginal     = marg_ec_eaa,
    contour      = contour_dt,
    sobol        = sobol_cons
  ))
}


if (sys.nframe() == 0) {
  run_conservative_analysis()
}
