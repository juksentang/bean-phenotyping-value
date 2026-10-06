# ============================================================================
# Layer 2 — Economic Parameter Sweep Driver
# ============================================================================
# Builds the 8-dimensional (5 L1 × 4 L2 × cell_id) NPV lookup table by
# combining every L1 cell from `cell_summaries.csv` with every point in
# the L2 economic grid defined in `config/econ_sweep_grid.yaml`.
#
# Output: layer2_economics/outputs/lookup_table.parquet  (~1.87M rows)

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/03_npv_irr.R"))

# --- Load grid config ---
load_econ_grid <- function(path = NULL) {
  if (is.null(path)) path <- file.path(getwd(), "config/econ_sweep_grid.yaml")
  if (!file.exists(path)) stop("econ_sweep_grid.yaml not found at: ", path)
  yaml::read_yaml(path)
}


# --- Build lookup table ---
build_lookup_table <- function(cell_summaries_path = NULL,
                                econ_grid = NULL,
                                econ_cfg = NULL,
                                breed_cfg = NULL,
                                output_path = NULL) {
  if (is.null(cell_summaries_path)) {
    cell_summaries_path <- l1_summary_path()
  }
  if (!all(file.exists(cell_summaries_path))) {
    stop("cell_summaries.csv not found at: ", cell_summaries_path)
  }
  if (is.null(econ_grid))  econ_grid  <- load_econ_grid()$econ_grid
  if (is.null(econ_cfg))   econ_cfg   <- load_config("economic_params")
  if (is.null(breed_cfg))  breed_cfg  <- load_config("breeding_params")

  cell_dt <- load_l1_summaries(cell_summaries_path)
  cat(sprintf("Loaded %d L1 cells from %s\n", nrow(cell_dt), cell_summaries_path))

  # --- L2 Cartesian product (one row per (price, r, A, area)) ---
  econ_cj <- CJ(
    bean_price_usd_per_kg = econ_grid$bean_price_usd_per_kg,
    discount_rate         = econ_grid$discount_rate,
    adoption_ceiling      = econ_grid$adoption_ceiling,
    total_area_ha         = econ_grid$total_area_ha
  )
  cat(sprintf("L2 grid points: %d\n", nrow(econ_cj)))

  # --- Cross join L1 cells × L2 grid via index expansion ---
  n_cells <- nrow(cell_dt)
  n_econ  <- nrow(econ_cj)
  n_total <- n_cells * n_econ
  cat(sprintf("Total lookup rows: %d\n", n_total))

  # Expand via repetition (fast, stays in data.table land)
  cell_idx <- rep(seq_len(n_cells), each = n_econ)
  econ_idx <- rep(seq_len(n_econ), times = n_cells)

  lookup <- data.table(
    cell_id           = cell_dt$cell_id[cell_idx],
    error_reduction   = cell_dt$error_reduction[cell_idx],
    cost_reduction    = cell_dt$cost_reduction[cell_idx],
    cycle_compression = cell_dt$cycle_compression[cell_idx],
    tool_fixed_cost    = cell_dt$tool_fixed_cost[cell_idx],
    total_budget      = cell_dt$total_budget[cell_idx],
    mean_dg_trad      = cell_dt$mean_dg_trad[cell_idx],
    mean_dg_tool       = cell_dt$mean_dg_tool[cell_idx],
    bean_price        = econ_cj$bean_price_usd_per_kg[econ_idx],
    discount_rate     = econ_cj$discount_rate[econ_idx],
    adoption_ceiling  = econ_cj$adoption_ceiling[econ_idx],
    total_area_ha     = econ_cj$total_area_ha[econ_idx]
  )

  # Derived cycle lengths
  lookup[, cycle_trad := breed_cfg$cycle_years]
  lookup[, cycle_tool  := breed_cfg$cycle_years - cycle_compression]

  # Constants from base econ_cfg (not varied in the sweep)
  k_val    <- econ_cfg$diffusion_rate_k
  tmid_val <- econ_cfg$diffusion_midpoint_years
  horizon  <- econ_cfg$evaluation_horizon

  # --- Row-wise NPV evaluation ---
  cat(sprintf("Computing NPV for %d rows (compute_irr = FALSE)...\n", n_total))
  t0 <- Sys.time()

  # Hoist columns into plain R vectors for fast indexing in the loop
  v_dg_trad       <- lookup$mean_dg_trad
  v_dg_tool        <- lookup$mean_dg_tool
  v_cycle_trad    <- lookup$cycle_trad
  v_cycle_tool     <- lookup$cycle_tool
  v_price         <- lookup$bean_price
  v_r             <- lookup$discount_rate
  v_area          <- lookup$total_area_ha
  v_A_max         <- lookup$adoption_ceiling
  v_total_budget  <- lookup$total_budget

  npv_out        <- numeric(n_total)
  bcr_out        <- numeric(n_total)
  cost_per_dg    <- numeric(n_total)
  incr_dg_out    <- numeric(n_total)

  progress_step <- max(1, n_total %/% 20)
  for (i in seq_len(n_total)) {
    ev <- eval_npv_core(
      dg_trad      = v_dg_trad[i],
      dg_tool       = v_dg_tool[i],
      cycle_trad   = v_cycle_trad[i],
      cycle_tool    = v_cycle_tool[i],
      price        = v_price[i],
      r            = v_r[i],
      area         = v_area[i],
      A_max        = v_A_max[i],
      k            = k_val,
      t_mid        = tmid_val,
      horizon      = horizon,
      total_budget = v_total_budget[i],
      compute_irr  = FALSE
    )
    npv_out[i]     <- ev$npv
    bcr_out[i]     <- ev$bcr
    cost_per_dg[i] <- ev$cost_per_unit_dg
    incr_dg_out[i] <- ev$incr_ann_dg
    if (i %% progress_step == 0) {
      cat(sprintf("  %d / %d (%.0f%%)\n", i, n_total, 100 * i / n_total))
    }
  }

  lookup[, `:=`(
    npv         = npv_out,
    bcr         = bcr_out,
    cost_per_dg = cost_per_dg,
    incr_ann_dg = incr_dg_out
  )]

  # Per-ha per-year equivalent annuity (interpretable for the reader).
  # See 01_gain_to_value.R:eaa_per_ha_per_yr.
  lookup[, eaa_per_ha_yr := eaa_per_ha_per_yr(npv, discount_rate, horizon, total_area_ha)]

  elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("NPV compute done in %.2f min\n", elapsed_min))

  # --- Write output ---
  if (is.null(output_path)) {
    output_dir <- l2_out_dir()
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    output_path <- file.path(output_dir, "lookup_table.parquet")
  }
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(lookup, output_path)
  } else {
    output_path <- sub("\\.parquet$", ".csv.gz", output_path)
    fwrite(lookup, output_path)
  }
  cat("Lookup table written:", output_path, "\n")
  cat(sprintf("Size: %.1f MB\n", file.info(output_path)$size / 1e6))

  invisible(lookup)
}


# --- CLI entry ---
if (sys.nframe() == 0) {
  build_lookup_table()
}
