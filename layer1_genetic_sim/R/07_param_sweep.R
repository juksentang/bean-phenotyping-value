# ============================================================================
# Layer 1 — Parameter Sweep Driver (sensitivity sweep, HPC array-job ready)
# ============================================================================
# Runs ONE cell of the 5D tool sensitivity grid defined in
# config/sensitivity_grid.yaml, calling run_monte_carlo() with a scenario
# override and writing one parquet file per cell.
#
# Usage:
#   Rscript layer1_genetic_sim/R/07_param_sweep.R <cell_id>          # run cell N
#   Rscript layer1_genetic_sim/R/07_param_sweep.R <cell_id> <reps>   # override reps
#   Rscript layer1_genetic_sim/R/07_param_sweep.R --print-grid       # list cells
#   Rscript layer1_genetic_sim/R/07_param_sweep.R --count            # print cell count
#
# Designed to be called from a SLURM array job with
#   Rscript layer1_genetic_sim/R/07_param_sweep.R $SLURM_ARRAY_TASK_ID
#
# Model version 2 (model_version column): single-draw phenotyping at the
# genotype-mean h2, gain measured vs the check varieties (delta_g) with the
# legacy founder-mean gain kept as delta_g_vs_base, and common random numbers
# across scenarios. Output goes to outputs/<output_subdir> (sweep_v2 in the
# config); override it with the TOOL_SWEEP_SUBDIR env var, e.g. for smoke tests:
#   TOOL_SWEEP_SUBDIR=sweep_v2_smoke Rscript layer1_genetic_sim/R/07_param_sweep.R 0 3
# The v1 results directory (outputs/sweep) is never written by this script.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

source(file.path(getwd(), "layer1_genetic_sim/R/05_monte_carlo.R"))

# Small helper for NULL coalesce — defined before run_cell uses it
`%||%` <- function(a, b) if (is.null(a)) b else a

# Model version stamped on every per-cell output (1 = original HPC sweep)
MODEL_VERSION <- 2L

# --- Resolve + guard the output directory ---
# Precedence: explicit output_dir arg > TOOL_SWEEP_SUBDIR env var > config
# output_subdir > "sweep_v2". Refuses the v1 results directory
# (outputs/sweep) AND anything inside it, however the path is spelled
# ("sweep", "./sweep", "../outputs/sweep", "sweep/summaries", ...).
resolve_output_dir <- function(sweep_cfg, output_dir = NULL) {
  if (is.null(output_dir)) {
    subdir <- Sys.getenv("TOOL_SWEEP_SUBDIR", unset = "")
    if (!nzchar(subdir)) subdir <- sweep_cfg$output_subdir %||% "sweep_v2"
    output_dir <- file.path(getwd(), "layer1_genetic_sim/outputs", subdir)
  }
  # Collapse "." / ".." lexically first (normalizePath leaves a not-yet-existing
  # path untouched), then resolve symlinks of whatever exists
  abs_dir <- if (startsWith(output_dir, "/")) output_dir else file.path(getwd(), output_dir)
  parts <- strsplit(abs_dir, "/", fixed = TRUE)[[1]]
  keep <- character()
  for (p in parts[nzchar(parts) & parts != "."]) {
    keep <- if (p == "..") head(keep, -1L) else c(keep, p)
  }
  target <- normalizePath(paste0("/", paste(keep, collapse = "/")), mustWork = FALSE)
  v1_dir <- normalizePath(file.path(getwd(), "layer1_genetic_sim/outputs/sweep"),
                          mustWork = FALSE)
  if (target == v1_dir || startsWith(target, paste0(v1_dir, "/")) ||
      basename(target) == "sweep") {
    stop("Refusing to use the v1 results directory 'sweep' (or a path inside it): ",
         output_dir,
         "\nSet output_subdir (config) or TOOL_SWEEP_SUBDIR to a v2 directory.")
  }
  output_dir
}

# --- Load sensitivity grid config ---
# Default grid: config/sensitivity_grid.yaml; TOOL_SWEEP_GRID=<path> selects
# another grid file (extension blocks, robustness grids).
load_sweep_grid <- function(path = NULL) {
  if (is.null(path)) {
    path <- Sys.getenv("TOOL_SWEEP_GRID", unset = "")
    if (!nzchar(path)) path <- file.path(getwd(), "config/sensitivity_grid.yaml")
  }
  if (!file.exists(path)) stop("sensitivity_grid.yaml not found at: ", path)
  yaml::read_yaml(path)
}

# --- Expand the 5D grid into a deterministic data.table ---
# Column order and row order are fixed; cell_id = row index - 1 (0-based).
# An optional `genetic_grid` block (h2_yield, gxe_cor_onstation) adds columns
# that override genetic_params.yaml per cell (robustness sweeps); without it
# the grid and cell ids are unchanged.
expand_sweep_grid <- function(sweep_cfg) {
  g <- sweep_cfg$tool_grid
  gg <- sweep_cfg$genetic_grid
  grid <- do.call(expand.grid, c(list(
    error_reduction   = g$error_reduction,
    cost_reduction    = g$cost_reduction,
    cycle_compression = g$cycle_compression,
    tool_fixed_cost    = g$tool_fixed_cost,
    total_budget      = g$total_budget),
    gg[intersect(names(gg), c("h2_yield", "gxe_cor_onstation"))],
    list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ))
  dt <- as.data.table(grid)
  dt[, cell_id := seq_len(.N) - 1L]
  setcolorder(dt, "cell_id")
  dt[]
}

# --- Run one cell ---
run_cell <- function(cell_id, n_reps_override = NULL,
                     output_dir = NULL, sweep_cfg = NULL,
                     master_cfg = NULL, gen_cfg = NULL, breed_cfg = NULL,
                     verbose_mc = FALSE) {
  if (is.null(sweep_cfg))  sweep_cfg  <- load_sweep_grid()
  if (is.null(master_cfg)) master_cfg <- load_config("master_config")
  if (is.null(gen_cfg))    gen_cfg    <- load_config("genetic_params")
  if (is.null(breed_cfg))  breed_cfg  <- load_config("breeding_params")

  # Resolve (and guard) the output dir up front so a bad target fails fast
  output_dir <- resolve_output_dir(sweep_cfg, output_dir)

  grid_dt <- expand_sweep_grid(sweep_cfg)
  if (cell_id < 0 || cell_id >= nrow(grid_dt)) {
    stop(sprintf("cell_id %d out of range [0, %d)", cell_id, nrow(grid_dt)))
  }
  # Avoid data.table i-scope shadowing: `cell_id` is also a column name.
  .row_idx <- cell_id + 1L
  row <- grid_dt[.row_idx]

  # Per-cell genetic overrides (robustness grids only)
  for (nm in intersect(names(row), c("h2_yield", "gxe_cor_onstation"))) {
    gen_cfg[[nm]] <- row[[nm]]
  }

  override <- list(
    total_budget      = row$total_budget,
    error_reduction   = row$error_reduction,
    cost_reduction    = row$cost_reduction,
    cycle_compression = row$cycle_compression,
    tool_fixed_cost    = row$tool_fixed_cost
  )

  n_reps <- if (!is.null(n_reps_override)) n_reps_override else sweep_cfg$mc_reps_per_cell
  # Stable per-cell seed base so cells never share seeds. This is NOT a per-rep
  # replay: runMacs() ignores set.seed(), so a rep's founder haplotypes depend
  # on the process's call history (run order, mc.cores layout) and a single rep
  # cannot be re-run in isolation. Common random numbers hold WITHIN a rep
  # (both scenarios share the founder and the pipeline seed), which is what the
  # paired tool - traditional difference relies on.
  seed_base <- sweep_cfg$mc_seed_base %||% 42L
  master_cfg$mc_seed <- seed_base + 10000L * cell_id

  cat(sprintf("[cell %d/%d] err=%.2f cost_red=%.2f cyc_comp=%d fixed=$%d budget=$%d h2=%.2f gxe=%.2f reps=%d\n",
              cell_id, nrow(grid_dt) - 1L,
              row$error_reduction, row$cost_reduction,
              row$cycle_compression, row$tool_fixed_cost,
              row$total_budget, gen_cfg$h2_yield, gen_cfg$gxe_cor_onstation, n_reps))

  t0 <- Sys.time()
  mc <- run_monte_carlo(
    master_cfg        = master_cfg,
    gen_cfg           = gen_cfg,
    breed_cfg         = breed_cfg,
    n_reps            = n_reps,
    use_parallel      = TRUE,
    scenario_override = override,
    verbose           = verbose_mc
  )
  elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Tag results with cell parameters
  res <- copy(mc$results)
  res[, `:=`(
    cell_id           = cell_id,
    model_version     = MODEL_VERSION,
    error_reduction   = row$error_reduction,
    cost_reduction    = row$cost_reduction,
    cycle_compression = row$cycle_compression,
    tool_fixed_cost    = row$tool_fixed_cost,
    total_budget      = row$total_budget,
    h2_yield          = gen_cfg$h2_yield,
    gxe_cor_onstation = gen_cfg$gxe_cor_onstation
  )]
  setcolorder(res, c("cell_id", "model_version", "error_reduction", "cost_reduction",
                     "cycle_compression", "tool_fixed_cost", "total_budget",
                     "rep", "scenario"))

  # Per-cell summary (paired ΔG diff)
  # delta_g = gain vs check varieties (primary); delta_g_vs_base = gain vs
  # founder mean (v1 definition). The paired diff is the same under both
  # baselines as long as both scenarios release the same number of lines.
  paired <- merge(
    res[scenario == "traditional", .(rep, dg_trad = delta_g)],
    res[scenario == "tool",         .(rep, dg_tool  = delta_g)],
    by = "rep"
  )
  paired[, dg_diff := dg_tool - dg_trad]
  cell_summary <- data.table(
    cell_id          = cell_id,
    model_version    = MODEL_VERSION,
    error_reduction  = row$error_reduction,
    cost_reduction   = row$cost_reduction,
    cycle_compression = row$cycle_compression,
    tool_fixed_cost   = row$tool_fixed_cost,
    total_budget     = row$total_budget,
    h2_yield         = gen_cfg$h2_yield,
    gxe_cor_onstation = gen_cfg$gxe_cor_onstation,
    n_reps           = n_reps,
    mean_dg_trad     = mean(res[scenario == "traditional", delta_g]),
    mean_dg_tool      = mean(res[scenario == "tool", delta_g]),
    mean_dg_diff     = mean(paired$dg_diff),
    sd_dg_diff       = sd(paired$dg_diff),
    p_tool_gt_trad    = mean(paired$dg_diff > 0),
    mean_dg_vs_base_trad = mean(res[scenario == "traditional", delta_g_vs_base]),
    mean_dg_vs_base_tool  = mean(res[scenario == "tool", delta_g_vs_base]),
    h2_plot_trad     = mean(res[scenario == "traditional", h2_plot]),
    h2_plot_tool      = mean(res[scenario == "tool", h2_plot]),
    elapsed_s        = elapsed_s
  )

  # Write parquet
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fname <- sprintf(sweep_cfg$cell_filename_pattern %||% "cell_%05d.parquet", cell_id)
  fpath <- file.path(output_dir, fname)

  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(res, fpath)
  } else {
    fpath <- sub("\\.parquet$", ".csv.gz", fpath)
    fwrite(res, fpath)
  }
  cat(sprintf("[cell %d] wrote %s (%d rows, %.1f s)\n",
              cell_id, fpath, nrow(res), elapsed_s))

  # Per-cell summary CSV in its own subdir — one file per cell so that
  # parallel SLURM array tasks never write to the same file. Aggregate
  # later with aggregate_sweep().
  summary_dir <- file.path(output_dir, "summaries")
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(cell_summary,
         file.path(summary_dir, sprintf("cell_%05d.csv", cell_id)))

  invisible(list(results = res, summary = cell_summary, path = fpath))
}

# --- Aggregate all per-cell summaries into one table ---
# Run once after the SLURM array completes.
aggregate_sweep <- function(output_dir = NULL, sweep_cfg = NULL) {
  if (is.null(sweep_cfg)) sweep_cfg <- load_sweep_grid()
  output_dir <- resolve_output_dir(sweep_cfg, output_dir)
  summary_dir <- file.path(output_dir, "summaries")
  files <- list.files(summary_dir, pattern = "^cell_\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No per-cell summaries found in ", summary_dir)
  # fill = TRUE tolerates column differences between model versions; summaries
  # without a model_version column are v1. Never mix versions in one table.
  all_summaries <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  if (!"model_version" %in% names(all_summaries)) all_summaries[, model_version := 1L]
  all_summaries[is.na(model_version), model_version := 1L]
  if (uniqueN(all_summaries$model_version) > 1) {
    stop("Mixed model_version in ", summary_dir, ": ",
         paste(sort(unique(all_summaries$model_version)), collapse = ", "))
  }
  setorder(all_summaries, cell_id)
  out_path <- file.path(output_dir, "cell_summaries.csv")
  fwrite(all_summaries, out_path)
  cat(sprintf("Aggregated %d cells -> %s\n", nrow(all_summaries), out_path))
  invisible(all_summaries)
}

# --- CLI entry ---
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("Usage: Rscript 07_param_sweep.R <cell_id> [n_reps] | --print-grid | --count")
  }
  if (args[1] == "--count") {
    sweep_cfg <- load_sweep_grid()
    cat(nrow(expand_sweep_grid(sweep_cfg)), "\n")
    quit(save = "no")
  }
  if (args[1] == "--print-grid") {
    sweep_cfg <- load_sweep_grid()
    print(expand_sweep_grid(sweep_cfg))
    quit(save = "no")
  }
  if (args[1] == "--aggregate") {
    aggregate_sweep()
    quit(save = "no")
  }
  cell_id <- as.integer(args[1])
  n_reps_override <- if (length(args) >= 2) as.integer(args[2]) else NULL
  run_cell(cell_id, n_reps_override = n_reps_override)
}
