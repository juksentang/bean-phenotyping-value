# ============================================================================
# Layer 2 — input/output locations (model version switch)
# ============================================================================
# Every Layer-2 script resolves its Layer-1 input and its output directory
# here, so switching model versions is one environment variable, not 20 edits.
#
#   TOOL_L1_SUBDIR  Layer-1 sweep subdir(s) under layer1_genetic_sim/outputs, comma-separated;
#                  blocks are stacked into one grid (default: v2 main grid + extension blocks)
#   TOOL_L2_OUT     Layer-2 output dir, relative to the project root       (default "layer2_economics/outputs_v2")
#
# The v1 results (layer2_economics/outputs, report/figures) are frozen: writing
# there is refused. Reading v1 Layer 1 (TOOL_L1_SUBDIR=sweep) stays allowed.

L1_DEFAULT_SUBDIRS <- "sweep_v2,sweep_v2_extA,sweep_v2_extB"

l1_summary_path <- function() {
  subs <- strsplit(Sys.getenv("TOOL_L1_SUBDIR", unset = L1_DEFAULT_SUBDIRS), ",")[[1]]
  file.path(getwd(), "layer1_genetic_sim/outputs", trimws(subs), "cell_summaries.csv")
}

l2_out_dir <- function() {
  rel <- sub("/+$", "", Sys.getenv("TOOL_L2_OUT", unset = "layer2_economics/outputs_v2"))
  if (normalizePath(file.path(getwd(), rel), mustWork = FALSE) ==
      normalizePath(file.path(getwd(), "layer2_economics/outputs"), mustWork = FALSE))
    stop("refusing to write Layer-2 results into the frozen v1 directory layer2_economics/outputs")
  d <- file.path(getwd(), rel)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

l2_fig_dir <- function() {
  d <- file.path(l2_out_dir(), "figures")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# --- Layer-1 cell summaries, calibrated ---
# Reads cell_summaries.csv and rescales mean_dg_trad / mean_dg_tool (and the
# paired diff) by k = target gain per cycle / mean simulated traditional gain,
# with the target from economic_params.yaml (gain_target_pct_per_yr). The
# factor is attached as attr(, "gain_scale"). Every Layer-2 consumer of
# Layer-1 gains must go through here.
load_l1_summaries <- function(path = l1_summary_path()) {
  if (!exists("load_config", mode = "function"))
    source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))
  # Stack sweep blocks; cell ids are only unique within a block, so renumber
  # and keep the origin in `block`
  cell_dt <- data.table::rbindlist(lapply(path, function(p)
    data.table::fread(p)[, block := basename(dirname(p))]), use.names = TRUE, fill = TRUE)
  key_cols <- c("error_reduction", "cost_reduction", "cycle_compression",
                "tool_fixed_cost", "total_budget")
  if (anyDuplicated(cell_dt[, ..key_cols]))
    stop("duplicate grid points across L1 blocks: ", paste(basename(dirname(path)), collapse = ", "))
  cell_dt[, cell_id := seq_len(.N) - 1L]
  econ  <- load_config("economic_params")
  breed <- load_config("breeding_params")
  k <- 1
  if (!is.null(econ$gain_target_pct_per_yr)) {
    target_per_cycle <- econ$gain_target_pct_per_yr / 100 *
                        econ$mean_yield_kgha * breed$cycle_years
    k <- target_per_cycle / mean(cell_dt$mean_dg_trad)
    for (col in intersect(c("mean_dg_trad", "mean_dg_tool", "mean_dg_diff", "sd_dg_diff"),
                          names(cell_dt)))
      data.table::set(cell_dt, j = col, value = cell_dt[[col]] * k)
  }
  attr(cell_dt, "gain_scale") <- k
  cell_dt
}
