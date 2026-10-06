# ============================================================================
# Layer 2 — Break-Even Surface Extraction
# ============================================================================
# Extracts NPV = 0 iso-surface (break-even frontier) from the lookup table
# produced by 07_sweep_economics.R. Holds economic parameters at their
# baseline values and slices the 5D L1 technical space.
#
# Output:
#   layer2_economics/outputs/break_even_surface.csv — (x, y, npv, contour_flag)
#   for 2D slices; used by plot_break_even_contour in 05_visualize.R.

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/01_gain_to_value.R"))

load_lookup <- function(path = NULL) {
  if (is.null(path)) {
    p1 <- file.path(l2_out_dir(), "lookup_table.parquet")
    p2 <- file.path(l2_out_dir(), "lookup_table.csv.gz")
    path <- if (file.exists(p1)) p1 else p2
  }
  if (!file.exists(path)) stop("lookup table not found: ", path)
  if (grepl("\\.parquet$", path)) {
    as.data.table(arrow::read_parquet(path))
  } else {
    fread(path)
  }
}


# --- Extract the baseline economic slice ---
# Fixes all 4 L2 parameters at their baseline values. Remaining variation
# is purely across the 5D L1 cell grid → 2993 rows.
baseline_slice <- function(lookup, econ_cfg = NULL) {
  if (is.null(econ_cfg)) econ_cfg <- load_config("economic_params")
  base_price <- econ_cfg$bean_price_usd_per_kg
  base_r     <- econ_cfg$discount_rate
  base_A     <- econ_cfg$adoption_ceiling
  base_area  <- econ_cfg$total_area_ha

  # Use nearest value in the grid (lookup may not contain the exact base)
  nearest <- function(x, grid) grid[which.min(abs(grid - x))]

  price_grid <- sort(unique(lookup$bean_price))
  r_grid     <- sort(unique(lookup$discount_rate))
  A_grid     <- sort(unique(lookup$adoption_ceiling))
  area_grid  <- sort(unique(lookup$total_area_ha))

  bp <- nearest(base_price, price_grid)
  br <- nearest(base_r,     r_grid)
  bA <- nearest(base_A,     A_grid)
  bAr<- nearest(base_area,  area_grid)

  cat(sprintf("Baseline slice: price=%.2f, r=%.2f, A_max=%.2f, area=%.0f\n",
              bp, br, bA, bAr))
  slice <- lookup[bean_price == bp & discount_rate == br &
                   adoption_ceiling == bA & total_area_ha == bAr]
  cat(sprintf("  → %d cells in baseline slice\n", nrow(slice)))
  slice
}


# --- 2D marginal: collapse 3 dimensions by median, return xy × NPV grid ---
# Returns a data.table with (x, y, npv) suitable for geom_tile / contour.
marginal_2d <- function(slice, x = "error_reduction", y = "cost_reduction",
                         collapse = c("cycle_compression", "tool_fixed_cost",
                                      "total_budget"),
                         stat = "median") {
  stat_fn <- switch(stat,
                    median = function(x) median(x, na.rm = TRUE),
                    mean   = function(x) mean(x, na.rm = TRUE),
                    max    = function(x) max(x, na.rm = TRUE),
                    min    = function(x) min(x, na.rm = TRUE))
  agg <- slice[, .(npv = stat_fn(npv)), by = c(x, y)]
  setnames(agg, c(x, y, "npv"))
  setorderv(agg, c(x, y))
  agg
}


# --- Extract NPV = 0 contour line from a 2D marginal ---
extract_contour <- function(marginal_dt, x_col, y_col, level = 0) {
  x_vals <- sort(unique(marginal_dt[[x_col]]))
  y_vals <- sort(unique(marginal_dt[[y_col]]))
  z_mat <- matrix(NA_real_, nrow = length(x_vals), ncol = length(y_vals))
  for (i in seq_along(x_vals)) {
    for (j in seq_along(y_vals)) {
      z_mat[i, j] <- marginal_dt[get(x_col) == x_vals[i] & get(y_col) == y_vals[j], npv][1]
    }
  }
  if (all(is.na(z_mat))) return(NULL)
  cl <- grDevices::contourLines(x_vals, y_vals, z_mat, levels = level)
  if (length(cl) == 0) return(data.table(x = numeric(), y = numeric(), segment = integer()))
  do.call(rbind, lapply(seq_along(cl), function(s) {
    data.table(x = cl[[s]]$x, y = cl[[s]]$y, segment = s)
  }))
}


# --- Stressed economic slice (pessimistic corner of the L2 grid) ---
stressed_slice <- function(lookup) {
  price_lo <- min(lookup$bean_price)
  r_hi     <- max(lookup$discount_rate)
  A_lo     <- min(lookup$adoption_ceiling)
  area_lo  <- min(lookup$total_area_ha)
  cat(sprintf("Stressed slice: price=%.2f, r=%.2f, A_max=%.2f, area=%.0f\n",
              price_lo, r_hi, A_lo, area_lo))
  slice <- lookup[bean_price == price_lo & discount_rate == r_hi &
                   adoption_ceiling == A_lo & total_area_ha == area_lo]
  cat(sprintf("  → %d cells in stressed slice\n", nrow(slice)))
  slice
}


# --- Main extractor ---
extract_break_even <- function(lookup = NULL, econ_cfg = NULL,
                                output_dir = NULL) {
  if (is.null(lookup)) lookup <- load_lookup()
  if (is.null(econ_cfg)) econ_cfg <- load_config("economic_params")
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  slice <- baseline_slice(lookup, econ_cfg)

  # Three canonical 2D slices (x, y) for the story:
  slice_specs <- list(
    list(x = "error_reduction",   y = "cost_reduction",    label = "err_vs_cost"),
    list(x = "cycle_compression", y = "total_budget",      label = "cyc_vs_budget"),
    list(x = "error_reduction",   y = "cycle_compression", label = "err_vs_cyc")
  )

  all_marginals <- list()
  all_contours  <- list()
  for (sp in slice_specs) {
    collapse_cols <- setdiff(
      c("error_reduction", "cost_reduction", "cycle_compression",
        "tool_fixed_cost", "total_budget"),
      c(sp$x, sp$y)
    )
    m <- marginal_2d(slice, x = sp$x, y = sp$y, collapse = collapse_cols)
    m[, slice := sp$label]
    setnames(m, c("x_val", "y_val", "npv", "slice"))
    m[, x_name := sp$x][, y_name := sp$y]
    all_marginals[[sp$label]] <- m

    c <- extract_contour(m[, .(x_val, y_val, npv)], "x_val", "y_val", level = 0)
    if (!is.null(c) && nrow(c) > 0) {
      c[, slice := sp$label]
      c[, x_name := sp$x][, y_name := sp$y]
      all_contours[[sp$label]] <- c
    }
  }

  marginals_dt <- rbindlist(all_marginals)
  contours_dt  <- if (length(all_contours) > 0) rbindlist(all_contours) else
                  data.table(x = numeric(), y = numeric(), segment = integer(),
                             slice = character(), x_name = character(), y_name = character())

  fwrite(marginals_dt, file.path(output_dir, "break_even_marginals.csv"))
  fwrite(contours_dt,  file.path(output_dir, "break_even_contours.csv"))
  cat("break_even_marginals.csv:", nrow(marginals_dt), "rows\n")
  cat("break_even_contours.csv: ", nrow(contours_dt),  "rows\n")

  # --- Stressed slice: where break-even actually exists ---
  sslice <- stressed_slice(lookup)
  stressed_marginals <- list()
  stressed_contours  <- list()
  for (sp in slice_specs) {
    m <- marginal_2d(sslice, x = sp$x, y = sp$y)
    m[, slice := sp$label]
    setnames(m, c("x_val", "y_val", "npv", "slice"))
    m[, x_name := sp$x][, y_name := sp$y]
    stressed_marginals[[sp$label]] <- m
    c <- extract_contour(m[, .(x_val, y_val, npv)], "x_val", "y_val", level = 0)
    if (!is.null(c) && nrow(c) > 0) {
      c[, slice := sp$label]
      c[, x_name := sp$x][, y_name := sp$y]
      stressed_contours[[sp$label]] <- c
    }
  }
  stressed_marginals_dt <- rbindlist(stressed_marginals)
  stressed_contours_dt  <- if (length(stressed_contours) > 0) rbindlist(stressed_contours) else
    data.table(x = numeric(), y = numeric(), segment = integer(),
               slice = character(), x_name = character(), y_name = character())
  fwrite(stressed_marginals_dt, file.path(output_dir, "break_even_marginals_stressed.csv"))
  fwrite(stressed_contours_dt,  file.path(output_dir, "break_even_contours_stressed.csv"))
  cat("break_even_marginals_stressed.csv:", nrow(stressed_marginals_dt), "rows\n")
  cat("break_even_contours_stressed.csv: ", nrow(stressed_contours_dt),  "rows\n")

  # --- Profitability surface across L2 economic space ---
  # For every (price, r) combination in the grid, compute fraction of L1
  # cells with NPV > 0. This shows how the "break-even exists" frontier
  # moves with economic stress.
  profit_frac <- lookup[, .(
    frac_positive = mean(npv > 0),
    median_npv    = median(npv)
  ), by = .(bean_price, discount_rate, adoption_ceiling, total_area_ha)]
  fwrite(profit_frac, file.path(output_dir, "profit_surface.csv"))
  cat("profit_surface.csv:", nrow(profit_frac), "rows\n")

  invisible(list(
    slice              = slice,
    marginals          = marginals_dt,
    contours           = contours_dt,
    stressed_slice     = sslice,
    stressed_marginals = stressed_marginals_dt,
    stressed_contours  = stressed_contours_dt,
    profit_frac        = profit_frac
  ))
}


if (sys.nframe() == 0) {
  extract_break_even()
}
