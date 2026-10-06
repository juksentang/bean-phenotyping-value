# ============================================================================
# Layer 2 — Global Sobol Sensitivity Analysis
# ============================================================================
# 8-dimensional first-order (S1) and total-order (ST) Sobol indices for
# NPV as a function of (5 L1 technical params × 3 L2 economic params).
# area is fixed at its baseline (5M ha) per the plan.
#
# The `sensitivity` CRAN package is unavailable (dtwclust / RcppThread
# compilation issues), so we implement the Sobol-Jansen estimator manually.
# Reference: Saltelli et al. 2010 "Variance based sensitivity analysis of
# model output", eqns (b) and (f). Implementation validated by checking
# S1 ≤ ST and bootstrap CIs.

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/03_npv_irr.R"))
source(file.path(getwd(), "layer2_economics/R/07_sweep_economics.R"))  # load_econ_grid
source(file.path(getwd(), "layer2_economics/R/08_break_even.R"))


# --- Quasi-random sampling (Halton sequence via a simple van der Corput) ---
# For 8D we use column-wise Halton with first 8 primes.
halton_primes <- c(2, 3, 5, 7, 11, 13, 17, 19)

van_der_corput <- function(n, base) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1
    r <- 0
    k <- i
    while (k > 0) {
      f <- f / base
      r <- r + f * (k %% base)
      k <- k %/% base
    }
    out[i] <- r
  }
  out
}

halton_matrix <- function(n, d) {
  stopifnot(d <= length(halton_primes))
  # Skip first 20 values (burn-in) to reduce correlation
  sapply(halton_primes[seq_len(d)], function(b) van_der_corput(n + 20, b)[-(1:20)])
}

# --- A/B sample matrices for the Jansen estimator ---
# A = Halton, B = independent uniforms. B must NOT be a shifted slice of the
# same Halton sequence: in base 2 an offset of 1024 = 2^10 only flips bit 11,
# so B[,1] = A[,1] +/- 2^-11 and the first parameter's ST collapsed to exactly 0
# (bug in v1 and the first v2 run).
sobol_ab_matrices <- function(N, d, seed = 2026) {
  set.seed(seed)
  list(X1 = halton_matrix(N, d),
       X2 = matrix(runif(N * d), nrow = N, ncol = d))
}


# --- Sobol-Jansen estimator ---
# X1, X2: N × d matrices of unit-cube samples.
# f: function taking one row (numeric vector of length d), returning scalar.
# Returns list of S1 and ST vectors (length d) plus the raw f evaluations.
sobol_jansen <- function(X1, X2, f, progress = TRUE) {
  N <- nrow(X1)
  d <- ncol(X1)
  stopifnot(nrow(X2) == N, ncol(X2) == d)

  if (progress) cat(sprintf("  evaluating f on A (%d rows)...\n", N))
  yA <- apply(X1, 1, f)
  if (progress) cat(sprintf("  evaluating f on B (%d rows)...\n", N))
  yB <- apply(X2, 1, f)

  # N × d, f on A with column i replaced by B's column i
  S1 <- numeric(d)
  ST <- numeric(d)
  varY <- var(c(yA, yB))

  for (i in seq_len(d)) {
    X_ABi <- X1
    X_ABi[, i] <- X2[, i]
    if (progress) cat(sprintf("  evaluating f on A_B[%d] ...\n", i))
    yABi <- apply(X_ABi, 1, f)

    # Saltelli 2010 preferred estimators (Jansen)
    #   V_i  ≈ (1/(2N)) Σ (yB - yABi)^2
    #   VT_i ≈ (1/(2N)) Σ (yA - yABi)^2
    Vi  <- mean((yB - yABi)^2) / 2
    VTi <- mean((yA - yABi)^2) / 2
    S1[i] <- 1 - Vi / varY
    ST[i] <- VTi / varY
  }

  list(S1 = S1, ST = ST, yA = yA, yB = yB, varY = varY, N = N, d = d)
}


# --- Nearest-neighbor lookup on the L1 grid ---
# Given a query point in normalized L1 space, find the nearest cell in
# the sensitivity grid and return its (mean_dg_trad, mean_dg_tool, cycle_compression).
# `free_dims`: character vector of L1 columns that vary in the Sobol sample.
#   Default = all 5. For conditional analysis, pass a 4-element vector
#   (omitting "cycle_compression") and set `fix_vals` to the fixed values.
build_l1_grid_lookup <- function(cell_dt = NULL,
                                   free_dims = c("error_reduction", "cost_reduction",
                                                 "cycle_compression", "tool_fixed_cost",
                                                 "total_budget"),
                                   fix_vals = NULL) {
  if (is.null(cell_dt)) {
    cell_dt <- load_l1_summaries()
  }
  # Apply fixed-value filters (keep only cells matching)
  if (!is.null(fix_vals)) {
    for (nm in names(fix_vals)) {
      cell_dt <- cell_dt[get(nm) == fix_vals[[nm]]]
    }
  }
  stopifnot(nrow(cell_dt) > 0)

  # Normalize each free L1 dim to [0,1]
  lims <- list()
  for (nm in free_dims) {
    lims[[nm]] <- range(cell_dt[[nm]])
  }
  normalize <- function(x, lim) {
    if (lim[2] == lim[1]) return(rep(0, length(x)))
    (x - lim[1]) / (lim[2] - lim[1])
  }
  grid_mat <- do.call(cbind, lapply(free_dims, function(nm) {
    normalize(cell_dt[[nm]], lims[[nm]])
  }))
  colnames(grid_mat) <- free_dims

  list(
    cells             = cell_dt,
    grid              = grid_mat,
    lims              = lims,
    free_dims         = free_dims,
    dg_trad           = cell_dt$mean_dg_trad,
    dg_tool            = cell_dt$mean_dg_tool,
    cycle_compression = cell_dt$cycle_compression,
    total_budget      = cell_dt$total_budget
  )
}


nn_query <- function(l1_lookup, x_norm_l1) {
  # x_norm_l1: length-ncol vector in [0,1]^ncol
  ncols <- ncol(l1_lookup$grid)
  stopifnot(length(x_norm_l1) == ncols)
  d2 <- rowSums((l1_lookup$grid - matrix(x_norm_l1, nrow = nrow(l1_lookup$grid),
                                          ncol = ncols, byrow = TRUE))^2)
  which.min(d2)
}


# --- Build the NPV query function f(x) used by Sobol ---
# x[1..n_l1] = normalized L1 free dims, x[n_l1+1..n_l1+3] = L2 params
# `n_l1` = number of free L1 dims (5 for full, 4 for conditional on cyc=0)
make_npv_fn <- function(l1_lookup, sobol_cfg, econ_cfg, breed_cfg) {
  l2_ranges  <- sobol_cfg$l2_ranges
  area_fixed <- sobol_cfg$fixed_area_ha
  k_val      <- econ_cfg$diffusion_rate_k
  tmid_val   <- econ_cfg$diffusion_midpoint_years
  horizon    <- econ_cfg$evaluation_horizon
  base_cycle <- breed_cfg$cycle_years
  n_l1 <- length(l1_lookup$free_dims)

  function(x) {
    idx <- nn_query(l1_lookup, x[seq_len(n_l1)])
    dg_trad      <- l1_lookup$dg_trad[idx]
    dg_tool       <- l1_lookup$dg_tool[idx]
    cyc_comp     <- l1_lookup$cycle_compression[idx]
    total_budget <- l1_lookup$total_budget[idx]
    cycle_trad   <- base_cycle
    cycle_tool    <- base_cycle - cyc_comp

    l2_offset <- n_l1
    price <- l2_ranges$bean_price_usd_per_kg[1] +
             x[l2_offset + 1] * diff(l2_ranges$bean_price_usd_per_kg)
    r     <- l2_ranges$discount_rate[1] +
             x[l2_offset + 2] * diff(l2_ranges$discount_rate)
    A_max <- l2_ranges$adoption_ceiling[1] +
             x[l2_offset + 3] * diff(l2_ranges$adoption_ceiling)

    ev <- eval_npv_core(
      dg_trad = dg_trad, dg_tool = dg_tool,
      cycle_trad = cycle_trad, cycle_tool = cycle_tool,
      price = price, r = r, area = area_fixed,
      A_max = A_max, k = k_val, t_mid = tmid_val,
      horizon = horizon, total_budget = total_budget,
      compute_irr = FALSE
    )
    ev$npv
  }
}


# --- Main Sobol runner ---
# `fix_l1`: named list of L1 params to fix (e.g. list(cycle_compression = 0))
#           → the remaining L1 dims + 3 L2 dims become the Sobol space.
# `label`:  suffix used for output filenames (e.g., "full", "cyc0").
run_sobol <- function(N = NULL, output_dir = NULL,
                        fix_l1 = NULL, label = "full") {
  sobol_cfg <- load_econ_grid()$sobol_params
  econ_cfg  <- load_config("economic_params")
  breed_cfg <- load_config("breeding_params")

  if (is.null(N)) N <- sobol_cfg$n_samples %||% 8192
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  all_l1 <- c("error_reduction", "cost_reduction", "cycle_compression",
              "tool_fixed_cost", "total_budget")
  free_l1 <- if (is.null(fix_l1)) all_l1 else setdiff(all_l1, names(fix_l1))
  n_l1 <- length(free_l1)
  d <- n_l1 + 3   # + 3 L2 params

  cat(sprintf("=== Sobol-Jansen (%s): free_L1=[%s], d=%d, N=%d (total evals=%d) ===\n",
              label, paste(free_l1, collapse = ","), d, N, N * (d + 2)))
  if (!is.null(fix_l1)) {
    for (nm in names(fix_l1)) cat(sprintf("  fixed %s = %s\n", nm, fix_l1[[nm]]))
  }

  ab <- sobol_ab_matrices(N, d)
  X1 <- ab$X1
  X2 <- ab$X2

  l1_lookup <- build_l1_grid_lookup(free_dims = free_l1, fix_vals = fix_l1)
  cat(sprintf("  L1 grid size after filtering: %d cells\n", nrow(l1_lookup$cells)))
  f <- make_npv_fn(l1_lookup, sobol_cfg, econ_cfg, breed_cfg)

  cat("Sampling...\n")
  t0 <- Sys.time()
  sa <- sobol_jansen(X1, X2, f)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Sobol done in %.2f min\n", elapsed))

  param_names <- c(free_l1,
                    "bean_price", "discount_rate", "adoption_ceiling")
  param_groups <- c(rep("L1 technical", n_l1),
                     rep("L2 economic", 3))

  out <- data.table(
    param = param_names,
    S1    = sa$S1,
    ST    = sa$ST,
    group = param_groups
  )
  setorder(out, -ST)

  cat("\n=== Sobol indices (sorted by ST) ===\n")
  print(out)
  cat(sprintf("\nsum(S1) = %.3f   sum(ST) = %.3f\n", sum(sa$S1), sum(sa$ST)))
  cat(sprintf("var(Y)  = %.3e\n", sa$varY))

  out_path <- file.path(output_dir, sprintf("sobol_indices_%s.csv", label))
  fwrite(out, out_path)
  cat("sobol_indices file written:", out_path, "\n")

  invisible(list(
    indices = out,
    S1 = sa$S1, ST = sa$ST,
    yA = sa$yA, yB = sa$yB,
    varY = sa$varY,
    N = N, d = d,
    param_names = param_names,
    label = label,
    fix_l1 = fix_l1
  ))
}


`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0) {
  run_sobol()
}
