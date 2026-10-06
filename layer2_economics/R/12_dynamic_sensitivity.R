# ============================================================================
# Layer 2 — Dynamic (Time-varying) Sobol Sensitivity
# ============================================================================
# Re-runs the 8-D Sobol decomposition on the *cumulative* NPV truncated at
# each year t ∈ {1, ..., horizon}. Reveals which parameters drive the tool value
# at different investment horizons:
#   early years  → breeding cost, fixed cost  (investment phase)
#   mid years    → adoption ramp
#   late years   → cycle compression (compound effect)
#
# Implementation: we share the Halton design with 09_sobol.R and evaluate a
# vector-valued f(x) that returns cumulative NPV over all years. One pass
# through the sample amortizes across all horizons.
#
# Output:
#   outputs/dynamic_sensitivity.csv     (long: year × param × S1 × ST)
#   outputs/dynamic_sensitivity_wide.csv
#   figures/dynamic_sensitivity.png

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))
suppressPackageStartupMessages(library(ggplot2))


# --- NPV core returning cumulative NPV at each year ----------------------
# Same math as eval_npv_core but returns a length-(horizon+1) vector of
# cumulative discounted incremental cash flow. Index t corresponds to NPV
# evaluated with horizon = t.
eval_cumnpv_core <- function(dg_trad, dg_tool, cycle_trad, cycle_tool,
                              price, r, area, A_max, k, t_mid,
                              horizon, total_budget) {
  ann_dg_trad <- dg_trad / cycle_trad
  ann_dg_tool  <- dg_tool  / cycle_tool
  release_trad <- cycle_trad
  release_tool  <- cycle_tool

  years <- 0:horizon
  benefits_trad <- numeric(horizon + 1)
  benefits_tool  <- numeric(horizon + 1)
  costs_trad    <- numeric(horizon + 1)
  costs_tool     <- numeric(horizon + 1)

  cost_per_yr_trad <- total_budget / cycle_trad
  cost_per_yr_tool  <- total_budget / cycle_tool

  for (i in seq_along(years)) {
    t <- years[i]
    if (t < cycle_trad) costs_trad[i] <- cost_per_yr_trad
    if (t < cycle_tool)  costs_tool[i]  <- cost_per_yr_tool
    if (t >= release_trad) {
      t_since <- t - release_trad
      adopt <- adoption_rate(t_since, A_max, k, t_mid)
      cum_dg <- ann_dg_trad * min(t_since + 1, cycle_trad)
      benefits_trad[i] <- cum_dg * price * adopt * area
    }
    if (t >= release_tool) {
      t_since <- t - release_tool
      adopt <- adoption_rate(t_since, A_max, k, t_mid)
      cum_dg <- ann_dg_tool * min(t_since + 1, cycle_tool)
      benefits_tool[i] <- cum_dg * price * adopt * area
    }
  }

  incr <- (benefits_tool - benefits_trad) - (costs_tool - costs_trad)
  disc <- 1 / (1 + r)^years
  cumsum(incr * disc)
}


# --- Vector-valued NPV query for Sobol -----------------------------------
make_cumnpv_fn <- function(l1_lookup, sobol_cfg, econ_cfg, breed_cfg) {
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

    eval_cumnpv_core(
      dg_trad = dg_trad, dg_tool = dg_tool,
      cycle_trad = cycle_trad, cycle_tool = cycle_tool,
      price = price, r = r, area = area_fixed,
      A_max = A_max, k = k_val, t_mid = tmid_val,
      horizon = horizon, total_budget = total_budget
    )
  }
}


# --- Vector-valued Sobol-Jansen ------------------------------------------
# f: x -> numeric vector of length L (L = horizon+1).
# Returns arrays S1 and ST of shape (L, d).
sobol_jansen_vec <- function(X1, X2, f, progress = TRUE) {
  N <- nrow(X1)
  d <- ncol(X1)

  if (progress) cat(sprintf("  f on A (%d rows, vector output)...\n", N))
  # First evaluation gives us the output dimension
  y1 <- f(X1[1, ])
  L <- length(y1)
  yA <- matrix(0, nrow = N, ncol = L)
  yA[1, ] <- y1
  for (i in 2:N) yA[i, ] <- f(X1[i, ])

  if (progress) cat(sprintf("  f on B (%d rows)...\n", N))
  yB <- matrix(0, nrow = N, ncol = L)
  for (i in seq_len(N)) yB[i, ] <- f(X2[i, ])

  varY <- matrix(0, nrow = 1, ncol = L)
  for (l in seq_len(L)) varY[l] <- var(c(yA[, l], yB[, l]))

  S1 <- matrix(0, nrow = d, ncol = L)
  ST <- matrix(0, nrow = d, ncol = L)

  for (i in seq_len(d)) {
    if (progress) cat(sprintf("  f on A_B[%d] ...\n", i))
    X_ABi <- X1
    X_ABi[, i] <- X2[, i]
    yABi <- matrix(0, nrow = N, ncol = L)
    for (row in seq_len(N)) yABi[row, ] <- f(X_ABi[row, ])

    for (l in seq_len(L)) {
      Vi  <- mean((yB[, l] - yABi[, l])^2) / 2
      VTi <- mean((yA[, l] - yABi[, l])^2) / 2
      vy <- varY[l]
      S1[i, l] <- if (vy > 0) 1 - Vi / vy else 0
      ST[i, l] <- if (vy > 0) VTi / vy else 0
    }
  }

  list(S1 = S1, ST = ST, varY = drop(varY), yA = yA, yB = yB, N = N, d = d, L = L)
}


# --- Main runner ---------------------------------------------------------
run_dynamic_sensitivity <- function(N = NULL, output_dir = NULL,
                                       fig_dir = NULL) {
  sobol_cfg <- load_econ_grid()$sobol_params
  econ_cfg  <- load_config("economic_params")
  breed_cfg <- load_config("breeding_params")

  if (is.null(N)) N <- sobol_cfg$n_samples %||% 8192
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  if (is.null(fig_dir)) {
    fig_dir <- l2_fig_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir,    showWarnings = FALSE, recursive = TRUE)

  all_l1 <- c("error_reduction", "cost_reduction", "cycle_compression",
              "tool_fixed_cost", "total_budget")
  d <- length(all_l1) + 3
  horizon <- econ_cfg$evaluation_horizon
  L <- horizon + 1

  cat(sprintf("=== Dynamic Sobol  (d=%d, N=%d, years=0..%d) ===\n",
              d, N, horizon))

  ab <- sobol_ab_matrices(N, d)
  X1 <- ab$X1
  X2 <- ab$X2

  l1_lookup <- build_l1_grid_lookup(free_dims = all_l1)
  f <- make_cumnpv_fn(l1_lookup, sobol_cfg, econ_cfg, breed_cfg)

  t0 <- Sys.time()
  sa <- sobol_jansen_vec(X1, X2, f)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Dynamic Sobol done in %.2f min\n", elapsed))

  param_names <- c(all_l1, "bean_price", "discount_rate", "adoption_ceiling")
  group <- c(rep("L1 technical", 5), rep("L2 economic", 3))

  # Long-format dataframe
  long_rows <- list()
  for (i in seq_len(d)) {
    for (l in seq_len(L)) {
      long_rows[[length(long_rows) + 1]] <- data.table(
        year  = l - 1,
        param = param_names[i],
        group = group[i],
        S1    = sa$S1[i, l],
        ST    = sa$ST[i, l],
        varY  = sa$varY[l]
      )
    }
  }
  long_dt <- rbindlist(long_rows)

  wide_ST <- dcast(long_dt, year ~ param, value.var = "ST")
  wide_S1 <- dcast(long_dt, year ~ param, value.var = "S1")

  fwrite(long_dt, file.path(output_dir, "dynamic_sensitivity.csv"))
  fwrite(wide_ST, file.path(output_dir, "dynamic_sensitivity_ST_wide.csv"))
  fwrite(wide_S1, file.path(output_dir, "dynamic_sensitivity_S1_wide.csv"))

  # --- Plot: ST over time --------------------------------------------------
  plot_dt <- long_dt[year > 0]   # year 0 has zero variance (all flows = 0)
  plot_dt[, param := factor(param, levels = param_names)]
  p <- ggplot(plot_dt, aes(x = year, y = ST, colour = param, linetype = group)) +
    geom_line(size = 0.9) +
    geom_point(size = 1.2) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      title = "Dynamic Sensitivity: Total-order Sobol indices over horizon",
      subtitle = "Cumulative NPV truncated at year t",
      x = "Evaluation horizon (years)",
      y = "Total-order Sobol index (ST)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "vertical")

  fig_path <- file.path(fig_dir, "dynamic_sensitivity_ST.png")
  ggsave(fig_path, p, width = 9, height = 6, dpi = 150)
  cat("Figure saved to:", fig_path, "\n")

  # --- Summary: top-3 drivers at year 5, 10, 15, 25 -----------------------
  cat("\n=== Top-3 drivers by year (ST) ===\n")
  for (y in c(5, 10, 15, 20, 25)) {
    if (y > horizon) next
    yd <- long_dt[year == y][order(-ST)][1:3]
    cat(sprintf("  year %d: %s  |  varY=%.2e\n",
                y, paste(sprintf("%s=%.2f", yd$param, yd$ST), collapse = ", "),
                yd$varY[1]))
  }

  invisible(list(long = long_dt, wide_ST = wide_ST, wide_S1 = wide_S1,
                  fig = fig_path, elapsed_min = elapsed))
}


if (sys.nframe() == 0) {
  run_dynamic_sensitivity()
}
