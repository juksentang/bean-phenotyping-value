# ============================================================================
# Layer 2 — PAWN (distribution-based) sensitivity
# ============================================================================
# PAWN (Pianosi & Wagener 2015, EMS 67:1–11; Pianosi et al. 2017,
# https://doi.org/10.1016/j.envsoft.2016.02.008) measures sensitivity via the
# Kolmogorov–Smirnov distance between the unconditional CDF of Y and the
# conditional CDF of Y given X_i lies in one of M strata. Unlike Sobol it is
# not driven by moments, so it captures tail and threshold effects that
# variance misses — e.g. "NPV > threshold" probability shifts.
#
# We use the *generic* (single-sample) estimator: one MC sample of size N is
# used for every parameter; for each X_i we bin the sample into M strata of
# equal width on [0,1] and compute KS(F_{Y|X_i∈bin_m}, F_Y). The summary
# statistic is the median KS across bins (Pianosi's recommendation; median is
# more robust than max to sparse bins).
#
# Output:
#   outputs/pawn_indices.csv         (param, PAWN_median, PAWN_max, CI_lo, CI_hi)
#   figures/pawn_kde.png             KDEs of unconditional vs per-bin conditional Y

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))
suppressPackageStartupMessages(library(ggplot2))


# --- PAWN core -----------------------------------------------------------
# f: unit-cube function → scalar
# N: sample size; M: number of strata per parameter
# Returns list of PAWN statistics per parameter + the raw sample (X, Y).
pawn_core <- function(f, d, N = 10000, M = 10, n_boot = 200,
                        progress = TRUE) {
  set.seed(2026)
  if (progress) cat(sprintf("  PAWN: sampling %d points in [0,1]^%d ...\n", N, d))
  X <- matrix(runif(N * d), nrow = N, ncol = d)
  Y <- numeric(N)
  for (i in seq_len(N)) Y[i] <- f(X[i, ])
  Y_sorted <- sort(Y)

  # Unconditional empirical CDF evaluator
  ecdf_uncond <- ecdf(Y)

  pawn_med <- numeric(d)
  pawn_max <- numeric(d)
  ci_lo    <- numeric(d)
  ci_hi    <- numeric(d)

  # Per-bin KS distances (d × M)
  ks_all <- matrix(0, nrow = d, ncol = M)

  edges <- seq(0, 1, length.out = M + 1)

  for (i in seq_len(d)) {
    bin_idx <- pmin(M, pmax(1, findInterval(X[, i], edges, rightmost.closed = TRUE)))
    for (m in seq_len(M)) {
      sel <- bin_idx == m
      n_bin <- sum(sel)
      if (n_bin < 5) { ks_all[i, m] <- NA_real_; next }
      # KS on pooled samples Y[sel] vs Y
      Y_bin <- Y[sel]
      # Unified eval grid
      grid <- sort(unique(c(Y_bin, Y)))
      ks_all[i, m] <- max(abs(ecdf(Y_bin)(grid) - ecdf_uncond(grid)))
    }
    pawn_med[i] <- median(ks_all[i, ], na.rm = TRUE)
    pawn_max[i] <- max(ks_all[i, ], na.rm = TRUE)

    # Bootstrap CIs on median KS (resample (X,Y) pairs)
    if (n_boot > 0) {
      boot_med <- numeric(n_boot)
      for (b in seq_len(n_boot)) {
        ii <- sample.int(N, N, replace = TRUE)
        Xb_i <- X[ii, i]
        Yb   <- Y[ii]
        bin_b <- pmin(M, pmax(1, findInterval(Xb_i, edges,
                                                 rightmost.closed = TRUE)))
        ks_boot <- numeric(M)
        for (m in seq_len(M)) {
          sel <- bin_b == m
          if (sum(sel) < 5) { ks_boot[m] <- NA_real_; next }
          Y_bin <- Yb[sel]
          grid <- sort(unique(c(Y_bin, Yb)))
          ks_boot[m] <- max(abs(ecdf(Y_bin)(grid) - ecdf(Yb)(grid)))
        }
        boot_med[b] <- median(ks_boot, na.rm = TRUE)
      }
      ci_lo[i] <- quantile(boot_med, 0.025, na.rm = TRUE)
      ci_hi[i] <- quantile(boot_med, 0.975, na.rm = TRUE)
    }
    if (progress) cat(sprintf("    [%d/%d] done\n", i, d))
  }

  list(
    pawn_median = pawn_med,
    pawn_max    = pawn_max,
    ci_lo       = ci_lo,
    ci_hi       = ci_hi,
    ks_matrix   = ks_all,
    X = X, Y = Y, edges = edges
  )
}


run_pawn <- function(N = 10000, M = 10, n_boot = 200,
                      output_dir = NULL, fig_dir = NULL) {
  sobol_cfg <- load_econ_grid()$sobol_params
  econ_cfg  <- load_config("economic_params")
  breed_cfg <- load_config("breeding_params")

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
  param_names <- c(all_l1, "bean_price", "discount_rate", "adoption_ceiling")
  group <- c(rep("L1 technical", 5), rep("L2 economic", 3))
  d <- length(param_names)

  l1_lookup <- build_l1_grid_lookup(free_dims = all_l1)
  f <- make_npv_fn(l1_lookup, sobol_cfg, econ_cfg, breed_cfg)

  cat(sprintf("=== PAWN  (d=%d, N=%d, M=%d strata, bootstrap=%d) ===\n",
              d, N, M, n_boot))
  t0 <- Sys.time()
  res <- pawn_core(f, d = d, N = N, M = M, n_boot = n_boot)
  cat(sprintf("elapsed: %.2f min\n",
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  out <- data.table(
    param       = param_names,
    PAWN_median = res$pawn_median,
    PAWN_max    = res$pawn_max,
    CI_lo       = res$ci_lo,
    CI_hi       = res$ci_hi,
    group       = group
  )
  setorder(out, -PAWN_median)

  cat("\n=== PAWN indices (sorted by median KS) ===\n"); print(out)

  fwrite(out, file.path(output_dir, "pawn_indices.csv"))

  # --- KDE plot: unconditional vs per-bin conditional for top-4 drivers ---
  top4 <- head(out$param, 4)
  dens_rows <- list()
  # Unconditional
  Y <- res$Y
  for (p in top4) {
    i <- which(param_names == p)
    bin_idx <- pmin(M, pmax(1, findInterval(res$X[, i], res$edges,
                                               rightmost.closed = TRUE)))
    for (m in seq_len(M)) {
      sel <- bin_idx == m
      if (sum(sel) < 20) next
      dens <- density(Y[sel])
      dens_rows[[length(dens_rows) + 1]] <- data.table(
        param = p, stratum = sprintf("bin %d/%d", m, M),
        x = dens$x, y = dens$y, kind = "conditional")
    }
    dens <- density(Y)
    dens_rows[[length(dens_rows) + 1]] <- data.table(
      param = p, stratum = "unconditional",
      x = dens$x, y = dens$y, kind = "unconditional")
  }
  dens_dt <- rbindlist(dens_rows)
  dens_dt[, param := factor(param, levels = top4)]

  p_kde <- ggplot(dens_dt, aes(x = x, y = y, group = stratum)) +
    geom_line(data = dens_dt[kind == "conditional"],
              aes(colour = stratum), alpha = 0.5, size = 0.5) +
    geom_line(data = dens_dt[kind == "unconditional"],
              colour = "black", size = 1) +
    facet_wrap(~param, scales = "free") +
    labs(
      title = "PAWN: Unconditional vs conditional NPV densities (top-4 drivers)",
      subtitle = "Black: unconditional.  Coloured: conditional on X_i bin.",
      x = "NPV (USD)", y = "density"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  fig_path <- file.path(fig_dir, "pawn_kde.png")
  ggsave(fig_path, p_kde, width = 10, height = 7, dpi = 150)
  cat("Figure saved:", fig_path, "\n")

  invisible(list(summary = out, raw = res, fig = fig_path))
}


if (sys.nframe() == 0) {
  run_pawn()
}
