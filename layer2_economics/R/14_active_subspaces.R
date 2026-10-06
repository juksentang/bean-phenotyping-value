# ============================================================================
# Layer 2 — Active Subspaces (Constantine 2015)
# ============================================================================
# Constantine PG (2015) "Active Subspaces: Emerging Ideas for Dimension
# Reduction in Parameter Studies", SIAM. Identifies the directions in
# parameter space along which the NPV response varies most on average:
#
#   C = E[∇f(x) ∇f(x)^T]   (d × d, symmetric PSD)
#   C = W Λ W^T             (eigendecomposition)
#
# Eigenvectors of C with the largest eigenvalues span the "active subspace".
# We compute ∇f via central finite differences; to tame the discrete L1
# look-up, we use a step size h = 0.04 (≈ 2 grid cells) so the FD average
# over a smoothed region rather than hitting a single piecewise-constant
# plateau.
#
# Output:
#   outputs/active_subspace_eigenvalues.csv
#   outputs/active_subspace_W.csv            (d × d: each column = eigenvector)
#   outputs/active_subspace_scores.csv        (N × first-k projected scores + y)
#   figures/active_subspace_spectrum.png
#   figures/active_subspace_1d.png            Sufficient-summary plot
#   figures/active_subspace_2d.png            2-D projection coloured by NPV

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))
suppressPackageStartupMessages(library(ggplot2))


# --- Central finite-difference gradient in [0,1]^d ----------------------
fd_grad <- function(f, x, h = 0.04) {
  d <- length(x)
  g <- numeric(d)
  for (i in seq_len(d)) {
    xp <- x; xp[i] <- min(1, xp[i] + h)
    xm <- x; xm[i] <- max(0, xm[i] - h)
    step <- xp[i] - xm[i]
    if (step <= 0) { g[i] <- 0; next }
    g[i] <- (f(xp) - f(xm)) / step
  }
  g
}


# --- Active subspace computation ----------------------------------------
active_subspace <- function(f, d, N = 500, h = 0.04, progress = TRUE) {
  set.seed(2026)
  X <- matrix(runif(N * d), nrow = N, ncol = d)
  Y <- numeric(N)
  G <- matrix(0, nrow = N, ncol = d)

  if (progress) cat(sprintf("  AS: gradient-sampling N=%d ...\n", N))
  for (n in seq_len(N)) {
    Y[n]   <- f(X[n, ])
    G[n, ] <- fd_grad(f, X[n, ], h = h)
    if (progress && n %% max(1, N %/% 10) == 0) {
      cat(sprintf("    %d/%d\n", n, N))
    }
  }

  # Normalize gradient scale — otherwise the huge spread of NPV values
  # makes the spectrum numerically unstable. We scale by sd(Y) so entries
  # of C are in a sensible range.
  sy <- sd(Y)
  if (sy > 0) G <- G / sy

  C <- t(G) %*% G / N
  # Symmetrize for numerical safety
  C <- (C + t(C)) / 2

  ee <- eigen(C, symmetric = TRUE)
  # Sort descending (eigen returns in decreasing order already, but ensure)
  ord <- order(-ee$values)
  list(
    eigenvalues  = ee$values[ord],
    W            = ee$vectors[, ord, drop = FALSE],
    X = X, Y = Y, G = G, C = C, scale_y = sy
  )
}


run_active_subspaces <- function(N = 500, h = 0.04,
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
  d <- length(param_names)

  l1_lookup <- build_l1_grid_lookup(free_dims = all_l1)
  f <- make_npv_fn(l1_lookup, sobol_cfg, econ_cfg, breed_cfg)

  cat(sprintf("=== Active Subspaces  (d=%d, N=%d gradients, h=%.3f) ===\n",
              d, N, h))
  t0 <- Sys.time()
  res <- active_subspace(f, d = d, N = N, h = h)
  cat(sprintf("elapsed: %.2f min\n",
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # --- Eigenvalue spectrum -----------------------------------------------
  lam <- res$eigenvalues
  tot <- sum(abs(lam))
  spec_dt <- data.table(
    k           = seq_len(d),
    eigenvalue  = lam,
    var_share   = abs(lam) / tot,
    cum_share   = cumsum(abs(lam)) / tot
  )
  cat("\n=== Eigenvalue spectrum ===\n"); print(spec_dt)
  fwrite(spec_dt, file.path(output_dir, "active_subspace_eigenvalues.csv"))

  # --- Eigenvectors (loadings) -------------------------------------------
  W_dt <- as.data.table(res$W)
  setnames(W_dt, paste0("w", seq_len(d)))
  W_dt[, param := param_names]
  setcolorder(W_dt, c("param", paste0("w", seq_len(d))))
  cat("\n=== First 3 eigenvectors (loadings) ===\n")
  print(W_dt[, c("param", "w1", "w2", "w3"), with = FALSE])
  fwrite(W_dt, file.path(output_dir, "active_subspace_W.csv"))

  # --- Scores (projections) ----------------------------------------------
  scores <- res$X %*% res$W[, 1:3, drop = FALSE]
  scores_dt <- data.table(
    u1 = scores[, 1], u2 = scores[, 2], u3 = scores[, 3],
    y  = res$Y
  )
  fwrite(scores_dt, file.path(output_dir, "active_subspace_scores.csv"))

  # --- Plots --------------------------------------------------------------
  spec_plot <- ggplot(spec_dt, aes(x = k, y = eigenvalue)) +
    geom_col(fill = "#2c7fb8") +
    scale_y_log10() +
    scale_x_continuous(breaks = 1:d) +
    labs(
      title = "Active Subspaces: eigenvalue spectrum of C",
      subtitle = sprintf("k=1 explains %.1f%%; k<=2 explains %.1f%%",
                          100 * spec_dt$cum_share[1],
                          100 * spec_dt$cum_share[2]),
      x = "component k", y = "eigenvalue (log)"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "active_subspace_spectrum.png"),
         spec_plot, width = 7, height = 5, dpi = 150)

  summary_1d <- ggplot(scores_dt, aes(x = u1, y = y)) +
    geom_point(alpha = 0.5, size = 1, colour = "#1b9e77") +
    geom_smooth(method = "loess", se = FALSE, colour = "#d95f02", size = 1) +
    labs(
      title = "Sufficient-summary plot along first active direction",
      subtitle = sprintf("w1 = %s",
                          paste(sprintf("%+.2f*%s",
                                          res$W[, 1], param_names),
                                 collapse = "  ")),
      x = "u1 = w1^T x", y = "NPV (USD)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8))
  ggsave(file.path(fig_dir, "active_subspace_1d.png"),
         summary_1d, width = 9, height = 6, dpi = 150)

  summary_2d <- ggplot(scores_dt, aes(x = u1, y = u2, colour = y)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_colour_viridis_c(option = "C", name = "NPV") +
    labs(
      title = "NPV in the 2-D active subspace",
      subtitle = "Colour = NPV. Horizontal axis = dominant direction.",
      x = "u1 = w1^T x", y = "u2 = w2^T x"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "active_subspace_2d.png"),
         summary_2d, width = 9, height = 6, dpi = 150)

  cat("Figures saved to:", fig_dir, "\n")
  invisible(list(
    eigenvalues = spec_dt,
    W           = W_dt,
    scores      = scores_dt,
    raw         = res
  ))
}


if (sys.nframe() == 0) {
  run_active_subspaces()
}
