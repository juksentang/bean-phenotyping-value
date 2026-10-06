# ============================================================================
# Layer 2 — Shapley Effects (strict, game-theoretic)
# ============================================================================
# Implements the random-permutation Shapley-effects estimator
# (Song, Nelson & Staum 2016, Owen 2014). Unlike Sobol, Shapley Effects
# remain well-defined under input correlation and always sum to 1.
#
# We compute TWO variants:
#   (1) independent inputs  — direct comparison with Sobol_full.
#   (2) Gaussian-copula correlated inputs — with a plausible correlation
#       structure reflecting:
#         discount_rate ↔ bean_price       +0.40  (inflation expectations)
#         adoption_ceiling ↔ total_budget  +0.30  (market structure)
#         cycle_compression ↔ error_reduction +0.30 (tool measurement quality)
#         cost_reduction ↔ tool_fixed_cost  -0.20  (amortization)
#
# Output:
#   outputs/shapley_indices_indep.csv
#   outputs/shapley_indices_copula.csv
#
# References:
#   Owen AB (2014). SIAM/ASA J. Uncertainty Quantification 2:245–251.
#   Song E, Nelson BL, Staum J (2016). SIAM/ASA JUQ 4:1060–1083.
#   Iooss B, Prieur C (2019). Int J Uncertainty Quantification 9:493–514.

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

source(file.path(getwd(), "layer2_economics/R/09_sobol.R"))
suppressPackageStartupMessages(library(MASS))


# --- Gaussian-copula helpers ---------------------------------------------
# Unit cube [0,1]^d  <-(Φ)->  R^d Gaussian  <-(Cholesky L)->  correlated Gaussian
# <-(Φ^{-1})-> correlated uniforms.
#
# Conditional sampling under Gaussian copula at the Gaussian level:
# given Z_S = z_S*, the conditional distribution of Z_{-S} is N(μ_cond, Σ_cond)
# with closed-form μ_cond, Σ_cond. We draw and map back to [0,1]^d.

# Build symmetric PSD correlation matrix from off-diagonal entries.
build_correlation <- function(d, pairs) {
  # pairs: list of list(i, j, rho)
  R <- diag(d)
  for (p in pairs) {
    R[p$i, p$j] <- p$rho
    R[p$j, p$i] <- p$rho
  }
  # Verify PSD — if not, nudge via eigenvalue flooring.
  ee <- eigen(R, symmetric = TRUE)
  if (min(ee$values) < 1e-8) {
    ee$values[ee$values < 1e-8] <- 1e-8
    R <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
    # Rescale to unit diagonal
    s <- sqrt(diag(R))
    R <- R / tcrossprod(s)
  }
  R
}


# Draw N joint samples from copula-correlated uniforms.
rcopula_unit <- function(N, R) {
  d <- nrow(R)
  Z <- MASS::mvrnorm(N, mu = rep(0, d), Sigma = R)
  pnorm(Z)
}


# Sample X_{-S} | X_S = u_S under a Gaussian copula with correlation R.
# u_S is a length-|S| vector in [0,1]. Returns a length-(d-|S|) sample in [0,1].
rcopula_conditional <- function(u_S, S, R, n = 1) {
  d <- nrow(R)
  free <- setdiff(seq_len(d), S)
  if (length(S) == 0) {
    # Unconditional draw of the free coords
    Z <- MASS::mvrnorm(n, mu = rep(0, length(free)),
                        Sigma = R[free, free, drop = FALSE])
    if (n == 1) Z <- matrix(Z, nrow = 1)
    return(pnorm(Z))
  }
  # Condition at the Gaussian level
  z_S <- qnorm(pmin(pmax(u_S, 1e-8), 1 - 1e-8))
  R_SS <- R[S, S, drop = FALSE]
  R_FS <- R[free, S, drop = FALSE]
  R_FF <- R[free, free, drop = FALSE]
  # Solve once for this conditioning set
  W <- R_FS %*% solve(R_SS)
  mu_F <- drop(W %*% z_S)
  Sigma_F <- R_FF - W %*% t(R_FS)
  # Symmetrize numerically
  Sigma_F <- (Sigma_F + t(Sigma_F)) / 2
  Z_F <- MASS::mvrnorm(n, mu = mu_F, Sigma = Sigma_F)
  if (n == 1) Z_F <- matrix(Z_F, nrow = 1)
  pnorm(Z_F)
}


# --- Random-permutation Shapley-effects estimator ------------------------
# f: function taking a length-d unit-cube vector, returning scalar.
# R: d×d Gaussian-copula correlation matrix (identity = independent).
# m: number of random permutations.
# N_o, N_i: outer and inner Monte Carlo sample sizes for conditional variance.
# N_v: MC sample size for Var(Y).
# Returns: length-d vector of Shapley Effects (normalized to sum to 1).
shapley_perm <- function(f, d, R, m = 300, N_o = 10, N_i = 3, N_v = 5000,
                          progress = TRUE) {
  # Marginal variance estimate
  U_v <- rcopula_unit(N_v, R)
  Y_v <- apply(U_v, 1, f)
  varY <- var(Y_v)
  EY   <- mean(Y_v)
  if (progress) cat(sprintf("  Var(Y) = %.3e  (n=%d)\n", varY, N_v))

  Sh <- numeric(d)

  for (k in seq_len(m)) {
    pi_k <- sample.int(d)
    # Accumulator for V[E[Y|X_{π[1:i]}]] — the "cost function" c(·)
    prev_c <- 0
    # The full set needed at step i is π[1:i]; marginal draws for outer loop.
    for (i in seq_len(d)) {
      S <- pi_k[seq_len(i)]
      free_i <- setdiff(seq_len(d), S)

      # Outer loop: N_o independent draws of X_S (from marginal = uniform)
      # BUT we need joint marginal consistent with copula — sample jointly
      # and keep only X_S.
      U_o <- rcopula_unit(N_o, R)[, S, drop = FALSE]

      EY_o <- numeric(N_o)
      for (o in seq_len(N_o)) {
        if (length(free_i) == 0) {
          # Fully conditioned — y is deterministic in f
          x_full <- numeric(d)
          x_full[S] <- U_o[o, ]
          EY_o[o] <- f(x_full)
        } else {
          # Inner loop: draw X_{free} | X_S = U_o[o, ] under the copula
          U_in <- rcopula_conditional(U_o[o, ], S, R, n = N_i)
          Y_in <- numeric(N_i)
          for (inn in seq_len(N_i)) {
            x_full <- numeric(d)
            x_full[S]      <- U_o[o, ]
            x_full[free_i] <- U_in[inn, ]
            Y_in[inn] <- f(x_full)
          }
          EY_o[o] <- mean(Y_in)
        }
      }
      c_i <- var(EY_o)

      # Last step: by construction V[E[Y|X_{all}]] = Var(Y)
      if (i == d) c_i <- varY

      j <- pi_k[i]
      Sh[j] <- Sh[j] + (c_i - prev_c)
      prev_c <- c_i
    }
    if (progress && (k %% max(1, m %/% 20) == 0)) {
      cat(sprintf("    perm %d/%d  running est: [%s]\n", k, m,
                   paste(sprintf("%.3f", Sh / (k * varY)), collapse = ",")))
    }
  }

  Sh / (m * varY)
}


# --- Correlation structure for our 8-D problem ---------------------------
# Parameter order matches `free_l1 ++ l2`:
#   1: error_reduction   2: cost_reduction   3: cycle_compression
#   4: tool_fixed_cost    5: total_budget     6: bean_price
#   7: discount_rate     8: adoption_ceiling
shapley_param_names <- c("error_reduction", "cost_reduction", "cycle_compression",
                          "tool_fixed_cost", "total_budget",
                          "bean_price", "discount_rate", "adoption_ceiling")

default_correlation_pairs <- list(
  list(i = 7, j = 6, rho =  0.40),  # discount_rate ↔ bean_price
  list(i = 8, j = 5, rho =  0.30),  # adoption_ceiling ↔ total_budget
  list(i = 3, j = 1, rho =  0.30),  # cycle_compression ↔ error_reduction
  list(i = 2, j = 4, rho = -0.20)   # cost_reduction ↔ tool_fixed_cost
)


# --- Main runner ---------------------------------------------------------
run_shapley <- function(m = 300, N_o = 10, N_i = 3, N_v = 5000,
                         output_dir = NULL,
                         correlation_pairs = default_correlation_pairs) {
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  sobol_cfg <- load_econ_grid()$sobol_params
  econ_cfg  <- load_config("economic_params")
  breed_cfg <- load_config("breeding_params")
  all_l1    <- c("error_reduction", "cost_reduction", "cycle_compression",
                  "tool_fixed_cost", "total_budget")

  l1_lookup <- build_l1_grid_lookup(free_dims = all_l1)
  f <- make_npv_fn(l1_lookup, sobol_cfg, econ_cfg, breed_cfg)

  d <- 8
  R_indep  <- diag(d)
  R_copula <- build_correlation(d, correlation_pairs)

  cat(sprintf("=== Shapley Effects  (d=%d, m=%d perms, N_o=%d, N_i=%d, N_v=%d) ===\n",
              d, m, N_o, N_i, N_v))
  cat("Correlation (copula) matrix:\n")
  rownames(R_copula) <- colnames(R_copula) <- shapley_param_names
  print(round(R_copula, 2))

  set.seed(2026)
  cat("\n--- [1/2] independent inputs ---\n")
  t0 <- Sys.time()
  sh_indep <- shapley_perm(f, d, R_indep, m = m, N_o = N_o, N_i = N_i, N_v = N_v)
  cat(sprintf("elapsed: %.2f min\n",
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  set.seed(2027)
  cat("\n--- [2/2] Gaussian-copula correlated inputs ---\n")
  t0 <- Sys.time()
  sh_cop <- shapley_perm(f, d, R_copula, m = m, N_o = N_o, N_i = N_i, N_v = N_v)
  cat(sprintf("elapsed: %.2f min\n",
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  group <- c(rep("L1 technical", 5), rep("L2 economic", 3))

  out_indep <- data.table(param = shapley_param_names,
                            shapley_effect = sh_indep, group = group)
  out_cop   <- data.table(param = shapley_param_names,
                            shapley_effect = sh_cop,   group = group)
  setorder(out_indep, -shapley_effect)
  setorder(out_cop,   -shapley_effect)

  cat("\n=== Shapley Effects — independent ===\n"); print(out_indep)
  cat(sprintf("sum = %.3f\n", sum(sh_indep)))
  cat("\n=== Shapley Effects — Gaussian copula ===\n"); print(out_cop)
  cat(sprintf("sum = %.3f\n", sum(sh_cop)))

  fwrite(out_indep, file.path(output_dir, "shapley_indices_indep.csv"))
  fwrite(out_cop,   file.path(output_dir, "shapley_indices_copula.csv"))

  # Also save the correlation matrix used so results are reproducible
  R_out <- as.data.table(R_copula, keep.rownames = "param")
  fwrite(R_out, file.path(output_dir, "shapley_correlation_matrix.csv"))

  invisible(list(indep = out_indep, copula = out_cop, R = R_copula))
}


if (sys.nframe() == 0) {
  run_shapley()
}
