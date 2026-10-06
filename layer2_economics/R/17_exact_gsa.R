# ============================================================================
# Layer 2 — Exact global sensitivity and value-of-information analysis on the
#           simulated grid (replaces the Monte Carlo estimators of 09/11/12/13/15)
# ============================================================================
# Why this script exists (revision of the manuscript, review findings on the
# attribution and decision analyses):
#   * The Layer-1 profiles form a FULL FACTORIAL grid (6 error x 8 cost x 4
#     compression x 6 fixed cost x 5 budget = 5,760 profiles) and the economic
#     inputs have five levels each. With area fixed at its baseline the whole
#     model is a table of 6*8*4*6*5*5*5*5 = 720,000 NPV values, so every
#     conditional mean E[Y | X_S] is computed EXACTLY by marginalisation.
#   * Sampling uniformly on the explicit grid levels gives one prior shared by
#     the attribution and the decision analysis (the old nearest-grid mapping of
#     continuous uniforms gave the end levels half the weight of the interior
#     levels, so the Sobol prior differed from the equal weights of the EVPPI).
#   * Sobol first/total-order indices, Shapley effects (independent inputs and a
#     Gaussian copula discretised on the grid), PAWN and horizon-dependent total
#     indices are therefore free of Monte Carlo error and satisfy
#     S1 <= Shapley <= ST exactly for independent inputs. The only remaining
#     uncertainty is the Layer-1 simulation noise in the profile gains, which a
#     parametric bootstrap (bootstrap_layer1_noise) propagates for the Sobol indices.
#   * Area is fixed at its baseline in the attribution AND in the decision
#     analysis (the old EVPPI averaged over area).
#
# NPV is separable: NPV = p * A * area * (G_o S(L_o, r) - G_t S(L_t, r))
#                         - B * (C(L_o, r) - C(L_t, r)),
# with S the discounted adoption-and-ramp factor and C the discounted spending
# factor (npv_factors). npv_factors reproduces eval_npv_core (03_npv_irr.R) to
# machine precision; check_against_lookup() verifies this against the stored
# lookup table.
#
# Outputs (layer2_economics/outputs_v2): sobol_indices_full.csv,
# sobol_indices_cyc0.csv, shapley_indices_indep.csv, shapley_indices_copula.csv,
# shapley_correlation_matrix.csv, pawn_indices.csv, dynamic_sensitivity_*.csv,
# evppi.csv, sobol_bootstrap_layer1.csv

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(mvtnorm)
})
if (!exists("load_config", mode = "function"))
  source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

# ---------------------------------------------------------------------------
# NPV building blocks
# ---------------------------------------------------------------------------
# Discounted adoption-and-ramp factor of one program (per unit of per-cycle gain,
# price, ceiling and area) and discounted spending factor (per unit of budget).
#   ramp = "paper"     : min(t - L + 1, L)/L, the published model (gain builds over L years, then holds)
#   ramp = "step"      : full per-cycle gain from release
#   ramp = "recurrent" : (t - L + 1)/L without the cap (cycles recur; spending continues every year)
# dt is the time step in years (1 = published model; 0.5 = seasonal steps).
npv_factors <- function(L, r, k, mid, H, ramp = "paper", dt = 1) {
  tt <- seq(0, H, by = dt)
  tt <- tt[tt <= H + 1e-9]
  disc <- (1 + r)^(-tt)
  after <- tt >= L - 1e-9
  tsince <- tt - L
  rampf <- switch(ramp,
    paper     = pmin(tsince + 1, L) / L,
    step      = rep(1, length(tt)),
    recurrent = (tsince + 1) / L,
    stop("unknown ramp ", ramp))
  logis <- 1 / (1 + exp(-k * (tsince - mid)))
  S <- sum(dt * after * rampf * logis * disc)
  cost_t <- if (ramp == "recurrent") rep(1, length(tt)) else as.numeric(tt < L - 1e-9)
  C <- sum(dt * cost_t * disc) / L
  c(S = S, C = C)
}

# Vectorised over discount rates (r) and cycle lengths (L); returns arrays S[L, r], C[L, r]
npv_factor_tables <- function(Ls, rs, k, mid, H, ramp = "paper", dt = 1) {
  S <- C <- matrix(NA_real_, length(Ls), length(rs), dimnames = list(as.character(Ls), as.character(rs)))
  for (i in seq_along(Ls)) for (j in seq_along(rs)) {
    f <- npv_factors(Ls[i], rs[j], k, mid, H, ramp, dt)
    S[i, j] <- f["S"]; C[i, j] <- f["C"]
  }
  list(S = S, C = C)
}

# ---------------------------------------------------------------------------
# Profile table (Layer 1) and economic grid
# ---------------------------------------------------------------------------
TECH <- c("error_reduction", "cost_reduction", "cycle_compression", "tool_fixed_cost", "total_budget")
ECON <- c("bean_price", "discount_rate", "adoption_ceiling")
PARAMS8 <- c(TECH, ECON)

# Layer-1 cells with the common calibration factor for a given target (% of mean yield per year)
load_cells <- function(target = NULL, extra_paths = NULL) {
  econ  <- load_config("economic_params"); breed <- load_config("breeding_params")
  if (is.null(target)) target <- econ$gain_target_pct_per_yr
  cells <- rbindlist(lapply(l1_summary_path(), fread), use.names = TRUE, fill = TRUE)
  k <- if (is.null(target)) 1 else
    target / 100 * econ$mean_yield_kgha * breed$cycle_years / mean(cells$mean_dg_trad)
  for (col in c("mean_dg_trad", "mean_dg_tool", "mean_dg_diff", "sd_dg_diff")) set(cells, j = col, value = cells[[col]] * k)
  attr(cells, "gain_scale") <- k
  cells
}

grid_levels <- function() {
  g <- yaml::read_yaml(file.path(getwd(), "config/econ_sweep_grid.yaml"))
  list(price = g$econ_grid$bean_price_usd_per_kg, disc = g$econ_grid$discount_rate,
       adopt = g$econ_grid$adoption_ceiling, area = g$econ_grid$total_area_ha,
       base_area = g$sobol_params$fixed_area_ha)
}

# Exact NPV array over the 8 inputs (area fixed). dims follow PARAMS8 order.
# opts: k, mid, H, ramp, dt, cycle_years; comp_levels optionally overrides the compression
# levels (then cells must supply gains for every level: gains are pooled over the compression
# cells of a profile, valid because compression leaves the simulated selection unchanged).
build_npv_array <- function(cells, opts = list(), comp_levels = NULL, area = NULL,
                            running_cost = 0, pool_comp = FALSE) {
  econ <- load_config("economic_params"); breed <- load_config("breeding_params")
  gl <- grid_levels()
  o <- modifyList(list(k = econ$diffusion_rate_k, mid = econ$diffusion_midpoint_years,
                       H = econ$evaluation_horizon, ramp = "paper", dt = 1,
                       L = breed$cycle_years), opts)
  if (is.null(area)) area <- gl$base_area
  lev <- lapply(TECH, function(p) sort(unique(cells[[p]]))); names(lev) <- TECH
  if (pool_comp || !is.null(comp_levels)) {
    # pool gains over the compression cells of each profile
    key <- c("error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget")
    pooled <- cells[, .(mean_dg_trad = mean(mean_dg_trad), mean_dg_tool = mean(mean_dg_tool)), by = key]
    if (!is.null(comp_levels)) { lev$cycle_compression <- sort(comp_levels) }
    cc <- CJ(error_reduction = lev$error_reduction, cost_reduction = lev$cost_reduction,
             cycle_compression = lev$cycle_compression, tool_fixed_cost = lev$tool_fixed_cost,
             total_budget = lev$total_budget)
    cc <- merge(cc, pooled, by = key, all.x = TRUE, sort = FALSE)
  } else {
    cc <- CJ(error_reduction = lev$error_reduction, cost_reduction = lev$cost_reduction,
             cycle_compression = lev$cycle_compression, tool_fixed_cost = lev$tool_fixed_cost,
             total_budget = lev$total_budget)
    cc <- merge(cc, cells[, c(TECH, "mean_dg_trad", "mean_dg_tool"), with = FALSE], by = TECH, all.x = TRUE, sort = FALSE)
  }
  stopifnot(!anyNA(cc$mean_dg_trad))
  setorderv(cc, rev(TECH))   # first dim (error_reduction) varies fastest, matching array() order
  dimsT <- vapply(lev, length, 1L)
  Lo <- o$L - lev$cycle_compression
  ft <- npv_factor_tables(c(o$L, Lo), gl$disc, o$k, o$mid, o$H, o$ramp, o$dt)
  St <- ft$S[as.character(o$L), ]; Ct <- ft$C[as.character(o$L), ]
  So <- ft$S[as.character(Lo), , drop = FALSE]; Co <- ft$C[as.character(Lo), , drop = FALSE]
  # arrays over the tech dims
  Gt <- array(cc$mean_dg_trad, dimsT); Go <- array(cc$mean_dg_tool, dimsT)
  Bu <- array(rep(lev$total_budget, each = prod(dimsT[1:4])), dimsT)
  ci <- array(rep(rep(seq_along(lev$cycle_compression), each = prod(dimsT[1:2])), times = prod(dimsT[4:5])), dimsT)
  nP <- length(gl$price); nR <- length(gl$disc); nA <- length(gl$adopt)
  dims <- c(dimsT, nP, nR, nA)
  Y <- array(NA_real_, dims)
  # AF: annuity factor of a running cost paid in years 0..H
  for (ir in seq_len(nR)) {
    r <- gl$disc[ir]
    So_c <- array(So[ci, ir], dimsT); Co_c <- array(Co[ci, ir], dimsT)
    ben <- Go * So_c - Gt * St[ir]             # per unit of price x ceiling x area
    cst <- Bu * (Co_c - Ct[ir])                # incremental spending PV
    runPV <- running_cost * sum((1 + r)^(-seq(0, o$H, by = o$dt)) * o$dt)
    for (ip in seq_len(nP)) for (ia in seq_len(nA)) {
      Y[, , , , , ip, ir, ia] <- gl$price[ip] * gl$adopt[ia] * area * ben - cst - runPV
    }
  }
  dimnames(Y) <- c(lev, list(bean_price = gl$price, discount_rate = gl$disc, adoption_ceiling = gl$adopt))
  names(dimnames(Y)) <- PARAMS8
  attr(Y, "opts") <- o
  Y
}

# ---------------------------------------------------------------------------
# Exact variance-based measures on a (weighted) array
# ---------------------------------------------------------------------------
# W: array of probabilities (sum 1), same dims as Y. Returns the vector c(S) = Var(E[Y | X_S])
# for every subset S, indexed by bitmask (S = {i : bit i-1 set}); c[1] is the empty set.
all_subset_variances <- function(Y, W = NULL) {
  d <- length(dim(Y)); full <- 2^d - 1
  if (is.null(W)) W <- array(1 / length(Y), dim(Y))
  mu <- sum(W * Y)
  N <- vector("list", full + 1); P <- vector("list", full + 1)
  N[[full + 1]] <- W * Y; P[[full + 1]] <- W
  inc <- function(m) which(bitwAnd(m, 2^(0:(d - 1))) > 0)
  cvec <- numeric(full + 1)
  for (m in full:0) {
    if (m < full) {
      miss <- setdiff(seq_len(d), inc(m)); i <- miss[1]
      parent <- bitwOr(m, 2^(i - 1))
      pdims <- inc(parent); pos <- match(i, pdims)
      keep <- setdiff(seq_along(pdims), pos)
      if (length(keep) == 0) {
        N[[m + 1]] <- sum(N[[parent + 1]]); P[[m + 1]] <- sum(P[[parent + 1]])
      } else {
        N[[m + 1]] <- apply(N[[parent + 1]], keep, sum); P[[m + 1]] <- apply(P[[parent + 1]], keep, sum)
      }
    }
    Pm <- P[[m + 1]]; Nm <- N[[m + 1]]
    cvec[m + 1] <- if (m == 0) 0 else sum(Pm * (Nm / Pm - mu)^2)
  }
  attr(cvec, "d") <- d; attr(cvec, "varY") <- cvec[full + 1]
  cvec
}

sobol_from_c <- function(cvec) {
  d <- attr(cvec, "d"); full <- 2^d - 1; V <- cvec[full + 1]
  S1 <- vapply(seq_len(d), function(i) cvec[2^(i - 1) + 1] / V, 0)
  ST <- vapply(seq_len(d), function(i) 1 - cvec[full - 2^(i - 1) + 1] / V, 0)
  list(S1 = S1, ST = ST, varY = V)
}

shapley_from_c <- function(cvec) {
  d <- attr(cvec, "d"); full <- 2^d - 1; V <- cvec[full + 1]
  pc <- vapply(0:full, function(m) sum(bitwAnd(m, 2^(0:(d - 1))) > 0), 0)
  Sh <- numeric(d)
  for (i in seq_len(d)) {
    bit <- 2^(i - 1)
    for (m in 0:full) {
      if (bitwAnd(m, bit) > 0) next
      s <- pc[m + 1]
      Sh[i] <- Sh[i] + factorial(s) * factorial(d - s - 1) / factorial(d) * (cvec[m + bit + 1] - cvec[m + 1])
    }
  }
  Sh / V
}

# Discretised Gaussian copula: probabilities of the grid cells under a Gaussian copula with
# the given pair correlations (disjoint pairs, so the joint law is a product of bivariate laws).
# Marginals are uniform over the levels (equal-width bins in the unit cube).
pair_pmf <- function(n1, n2, rho) {
  c1 <- qnorm(seq(0, 1, length.out = n1 + 1)); c2 <- qnorm(seq(0, 1, length.out = n2 + 1))
  sig <- matrix(c(1, rho, rho, 1), 2)
  cdf <- function(a, b) if (is.infinite(a) && a < 0 || is.infinite(b) && b < 0) 0 else
    if (is.infinite(a) && is.infinite(b)) 1 else if (is.infinite(a)) pnorm(b) else if (is.infinite(b)) pnorm(a) else
      pmvnorm(upper = c(a, b), corr = sig)[1]
  F <- outer(seq_along(c1), seq_along(c2), Vectorize(function(i, j) cdf(c1[i], c2[j])))
  p <- F[-1, -1] - F[-nrow(F), -1] - F[-1, -ncol(F)] + F[-nrow(F), -ncol(F)]
  p / sum(p)
}

copula_weights <- function(dims, pairs) {
  # dims: named integer vector over PARAMS8; pairs: list(list(a=, b=, rho=)) on disjoint inputs
  W <- array(1, dims)
  idx <- as.matrix(expand.grid(lapply(dims, seq_len)))
  used <- character()
  for (p in pairs) {
    ia <- match(p$a, names(dims)); ib <- match(p$b, names(dims))
    pm <- pair_pmf(dims[ia], dims[ib], p$rho)
    W <- W * array(pm[cbind(idx[, ia], idx[, ib])], dims)
    used <- c(used, p$a, p$b)
  }
  for (nm in setdiff(names(dims), used)) W <- W / dims[[nm]]
  W / sum(W)
}

default_copula_pairs <- list(
  list(a = "discount_rate",     b = "bean_price",   rho =  0.40),
  list(a = "adoption_ceiling",  b = "total_budget", rho =  0.30),
  list(a = "cycle_compression", b = "error_reduction", rho = 0.30),
  list(a = "cost_reduction",    b = "tool_fixed_cost", rho = -0.20))

# PAWN on discrete inputs: median (and max) over input levels of the Kolmogorov-Smirnov distance
# between the conditional and the unconditional distribution of Y (independent equal-weight prior)
pawn_exact <- function(Y) {
  y <- as.vector(Y); n <- length(y); o <- order(y); ys <- y[o]
  dims <- dim(Y); res <- list()
  Fall <- seq_len(n) / n
  # collapse ties so the CDF is evaluated at distinct values
  last <- c(ys[-1] != ys[-n], TRUE)
  idxarr <- array(seq_len(n), dims)
  for (i in seq_along(dims)) {
    ks <- numeric(dims[i])
    for (l in seq_len(dims[i])) {
      sel <- rep(FALSE, n)
      ind <- slice.index(Y, i) == l
      sel <- as.vector(ind)[o]
      Fc <- cumsum(sel) / sum(sel)
      ks[l] <- max(abs(Fc[last] - Fall[last]))
    }
    res[[i]] <- c(median = median(ks), max = max(ks), mean = mean(ks))
  }
  names(res) <- names(dimnames(Y))
  res
}

# ---------------------------------------------------------------------------
# Horizon-dependent total-order indices (NPV truncated at year h)
# ---------------------------------------------------------------------------
dynamic_indices <- function(cells, opts = list()) {
  econ <- load_config("economic_params")
  H <- econ$evaluation_horizon
  out_ST <- out_S1 <- list()
  for (h in 0:H) {
    Yh <- build_npv_array(cells, modifyList(opts, list(H = h)))
    cv <- all_subset_variances(Yh)
    sb <- sobol_from_c(cv)
    out_ST[[length(out_ST) + 1]] <- data.table(year = h, t(setNames(sb$ST, PARAMS8)))
    out_S1[[length(out_S1) + 1]] <- data.table(year = h, t(setNames(sb$S1, PARAMS8)))
  }
  list(ST = rbindlist(out_ST), S1 = rbindlist(out_S1))
}

# ---------------------------------------------------------------------------
# EVPPI on the grid (decision: invest in the tool, payoff NPV, versus a fixed payoff tau)
# ---------------------------------------------------------------------------
evppi_array <- function(Y, tau, W = NULL) {
  d <- length(dim(Y)); nm <- names(dimnames(Y))
  if (is.null(W)) W <- array(1 / length(Y), dim(Y))
  prior <- max(sum(W * Y), tau)
  rbindlist(lapply(seq_len(d), function(i) {
    Nm <- apply(W * Y, i, sum); Pm <- apply(W, i, sum)
    ev <- sum(Pm * pmax(Nm / Pm, tau))
    data.table(param = nm[i], threshold = tau, EV_prior = prior, EV_perfect = ev, EVPPI = ev - prior)
  }))
}

# Layer-1 noise bootstrap for the Sobol indices: each profile's tool gain is perturbed by
# N(0, SE_diff^2) (SE of the paired per-cycle difference; manual gain held fixed) and the exact
# indices are recomputed.
bootstrap_layer1_noise <- function(cells, B = 200, seed = 2026) {
  set.seed(seed)
  res <- vector("list", B)
  se <- cells$sd_dg_diff / sqrt(cells$n_reps)
  base <- copy(cells)
  for (b in seq_len(B)) {
    cb <- copy(base)
    cb[, mean_dg_tool := mean_dg_tool + rnorm(.N) * se]
    Y <- build_npv_array(cb)
    sb <- sobol_from_c(all_subset_variances(Y))
    res[[b]] <- data.table(b = b, param = PARAMS8, S1 = sb$S1, ST = sb$ST)
  }
  rbindlist(res)
}

check_against_lookup <- function(cells) {
  lk <- as.data.table(arrow::read_parquet(file.path(l2_out_dir(), "lookup_table.parquet")))
  gl <- grid_levels(); econ <- load_config("economic_params")
  Y <- build_npv_array(cells, area = gl$base_area)
  sub <- lk[total_area_ha == gl$base_area]
  sub <- sub[sample(.N, 20000)]
  idx <- cbind(match(sub$error_reduction, dimnames(Y)$error_reduction), match(sub$cost_reduction, dimnames(Y)$cost_reduction),
               match(sub$cycle_compression, dimnames(Y)$cycle_compression), match(sub$tool_fixed_cost, dimnames(Y)$tool_fixed_cost),
               match(sub$total_budget, dimnames(Y)$total_budget), match(sub$bean_price, gl$price),
               match(sub$discount_rate, gl$disc), match(sub$adoption_ceiling, gl$adopt))
  stopifnot(!anyNA(idx))
  max(abs(Y[idx] - sub$npv) / pmax(1, abs(sub$npv)))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_exact_gsa <- function(B = 200) {
  out_dir <- l2_out_dir()
  cells <- load_cells()
  cat(sprintf("Cells: %d, gain scale %.4f\n", nrow(cells), attr(cells, "gain_scale")))
  set.seed(1)
  cat(sprintf("Max relative deviation from stored lookup table: %.2e\n", check_against_lookup(cells)))
  Y <- build_npv_array(cells)
  cat("NPV array dims:", dim(Y), " mean NPV (M USD):", mean(Y) / 1e6, "\n")
  grp <- c(rep("L1 technical", 5), rep("L2 economic", 3))

  # --- Sobol (exact, independent grid prior) ---
  cv <- all_subset_variances(Y); sb <- sobol_from_c(cv)
  sob <- data.table(param = PARAMS8, S1 = sb$S1, ST = sb$ST, group = grp); setorder(sob, -ST)
  fwrite(sob, file.path(out_dir, "sobol_indices_full.csv")); print(sob)
  # --- Sobol with compression fixed at zero ---
  Y0 <- Y[, , 1, , , , , , drop = FALSE]
  dim(Y0) <- dim(Y)[-3]; dimnames(Y0) <- dimnames(Y)[-3]
  cv0 <- all_subset_variances(Y0); sb0 <- sobol_from_c(cv0)
  p0 <- setdiff(PARAMS8, "cycle_compression")
  sob0 <- data.table(param = p0, S1 = sb0$S1, ST = sb0$ST, group = c(rep("L1 technical", 4), rep("L2 economic", 3))); setorder(sob0, -ST)
  fwrite(sob0, file.path(out_dir, "sobol_indices_cyc0.csv")); print(sob0)

  # --- Shapley (exact): independent and discretised Gaussian copula ---
  sh_i <- shapley_from_c(cv)
  shi <- data.table(param = PARAMS8, shapley_effect = sh_i, group = grp,
                    S1 = sb$S1, ST = sb$ST, bracket_ok = sh_i >= sb$S1 - 1e-12 & sh_i <= sb$ST + 1e-12)
  setorder(shi, -shapley_effect); fwrite(shi, file.path(out_dir, "shapley_indices_indep.csv")); print(shi)
  dims <- setNames(dim(Y), PARAMS8)
  W <- copula_weights(dims, default_copula_pairs)
  cvc <- all_subset_variances(Y, W); sh_c <- shapley_from_c(cvc)
  shc <- data.table(param = PARAMS8, shapley_effect = sh_c, group = grp); setorder(shc, -shapley_effect)
  fwrite(shc, file.path(out_dir, "shapley_indices_copula.csv")); print(shc)
  Rm <- diag(8); dimnames(Rm) <- list(PARAMS8, PARAMS8)
  for (p in default_copula_pairs) { Rm[p$a, p$b] <- p$rho; Rm[p$b, p$a] <- p$rho }
  fwrite(as.data.table(Rm, keep.rownames = "param"), file.path(out_dir, "shapley_correlation_matrix.csv"))

  # --- PAWN (exact) ---
  pw <- pawn_exact(Y)
  pawn <- data.table(param = names(pw), PAWN_median = vapply(pw, `[`, 0, "median"),
                     PAWN_max = vapply(pw, `[`, 0, "max"), PAWN_mean = vapply(pw, `[`, 0, "mean"))
  setorder(pawn, -PAWN_median); fwrite(pawn, file.path(out_dir, "pawn_indices.csv")); print(pawn)

  # --- Dynamic (horizon-dependent) indices ---
  dyn <- dynamic_indices(cells)
  setnames(dyn$ST, c("year", PARAMS8)); setnames(dyn$S1, c("year", PARAMS8))
  fwrite(dyn$ST, file.path(out_dir, "dynamic_sensitivity_ST_wide.csv"))
  fwrite(dyn$S1, file.path(out_dir, "dynamic_sensitivity_S1_wide.csv"))
  print(dyn$ST[year %in% c(0, 5, 7, 10, 25)])

  # --- EVPPI ---
  mu <- mean(Y)
  hurdle_mult <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3)
  ev <- rbindlist(lapply(hurdle_mult, function(hm) { r <- evppi_array(Y, mu * hm); r[, hurdle_mult := hm]; r }))
  ev[, `:=`(EVPPI_M = EVPPI / 1e6, threshold_M = threshold / 1e6, at_indifference = abs(hurdle_mult - 1) < 1e-9)]
  setorder(ev, -EVPPI); fwrite(ev, file.path(out_dir, "evppi.csv")); print(head(ev[, .(param, threshold_M, EVPPI_M)], 12))

  # --- Layer-1 noise bootstrap of the Sobol indices ---
  cat(sprintf("Layer-1 noise bootstrap (B = %d)...\n", B))
  bs <- bootstrap_layer1_noise(cells, B)
  bsum <- bs[, .(S1_lo = quantile(S1, .025), S1_hi = quantile(S1, .975),
                 ST_lo = quantile(ST, .025), ST_hi = quantile(ST, .975)), by = param]
  fwrite(bsum, file.path(out_dir, "sobol_bootstrap_layer1.csv")); print(bsum)
  writeLines(sprintf("exact grid analysis: levels per input = %s; n = %d; layer1_bootstrap B=%d; copula pairs disjoint",
                     paste(dim(Y), collapse = "x"), length(Y), B), file.path(out_dir, "exact_gsa_run.log"))
  invisible(NULL)
}

if (sys.nframe() == 0) run_exact_gsa()
