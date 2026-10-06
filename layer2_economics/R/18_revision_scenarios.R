# ============================================================================
# Layer 2 — Robustness of the attribution and decision results to the modeling
#           choices that favor cycle compression (revision analyses)
# ============================================================================
# Every scenario re-evaluates the exact NPV table of 17_exact_gsa.R (no new Layer 1
# simulation) under one changed assumption and reports the same summary measures:
#   * compression priors  : which levels of cycle compression the prior covers
#                           (grid levels 0-3; levels 0 and 1; mass at zero; seasonal steps 0, 0.5, 1)
#   * structural variants : gain ramp (published ramp, no ramp, recurrent cycles),
#                           horizon of 40 years, adoption rate k = 0.10 and 0.15
#   * calibration targets : 0.25 and 1.0 percent of mean yield per year (published: 0.5)
#   * tool running cost   : an incremental annual cost over the horizon (net value of the tool)
#   * stress block        : a tool noisier than manual scoring (error reduction below zero)
#
# Outputs: layer2_economics/outputs_v2/scenario_summary.csv (long format),
#          stress_block_summary.csv

source(file.path(getwd(), "layer2_economics/R/17_exact_gsa.R"))

econ0  <- load_config("economic_params")
scn    <- load_config("revision_scenarios")
gl     <- grid_levels()

# --- summary measures of one NPV array (optionally with a weight array) ---
scenario_measures <- function(Y, W = NULL, hurdle0 = FALSE) {
  d <- length(dim(Y)); nm <- names(dimnames(Y))
  if (is.null(W)) W <- array(1 / length(Y), dim(Y))
  cv <- all_subset_variances(Y, W); sb <- sobol_from_c(cv)
  mu <- sum(W * Y)
  ev <- evppi_array(Y, mu, W)                         # indifference hurdle E[NPV]
  out <- c(EV = mu / 1e6, PrPos = 100 * sum(W * (Y > 0)),
           setNames(sb$ST, paste0("ST_", nm)), setNames(sb$S1, paste0("S1_", nm)),
           setNames(ev$EVPPI / 1e6, paste0("EVPPIind_", nm)))
  if (hurdle0) {
    ev0 <- evppi_array(Y, 0, W)
    out <- c(out, setNames(ev0$EVPPI / 1e6, paste0("EVPPIzero_", nm)))
  }
  out
}

# median NPV over the technical inputs at the baseline economic state, by compression level
baseline_state_idx <- function(Y) {
  econ <- econ0; dn <- dimnames(Y)
  nearest <- function(x, g) which.min(abs(g - x))
  c(price = nearest(econ$bean_price_usd_per_kg, gl$price), disc = nearest(econ$discount_rate, gl$disc),
    adopt = nearest(econ$adoption_ceiling, gl$adopt))
}
baseline_medians <- function(Y) {
  b <- baseline_state_idx(Y); levc <- as.numeric(dimnames(Y)$cycle_compression)
  setNames(vapply(seq_along(levc), function(i) median(Y[, , i, , , b["price"], b["disc"], b["adopt"]]) / 1e6, 0),
           paste0("NPVbase_c", levc))
}

res <- list()
add <- function(id, label, vals, extra = NULL) {
  res[[length(res) + 1]] <<- data.table(scenario = id, label = label, metric = names(vals), value = as.numeric(vals))
  if (!is.null(extra)) res[[length(res) + 1]] <<- data.table(scenario = id, label = label, metric = names(extra), value = as.numeric(extra))
}

cells <- load_cells()
Ybase <- build_npv_array(cells)
dimsB <- setNames(dim(Ybase), PARAMS8)

# --- 0. published structure, grid prior -------------------------------------
bm <- baseline_medians(Ybase)
add("Base", "Published structure, equal weight on the four compression levels",
    scenario_measures(Ybase, hurdle0 = TRUE),
    bm)

# --- 1. priors on cycle compression ------------------------------------------
prior_array <- function(Y, p) {                      # p: probabilities over the compression levels of Y
  keep <- which(p > 0); Yk <- Y[, , keep, , , , , , drop = FALSE]
  Wk <- array(1, dim(Yk)); nd <- length(dim(Yk))
  for (j in seq_along(keep)) Wk[, , j, , , , , ] <- p[keep[j]]
  Wk <- Wk / prod(dim(Yk)[-3]); Wk <- Wk / sum(Wk)
  dimnames(Wk) <- dimnames(Yk)
  list(Y = Yk, W = Wk)
}
pr <- prior_array(Ybase, c(.5, .5, 0, 0));       add("PriorZeroOne", "Compression 0 or 1 year, equal weight", scenario_measures(pr$Y, pr$W))
z1 <- scn$compression_zero_mass[[1]]; z2 <- scn$compression_zero_mass[[2]]
pr <- prior_array(Ybase, c(z1, rep((1 - z1) / 3, 3)));  add("PriorZeroMass", sprintf("Compression 0 with probability %.1f, 1 to 3 years share the rest", z1), scenario_measures(pr$Y, pr$W))
pr <- prior_array(Ybase, c(z2, 1 - z2, 0, 0));          add("PriorZeroHeavy", sprintf("Compression 0 with probability %.1f, 1 year with the rest", z2), scenario_measures(pr$Y, pr$W))
# seasonal steps: compression 0, 0.5 and 1 year with half-year time steps and gains pooled over compression cells
Ys <- build_npv_array(cells, opts = list(dt = scn$seasonal_time_step_years), comp_levels = unlist(scn$seasonal_compression_years))
add("PriorSeason", "Compression 0, 0.5 or 1 year (half-year steps), equal weight", scenario_measures(Ys))
Ys0 <- build_npv_array(cells, opts = list(dt = scn$seasonal_time_step_years), comp_levels = c(0, 1))
add("SeasonRef", "Compression 0 or 1 year with half-year steps (reference for the seasonal prior)", scenario_measures(Ys0))

# --- 2. structural variants ----------------------------------------------------
add("RampStep",  "No ramp: full per-cycle gain from release",       scenario_measures(build_npv_array(cells, list(ramp = "step"))),
    baseline_medians(build_npv_array(cells, list(ramp = "step"))))
add("Recurrent", "Recurrent cycles: gain accrues every cycle, spending continues", scenario_measures(build_npv_array(cells, list(ramp = "recurrent"))),
    baseline_medians(build_npv_array(cells, list(ramp = "recurrent"))))
add("HorizonForty", "Horizon of 40 years", scenario_measures(build_npv_array(cells, list(H = scn$horizon_long_years))),
    baseline_medians(build_npv_array(cells, list(H = scn$horizon_long_years))))
add("KTen",     "Adoption rate k = 0.10", scenario_measures(build_npv_array(cells, list(k = scn$adoption_rate_k[[1]]))),
    baseline_medians(build_npv_array(cells, list(k = scn$adoption_rate_k[[1]]))))
add("KFifteen", "Adoption rate k = 0.15", scenario_measures(build_npv_array(cells, list(k = scn$adoption_rate_k[[2]]))),
    baseline_medians(build_npv_array(cells, list(k = scn$adoption_rate_k[[2]]))))
Ycons <- build_npv_array(cells, list(ramp = "step", k = scn$adoption_rate_k[[1]], H = scn$horizon_long_years))
add("Conservative", "No ramp, k = 0.10, horizon 40 years", scenario_measures(Ycons), baseline_medians(Ycons))
pr <- prior_array(Ycons, c(.5, .5, 0, 0));        add("ConservativeZeroOne", "No ramp, k = 0.10, horizon 40 years, compression 0 or 1 year", scenario_measures(pr$Y, pr$W))

# --- 3. calibration targets ------------------------------------------------------
for (tg in unlist(scn$gain_targets_pct_per_yr)) {
  cl <- load_cells(target = tg); Yt <- build_npv_array(cl)
  id <- if (tg < econ0$gain_target_pct_per_yr) "TargetQuarter" else "TargetOne"
  add(id, sprintf("Calibration target %.2f percent per year (scale %.3f)", tg, attr(cl, "gain_scale")),
      scenario_measures(Yt), c(baseline_medians(Yt), GainScale = attr(cl, "gain_scale")))
}

# --- 4. tool running cost (net value) ----------------------------------------------
runcost <- econ0$tool_running_cost_scenarios_usd_per_yr
if (is.null(runcost)) runcost <- c(5e5, 1e6, 2e6)
for (cst in runcost) {
  Yn <- build_npv_array(cells, running_cost = cst)
  id <- paste0("RunCost", format(cst / 1e6, nsmall = 1))
  add(id, sprintf("Incremental running cost %.2f million USD per year", cst / 1e6),
      scenario_measures(Yn, hurdle0 = TRUE), c(RunCostM = cst / 1e6, baseline_medians(Yn)))
}

# --- 5. value of accuracy versus one year of compression in NPV terms ------------------
b <- baseline_state_idx(Ybase)
npv_err_c0 <- apply(Ybase[, , 1, , , b["price"], b["disc"], b["adopt"]], 1, median) / 1e6
npv_c1_err0 <- median(Ybase[1, , 2, , , b["price"], b["disc"], b["adopt"]]) / 1e6
add("ValueOfAccuracy", "Baseline state: median NPV by error-reduction level at zero compression, and at one year with no error reduction",
    c(setNames(npv_err_c0, paste0("NPV_err_", dimnames(Ybase)$error_reduction)), NPV_c1_err0 = npv_c1_err0))

dt <- rbindlist(res)
fwrite(dt, file.path(l2_out_dir(), "scenario_summary.csv"))
cat("Wrote scenario_summary.csv with", length(unique(dt$scenario)), "scenarios\n")
print(dcast(dt[metric %in% c("EV", "PrPos", "ST_cycle_compression", "ST_error_reduction", "S1_cycle_compression", "S1_error_reduction",
                             "EVPPIind_cycle_compression", "EVPPIind_error_reduction")],
            scenario ~ metric, value.var = "value")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

# --- 6. stress block (error reduction below zero) --------------------------------------
stress_path <- file.path(getwd(), "layer1_genetic_sim/outputs/sweep_v2_stress/cell_summaries.csv")
if (file.exists(stress_path)) {
  sc <- fread(stress_path); k <- attr(cells, "gain_scale")
  for (col in c("mean_dg_trad", "mean_dg_tool", "mean_dg_diff", "sd_dg_diff")) set(sc, j = col, value = sc[[col]] * k)
  o <- list(k = econ0$diffusion_rate_k, mid = econ0$diffusion_midpoint_years, H = econ0$evaluation_horizon,
            L = load_config("breeding_params")$cycle_years)
  rows <- list()
  for (ci in seq_len(nrow(sc))) for (cc in 0:3) {
    ft_t <- sapply(gl$disc, function(r) npv_factors(o$L, r, o$k, o$mid, o$H))
    ft_o <- sapply(gl$disc, function(r) npv_factors(o$L - cc, r, o$k, o$mid, o$H))
    for (ir in seq_along(gl$disc)) for (ip in seq_along(gl$price)) for (ia in seq_along(gl$adopt)) {
      npv <- gl$price[ip] * gl$adopt[ia] * gl$base_area * (sc$mean_dg_tool[ci] * ft_o["S", ir] - sc$mean_dg_trad[ci] * ft_t["S", ir]) -
        sc$total_budget[ci] * (ft_o["C", ir] - ft_t["C", ir])
      rows[[length(rows) + 1]] <- data.table(cell = ci, cost_reduction = sc$cost_reduction[ci], total_budget = sc$total_budget[ci],
                                              compression = cc, price = gl$price[ip], disc = gl$disc[ir], adopt = gl$adopt[ia], npv = npv)
    }
  }
  st <- rbindlist(rows)
  bsl <- st[price == gl$price[b["price"]] & disc == gl$disc[b["disc"]] & adopt == gl$adopt[b["adopt"]]]
  summ <- st[, .(medNPV = median(npv) / 1e6, PrPos = 100 * mean(npv > 0)), by = compression]
  summ_b <- bsl[, .(medNPV_base = median(npv) / 1e6), by = compression]
  upl <- sc[, .(uplift_pct = 100 * mean(mean_dg_diff) / mean(mean_dg_trad), diff_kgha = mean(mean_dg_diff), se_kgha = mean(sd_dg_diff / sqrt(n_reps)))]
  stress <- merge(summ, summ_b, by = "compression"); stress[, `:=`(uplift_pct = upl$uplift_pct, diff_kgha = upl$diff_kgha, se_kgha = upl$se_kgha)]
  fwrite(stress, file.path(l2_out_dir(), "stress_block_summary.csv")); print(stress)
} else cat("Stress block not available yet (", stress_path, ")\n")
