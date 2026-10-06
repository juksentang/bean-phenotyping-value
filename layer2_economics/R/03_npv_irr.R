# ============================================================================
# Layer 2 — NPV, IRR, and Cost-Benefit Analysis
# ============================================================================

source(file.path(getwd(), "layer2_economics/R/01_gain_to_value.R"))
source(file.path(getwd(), "layer2_economics/R/02_diffusion_model.R"))

# --- NPV calculation ---

calc_npv <- function(benefits, costs, discount_rate) {
  # benefits, costs: vectors of annual values (year 0, 1, 2, ...)
  horizon <- length(benefits) - 1
  years <- 0:horizon
  net <- benefits - costs
  sum(net / (1 + discount_rate)^years)
}


# --- IRR calculation ---

calc_irr <- function(benefits, costs, lower = -0.5, upper = 2.0) {
  npv_fn <- function(r) calc_npv(benefits, costs, r)
  # Check that there's a sign change
  if (npv_fn(lower) * npv_fn(upper) > 0) return(NA_real_)
  tryCatch(
    uniroot(npv_fn, interval = c(lower, upper), tol = 1e-6)$root,
    error = function(e) NA_real_
  )
}


# --- Full economic evaluation — scalar-argument core ---
# Pure function: all parameters passed in, no config side-effects.
# Returns a list (same keys as the old evaluate_one_rep). The year loop,
# NPV, IRR, BCR, and cost_per_dg logic below is unchanged — only the
# parameter source was hoisted from econ_cfg/master_cfg into function args
# so that the sweep driver can vectorize over the economic parameter grid.

eval_npv_core <- function(dg_trad, dg_tool, cycle_trad, cycle_tool,
                           price, r, area, A_max, k, t_mid,
                           horizon, total_budget,
                           compute_irr = TRUE) {
  # Annual genetic gain (kg/ha/yr)
  ann_dg_trad <- dg_trad / cycle_trad
  ann_dg_tool  <- dg_tool  / cycle_tool

  # Release years (traditional releases at cycle_trad, tool at cycle_tool)
  release_trad <- cycle_trad
  release_tool  <- cycle_tool

  years <- 0:horizon
  benefits_trad <- numeric(horizon + 1)
  benefits_tool  <- numeric(horizon + 1)
  costs_trad    <- numeric(horizon + 1)
  costs_tool     <- numeric(horizon + 1)

  # Breeding costs: incurred during the breeding cycle
  cost_per_yr_trad <- total_budget / cycle_trad
  cost_per_yr_tool  <- total_budget / cycle_tool

  for (i in seq_along(years)) {
    t <- years[i]

    # Costs: breeding program runs for the cycle duration
    if (t < cycle_trad) costs_trad[i] <- cost_per_yr_trad
    if (t < cycle_tool)  costs_tool[i]  <- cost_per_yr_tool

    # Benefits: after release, cumulative genetic gain accrues with adoption
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

  # Incremental analysis: tool vs Traditional
  incr_benefits <- benefits_tool - benefits_trad
  incr_costs    <- costs_tool - costs_trad

  npv_incr <- calc_npv(incr_benefits, incr_costs, r)
  irr_incr <- if (compute_irr) calc_irr(incr_benefits, incr_costs) else NA_real_

  # Benefit-cost ratio (on incremental flows)
  pv_benefits <- sum(pmax(incr_benefits, 0) / (1 + r)^years)
  pv_costs    <- sum(pmax(-pmin(incr_benefits - incr_costs, 0), 0) / (1 + r)^years)
  bcr <- if (pv_costs > 0) pv_benefits / pv_costs else Inf

  # Cost per unit genetic gain
  total_incr_cost <- sum(pmax(incr_costs, 0))
  incr_dg <- ann_dg_tool - ann_dg_trad
  cost_per_dg <- if (incr_dg > 0) total_incr_cost / incr_dg else NA_real_

  list(
    npv                 = npv_incr,
    irr                 = irr_incr,
    bcr                 = bcr,
    cost_per_unit_dg    = cost_per_dg,
    ann_dg_trad         = ann_dg_trad,
    ann_dg_tool          = ann_dg_tool,
    incr_ann_dg         = incr_dg,
    benefits_trad       = benefits_trad,
    benefits_tool        = benefits_tool,
    costs_trad          = costs_trad,
    costs_tool           = costs_tool
  )
}


# --- Backward-compatible wrapper: takes econ_cfg/master_cfg ---

evaluate_one_rep <- function(dg_trad, dg_tool, cycle_trad, cycle_tool,
                              econ_cfg, master_cfg) {
  eval_npv_core(
    dg_trad      = dg_trad,
    dg_tool       = dg_tool,
    cycle_trad   = cycle_trad,
    cycle_tool    = cycle_tool,
    price        = econ_cfg$bean_price_usd_per_kg,
    r            = econ_cfg$discount_rate,
    area         = econ_cfg$total_area_ha,
    A_max        = econ_cfg$adoption_ceiling,
    k            = econ_cfg$diffusion_rate_k,
    t_mid        = econ_cfg$diffusion_midpoint_years,
    horizon      = econ_cfg$evaluation_horizon,
    total_budget = master_cfg$total_budget
  )
}


# --- Evaluate all MC reps ---

evaluate_all_reps <- function(mc_results, econ_cfg = NULL, master_cfg = NULL) {
  if (is.null(econ_cfg))   econ_cfg   <- load_config("economic_params")
  if (is.null(master_cfg)) master_cfg <- load_config("master_config")

  reps <- unique(mc_results$rep)
  results <- data.table(
    rep = integer(), npv = numeric(), irr = numeric(),
    bcr = numeric(), cost_per_unit_dg = numeric(),
    ann_dg_trad = numeric(), ann_dg_tool = numeric(), incr_ann_dg = numeric()
  )

  for (r in reps) {
    trad <- mc_results[rep == r & scenario == "traditional"]
    tool  <- mc_results[rep == r & scenario == "tool"]
    if (nrow(trad) == 0 || nrow(tool) == 0) next

    ev <- evaluate_one_rep(
      dg_trad = trad$delta_g, dg_tool = tool$delta_g,
      cycle_trad = trad$cycle_years, cycle_tool = tool$cycle_years,
      econ_cfg = econ_cfg, master_cfg = master_cfg
    )

    results <- rbindlist(list(results, data.table(
      rep = r, npv = ev$npv, irr = ev$irr,
      bcr = ev$bcr, cost_per_unit_dg = ev$cost_per_unit_dg,
      ann_dg_trad = ev$ann_dg_trad, ann_dg_tool = ev$ann_dg_tool,
      incr_ann_dg = ev$incr_ann_dg
    )))
  }

  # Summary
  summary_dt <- results[, .(
    mean_npv     = mean(npv, na.rm = TRUE),
    sd_npv       = sd(npv, na.rm = TRUE),
    median_npv   = median(npv, na.rm = TRUE),
    q025_npv     = quantile(npv, 0.025, na.rm = TRUE),
    q975_npv     = quantile(npv, 0.975, na.rm = TRUE),
    mean_irr     = mean(irr, na.rm = TRUE),
    mean_bcr     = mean(bcr[is.finite(bcr)], na.rm = TRUE),
    prob_positive = mean(npv > 0, na.rm = TRUE),
    mean_incr_dg = mean(incr_ann_dg, na.rm = TRUE)
  )]

  list(results = results, summary = summary_dt)
}
