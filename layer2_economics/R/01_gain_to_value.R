# ============================================================================
# Layer 2 — Genetic Gain to Economic Value Conversion
# ============================================================================
# Converts genetic gain (kg/ha) to monetary value (USD).

source(file.path(getwd(), "layer1_genetic_sim/R/utils.R"))

# --- Convert genetic gain to USD ---

gain_to_value_per_ha <- function(delta_g_kgha, price_usd_per_kg) {
  # delta_g_kgha: genetic gain in kg/ha (from Layer 1)
  # Returns: incremental value in USD/ha
  delta_g_kgha * price_usd_per_kg
}

gain_to_total_value <- function(delta_g_kgha, price_usd_per_kg,
                                 adoption_rate, total_area_ha) {
  # Total value across all adopting area
  delta_g_kgha * price_usd_per_kg * adoption_rate * total_area_ha
}


# --- Annual genetic gain (per year, accounting for cycle length) ---

annual_gain <- function(delta_g_kgha, cycle_years) {
  delta_g_kgha / cycle_years
}


# --- Equivalent Annual Annuity (EAA) converter ---
# Converts a NPV (present value) back into a level annual stream of the
# same present value over `horizon` years at discount rate `r`.
#   EAA = NPV * r / (1 - (1+r)^(-horizon))
# Divide by `area` to get per-ha per-year incremental value, which is the
# only human-interpretable scale for this kind of infrastructure benefit.

eaa_from_npv <- function(npv, r, horizon) {
  # Vectorized in all args
  annuity_factor <- r / (1 - (1 + r)^(-horizon))
  npv * annuity_factor
}

eaa_per_ha_per_yr <- function(npv, r, horizon, area) {
  eaa_from_npv(npv, r, horizon) / area
}


# --- Process Layer 1 MC results into economic inputs ---

prepare_economic_inputs <- function(mc_results, econ_cfg = NULL) {
  if (is.null(econ_cfg)) econ_cfg <- load_config("economic_params")

  # mc_results: data.table with columns: rep, scenario, delta_g, cycle_years
  dt <- copy(mc_results)

  # Per-year genetic gain
  dt[, annual_dg := delta_g / cycle_years]

  # Per-ha value per year
  dt[, value_per_ha_yr := annual_dg * econ_cfg$bean_price_usd_per_kg]

  dt
}
