# ============================================================================
# Layer 2 — Logistic Variety Adoption Diffusion Model
# ============================================================================
# Models the adoption of improved varieties over time using a logistic curve.
# A(t) = A_max / (1 + exp(-k * (t - t_mid)))

# --- Logistic adoption function ---

adoption_rate <- function(t, A_max, k, t_mid) {
  # t: years since variety release
  # A_max: ceiling adoption rate (fraction of area)
  # k: diffusion speed (per year)
  # t_mid: years to 50% of ceiling adoption
  A_max / (1 + exp(-k * (t - t_mid)))
}


# --- Cumulative adoption over evaluation horizon ---

adoption_trajectory <- function(horizon, A_max, k, t_mid, release_year = 0) {
  # Returns a data.table with year and adoption rate
  years <- seq(0, horizon)
  t_since_release <- years - release_year
  # Before release, adoption = 0
  rates <- ifelse(t_since_release < 0, 0,
                   adoption_rate(t_since_release, A_max, k, t_mid))
  data.table(year = years, t_since_release = t_since_release,
             adoption = rates)
}


# --- Compare adoption trajectories for two scenarios ---

compare_adoption <- function(horizon, A_max, k, t_mid,
                              release_year_trad, release_year_tool) {
  trad <- adoption_trajectory(horizon, A_max, k, t_mid, release_year_trad)
  trad[, scenario := "traditional"]
  tool <- adoption_trajectory(horizon, A_max, k, t_mid, release_year_tool)
  tool[, scenario := "tool"]
  rbindlist(list(trad, tool))
}
