# ============================================================================
# Layer 2 — Decision-theoretic extension: Scenario Decomposition + VoI
# ============================================================================
# Computational prototype for the methods described in
#   report/Decision.md
#
# Two computations:
#   (1) Scenario decomposition of E[NPV] on the lookup table (pessimistic /
#       baseline / optimistic tercile combinations of the 4 L2 drivers).
#   (2) Expected Value of Perfect Partial Information (EVPPI) — per-parameter
#       "how much is it worth to learn the true value?" against a
#       "invest in the tool vs. status-quo" decision.
#
# Output:
#   outputs/scenario_decomposition.csv
#   outputs/evppi_lookup_all_areas.csv  (EVPPI averaged over all four economic inputs incl. area)

source(file.path(getwd(), "layer2_economics/R/00_paths.R"))

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

LOOKUP_PATH <- file.path(l2_out_dir(), "lookup_table.parquet")

load_lookup <- function() {
  stopifnot(file.exists(LOOKUP_PATH))
  dt <- setDT(as.data.frame(read_parquet(LOOKUP_PATH)))
  dt
}


# --- Scenario decomposition ---------------------------------------------
# Build 3 scenarios by tercile-filtering the 4 drivers:
#   pessimistic:  low bean_price, high discount_rate, low A_max, slow cycle_compression
#   baseline:     median of each
#   optimistic:   high price, low discount, high A_max, fast compression
scenario_decomposition <- function(dt) {
  # Tercile cut-points on each var
  q_price <- quantile(dt$bean_price, c(1/3, 2/3))
  q_rate  <- quantile(dt$discount_rate, c(1/3, 2/3))
  q_amax  <- quantile(dt$adoption_ceiling, c(1/3, 2/3))
  q_cyc   <- quantile(dt$cycle_compression, c(1/3, 2/3))

  label_scenario <- function(price, rate, amax, cyc) {
    pp <- ifelse(price <= q_price[1], "lo",
           ifelse(price >= q_price[2], "hi", "md"))
    rr <- ifelse(rate  <= q_rate[1],  "lo",
           ifelse(rate  >= q_rate[2],  "hi", "md"))
    aa <- ifelse(amax  <= q_amax[1],  "lo",
           ifelse(amax  >= q_amax[2],  "hi", "md"))
    cc <- ifelse(cyc   <= q_cyc[1],   "lo",
           ifelse(cyc   >= q_cyc[2],   "hi", "md"))

    fifelse(pp == "lo" & rr == "hi" & aa == "lo" & cc == "lo", "pessimistic",
    fifelse(pp == "hi" & rr == "lo" & aa == "hi" & cc == "hi", "optimistic",
    fifelse(pp == "md" & rr == "md" & aa == "md" & cc == "md", "baseline",
            "other")))
  }

  dt[, scen := label_scenario(bean_price, discount_rate,
                                adoption_ceiling, cycle_compression)]
  out <- dt[, .(
    n          = .N,
    prob       = .N / nrow(dt),
    mean_npv   = mean(npv),
    median_npv = median(npv),
    pr_neg     = mean(npv < 0)
  ), by = scen][order(mean_npv)]
  out
}


# --- EVPPI via per-parameter "conditional-on-truth" loop ---------------
# Decision D: choose action A1 = invest in the tool (payoff = NPV) vs
#   A0 = reference alternative (payoff = threshold, e.g. a competing
#   program or a required hurdle rate). The program's real-world alternatives
#   are discrete: $0 (do nothing) or ~$200M+ (alternative breeding
#   technology). We parameterize by `threshold`.
#
# EV under uncertainty:        max_A E[payoff_A] = max(E[NPV], threshold)
# EV with perfect info on X_i: E_{X_i}[max(E[NPV|X_i], threshold)]
# EVPPI_i(threshold) = second - first.
compute_evppi <- function(dt, params, threshold = 0) {
  baseline <- max(mean(dt$npv), threshold)
  out <- lapply(params, function(p) {
    cond <- dt[, .(E_npv = mean(npv), n = .N), by = c(p)]
    cond[, w := n / sum(n)]
    ev_perfect <- sum(cond$w * pmax(cond$E_npv, threshold))
    evppi      <- ev_perfect - baseline
    data.table(param = p, threshold = threshold,
                EV_prior = baseline, EV_perfect = ev_perfect,
                EVPPI = evppi)
  })
  rbindlist(out)
}


run_decision <- function(output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- l2_out_dir()
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("Loading lookup_table.parquet ...\n")
  dt <- load_lookup()
  cat(sprintf("  %d rows\n", nrow(dt)))

  # --- Scenario decomposition -----------------------------------------------
  cat("\n--- Scenario decomposition ---\n")
  scn <- scenario_decomposition(dt)
  print(scn)
  fwrite(scn, file.path(output_dir, "scenario_decomposition.csv"))

  # --- EVPPI at multiple decision thresholds ---------------------------------
  # threshold = payoff of the "reference" alternative the program could choose
  #   instead of investing in the tool. We sweep across plausible hurdle values
  #   (from 0 = "do nothing" to $500M = "deploy competing tech"). EVPPI is
  #   highest near the indifference threshold.
  cat("\n--- EVPPI across thresholds ---\n")
  params <- c("cycle_compression", "bean_price", "discount_rate",
              "adoption_ceiling", "total_budget", "tool_fixed_cost",
              "error_reduction", "cost_reduction")
  params <- intersect(params, colnames(dt))

  # Fixed grid plus the indifference point tau = E[NPV], where EVPPI peaks;
  # the fixed grid alone misses the peak whenever E[NPV] moves (v1 -> v2).
  # Hurdles as multiples of E[NPV] so the grid brackets the indifference point
  # whatever the calibration (the v1 fixed $ grid missed it once E[NPV] moved)
  hurdle_mult <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3)
  thresholds <- mean(dt$npv) * hurdle_mult
  all_ev <- rbindlist(lapply(thresholds, function(th) {
    compute_evppi(dt, params, threshold = th)
  }))
  all_ev[, EVPPI_M := EVPPI / 1e6]
  all_ev[, threshold_M := threshold / 1e6]
  all_ev[, hurdle_mult := threshold / mean(dt$npv)]
  all_ev[, at_indifference := abs(hurdle_mult - 1) < 1e-9]
  setorder(all_ev, -EVPPI)

  cat("Top rows of EVPPI × threshold:\n")
  print(head(all_ev[, .(param, threshold_M, EVPPI_M)], 15))
  # The manuscript's EVPPI (area fixed, same prior as the Sobol analysis) comes from 17_exact_gsa.R as evppi.csv;
  # this lookup-table version averages over bean area as well and is kept under another name.
  fwrite(all_ev, file.path(output_dir, "evppi_lookup_all_areas.csv"))

  # --- Interpretation: at which threshold does each param matter most? --
  cat("\n--- At each threshold, top driver of VoI ---\n")
  for (th in thresholds) {
    sub <- all_ev[threshold == th][order(-EVPPI)]
    cat(sprintf("  threshold $%5.0fM:  %s ($%.1fM)   2nd: %s ($%.1fM)\n",
                 th / 1e6,
                 sub$param[1], sub$EVPPI_M[1],
                 sub$param[2], sub$EVPPI_M[2]))
  }

  invisible(list(scenarios = scn, evppi = all_ev))
}


if (sys.nframe() == 0) {
  run_decision()
}
