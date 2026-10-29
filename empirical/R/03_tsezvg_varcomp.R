# ============================================================================
# Empirical calibration — Step 3: variance components for TSEZVG grain yield
# ============================================================================
# Source: Harvard Dataverse doi:10.7910/DVN/TSEZVG (CC BY 4.0); Saradadevi et al.
# 2021, Plant Genome 14:e20156, doi:10.1002/tpg2.20156. Plot-level phenotypes of
# 33 trials (Kawanda and Kachwekano, Uganda; Kitengule, Tanzania; 2015-2018).
# GID is a germplasm identifier that is consistent across trials.
# Input : TPG_Phenotypedata.tab
# Output: empirical/outputs/tsezvg_trial_h2.csv
#         empirical/outputs/tsezvg_varcomp.csv
#         empirical/outputs/tsezvg_varcomp_summary.md
#
# Usage (from the project root):
#   Rscript empirical/R/03_tsezvg_varcomp.R [n_bootstrap] [path/to/TPG_Phenotypedata.tab]

source("empirical/R/varcomp_utils.R")

args <- commandArgs(trailingOnly = TRUE)
NSIM <- if (length(args) >= 1) as.integer(args[1]) else 200L
src  <- if (length(args) >= 2) args[2] else
  "data/external/bean_public_2026-10-10/datasets/TSEZVG/x/TPG_Phenotypedata.tab"
NCPU <- max(1L, min(8L, parallel::detectCores() - 1L))
OUTLIER_Z <- 3.5
out_dir <- "empirical/outputs"

raw <- fread(src)
setnames(raw, make.names(names(raw)))
n_raw <- nrow(raw)
x <- raw[, .(trial_id = Trial.Code, trial_no = TrialNo, site = sub(",.*", "", LOCATION),
             year = Year, genotype = as.character(GID), designation = DESIGNATION,
             rep = REP_NO, block = BLOCK_NO, growth_habit = GHT, y = GY)]
x[, env := trial_id]

# --- Cleaning -----------------------------------------------------------------
x[, exclude_reason := NA_character_]
x[is.na(y) | y <= 0, exclude_reason := "missing_or_zero_yield"]
x[, outlier_z := NA_real_]
for (tid in unique(x$trial_id)) {
  idx <- x[, which(trial_id == tid & is.na(exclude_reason))]
  d <- x[idx]
  if (nrow(d) < 10 || uniqueN(d$rep) < 2) next
  fit <- suppressMessages(lmer(y ~ factor(rep) + (1 | genotype), data = d))
  x[idx, outlier_z := residuals(fit) / sigma(fit)]
}
x[is.na(exclude_reason) & !is.na(outlier_z) & abs(outlier_z) > OUTLIER_Z,
  exclude_reason := "residual_outlier"]
cl <- x[is.na(exclude_reason)]
habit_by_trial <- cl[, .(n_habit = uniqueN(growth_habit), habit = growth_habit[1]), by = trial_id]
stopifnot(all(habit_by_trial$n_habit == 1))

# --- 1. Per-trial h2_w --------------------------------------------------------------
trial_h2 <- rbindlist(lapply(split(cl, cl$trial_id), function(d) {
  cbind(data.table(trial_id = d$trial_id[1], trial_no = d$trial_no[1], site = d$site[1],
                   year = d$year[1], growth_habit = d$growth_habit[1]), fit_trial(d))
}))
setorder(trial_h2, trial_no)
fwrite(trial_h2, file.path(out_dir, "tsezvg_trial_h2.csv"))

# --- 2. Connected sets ----------------------------------------------------------
cs <- connected_sets(unique(cl[, .(env, genotype)]))
cs <- merge(cs, unique(cl[, .(env, trial_no, growth_habit)]), by = "env")
setorder(cs, component, trial_no)

# Balanced core of a component: start from its largest trial and add trials
# that contain at least 90% of that trial's genotypes.
core_of <- function(d) {
  g <- lapply(split(d$genotype, d$env), unique)
  ref <- names(which.max(lengths(g)))
  keep <- names(g)[sapply(g, function(z) mean(g[[ref]] %in% z) >= 0.9)]
  keep
}

met <- list(); set_rows <- list()
for (k in sort(unique(cs$component))) {
  envs <- cs[component == k, env]
  d <- cl[env %in% envs]
  des <- set_design(d)
  hab <- paste(unique(cs[component == k, growth_habit]), collapse = "+")
  nm  <- sprintf("TSEZVG_%s_%s", hab, paste(cs[component == k, trial_no], collapse = "-"))
  set_rows[[nm]] <- data.table(set = nm, component = k, trials = paste(cs[component == k, trial_no], collapse = ","),
                               n_env = des$n_env, n_geno = des$n_geno, n_geno_common = des$n_geno_common)
  if (des$n_env < 3 || des$n_geno_common < 10) next
  cat("Fitting", nm, "\n")
  met[[nm]] <- cbind(data.table(dataset = "TSEZVG", growth_habit = hab, set_type = "connected_set",
                                trials = paste(cs[component == k, trial_no], collapse = ",")),
                     analyze_set(d, nm, nsim = NSIM, ncpus = NCPU))
  core <- core_of(d)
  if (length(core) >= 3 && length(core) < length(envs)) {
    dc  <- d[env %in% core]
    tn  <- sort(unique(dc$trial_no))
    nmc <- sprintf("TSEZVG_%s_%s_core", hab, paste(tn, collapse = "-"))
    cat("Fitting", nmc, "\n")
    met[[nmc]] <- cbind(data.table(dataset = "TSEZVG", growth_habit = hab, set_type = "balanced_core",
                                   trials = paste(tn, collapse = ",")),
                        analyze_set(dc, nmc, nsim = NSIM, ncpus = NCPU))
  }
}
set_rows <- rbindlist(set_rows)

# --- 3. Pooled models (variance components shared across sets) -----------------------
pooled <- list(Bush = cl[growth_habit == "Bush"], Climbing = cl[growth_habit == "Climbing"], All = cl)
for (h in names(pooled)) {
  cat("Fitting pooled", h, "\n")
  met[[paste0("pooled_", h)]] <- cbind(
    data.table(dataset = "TSEZVG", growth_habit = h, set_type = "pooled",
               trials = paste(sort(unique(pooled[[h]]$trial_no)), collapse = ",")),
    analyze_set(pooled[[h]], paste0("TSEZVG_pooled_", h), boot = FALSE))
}
met <- rbindlist(met, fill = TRUE)
fwrite(met, file.path(out_dir, "tsezvg_varcomp.csv"))

# --- 4. Summary ----------------------------------------------------------------
L <- character(); w <- function(...) L <<- c(L, paste0(...))
w("# TSEZVG grain-yield variance components")
w("")
w("Generated by `empirical/R/03_tsezvg_varcomp.R` (", format(Sys.Date()), "). Source: Harvard ",
  "Dataverse doi:10.7910/DVN/TSEZVG (CC BY 4.0), `TPG_Phenotypedata.tab`; Saradadevi et al. 2021, ",
  "Plant Genome, doi:10.1002/tpg2.20156. Definitions as in `zudpep_varcomp_summary.md`.")
w("")
w("## Data and cleaning")
w("")
w("- Rows read: ", n_raw, " plots in ", uniqueN(raw$TrialNo), " trials. Genotype = GID (one GID ",
  "per DESIGNATION and vice versa; no name harmonization needed). Environment = trial ",
  "(Trial Code; the file gives the year but not the season).")
w("- GY is analyzed as recorded (kg/ha; trial means ", fmt0(min(trial_h2$mean_kg_ha)), "-",
  fmt0(max(trial_h2$mean_kg_ha)), "). How GY was scaled from plot weight is not documented in the file.")
w("- Excluded: missing GY ", x[exclude_reason == "missing_or_zero_yield", .N], " plots (",
  "including all ", raw[TrialNo == 25, .N], " plots of trial 25, which has no GY data); ",
  "residual outliers (|standardized conditional residual| > ", OUTLIER_Z,
  " from yield ~ rep + (1|genotype) per trial) ", x[exclude_reason == "residual_outlier", .N], " plots.")
w("- Growth habit (GHT) is constant within trials: ", sum(habit_by_trial$habit == "Bush"), " bush and ",
  sum(habit_by_trial$habit == "Climbing"), " climbing trials analyzed. Incomplete blocks are ",
  "numbered within replicate (rep:block).")
w("")
w("## Per-trial within-environment plot heritability")
w("")
for (h in c("Bush", "Climbing", "All")) {
  tz <- if (h == "All") trial_h2 else trial_h2[growth_habit == h]
  ds <- dist_summary(tz$h2_w)
  w("- ", h, ": ", nrow(tz), " trials; h2_w median ", fmt(ds$median), " (IQR ", fmt(ds$q25), "-",
    fmt(ds$q75), "; range ", fmt(ds$min), "-", fmt(ds$max), "); Vg at boundary in ",
    sum(tz$vg_boundary), " trials.")
}
w("")
w("| trial | site | year | habit | genotypes | reps | mean kg/ha | Vg | Ve | h2_w | boundary |")
w("|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(trial_h2))) with(trial_h2[i], w("| ", trial_id, " | ", site, " | ", year, " | ",
  growth_habit, " | ", n_geno, " | ", n_rep, " | ", fmt0(mean_kg_ha), " | ", fmt0(Vg), " | ", fmt0(Ve),
  " | ", fmt(h2_w), " | ", ifelse(vg_boundary, "yes", ""), " |"))
w("")
w("## Connected sets")
w("")
w("Trials are linked when they share >= 10 GIDs and >= 50% of the smaller trial's GIDs; sets are ",
  "connected components. Eligible for a multi-environment fit: >= 3 trials and >= 10 GIDs common ",
  "to all trials of the set. Where a set is unbalanced, its balanced core (trials containing ",
  ">= 90% of the GIDs of the largest trial) is also fitted.")
w("")
w("| set | trials | n_env | n_geno | common to all |")
w("|---|---|---|---|---|")
for (i in seq_len(nrow(set_rows))) with(set_rows[i], w("| ", set, " | ", trials, " | ", n_env, " | ",
  n_geno, " | ", n_geno_common, " |"))
w("")
w("## Multi-environment estimates")
w("")
w("Model: yield ~ env + (1|genotype) + (1|genotype:env) + (1|env:rep) + (1|env:rep:block), REML, ",
  "residual variance per environment (glmmTMB) where it converged; Ve = mean environment residual ",
  "variance. 95% intervals: parametric bootstrap (", NSIM, " simulations) of the same model; ",
  "kappa interval derived from that of rho_env. Pooled rows share one set of variance components ",
  "across all trials of a growth habit (genotypes are connected across sets only through a few checks).")
w("")
w("| set | type | n_env | n_geno (common) | residual | Vg | Vge | Ve | h2_w [95% CI] | kappa [95% CI] | rho_env | h2_total | H2_mean |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(met))) with(met[i], w("| ", set, " | ", set_type, " | ", n_env, " | ", n_geno, " (", n_geno_common, ") | ",
  sub("_by_env", "", resid_model), " | ", fmt0(Vg), " | ", fmt0(Vge), " | ", fmt0(Ve), " | ",
  fmt(h2_w), if (!is.na(h2_w_lo)) paste0(" [", fmt(h2_w_lo), ", ", fmt(h2_w_hi), "]") else "",
  " | ", fmt(kappa), if (!is.na(kappa_lo)) paste0(" [", fmt(kappa_lo), ", ", fmt_k(kappa_hi), "]") else "",
  " | ", fmt(rho_env), " | ", fmt(h2_total), " | ", fmt(H2_mean), " |"))
writeLines(L, file.path(out_dir, "tsezvg_varcomp_summary.md"))
cat("Done.\n")
