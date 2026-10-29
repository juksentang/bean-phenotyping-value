# ============================================================================
# Empirical calibration — Step 5: plot-stand covariate and plot error variance
# ============================================================================
# Question: by how much does adjusting area-based plot yield for the number of
# plants in the plot reduce the plot residual variance (delta_e), and how much
# genetic signal does the adjustment remove?
# Input : empirical/outputs/zudpep_plots.csv (yield_area_kg_ha, plants_harvested)
# Output: empirical/outputs/stand_covariate.csv          per-trial results
#         empirical/outputs/stand_covariate_pooled.csv   per-series results
#         empirical/outputs/stand_covariate_summary.md
#
# Definitions (per trial, REML):
#   unadjusted  y ~ rep + (1|rep:block) + (1|genotype)
#   adjusted    the same plus a stand term f(s); f chosen by ML AIC among
#               linear, linear + quadratic and log(s)
#   delta_e   = 1 - Ve_adj / Ve_unadj
#   delta_g   = 1 - Vg_adj / Vg_unadj   (genetic variance removed)
#   "within"  variant: stand centered on the genotype mean within the trial,
#             so only non-genetic stand deviations enter the adjustment.
#
# Usage (from the project root):
#   Rscript empirical/R/05_stand_covariate.R

source("empirical/R/varcomp_utils.R")
out_dir <- "empirical/outputs"
OUTLIER_Z <- 3.5

p <- fread(file.path(out_dir, "zudpep_plots.csv"))
n_all <- nrow(p)
# Stand-related exclusions from Step 1 are undone here: low-stand plots are the
# plots a stand adjustment is meant to correct. Plots need an area yield and a
# positive plant count.
p <- p[!is.na(yield_area_kg_ha) & yield_area_kg_ha > 0 & !is.na(plants_harvested) & plants_harvested > 0]
p[, `:=`(y = yield_area_kg_ha, s = plants_harvested, env = trial_id)]
n_use0 <- nrow(p)

# Outlier screen on the area-based yield (unadjusted model), applied to both
# the adjusted and the unadjusted fits so they use identical plots.
p[, z := NA_real_]
for (tid in unique(p$trial_id)) {
  idx <- p[, which(trial_id == tid)]
  f <- suppressMessages(lmer(y ~ factor(rep) + (1 | genotype), data = p[idx]))
  p[idx, z := residuals(f) / sigma(f)]
}
n_out <- p[abs(z) > OUTLIER_Z, .N]
p <- p[abs(z) <= OUTLIER_Z]

# --- Model helpers ------------------------------------------------------------
rand_part <- function(blk) if (blk) "(1 | genotype) + (1 | repblock)" else "(1 | genotype)"
stand_terms <- list(none = "", linear = "+ sc", quadratic = "+ sc + I(sc^2)", log = "+ ls")

fit_one <- function(d, term, blk, reml = TRUE) {
  f <- as.formula(paste("y ~ rep", term, "+", rand_part(blk)))
  suppressMessages(lmer(f, data = d, REML = reml))
}
vg_ve <- function(m) {
  vc <- as.data.table(VarCorr(m))
  c(Vg = vc[grp == "genotype", vcov], Ve = vc[grp == "Residual", vcov])
}
# Univariate REML genotype and residual variance of an arbitrary response
uni_vc <- function(d, resp, blk) {
  dd <- copy(d); dd[, yy := get(resp)]
  f <- as.formula(paste("yy ~ rep +", rand_part(blk)))
  vg_ve(suppressMessages(lmer(f, data = dd)))
}
# Genetic and residual covariance by the variance-of-sum identity
# Cov(a,b) = [V(a+b) - V(a) - V(b)] / 2, applied to the genotype and residual
# components of separate univariate REML fits (b rescaled to the SD of a).
cor_by_sum <- function(d, blk) {
  dd <- copy(d)
  k  <- sd(dd$y) / sd(dd$s)
  dd[, s_k := s * k][, ysum := y + s_k]
  a <- uni_vc(dd, "y", blk); b <- uni_vc(dd, "s_k", blk); ab <- uni_vc(dd, "ysum", blk)
  cg <- (ab[["Vg"]] - a[["Vg"]] - b[["Vg"]]) / 2
  ce <- (ab[["Ve"]] - a[["Ve"]] - b[["Ve"]]) / 2
  rg <- if (a[["Vg"]] > 1e-6 * a[["Ve"]] && b[["Vg"]] > 1e-6 * b[["Ve"]]) cg / sqrt(a[["Vg"]] * b[["Vg"]]) else NA_real_
  re <- ce / sqrt(a[["Ve"]] * b[["Ve"]])
  c(r_g = rg, r_e = re)
}

# --- 1-2. Per trial ---------------------------------------------------------------
res <- list()
for (tid in sort(unique(p$trial_id))) {
  d <- copy(p[trial_id == tid])
  d[, `:=`(rep = factor(rep), genotype = factor(genotype))]
  blk <- d[, uniqueN(na.omit(block)) > 1 && all(!is.na(block))]
  if (blk) d[, repblock := factor(paste(rep, block))]
  d[, sc := s - mean(s)]
  d[, ls := log(s) - mean(log(s))]
  d[, sw := s - mean(s), by = genotype]          # within-genotype stand deviation

  aic <- sapply(stand_terms, function(t) AIC(fit_one(d, t, blk, reml = FALSE)))
  best <- names(which.min(aic[-1]))
  m0 <- fit_one(d, "", blk); m1 <- fit_one(d, stand_terms[[best]], blk)
  mw <- fit_one(d, "+ sw", blk)
  v0 <- vg_ve(m0); v1 <- vg_ve(m1); vw <- vg_ve(mw)
  slope <- fixef(fit_one(d, "+ sc", blk))[["sc"]]

  vs <- uni_vc(d, "s", blk)
  cr <- tryCatch(cor_by_sum(d, blk), error = function(e) c(r_g = NA_real_, r_e = NA_real_))

  res[[tid]] <- data.table(
    trial_id = tid, set = d$set[1], n_plots = nrow(d), n_geno = nlevels(d$genotype),
    mean_yield_area = mean(d$y), mean_stand = mean(d$s), cv_stand_pct = 100 * sd(d$s) / mean(d$s),
    best_stand_term = best, aic_none_minus_best = aic[["none"]] - aic[[best]],
    slope_kg_ha_per_plant = slope,
    Vg_unadj = v0[["Vg"]], Ve_unadj = v0[["Ve"]], Vg_adj = v1[["Vg"]], Ve_adj = v1[["Ve"]],
    h2_w_unadj = v0[["Vg"]] / sum(v0), h2_w_adj = v1[["Vg"]] / sum(v1),
    delta_e = 1 - v1[["Ve"]] / v0[["Ve"]],
    delta_g = if (v0[["Vg"]] > 0) 1 - v1[["Vg"]] / v0[["Vg"]] else NA_real_,
    Vg_within = vw[["Vg"]], Ve_within = vw[["Ve"]],
    delta_e_within = 1 - vw[["Ve"]] / v0[["Ve"]],
    delta_g_within = if (v0[["Vg"]] > 0) 1 - vw[["Vg"]] / v0[["Vg"]] else NA_real_,
    h2_w_within = vw[["Vg"]] / sum(vw),
    Vg_stand = vs[["Vg"]], Ve_stand = vs[["Ve"]], h2_stand = vs[["Vg"]] / sum(vs),
    r_g_stand_yield = cr[["r_g"]], r_e_stand_yield = cr[["r_e"]],
    vg_unadj_boundary = v0[["Vg"]] < 1e-6 * v0[["Ve"]])
}
res <- rbindlist(res)
res[, delta_h2_w := h2_w_adj - h2_w_unadj]
fwrite(res, file.path(out_dir, "stand_covariate.csv"))

# --- Pooled multi-environment version per series -------------------------------------
# Homogeneous-residual lme4 model of Step 2 with environment-specific stand
# slopes (linear, stand centered within environment); same plots with and
# without the covariate.
pool <- list()
for (s_ in sort(unique(p$set))) {
  d <- prep_met(copy(p[set == s_]))
  d[, sc := s - mean(s), by = env]
  d[, sw := s - mean(s), by = .(env, genotype)]
  blk <- d[, any(blk_ok)]
  f0 <- met_formula(blk)
  f1 <- update(f0, . ~ . + env:sc)
  fw <- update(f0, . ~ . + env:sw)
  m0 <- suppressMessages(lmer(f0, data = d)); m1 <- suppressMessages(lmer(f1, data = d))
  mw <- suppressMessages(lmer(fw, data = d))
  a <- vc_hom(m0); b <- vc_hom(m1); w <- vc_hom(mw)
  pool[[s_]] <- data.table(
    set = s_, n_env = nlevels(d$env), n_plots = nrow(d),
    Vg_unadj = a[["Vg"]], Vge_unadj = a[["Vge"]], Ve_unadj = a[["Ve"]],
    Vg_adj = b[["Vg"]], Vge_adj = b[["Vge"]], Ve_adj = b[["Ve"]],
    delta_e = 1 - b[["Ve"]] / a[["Ve"]],
    delta_g = if (a[["Vg"]] > 0) 1 - b[["Vg"]] / a[["Vg"]] else NA_real_,
    delta_ge = if (a[["Vge"]] > 0) 1 - b[["Vge"]] / a[["Vge"]] else NA_real_,
    h2_w_unadj = a[["Vg"]] / (a[["Vg"]] + a[["Ve"]]), h2_w_adj = b[["Vg"]] / (b[["Vg"]] + b[["Ve"]]),
    delta_e_within = 1 - w[["Ve"]] / a[["Ve"]],
    delta_g_within = if (a[["Vg"]] > 0) 1 - w[["Vg"]] / a[["Vg"]] else NA_real_)
}
pool <- rbindlist(pool)
fwrite(pool, file.path(out_dir, "stand_covariate_pooled.csv"))

# --- TSEZVG check for a stand column ---------------------------------------------------
ts_cols <- names(fread("data/external/bean_public_2026-10-10/datasets/TSEZVG/x/TPG_Phenotypedata.tab", nrows = 1))

# --- Summary --------------------------------------------------------------------------
L <- character(); w <- function(...) L <<- c(L, paste0(...))
ds_line <- function(x, d = 2) { z <- dist_summary(x)
  paste0("median ", fmt(z$median, d), " (IQR ", fmt(z$q25, d), " to ", fmt(z$q75, d), "; range ",
         fmt(z$min, d), " to ", fmt(z$max, d), "; n = ", z$n, ")") }
w("# Stand covariate and plot error variance (ZUDPEP)")
w("")
w("Generated by `empirical/R/05_stand_covariate.R` (", format(Sys.Date()), ").")
w("")
w("## Data")
w("")
w("- Response: area-based plot yield `yield_area_kg_ha` = plot seed weight (g) x 10 / 4.5 m2. ",
  "The per-plant-scaled YDHA is not used because it already divides by stand.")
w("- Stand: `plants_harvested` (PLPLOT, plants per plot at harvest).")
w("- Plots: ", n_all, " in the cleaned file; ", n_use0, " with area yield > 0 and stand > 0 ",
  "(plots excluded in Step 1 for low stand are retained here, since they are what a stand ",
  "adjustment corrects); ", n_out, " residual outliers on the area-yield scale removed ",
  "(|standardized residual| > ", OUTLIER_Z, "); ", nrow(p), " plots in 33 trials analyzed. ",
  "Adjusted and unadjusted fits use the same plots.")
w("- TSEZVG has no stand-like column (columns: ", paste(ts_cols, collapse = ", "),
  "; GHT and GHP are growth habit), so this analysis is ZUDPEP-only.")
w("")
w("## Per-trial results")
w("")
w("Model: y ~ rep + (1|rep:block) + (1|genotype), REML, with and without a stand term. The ",
  "stand term (linear, linear + quadratic, or log) is chosen per trial by ML AIC. ",
  "\"Within\" uses stand centered on each genotype's mean in that trial, so only the ",
  "non-genetic part of stand variation is regressed out.")
w("")
w("- Best stand term: ", paste(names(table(res$best_stand_term)), table(res$best_stand_term), sep = " ", collapse = ", "),
  "; the stand term lowers AIC by more than 2 in ", res[aic_none_minus_best > 2, .N], " of ", nrow(res), " trials.")
w("- Stand CV within trials: ", ds_line(res$cv_stand_pct, 1), " %.")
w("- delta_e (phenotypic stand adjustment): ", ds_line(res$delta_e), ".")
w("- delta_g (share of Vg removed), trials with Vg > 0: ", ds_line(res[vg_unadj_boundary == FALSE, delta_g]), ".")
w("- h2_w unadjusted ", ds_line(res$h2_w_unadj), "; adjusted ", ds_line(res$h2_w_adj),
  "; change ", ds_line(res$delta_h2_w), ".")
w("- delta_e (within-genotype stand adjustment): ", ds_line(res$delta_e_within),
  "; delta_g within: ", ds_line(res[vg_unadj_boundary == FALSE, delta_g_within]), ".")
w("")
w("## Is stand genetic?")
w("")
w("- Stand heritability per trial, h2 = Vg/(Vg+Ve) from stand ~ rep + (1|rep:block) + (1|genotype): ",
  ds_line(res$h2_stand), ".")
w("- Genetic correlation between stand and yield (r_g): ", ds_line(res$r_g_stand_yield),
  "; residual (plot) correlation (r_e): ", ds_line(res$r_e_stand_yield), ".")
w("- Method: variance-of-sum approximation, not a bivariate REML fit. Separate univariate REML ",
  "fits of yield, rescaled stand and their sum give Cov = [V(a+b) - V(a) - V(b)] / 2 for the ",
  "genotype and residual components. r_g is not estimated when either Vg is at zero, and can ",
  "fall outside [-1, 1] when genetic variances are small.")
w("")
w("| trial | stand mean | stand CV % | term | slope (kg/ha per plant) | Ve unadj | Ve adj | delta_e | delta_g | h2_w unadj | h2_w adj | delta_e within | h2 stand | r_g | r_e |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(res))) with(res[i], w("| ", trial_id, " | ", fmt(mean_stand, 1), " | ", fmt(cv_stand_pct, 1), " | ",
  best_stand_term, " | ", fmt(slope_kg_ha_per_plant, 1), " | ", fmt0(Ve_unadj), " | ", fmt0(Ve_adj), " | ",
  fmt(delta_e), " | ", ifelse(is.na(delta_g), "", fmt(delta_g)), " | ", fmt(h2_w_unadj), " | ", fmt(h2_w_adj), " | ",
  fmt(delta_e_within), " | ", fmt(h2_stand), " | ", ifelse(is.na(r_g_stand_yield), "", fmt(r_g_stand_yield)), " | ",
  fmt(r_e_stand_yield), " |"))
w("")
w("## Pooled per series (multi-environment)")
w("")
w("Model of Step 2 (homogeneous residual, lme4) with environment-specific linear stand slopes.")
w("")
w("| set | n_env | Ve unadj | Ve adj | delta_e | delta_g | delta_ge | h2_w unadj | h2_w adj | delta_e within | delta_g within |")
w("|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(pool))) with(pool[i], w("| ", set, " | ", n_env, " | ", fmt0(Ve_unadj), " | ", fmt0(Ve_adj), " | ",
  fmt(delta_e), " | ", ifelse(is.na(delta_g), "", fmt(delta_g)), " | ", ifelse(is.na(delta_ge), "", fmt(delta_ge)), " | ",
  fmt(h2_w_unadj), " | ", fmt(h2_w_adj), " | ", fmt(delta_e_within), " | ", ifelse(is.na(delta_g_within), "", fmt(delta_g_within)), " |"))
w("")
w("## Interpretation and caveats")
w("")
w("- When stand is partly genetic (h2_stand > 0, r_g != 0), the phenotypic adjustment also ",
  "removes or adds genetic variance (delta_g). delta_e then overstates the gain in selection ",
  "accuracy; the relevant quantity is the change in h2_w (or Vg/Ve). The within-genotype ",
  "adjustment removes only non-genetic stand variation and is the cleaner estimate of an ",
  "error-variance reduction that leaves Vg intact.")
w("- `plants_harvested` is a harvest count, not an emergence count. It reflects emergence plus ",
  "later plant loss to disease, drought, lodging and damage, which are partly genotype ",
  "dependent (vigor, disease resistance). An emergence count taken by smartphone early in the ",
  "season would capture the seed and establishment component and would differ in both its ",
  "error-reducing power and its genetic content.")
w("- Plot size is 4.5 m2 with about 20 to 75 harvested plants; results apply to similar ",
  "small-plot preliminary yield trials.")
writeLines(L, file.path(out_dir, "stand_covariate_summary.md"))
cat("Done.\n")
