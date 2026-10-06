#!/usr/bin/env Rscript
# =============================================================================
# meta_analysis.R  --  literature meta-analyses for the paper
#   "Faster, not more accurate? Where the economic value of AI phenotyping lies
#    in smallholder bean breeding"
#
# Run from the project root:   Rscript paper/lit/meta_analysis.R
# Needs: metafor (>= 4), base graphics only for plots.
#
# Inputs (data extraction sheets; every value was read from the full text of
# the cited article, with its page/table location recorded):
#   paper/lit/extract/meta_genetic_gain/extraction.csv  (+ candidates_status.csv)
#   paper/lit/extract/meta_pheno_accuracy/extraction.csv
#   paper/lit/extract/meta_pheno_cost/extraction.csv    (+ studies_assessed.csv)
#   paper/lit/extract/synthesis/crossref_verification.tsv  (Crossref metadata
#                                                          for every included DOI)
# Outputs:
#   paper/lit/meta_results.csv            one row per fitted model
#   paper/lit/meta_leave_one_out.csv      leave-one-cluster/study-out ranges
#   paper/lit/meta_inputs_used.csv        every effect that entered any model + flags
#   paper/lit/figures/*.png               forest / funnel / translation plots
#
# INCLUSION RULE (applied in code): an effect enters a model only if its study has
#   included = TRUE, verified_crossref = TRUE (also re-checked live), fulltext_read = TRUE.
# Preprints (Crossref type posted-content) are kept in the primary sets but are
# dropped in a dedicated sensitivity analysis.
#
# -----------------------------------------------------------------------------
# GENERAL STATISTICAL CHOICES
#  * Random effects, REML for tau^2. Confidence intervals use the Knapp-Hartung (t)
#    adjustment because every k here is small. 95% prediction intervals (PI) are
#    reported next to every pooled estimate and are the number that matters for
#    parameterising a simulation (what a NEW programme/tool might show).
#  * Dependence: several effects come from the same programme/paper. Primary
#    models for MA1/MA2 are three-level (cluster/effect) REML models (rma.mv);
#    naive two-level, cluster-robust and one-effect-per-cluster versions are
#    sensitivity analyses.
#  * Small-study effects (Egger-type regression, rank test, trim-and-fill, funnel)
#    are examined only where k >= 10 independent units exist (MA2 at study level).
#  * k < 3 usable effects  ->  descriptive only; the CSV says so.
# =============================================================================

suppressPackageStartupMessages(library(metafor))
options(stringsAsFactors = FALSE, warn = 1)
set.seed(20261005)

if (!file.exists("paper/lit/extract/meta_genetic_gain/extraction.csv"))
  stop("Run from the project root.")
LIT  <- "paper/lit"
FIG  <- file.path(LIT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Live Crossref verification table ------------------------------------
cr <- read.delim(file.path(LIT, "extract/synthesis/crossref_verification.tsv"),
                 quote = "", check.names = FALSE)
cr_ok <- function(doi) doi %in% cr$doi
cr_type <- function(doi) cr$type[match(doi, cr$doi)]
# Crossref marks preprints with an abstract-less type; the posted-content type is
# what we need. Re-derive from the JSON-based table (column "type").
is_preprint <- function(doi) {
  t <- cr$type[match(doi, cr$doi)]
  !is.na(t) & t == "posted-content"
}

# ---- 1. helpers ---------------------------------------------------------------
RES <- list()
add_row <- function(id, analysis, stratum, role, model, k, kcl, est, lb, ub, plb, pub,
                    tau2, I2, Q, Qp, units, ci_method, notes = "") {
  RES[[length(RES) + 1]] <<- data.frame(
    analysis_id = id, analysis = analysis, stratum = stratum, role = role, model = model,
    k_effects = k, k_clusters = kcl, estimate = est, ci_lb = lb, ci_ub = ub,
    pi_lb = plb, pi_ub = pub, tau2 = tau2, I2_pct = I2, Q = Q, Q_p = Qp,
    units = units, ci_method = ci_method, notes = notes)
}
I2_mv <- function(res) {            # Cheung/Viechtbauer I^2 for rma.mv, total over variance components
  W <- diag(1 / res$vi); X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k - res$p) / sum(diag(P)))
}
# generic summariser: returns list on the model scale
summ <- function(res, tf = identity) {
  pr <- predict(res)
  mv <- inherits(res, "rma.mv")
  list(est = tf(as.numeric(res$beta[1])), lb = tf(pr$ci.lb), ub = tf(pr$ci.ub),
       plb = tf(pr$pi.lb), pub = tf(pr$pi.ub),
       tau2 = if (mv) sum(res$sigma2) else res$tau2,
       I2 = if (mv) I2_mv(res) else res$I2,
       Q = res$QE, Qp = res$QEp, k = res$k)
}
fit_uni <- function(yi, vi) rma(yi, vi, method = "REML", test = "knha")
fit_mv  <- function(d) rma.mv(yi, vi, random = ~ 1 | cluster / eff, data = d,
                              method = "REML", test = "t")
log_row <- function(id, analysis, stratum, role, model, res, kcl, tf = identity,
                    units, ci_method, notes = "") {
  s <- summ(res, tf)
  add_row(id, analysis, stratum, role, model, s$k, kcl, s$est, s$lb, s$ub, s$plb, s$pub,
          s$tau2, s$I2, s$Q, s$Qp, units, ci_method, notes)
  invisible(s)
}
desc_row <- function(id, analysis, stratum, vals, units, notes, k_cl = NA) {
  vals <- vals[!is.na(vals)]
  add_row(id, analysis, stratum, "descriptive", "median (range) -- NOT pooled", length(vals), k_cl,
          median(vals), min(vals), max(vals), NA, NA, NA, NA, NA, NA, units,
          "none: estimate=median, ci_lb/ci_ub=min/max",
          paste0("k<3 usable pooled effects or no variances -> descriptive only. ", notes))
}
LOO <- list()
add_loo <- function(id, dropped, est, lb, ub) {
  LOO[[length(LOO) + 1]] <<- data.frame(analysis_id = id, left_out = dropped,
                                         estimate = est, ci_lb = lb, ci_ub = ub)
}
USED <- list()
add_used <- function(ma, d, cols) {
  d <- d[, cols, drop = FALSE]; d$analysis_family <- ma
  USED[[length(USED) + 1]] <<- d
}

# ---- plotting ----------------------------------------------------------------
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#d9d8d2"; BG <- "#fcfcfb"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
# Minimal custom forest plot (base graphics): rows is a data.frame with
#   label, est, lb, ub, type in {study, pooled, header, note}, annot (text at right),
#   plb/pub (prediction interval for pooled rows), size (relative square size)
forest_plot <- function(file, rows, xlim, xlab, title, subtitle = NULL, ref = 0,
                        xticks = NULL, xtick_labels = NULL, width = 2300, logx = FALSE,
                        label_cex = 0.82, secondary_axis = NULL) {
  n <- nrow(rows); H <- max(900, 120 * n + 560)
  png(file, width = width, height = H, res = 200, type = "cairo", bg = BG)
  on.exit(dev.off())
  par(mar = c(7.0, 17.5, 6.6, 11), xpd = NA, bg = BG, fg = INK, col.axis = INK2,
      col.lab = INK2, family = "sans")
  plot(NA, xlim = xlim, ylim = c(n + 0.6, 0.4), xlab = "", ylab = "", axes = FALSE,
       xaxs = "i", yaxs = "i", log = if (logx) "x" else "")
  usr <- par("usr")
  xl <- if (logx) 10^usr[1:2] else usr[1:2]
  xspan <- if (logx) NULL else diff(xl)
  xleft  <- if (logx) 10^(usr[1] - 0.025 * diff(usr[1:2])) else xl[1] - 0.025 * xspan
  xright <- if (logx) 10^(usr[2] + 0.025 * diff(usr[1:2])) else xl[2] + 0.025 * xspan
  # reference line
  if (!is.null(ref)) segments(ref, 0.5, ref, n + 0.5, col = GRID, lty = 2, lwd = 1.4)
  for (i in seq_len(n)) {
    r <- rows[i, ]
    y <- i
    if (r$type == "header") {
      text(xleft, y, r$label, adj = c(1, 0.5), font = 2, cex = label_cex + 0.04, col = INK)
    } else if (r$type == "note") {
      text(xleft, y, r$label, adj = c(1, 0.5), cex = label_cex - 0.08, col = INK2, font = 3)
      if (!is.na(r$annot)) text(xright, y, r$annot, adj = c(0, 0.5), cex = label_cex - 0.08, col = INK2, font = 3)
    } else if (r$type == "study") {
      text(xleft, y, r$label, adj = c(1, 0.5), cex = label_cex, col = INK)
      if (!is.na(r$lb) && !is.na(r$ub)) {
        lb <- max(r$lb, xl[1]); ub <- min(r$ub, xl[2])
        segments(lb, y, ub, y, col = INK2, lwd = 1.6)
        if (r$lb < xl[1]) arrows(xl[1] * 1.0001, y, xl[1], y, length = 0.06, col = INK2)
        if (r$ub > xl[2]) arrows(xl[2] * 0.9999, y, xl[2], y, length = 0.06, col = INK2)
      }
      sz <- if (is.na(r$size)) 1.2 else r$size
      points(r$est, y, pch = 15, cex = sz, col = INK2)
      text(xright, y, r$annot, adj = c(0, 0.5), cex = label_cex - 0.06, col = INK2)
    } else if (r$type == "pooled") {
      text(xleft, y, r$label, adj = c(1, 0.5), cex = label_cex, font = 2, col = INK)
      if (!is.na(r$plb)) {
        segments(max(r$plb, xl[1]), y + 0.32, min(r$pub, xl[2]), y + 0.32, col = ORANGE, lwd = 2.4)
        segments(c(max(r$plb, xl[1]), min(r$pub, xl[2])), y + 0.22,
                 c(max(r$plb, xl[1]), min(r$pub, xl[2])), y + 0.42, col = ORANGE, lwd = 2.4)
      }
      if (!is.na(r$lb))
        polygon(c(r$lb, r$est, r$ub, r$est), c(y, y - 0.28, y, y + 0.28),
                col = BLUE, border = BLUE)
      else points(r$est, y, pch = 18, cex = 2.2, col = BLUE)
      text(xright, y, r$annot, adj = c(0, 0.5), cex = label_cex - 0.06, font = 2, col = INK)
    }
  }
  # axis
  if (is.null(xticks)) xticks <- pretty(xlim)
  axis(1, at = xticks, labels = if (is.null(xtick_labels)) xticks else xtick_labels,
       col = GRID, col.ticks = GRID, cex.axis = 0.85, pos = n + 0.5)
  mtext(xlab, side = 1, line = 2.9, cex = 0.9, col = INK2)
  if (!is.null(secondary_axis)) mtext(secondary_axis, side = 1, line = 3.9, cex = 0.8, col = INK2)
  nsub <- length(subtitle)
  mtext(title, side = 3, line = 2.2 + 1.1 * nsub, adj = 0, at = xleft, cex = 1.05, font = 2, col = INK)
  if (!is.null(subtitle)) for (j in seq_len(nsub))
    mtext(subtitle[j], side = 3, line = 2.2 + 1.1 * (nsub - j), adj = 0, at = xleft, cex = 0.76, col = INK2)
  # legend (key)
  legend_y <- 0.0
  mtext("grey square = study estimate (area ~ weight) | blue diamond = pooled mean (95% CI) | orange bar = 95% prediction interval",
        side = 1, line = 5.6, cex = 0.7, col = INK2, adj = 0, at = xleft)
}
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# =============================================================================
# MA1  Realised genetic gain (% of mean yield per year)
# =============================================================================
g  <- read.csv(file.path(LIT, "extract/meta_genetic_gain/extraction.csv"))
gs <- read.csv(file.path(LIT, "extract/meta_genetic_gain/candidates_status.csv"))
g$verified_crossref <- tolower(gs$verified_crossref[match(g$doi, gs$doi)]) == "yes"
g$fulltext_read     <- gs$fulltext_read[match(g$doi, gs$doi)] == 1
g$live_crossref_ok  <- cr_ok(g$doi)
g$included          <- g$primary == "Y" & g$verified_crossref & g$fulltext_read & g$live_crossref_ok
# consistency of the extracted bibliographic metadata with the Crossref record
chk <- merge(unique(g[, c("doi", "crossref_first_author", "crossref_year")]),
             cr[, c("doi", "first_author_family", "issued_year")], by = "doi")
stopifnot(all(tolower(chk$crossref_first_author) == tolower(chk$first_author_family)),
          all(chk$crossref_year == chk$issued_year))
# effect size: gain as % of mean (or intercept) per year; SE on the same scale
g$yi_pct  <- ifelse(!is.na(g$pct_slope_over_ref), g$pct_slope_over_ref, g$pct_reported)
g$sei_pct <- ifelse(!is.na(g$se_pct_slope_over_ref), g$se_pct_slope_over_ref, g$se_pct)
g$stratum3 <- with(g, ifelse(stratum == "1_common_bean", "common bean",
                      ifelse(grepl("legume", stratum), "other grain legumes", "other crops")))
g$crop_group <- g$crop
# SE imputation rule (sensitivity analyses only; primary models use reported/derived SEs only):
#   SE_imputed = median (rule A) or 75th percentile (rule B, conservative) of the SEs of
#   all included rows (every stratum) that have an SE. Imputed rows are always flagged.
se_pool <- g$sei_pct[g$included & !is.na(g$sei_pct)]
SE_A <- median(se_pool); SE_B <- as.numeric(quantile(se_pool, 0.75))
cat(sprintf("MA1: SE imputation rule A (median of included SEs) = %.3f ; rule B (P75) = %.3f %%/yr\n", SE_A, SE_B))

# ---- 1a common bean ----------------------------------------------------------
cb <- subset(g, stratum3 == "common bean" & included & !is.na(sei_pct))
cb$vi <- cb$sei_pct^2; cb$yi <- cb$yi_pct; cb$cluster <- cb$cluster; cb$eff <- cb$effect_id
stopifnot(nrow(cb) >= 3)
add_used("MA1", transform(g, role = ifelse(included, "primary", "not in primary set")),
         c("effect_id", "stratum3", "cluster", "doi", "crossref_first_author", "crossref_year",
           "yi_pct", "sei_pct", "se_basis", "primary", "included", "verified_crossref",
           "fulltext_read", "live_crossref_ok", "location"))

r_cb3 <- fit_mv(cb)
s <- log_row("MA1.CB.primary", "Genetic gain, common bean: three-level REML (cluster/effect)", "common bean",
             "primary", "rma.mv REML, random=~1|programme/effect; t-based CI", r_cb3,
             length(unique(cb$cluster)), units = "% of mean yield per year", ci_method = "t (df = k-1)",
             notes = paste0("8 effects from 4 programmes (IAC Brazil x2 periods, Embrapa Brazil, Zeffa Brazil era trial, PABRA Africa x4 variety groups). ",
                            "Estimates are slope/mean of era or cycle regressions, SE from author tables or derived (see se_basis in meta_inputs_used.csv). ",
                            "CANDIDATE for gain_target_pct_per_yr."))
cat("MA1 common bean primary:", sprintf("%.3f [%.3f, %.3f] PI [%.3f, %.3f] tau2=%.3f I2=%.0f%%\n",
                                       s$est, s$lb, s$ub, s$plb, s$pub, s$tau2, s$I2))
# naive two-level (ignores dependence)
r_cb2 <- fit_uni(cb$yi, cb$vi)
log_row("MA1.CB.naive", "Genetic gain, common bean: two-level REML-KH, ignoring dependence", "common bean",
        "sensitivity", "rma REML, test=knha", r_cb2, length(unique(cb$cluster)),
        units = "% of mean yield per year", ci_method = "Knapp-Hartung t (df=k-1)")
# cluster-robust (CR1, df = clusters-1)
r_cb_rob <- robust(r_cb2, cluster = cb$cluster)
add_row("MA1.CB.robust", "Genetic gain, common bean: cluster-robust variance (CR1)", "common bean", "sensitivity",
        "rma REML + robust(cluster)", r_cb2$k, length(unique(cb$cluster)),
        as.numeric(r_cb_rob$beta), r_cb_rob$ci.lb, r_cb_rob$ci.ub, NA, NA, r_cb2$tau2, r_cb2$I2, r_cb2$QE, r_cb2$QEp,
        "% of mean yield per year", sprintf("sandwich t, df=%d (4 clusters: unreliable)", r_cb_rob$dfs),
        "With only 4 clusters the robust CI is indicative at best.")
# one effect per cluster (FE-pooled within cluster)
cl <- do.call(rbind, lapply(split(cb, cb$cluster), function(x) {
  w <- 1 / x$vi; data.frame(cluster = x$cluster[1], yi = sum(w * x$yi) / sum(w), vi = 1 / sum(w))
}))
r_cb_cl <- fit_uni(cl$yi, cl$vi)
log_row("MA1.CB.clustercollapsed", "Genetic gain, common bean: one estimate per programme", "common bean",
        "sensitivity", "FE-pool within programme, then rma REML-KH", r_cb_cl, nrow(cl),
        units = "% of mean yield per year", ci_method = "Knapp-Hartung t (df=3)",
        notes = "IAC periods 1989-96 (+1.07) and 1997-2007 (-0.25) are pooled by fixed effect here, which understates their disagreement.")
# Embrapa: swap direct (B04) for mixed-model-adjusted indirect estimate (B08)
b8 <- subset(g, effect_id == "B08")
cb_sw <- cb; i4 <- which(cb_sw$effect_id == "B04")
cb_sw$yi[i4] <- b8$yi_pct; cb_sw$vi[i4] <- b8$sei_pct^2
log_row("MA1.CB.embrapaB08", "Genetic gain, common bean: Embrapa direct estimate replaced by mixed-model indirect (B08)", "common bean",
        "sensitivity", "rma.mv REML", fit_mv(cb_sw), length(unique(cb_sw$cluster)),
        units = "% of mean yield per year", ci_method = "t",
        notes = "Embrapa 1993-2008: direct 1.34 vs 0.45 %/yr by method (de Faria 2018 Table 3); shows method dependence.")
# drop IAC 1997-2007 (regime change in breeding goals)
cb_b2 <- subset(cb, effect_id != "B02")
log_row("MA1.CB.exB02", "Genetic gain, common bean: excluding IAC 1997-2007 (null window after change of breeding goals)", "common bean",
        "sensitivity", "rma.mv REML", fit_mv(cb_b2), length(unique(cb_b2$cluster)),
        units = "% of mean yield per year", ci_method = "t")
# imputed SE: add B17 (Brazil black-bean cultivars, no SE; slope conditional on a 1988 breakpoint -> upward-biased)
for (rule in c("A", "B")) {
  cb_imp <- rbind(cb[, c("effect_id", "cluster", "yi", "vi", "eff")],
                  with(subset(g, effect_id == "B17"),
                       data.frame(effect_id = effect_id, cluster = cluster, yi = yi_pct,
                                  vi = (if (rule == "A") SE_A else SE_B)^2, eff = effect_id)))
  log_row(paste0("MA1.CB.imputed", rule), sprintf("Genetic gain, common bean: + B17 with imputed SE (rule %s = %.2f)", rule, if (rule == "A") SE_A else SE_B),
          "common bean", "sensitivity", "rma.mv REML", fit_mv(cb_imp), length(unique(cb_imp$cluster)),
          units = "% of mean yield per year", ci_method = "t",
          notes = "B17 = post-1988 segment slope of a bi-segmented model (upward bias); SE imputed.")
}
# leave-one-programme-out
for (cl_i in unique(cb$cluster)) {
  x <- subset(cb, cluster != cl_i); r <- fit_mv(x); p <- predict(r)
  add_loo("MA1.CB.primary", cl_i, as.numeric(r$beta), p$ci.lb, p$ci.ub)
}
# ---- 1b other crops (SSA maize/rice) ------------------------------------------
oc_all <- subset(g, stratum3 == "other crops")
oc <- subset(oc_all, included & !is.na(sei_pct)); oc$vi <- oc$sei_pct^2; oc$yi <- oc$yi_pct
stopifnot(nrow(oc) >= 3)
r_oc <- fit_uni(oc$yi, oc$vi)
s <- log_row("MA1.OC.primary", "Genetic gain, other crops (SSA maize and rice): REML-KH", "other crops", "primary",
             "rma REML, test=knha (one primary row per programme -> independent)", r_oc, nrow(oc),
             units = "% of mean yield per year", ci_method = "Knapp-Hartung t (df=k-1)",
             notes = "Primary rows: Zimbabwe national maize AVT (M02), CIMMYT ESA early hybrids random stress (M10), Uganda NPT maize (A01), AfricaRice irrigated lowland (M14). Maize dominates (3/4).")
cat("MA1 other crops primary:", sprintf("%.3f [%.3f, %.3f] PI [%.3f, %.3f]\n", s$est, s$lb, s$ub, s$plb, s$pub))
ocm <- subset(oc, !grepl("rice", crop))
log_row("MA1.OC.maizeOnly", "Genetic gain, SSA maize only", "other crops", "sensitivity", "rma REML-KH",
        fit_uni(ocm$yi, ocm$vi), nrow(ocm), units = "% of mean yield per year", ci_method = "Knapp-Hartung t")
for (i in seq_len(nrow(oc))) {
  r <- fit_uni(oc$yi[-i], oc$vi[-i]); p <- predict(r)
  add_loo("MA1.OC.primary", oc$cluster[i], as.numeric(r$beta), p$ci.lb, p$ci.ub)
}
# imputed: add A06 (IITA West Africa maize, no SE)
a6 <- subset(g, effect_id == "A06")
for (rule in c("A", "B")) {
  se <- if (rule == "A") SE_A else SE_B
  x <- rbind(oc[, c("effect_id", "yi", "vi")], data.frame(effect_id = "A06", yi = a6$yi_pct, vi = se^2))
  log_row(paste0("MA1.OC.imputed", rule), sprintf("Genetic gain, other crops: + A06 with imputed SE (rule %s = %.2f)", rule, se),
          "other crops", "sensitivity", "rma REML-KH", fit_uni(x$yi, x$vi), nrow(x),
          units = "% of mean yield per year", ci_method = "Knapp-Hartung t")
}
# all SE-bearing rows of the stratum, three-level (overlapping sub-analyses of the same programmes)
oca <- subset(oc_all, !is.na(sei_pct) & verified_crossref & fulltext_read & live_crossref_ok)
oca$vi <- oca$sei_pct^2; oca$yi <- oca$yi_pct; oca$eff <- oca$effect_id
log_row("MA1.OC.allrows", "Genetic gain, other crops: ALL rows with SE (incl. non-primary sub-analyses), three-level",
        "other crops", "sensitivity", "rma.mv REML, random=~1|programme/effect", fit_mv(oca),
        length(unique(oca$cluster)), units = "% of mean yield per year", ci_method = "t",
        notes = "Sub-rows overlap (e.g. IVT/AVT, management regimes); weights are not independent. Includes low-quality Michael 2020 rows.")
add_used("MA1-allrows", transform(oca, role = "sensitivity (all SE rows)"),
         c("effect_id", "stratum3", "cluster", "doi", "yi_pct", "sei_pct", "role"))
# ---- 1c other grain legumes: descriptive ------------------------------------------
lg <- subset(g, stratum3 == "other grain legumes")
lg_prim <- subset(lg, included)
desc_row("MA1.LG.primary", "Genetic gain, other grain legumes (primary rows: soybean Brazil L01 only; no SE)",
         "other grain legumes", lg_prim$yi_pct, "% of mean yield per year",
         "Only 1 usable primary effect (soybean South Brazil 1965-2011, 2.4 %/yr, no SE). SSA legume papers (cowpea, chickpea, lentil, Uganda beans) could not be obtained as legal full text.", 1)
lg_all <- subset(lg, verified_crossref & fulltext_read & live_crossref_ok)
desc_row("MA1.LG.sens_all", "Genetic gain, other grain legumes: soybean L01 + ICRISAT groundnut rows L02-L06 (no SE)",
         "other grain legumes", lg_all$yi_pct, "% of mean yield per year",
         "Groundnut rows are five trials of one ICRISAT programme (India), outside SSA/LatAm. Median across 6 rows shown.", 2)
x <- rbind(data.frame(cluster = lg_all$cluster, eff = lg_all$effect_id, yi = lg_all$yi_pct,
                      vi = SE_A^2))
r_lg <- fit_mv(x)
log_row("MA1.LG.sens_imputedA", sprintf("Genetic gain, other grain legumes: L01 + groundnut with imputed SE (rule A = %.2f), three-level", SE_A),
        "other grain legumes", "sensitivity", "rma.mv REML (all SEs imputed)", r_lg, length(unique(x$cluster)),
        units = "% of mean yield per year", ci_method = "t (df=k-1)",
        notes = "ALL SEs imputed and only 2 independent programmes -> indicative only; do not use as a parameter.")
# ---- 1d moderator test common bean vs other crops (SE-bearing primary rows)
mod <- rbind(data.frame(cluster = cb$cluster, eff = cb$effect_id, yi = cb$yi, vi = cb$vi, grp = "common bean"),
             data.frame(cluster = oc$cluster, eff = oc$effect_id, yi = oc$yi, vi = oc$vi, grp = "other crops"))
r_mod <- rma.mv(yi, vi, mods = ~ grp, random = ~ 1 | cluster / eff, data = mod, method = "REML", test = "t")
add_row("MA1.MOD", "Genetic gain: moderator test, other crops minus common bean (primary SE-bearing rows)",
        "common bean vs other crops", "sensitivity", "rma.mv REML, mods=~stratum",
        r_mod$k, length(unique(mod$cluster)), as.numeric(r_mod$beta[2]), r_mod$ci.lb[2], r_mod$ci.ub[2], NA, NA,
        sum(r_mod$sigma2), NA, r_mod$QM, r_mod$QMp, "difference in % of mean yield per year", "t",
        "Estimate = other crops - common bean; Q column = omnibus QM test of the moderator, Q_p its p-value.")

# ---- MA1 forest plot -------------------------------------------------------------
mk_rows <- function(d, grp_label) {
  w <- 1 / d$vi; sz <- 0.8 + 1.6 * (w / max(w))
  data.frame(label = paste0(d$crossref_first_author, " ", d$crossref_year, ifelse(grepl("^Amongi", d$cluster), "*", ""), " - ", d$effect_id),
             est = d$yi, lb = d$yi - 1.96 * sqrt(d$vi), ub = d$yi + 1.96 * sqrt(d$vi),
             type = "study", annot = sprintf("%s [%s, %s]", fmt(d$yi), fmt(d$yi - 1.96 * sqrt(d$vi)), fmt(d$yi + 1.96 * sqrt(d$vi))),
             size = sz, plb = NA, pub = NA)
}
lab_cb <- cb; lab_cb$crossref_first_author <- cb$crossref_first_author; lab_cb$crossref_year <- ifelse(grepl("^Amongi", cb$cluster), 2023, cb$crossref_year)
sm <- summ(r_cb3); so <- summ(r_oc)
rows <- rbind(
  data.frame(label = "Common bean (k = 8, 4 programmes)", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  mk_rows(lab_cb, ""),
  data.frame(label = "Pooled, three-level RE", est = sm$est, lb = sm$lb, ub = sm$ub, type = "pooled",
             annot = sprintf("%s [%s, %s]", fmt(sm$est), fmt(sm$lb), fmt(sm$ub)), size = NA, plb = sm$plb, pub = sm$pub),
  data.frame(label = "", est = NA, lb = NA, ub = NA, type = "note",
             annot = sprintf("PI [%s, %s]; tau=%s; I2=%s%%", fmt(sm$plb), fmt(sm$pub), fmt(sqrt(sm$tau2)), fmt(sm$I2, 0)), size = NA, plb = NA, pub = NA),
  data.frame(label = "Other crops: SSA maize and rice (k = 4)", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  mk_rows(transform(oc, crossref_first_author = crossref_first_author), ""),
  data.frame(label = "Pooled, RE (Knapp-Hartung)", est = so$est, lb = so$lb, ub = so$ub, type = "pooled",
             annot = sprintf("%s [%s, %s]", fmt(so$est), fmt(so$lb), fmt(so$ub)), size = NA, plb = so$plb, pub = so$pub),
  data.frame(label = "", est = NA, lb = NA, ub = NA, type = "note",
             annot = sprintf("PI [%s, %s]; tau=%s; I2=%s%%", fmt(so$plb), fmt(so$pub), fmt(sqrt(so$tau2)), fmt(so$I2, 0)), size = NA, plb = NA, pub = NA),
  data.frame(label = "Other grain legumes (k = 1: not pooled)", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  data.frame(label = "Todeschini 2019 - L01 soybean (no SE)", est = lg_prim$yi_pct, lb = NA, ub = NA, type = "study",
             annot = sprintf("%s [no SE]", fmt(lg_prim$yi_pct)), size = 1.2, plb = NA, pub = NA)
)
forest_plot(file.path(FIG, "fig_ma1_forest_genetic_gain.png"), rows, xlim = c(-2, 5),
            xlab = "Realised genetic gain in yield (% of mean yield per year)",
            title = "Realised yield gain per year, by stratum",
            subtitle = c("Historical / era-trial regressions of yield on release year or selection cycle; dashed line = no gain.",
                         "* Amongi et al. is printed 2023 (Crossref issued-year 1970 is an African Journals Online placeholder)."),
            ref = 0, xticks = -2:5)

# =============================================================================
# MA2  Phenotyping accuracy (correlation of image/CV tool with manual reference)
# =============================================================================
a <- read.csv(file.path(LIT, "extract/meta_pheno_accuracy/extraction.csv"))
a$verified_crossref <- as.logical(a$verified_crossref); a$fulltext_read <- as.logical(a$fulltext_read)
a$live_crossref_ok <- cr_ok(a$doi)
a$preprint <- is_preprint(a$doi)
a$included <- a$analysis_set == "primary" & a$verified_crossref & a$fulltext_read & a$live_crossref_ok
a$row <- seq_len(nrow(a)); a$study <- a$study_id
a$legume <- ifelse(grepl("legume", a$crop_group) & !grepl("non-legume", a$crop_group), "legume", "non-legume")
a$z <- ifelse(!is.na(a$r) & abs(a$r) < 1, atanh(a$r), NA)
chk <- merge(unique(a[, c("doi", "crossref_first_author", "crossref_year")]),
             cr[, c("doi", "first_author_family", "issued_year")], by = "doi")
stopifnot(all(tolower(chk$crossref_first_author) == tolower(chk$first_author_family)))
add_used("MA2", a, c("study_id", "doi", "crossref_first_author", "crossref_year", "label", "trait_class", "tool_class", "metric",
                      "r", "n", "analysis_set", "included", "verified_crossref", "fulltext_read", "live_crossref_ok", "preprint", "location"))
d2 <- subset(a, included & !is.na(n) & !is.na(z))
d2$vi <- 1 / (d2$n - 3); d2$yi <- d2$z; d2$cluster <- d2$study; d2$eff <- d2$row
cat(sprintf("\nMA2: %d primary effects from %d studies\n", nrow(d2), length(unique(d2$study))))

r2_mv <- fit_mv(d2)
s <- log_row("MA2.primary", "Phenotyping accuracy: tool-vs-manual correlation, three-level REML on Fisher z (back-transformed r)",
             "all crops/traits", "primary", "rma.mv REML on Fisher z, random=~1|study/effect, v=1/(n-3)",
             r2_mv, length(unique(d2$study)), tf = tanh, units = "Pearson-equivalent r (Spearman rho and sqrt(R2) treated as r)",
             ci_method = "t (df=k-1), back-transformed",
             notes = "Primary = analysis_set 'primary' with known n. Mix of Pearson, Spearman, sqrt(R2) (sign assumed positive); reference standards differ (ruler, manual count, visual score). I2 on the z scale.")
cat("MA2 primary pooled r:", sprintf("%.3f [%.3f, %.3f] PI [%.3f, %.3f] tau2(z)=%.3f I2=%.0f%%\n", s$est, s$lb, s$ub, s$plb, s$pub, s$tau2, s$I2))
MU2 <- as.numeric(r2_mv$beta); SEMU2 <- r2_mv$se
log_row("MA2.naive", "Phenotyping accuracy: two-level REML-KH on rows (ignores dependence)", "all crops/traits", "sensitivity",
        "rma REML, test=knha", fit_uni(d2$yi, d2$vi), length(unique(d2$study)), tf = tanh, units = "r", ci_method = "Knapp-Hartung t")
# study-level (FE-pool within study)
st <- do.call(rbind, lapply(split(d2, d2$study), function(x) {
  w <- 1 / x$vi
  data.frame(study = x$study[1], yi = sum(w * x$yi) / sum(w), vi = 1 / sum(w), n = sum(x$n), legume = x$legume[1],
             trait = x$trait_class[1], tool = x$tool_class[1], preprint = any(x$preprint),
             first_author = x$crossref_first_author[1], year = x$crossref_year[1], doi = x$doi[1], nrows = nrow(x))
}))
r2_st <- fit_uni(st$yi, st$vi)
s_st <- log_row("MA2.studylevel", "Phenotyping accuracy: one effect per study (FE-pooled within study), REML-KH", "all crops/traits",
                "sensitivity", "rma REML-KH on study-level z", r2_st, nrow(st), tf = tanh, units = "r", ci_method = "Knapp-Hartung t (df=13)")
# subgroups (three-level within each subgroup)
sub_fit <- function(dd, id, label, role = "sensitivity", note = "") {
  if (length(unique(dd$study)) < 3) {
    desc_row(id, label, "subgroup", tanh(dd$yi), "r", paste0("fewer than 3 independent studies. ", note), length(unique(dd$study)))
    return(invisible(NULL))
  }
  r <- if (length(unique(dd$study)) == nrow(dd)) fit_uni(dd$yi, dd$vi) else fit_mv(dd)
  log_row(id, label, "subgroup", role, if (inherits(r, "rma.mv")) "rma.mv REML" else "rma REML-KH", r,
          length(unique(dd$study)), tf = tanh, units = "r", ci_method = "t", notes = note)
}
sub_fit(subset(d2, legume == "legume"), "MA2.legume", "Accuracy: legumes (bean, soybean, faba)")
sub_fit(subset(d2, legume == "non-legume"), "MA2.nonlegume", "Accuracy: non-legumes (wheat, maize)")
r_leg <- rma.mv(yi, vi, mods = ~ legume, random = ~ 1 | cluster / eff, data = d2, method = "REML", test = "t")
add_row("MA2.MOD.legume", "Accuracy: moderator test legume vs non-legume (difference in Fisher z, non-legume - legume)", "subgroup",
        "sensitivity", "rma.mv REML, mods=~legume", r_leg$k, length(unique(d2$study)), as.numeric(r_leg$beta[2]), r_leg$ci.lb[2], r_leg$ci.ub[2],
        NA, NA, sum(r_leg$sigma2), NA, r_leg$QM, r_leg$QMp, "difference in Fisher z", "t", "Q = QM omnibus, Q_p its p-value.")
for (tc in c("plant count", "height", "disease severity", "seed size")) {
  sub_fit(subset(d2, trait_class == tc), paste0("MA2.trait.", gsub(" ", "_", tc)), paste0("Accuracy: trait = ", tc))
}
for (tc in c("smartphone", "UAV-RGB", "camera")) {
  dd <- subset(d2, grepl(tc, tool_class, ignore.case = TRUE))
  if (nrow(dd) > 0) sub_fit(dd, paste0("MA2.tool.", tc), paste0("Accuracy: tool class matching '", tc, "'"))
}
# bean-only (Volpato dry bean, Parrella common bean) + Lasdun (no n -> imputed)
bean <- subset(a, grepl("bean", crop, ignore.case = TRUE) & !grepl("faba", crop, ignore.case = TRUE) & !is.na(z) &
                 analysis_set %in% c("primary", "primary_no_n") & verified_crossref & fulltext_read)
bean$n_use <- ifelse(is.na(bean$n), 30, bean$n)   # imputed n = 30 where missing (conservative)
bean$vi <- 1 / (bean$n_use - 3); bean$yi <- bean$z; bean$cluster <- bean$study; bean$eff <- bean$row
r_bean <- fit_mv(bean)
log_row("MA2.beanonly_imputed", "Accuracy: common/dry bean only (Volpato, Parrella, Lasdun); missing n imputed = 30", "common bean",
        "sensitivity", "rma.mv REML", r_bean, length(unique(bean$study)), tf = tanh, units = "r", ci_method = "t",
        notes = "Lasdun 2024 reports no n (7 stand-count r values); n=30 assumed. Volpato is a preprint. 3 studies only: indicative.")
# sensitivity: impute n for rows without n
for (nv in c(30, round(median(d2$n)))) {
  imp <- subset(a, analysis_set %in% c("primary", "primary_no_n") & verified_crossref & fulltext_read & !is.na(z))
  imp$n_use <- ifelse(is.na(imp$n), nv, imp$n); imp$vi <- 1 / (imp$n_use - 3); imp$yi <- imp$z
  imp$cluster <- imp$study; imp$eff <- imp$row
  log_row(paste0("MA2.imputed_n", nv), sprintf("Accuracy: primary + 'primary_no_n' rows (n imputed = %d)", nv), "all crops/traits",
          "sensitivity", "rma.mv REML", fit_mv(imp), length(unique(imp$cluster)), tf = tanh, units = "r", ci_method = "t",
          notes = "Adds Lasdun (7 rows) and Volpato height; n imputed where missing.")
}
# excluding preprints / Spearman / sqrt(R2) conversions
dd <- subset(d2, !preprint)
log_row("MA2.nopreprint", "Accuracy: peer-reviewed studies only (preprints removed)", "all crops/traits", "sensitivity", "rma.mv REML",
        fit_mv(dd), length(unique(dd$study)), tf = tanh, units = "r", ci_method = "t",
        notes = "Removed: Volpato 2023, Zang 2022, Karisto 2017, Gijare 2026 (Crossref type posted-content).")
dd <- subset(d2, metric == "Pearson r")
log_row("MA2.pearsononly", "Accuracy: Pearson r rows only (no Spearman, no sqrt(R2))", "all crops/traits", "sensitivity", "rma.mv REML",
        fit_mv(dd), length(unique(dd$study)), tf = tanh, units = "r", ci_method = "t")
# leave-one-study-out
for (sid in unique(d2$study)) {
  x <- subset(d2, study != sid); r <- fit_mv(x); p <- predict(r)
  add_loo("MA2.primary", sid, tanh(as.numeric(r$beta)), tanh(p$ci.lb), tanh(p$ci.ub))
}
# small-study effects: study level, k >= 10
cat(sprintf("MA2 small-study checks on %d study-level effects\n", nrow(st)))
if (nrow(st) >= 10) {
  eg <- regtest(r2_st, model = "lm", predictor = "sei")
  rk <- ranktest(r2_st)
  tf <- trimfill(r2_st)
  add_row("MA2.egger", "Accuracy: small-study effects, Egger-type regression test (study level)", "all crops/traits", "sensitivity",
          "regtest(model='lm', predictor='sei')", nrow(st), nrow(st), tanh(as.numeric(eg$est)), NA, NA, NA, NA, NA, NA, NA, eg$pval, "r at se->0 (regression limit estimate); p of the asymmetry test in Q_p column", "-",
          sprintf("Egger-type test z=%.2f, p=%.3f, k=%d studies (low power; no evidence of asymmetry). Estimate column = limit estimate as se->0, back-transformed.", eg$zval, eg$pval, nrow(st)))
  add_row("MA2.ranktest", "Accuracy: small-study effects, Begg rank correlation (study level)", "all crops/traits", "sensitivity",
          "ranktest", nrow(st), nrow(st), as.numeric(rk$tau), NA, NA, NA, NA, NA, NA, NA, rk$pval, "Kendall tau; p in Q_p column", "-", "")
  r_tf <- tf
  add_row("MA2.trimfill", "Accuracy: trim-and-fill adjusted pooled r (study level, R0 estimator)", "all crops/traits", "sensitivity",
          "trimfill", r_tf$k, nrow(st), tanh(as.numeric(r_tf$beta)), tanh(r_tf$ci.lb), tanh(r_tf$ci.ub), NA, NA, r_tf$tau2, r_tf$I2, r_tf$QE, r_tf$QEp,
          "r", "t/z", sprintf("%d studies imputed on the %s side. Trim-and-fill is exploratory; extreme heterogeneity invalidates its assumptions.", r_tf$k - nrow(st), r_tf$side))
  png(file.path(FIG, "fig_ma2_funnel.png"), width = 1500, height = 1300, res = 200, type = "cairo", bg = BG)
  par(mar = c(4.5, 4.5, 3.8, 1), bg = BG, fg = INK, col.axis = INK2, col.lab = INK2)
  funnel(r2_st, xlab = "Fisher z (tool vs manual reference)", ylab = "Standard error",
         level = c(90, 95, 99), shade = c("#efeee9", "#e3e2da", "#d3d2c8"), refline = as.numeric(r2_st$beta),
         pch = 19, col = INK2, back = BG, hlines = NA)
  mtext("Funnel plot, one point per study (k = 14)", side = 3, line = 2.3, adj = 0, font = 2, cex = 1.0)
  mtext(sprintf("Egger-type test p = %.2f; very high heterogeneity (I2 = %.0f%%), so a funnel is not expected to be symmetric.", eg$pval, r2_st$I2),
        side = 3, line = 1.0, adj = 0, cex = 0.72, col = INK2)
  dev.off()
}
# ---- manual reliability (context for the translation to error reduction)
rel <- subset(a, analysis_set == "context_manual_reliability" & verified_crossref & fulltext_read & !is.na(z))
# keep rater-reliability style rows (rater vs rater / rater vs measured truth / human annotation vs field count);
# drop Dobbels visual date1-vs-date2 stability (r=0.80) because it conflates rater error with real change
rel <- subset(rel, !(study == "Dobbels2019" & r == 0.80))
rel$n_use <- rel$n
relst <- do.call(rbind, lapply(split(rel, rel$study), function(x) {
  data.frame(study = x$study[1], yi = mean(x$z), vi = 1 / (mean(x$n) - 3), nrows = nrow(x), r_mean = tanh(mean(x$z)))
}))
cat("\nManual-reliability inputs (study-level):\n"); print(relst)
r_rel <- fit_uni(relst$yi, relst$vi)
s_rel <- log_row("MA2.manualrel", "Manual-scoring reliability (rater-vs-rater or rater-vs-measured truth), RE on Fisher z", "context",
                 "primary", "rma REML-KH, 3 studies", r_rel, nrow(relst), tf = tanh, units = "r",
                 ci_method = "Knapp-Hartung t (df=2)",
                 notes = "Studies: Librelon 2015 (bean angular leaf spot; 7 rows averaged), Dobbels 2019 (soybean IDC inter-rater), Volpato 2023 preprint (human annotation vs field count). k=3: indicative. Bock 2021 benchmark (full text read): trained raters/standard area diagrams Lin's rho_c 0.85-0.95; inexperienced <=0.60.")
cat("Manual reliability pooled r:", sprintf("%.3f [%.3f, %.3f]\n", s_rel$est, s_rel$lb, s_rel$ub))
add_used("MA2-manualrel", rel, c("study_id", "doi", "label", "r", "n", "analysis_set", "location"))

# ---- implied measurement-error reduction -----------------------------------------
# Notation: rho_x = correlation of measurement x with the (latent) true value; error share E_x = 1 - rho_x^2
# (error variance as a fraction of measured variance). Model parameter error_reduction = proportional
# reduction of plot error variance, ER = 1 - E_tool / E_manual.
#  T1 (reference treated as error-free): rho_tool = r_pool.
#  T2 (reference is manual with error, independent errors): r_pool = rho_tool * rho_man => rho_tool = r_pool/rho_man.
#  Break-even (ER = 0):  T2: rho_man = sqrt(r_pool).
er_fun <- function(r, rho_m, T = "T2") {
  rho_t <- if (T == "T1") r else r / rho_m
  rho_t[rho_t > 1] <- NA            # T2 is only defined when rho_manual >= r (r = rho_tool * rho_manual <= rho_manual)
  1 - (1 - rho_t^2) / (1 - rho_m^2)
}
r_pool <- s$est
cat(sprintf("Break-even manual reliability (T2) = sqrt(r_pool) = %.3f\n", sqrt(r_pool)))
rho_grid <- c(0.60, 0.70, 0.80, 0.85, 0.90, 0.95)
tab_er <- expand.grid(rho_m = rho_grid, T = c("T1", "T2"))
tab_er$ER_at_pooled_r <- mapply(function(rm, T) er_fun(r_pool, rm, T), tab_er$rho_m, tab_er$T)
tab_er$ER_at_r_lb <- mapply(function(rm, T) er_fun(s$lb, rm, T), tab_er$rho_m, tab_er$T)
tab_er$ER_at_r_ub <- mapply(function(rm, T) er_fun(s$ub, rm, T), tab_er$rho_m, tab_er$T)
write.csv(tab_er, file.path(LIT, "meta_error_reduction_grid.csv"), row.names = FALSE)
# Monte-Carlo propagation of both pooled correlations (means, not individual studies)
N <- 1e5
zt <- rnorm(N, MU2, SEMU2)
zm <- rnorm(N, as.numeric(r_rel$beta), r_rel$se)
rt <- tanh(zt); rm_ <- tanh(zm)
mc <- function(T) {
  e <- er_fun(rt, rm_, T); infeasible <- mean(is.na(e)); e <- e[!is.na(e)]    # T2: discard draws with r > rho_manual
  c(med = median(e), lo = quantile(e, .025), hi = quantile(e, .975), p_pos = mean(e > 0), p_ge10 = mean(e >= 0.10), infeasible = infeasible)
}
mc_T1 <- mc("T1"); mc_T2 <- mc("T2")
for (nm in c("T1", "T2")) {
  m <- if (nm == "T1") mc_T1 else mc_T2
  add_row(paste0("MA2.ER.", nm), paste0("Implied error-variance reduction vs manual, translation ", nm),
          "translation", "sensitivity", "Monte-Carlo over pooled r(tool) and pooled r(manual reliability), 1e5 draws",
          NA, NA, m["med"], m["lo.2.5%"], m["hi.97.5%"], NA, NA, NA, NA, NA, NA, "fraction of plot error variance", "2.5-97.5% of MC draws",
          sprintf("P(ER>0)=%.2f, P(ER>=0.10)=%.2f, share of draws infeasible (r > rho_manual) and discarded = %.2f. %s", m["p_pos"], m["p_ge10"], m["infeasible"],
                  if (nm == "T1") "T1: reference treated as error-free -> tool error share = 1 - r^2 (pessimistic for tools)."
                  else "T2: reference manual with reliability from MA2.manualrel, independent errors (central)."))
}
# ---- REVISION: manual-reliability rows estimate three different quantities -----------------
# The 3-study pool above (MA2.manualrel) enters every row as a correlation with the true value (a validity
# coefficient rho_m). The rows are in fact of three kinds:
#   rater-vs-truth r      : correlation of a visual estimate with a digital measurement of actual severity
#                           -> estimates rho_m directly
#   rater-vs-rater r      : correlation of two raters (or of a human annotation and a manual count)
#                           -> for parallel measures with independent errors r_rr = rho_m^2, so rho_m = sqrt(r_rr)
#   rater-vs-rater R2     : squared correlation of two ratings (intra-rater repeatability, mean pairwise R2)
#                           -> r_rr = sqrt(R2), rho_m = sqrt(r_rr) = R2^(1/4)
# MA2.manualrel stays as the "as entered" pool. MA2.manualrel.validity converts every row to rho_m before pooling.
# (The new Monte-Carlo block restores the random-number stream afterwards, so no earlier result changes.)
rel$rel_type <- ifelse(rel$study == "Librelon2015" & grepl("intra-rater|reproducibility", rel$label), "rater-vs-rater R2",
                ifelse(rel$study == "Librelon2015" & grepl("actual severity", rel$label), "rater-vs-truth r", "rater-vs-rater r"))
rel$rho_m <- ifelse(rel$rel_type == "rater-vs-truth r", rel$r,
              ifelse(rel$rel_type == "rater-vs-rater r", sqrt(rel$r), rel$r^(1/4)))
rel$z_val <- atanh(rel$rho_m)
cat("\nManual-reliability rows by type (validity conversion):\n"); print(rel[, c("study", "rel_type", "r", "rho_m", "n")])
relv <- do.call(rbind, lapply(split(rel, rel$study), function(x) {
  data.frame(study = x$study[1], yi = mean(x$z_val), vi = 1 / (mean(x$n) - 3), nrows = nrow(x), r_mean = tanh(mean(x$z_val)))
}))
r_relv <- fit_uni(relv$yi, relv$vi)
s_relv <- log_row("MA2.manualrel.validity", "Manual-scoring validity coefficient after converting rater-vs-rater r (sqrt) and R2 (fourth root), RE on Fisher z", "context",
                  "sensitivity", "rma REML-KH, 3 studies", r_relv, nrow(relv), tf = tanh, units = "r (correlation with the true value)",
                  ci_method = "Knapp-Hartung t (df=2)",
                  notes = "Same studies as MA2.manualrel; each row converted to a correlation with the true value: rater-vs-truth r kept, rater-vs-rater r -> sqrt(r), rater-vs-rater R2 -> R2^(1/4). Assumes parallel raters with independent errors; intra-rater repeatability also contains persistent rater bias, so it overstates validity.")
cat("Manual validity pooled rho_m:", sprintf("%.3f [%.3f, %.3f]\n", s_relv$est, s_relv$lb, s_relv$ub))
tr <- subset(rel, rel_type == "rater-vs-truth r")
desc_row("MA2.manualrel.truthonly", "Manual-scoring correlation with a digital measure of the truth (rater-vs-truth rows only, one study)", "context",
         tr$r, "r", sprintf("Librelon 2015 only (%d rows: rater experience x diagrammatic scale); no pooling possible.", nrow(tr)), 1)
{  # Monte Carlo for T2 with the validity coefficient; RNG state restored afterwards
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(20261006)
  zm_v <- rnorm(N, as.numeric(r_relv$beta), r_relv$se)
  e_v <- er_fun(rt, tanh(zm_v), "T2"); inf_v <- mean(is.na(e_v)); e_v <- e_v[!is.na(e_v)]
  mc_T2v <- c(med = median(e_v), lo = unname(quantile(e_v, .025)), hi = unname(quantile(e_v, .975)),
              p_pos = mean(e_v > 0), p_ge10 = mean(e_v >= 0.10), infeasible = inf_v)
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
}
add_row("MA2.ER.T2v", "Implied error-variance reduction vs manual, translation T2 with the validity coefficient of manual scoring",
        "translation", "sensitivity", "Monte-Carlo over pooled r(tool) and pooled manual validity, 1e5 draws",
        NA, NA, mc_T2v["med"], mc_T2v["lo"], mc_T2v["hi"], NA, NA, NA, NA, NA, NA, "fraction of plot error variance", "2.5-97.5% of MC draws",
        sprintf("P(ER>0)=%.2f, P(ER>=0.10)=%.2f, share of draws infeasible (r > rho_manual) and discarded = %.2f. Break-even manual validity is sqrt(r_pool).",
                mc_T2v["p_pos"], mc_T2v["p_ge10"], mc_T2v["infeasible"]))
cat(sprintf("T2 with validity: median ER %.2f [%.2f, %.2f], P(ER>0)=%.2f, infeasible=%.3f\n", mc_T2v["med"], mc_T2v["lo"], mc_T2v["hi"], mc_T2v["p_pos"], mc_T2v["infeasible"]))
add_used("MA2-manualrel-validity", rel, c("study_id", "doi", "label", "r", "n", "analysis_set", "location", "rel_type", "rho_m"))

# direct estimate from heritability pairs (primary_h2 rows): error-to-signal ratio (1-h2)/h2
h2 <- subset(a, analysis_set == "primary_h2" & verified_crossref & fulltext_read & !is.na(errvar_ratio_image_over_manual))
h2$lr <- log(h2$errvar_ratio_image_over_manual)
h2s <- do.call(rbind, lapply(split(h2, h2$study), function(x) data.frame(study = x$study[1], yi = mean(x$lr), nrows = nrow(x),
                                                                      ratio = exp(mean(x$lr)))))
cat("\nH2-based noise/signal ratios (image/manual) by study:\n"); print(h2s)
add_used("MA2-h2", h2, c("study_id", "doi", "label", "h2_image", "h2_manual", "errvar_ratio_image_over_manual", "location"))
r_h2 <- rma(h2s$yi, vi = rep(1e-8, nrow(h2s)), method = "REML", test = "knha")
log_row("MA2.ER.h2", "Error-to-signal ratio (image/manual) from heritability pairs; reported as error reduction = 1 - ratio",
        "translation", "sensitivity", "unweighted RE on log ratio (no variances reported; nominal v=1e-8)", r_h2, nrow(h2s),
        tf = function(x) 1 - exp(x), units = "fraction (positive = image tool has LESS error per unit signal)",
        ci_method = "Knapp-Hartung t (df=3); PI/CI limits flipped by transform",
        notes = "Studies: Zang 2023 (wheat height, UAV), Makanza 2018 (maize senescence, UAV), Gijare 2026 preprint (soybean seed area, two environments averaged), Gage 2018 (maize tassel, image traits from 1 env vs manual 3 env). Burner 2025 excluded (proxy index). Median study ratio shown separately.")
# the transform 1 - exp(x) reverses the order of CI/PI limits: put them back in ascending order
sh <- summ(r_h2, function(x) 1 - exp(x))
RES[[length(RES)]]$ci_lb <- min(sh$lb, sh$ub); RES[[length(RES)]]$ci_ub <- max(sh$lb, sh$ub)
RES[[length(RES)]]$pi_lb <- min(sh$plb, sh$pub); RES[[length(RES)]]$pi_ub <- max(sh$plb, sh$pub)
desc_row("MA2.ER.h2.median", "Error reduction from H2 pairs: median across 4 studies (range)", "translation",
         1 - h2s$ratio, "fraction", "Median is the robust summary because Gage 2018 is a large negative outlier (image traits from one environment).", nrow(h2s))
# error-reduction figure
png(file.path(FIG, "fig_ma2_error_reduction.png"), width = 1900, height = 1300, res = 200, type = "cairo", bg = BG)
par(mar = c(4.8, 4.8, 4.8, 1.2), bg = BG, fg = INK, col.axis = INK2, col.lab = INK2)
rg <- seq(0.55, 0.97, by = 0.0025)
plot(NA, xlim = c(0.55, 0.97), ylim = c(-1.0, 1.0), xlab = "Reliability (correlation with true value) of the manual method, rho_manual",
     ylab = "Implied reduction in measurement error variance", axes = FALSE)
rect(0.55, 0.10, 0.97, 0.50, col = "#e8f0fb", border = NA)
# uncertainty band of T2 from the CI of pooled r (only where defined)
lo_c <- er_fun(s$lb, rg, "T2"); hi_c <- er_fun(s$ub, rg, "T2")
hi_c[is.na(hi_c) & rg >= s$ub] <- NA
ok <- !is.na(lo_c) & !is.na(hi_c)
polygon(c(rg[ok], rev(rg[ok])), c(hi_c[ok], rev(lo_c[ok])), col = "#d8e6f8", border = NA)
rect(0.55, -1, r_pool, 1, col = "#efeee9", border = NA)
text(0.555, 0.93, "region infeasible under T2:\nrho_manual must be >= r (= rho_tool x rho_manual)", adj = c(0, 1), cex = 0.68, col = INK2)
rect(0.85, -1, 0.95, 1, col = NA, border = GRID, lty = 3)
text(0.90, 1.07, "Bock 2021: trained raters, Lin's rho_c 0.85-0.95", cex = 0.68, col = INK2, adj = c(0.5, 0), xpd = NA)
text(0.552, 0.115, "model's current error_reduction range (10-50%)", adj = c(0, 0), cex = 0.68, col = INK2)
lines(rg, er_fun(r_pool, rg, "T2"), col = BLUE, lwd = 2.4)
lines(rg, er_fun(r_pool, rg, "T1"), col = ORANGE, lwd = 2.0, lty = 2)
abline(h = 0, col = INK2, lwd = 1)
axis(1, at = seq(0.55, 0.95, 0.05), col = GRID, col.ticks = GRID, cex.axis = 0.85)
axis(2, at = seq(-1, 1, 0.25), labels = paste0(seq(-100, 100, 25), "%"), col = GRID, col.ticks = GRID, cex.axis = 0.85, las = 1)
abline(v = sqrt(r_pool), col = INK2, lty = 3)
text(sqrt(r_pool) + 0.004, 0.62, sprintf("break-even\nrho_manual = %.2f", sqrt(r_pool)), cex = 0.7, col = INK2, adj = c(0, 0.5))
legend("bottomleft", bty = "n", cex = 0.75, text.col = INK2, inset = c(0.0, 0.0),
       legend = c(sprintf("T2: manual reference has error (pooled r = %.2f; band = 95%% CI of r)", r_pool), "T1: reference error-free (tool error share = 1 - r^2)"),
       col = c(BLUE, ORANGE), lty = c(1, 2), lwd = c(2.4, 2))
mtext(sprintf("Pooled agreement (r = %.2f) does not identify error reduction; it depends on how good manual scoring is", r_pool), side = 3, line = 3.2, adj = 0, font = 2, cex = 1)
mtext(sprintf("Pooled manual reliability (3 studies) = %.2f [%.2f, %.2f]; T2 median ER = %.0f%% [%.0f%%, %.0f%%] (feasible draws only).",
              s_rel$est, s_rel$lb, s_rel$ub, 100 * mc_T2["med"], 100 * mc_T2["lo.2.5%"], 100 * mc_T2["hi.97.5%"]),
      side = 3, line = 2.0, adj = 0, cex = 0.72, col = INK2)
dev.off()

# ---- MA2 forest plot ---------------------------------------------------------------
d2o <- d2[order(d2$legume, d2$trait_class, d2$study, d2$row), ]
mk2 <- function(d) {
  lo <- tanh(d$yi - 1.96 * sqrt(d$vi)); hi <- tanh(d$yi + 1.96 * sqrt(d$vi))
  w <- 1 / d$vi; sz <- 0.6 + 1.5 * (w / max(w))
  fa <- paste0(toupper(substr(d$crossref_first_author, 1, 1)), substring(d$crossref_first_author, 2))
  lab <- paste0(fa, " ", d$crossref_year, ifelse(d$preprint, " (preprint)", ""), " - ", d$trait_class)
  data.frame(label = lab, est = tanh(d$yi), lb = lo, ub = hi, type = "study",
             annot = sprintf("%.2f [%.2f, %.2f]  n=%d", tanh(d$yi), lo, hi, d$n), size = sz, plb = NA, pub = NA)
}
rows2 <- rbind(
  data.frame(label = "Legumes", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  mk2(subset(d2o, legume == "legume")),
  data.frame(label = "Non-legumes (wheat, maize)", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  mk2(subset(d2o, legume == "non-legume")),
  data.frame(label = "Pooled, three-level RE (all)", est = tanh(MU2), lb = s$lb, ub = s$ub, type = "pooled",
             annot = sprintf("%.2f [%.2f, %.2f]", s$est, s$lb, s$ub), size = NA, plb = s$plb, pub = s$pub),
  data.frame(label = "", est = NA, lb = NA, ub = NA, type = "note",
             annot = sprintf("PI [%.2f, %.2f]; I2=%.0f%%", s$plb, s$pub, s$I2), size = NA, plb = NA, pub = NA))
forest_plot(file.path(FIG, "fig_ma2_forest_accuracy.png"), rows2, xlim = c(0, 1), xlab = "Correlation of image/CV tool with manual reference (r)",
            title = sprintf("Tool-vs-manual agreement (k = %d effects, %d studies)", nrow(d2), length(unique(d2$study))),
            subtitle = c("Pearson r; Spearman rho and sqrt(R2) treated as r. Preprints labelled.",
                         "n = plots / plants / images, not independent genotypes."),
            ref = NULL, xticks = seq(0, 1, 0.2), width = 2400, label_cex = 0.78)
# sample-size caveat
# =============================================================================
# MA3  Time / cost savings of digital phenotyping vs manual
# =============================================================================
cst <- read.csv(file.path(LIT, "extract/meta_pheno_cost/extraction.csv"))
sa  <- read.csv(file.path(LIT, "extract/meta_pheno_cost/studies_assessed.csv"))
cst$crossref_verified <- as.logical(cst$crossref_verified)
cst$included_study <- as.logical(sa$included[match(cst$doi, sa$doi)])
cst$fulltext_read  <- as.logical(sa$fulltext_read[match(cst$doi, sa$doi)])
cst$live_crossref_ok <- cr_ok(cst$doi)
cst$ok <- cst$included_study & cst$fulltext_read & cst$crossref_verified & cst$live_crossref_ok
chk <- merge(unique(cst[, c("doi", "first_author_crossref", "year_crossref")]), cr[, c("doi", "first_author_family", "issued_year")], by = "doi")
stopifnot(all(tolower(chk$first_author_crossref) == tolower(chk$first_author_family)), all(chk$year_crossref == chk$issued_year))
num <- function(x) suppressWarnings(as.numeric(x))
cst$base <- num(cst$baseline_value); cst$new <- num(cst$new_value)
cst$lr <- ifelse(!is.na(cst$base) & !is.na(cst$new) & cst$base > 0 & cst$new > 0, log(cst$new / cst$base), NA)
# effect units = paper x digital tool (rows of the same unit are different outcomes of the same comparison)
unit_map <- data.frame(
  effect_id = c("E01", "E04", "E04b", "E05", "E06", "E07", "E08", "E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", "E20"),
  unit = c("Lasdun2024_phone", "Walter2019_handheld", "Walter2019_handheld_total", "Walter2019_boom", "Tanger2017_tractor", "Tanger2017_tractor",
           "DeBruin2025_robot", "DeBruin2025_robot", "DeBruin2025_robot", "Mwanje2026_ML", "Mwanje2026_ImageJ", "Komyshev2017_SeedCounter",
           "Komyshev2017_count", "Zu2025_DL", "Zu2025_IP", "Zu2025_DL", "Haghighattalab2016_UAS", "Tattaris2016_UAV", "Reynolds2019_model"),
  paper = c("Lasdun2024", "Walter2019", "Walter2019", "Walter2019", "Tanger2017", "Tanger2017", "DeBruin2025", "DeBruin2025", "DeBruin2025",
            "Mwanje2026", "Mwanje2026", "Komyshev2017", "Komyshev2017", "Zu2025", "Zu2025", "Zu2025", "Haghighattalab2016", "Tattaris2016", "Reynolds2019"),
  setting = c("field plot", "field plot", "field plot", "field plot", "field plot", "field plot", "field plot", "field plot", "field plot",
              "lab/sample", "lab/sample", "lab/sample", "lab/sample", "lab/sample", "lab/sample", "lab/sample", "field plot", "field plot", "field plot"),
  platform = c("smartphone/handheld", "smartphone/handheld", "smartphone/handheld", "vehicle/UAV/robot", "vehicle/UAV/robot", "vehicle/UAV/robot",
               "vehicle/UAV/robot", "vehicle/UAV/robot", "vehicle/UAV/robot", "smartphone/handheld", "smartphone/handheld", "smartphone/handheld",
               "smartphone/handheld", "smartphone/handheld", "smartphone/handheld", "smartphone/handheld", "vehicle/UAV/robot", "vehicle/UAV/robot", "vehicle/UAV/robot"))
cst <- merge(cst, unit_map, by = "effect_id", all.x = TRUE)
add_used("MA3", cst, c("effect_id", "unit", "doi", "first_author_crossref", "year_crossref", "outcome", "baseline_value", "baseline_unit",
                        "new_value", "new_unit", "pct_reduction", "evidence_class", "included_study", "fulltext_read", "crossref_verified",
                        "live_crossref_ok", "location"))
unit_level <- function(rows) {
  rows <- rows[!is.na(rows$lr) & rows$ok, ]
  do.call(rbind, lapply(split(rows, rows$unit), function(x)
    data.frame(unit = x$unit[1], paper = x$paper[1], setting = x$setting[1], platform = x$platform[1], yi = mean(x$lr), nrows = nrow(x),
               pct_red = 100 * (1 - exp(mean(x$lr))), classes = paste(sort(unique(x$evidence_class)), collapse = "/"))))
}
ma3_fit <- function(u) rma(u$yi, vi = rep(1e-8, nrow(u)), method = "REML", test = "knha")
ma3_row <- function(id, label, u, role = "sensitivity", note = "") {
  if (nrow(u) < 3) {
    desc_row(id, label, "time savings", u$pct_red, "% time reduction", paste0("k<3 units. ", note), length(unique(u$paper)))
    return(invisible(NULL))
  }
  r <- ma3_fit(u)
  sh <- summ(r, function(x) 100 * (1 - exp(x)))
  add_row(id, label, "time savings", role, "unweighted RE (REML) on log(new/baseline), nominal v=1e-8 (no variances reported)", nrow(u), length(unique(u$paper)),
          sh$est, min(sh$lb, sh$ub), max(sh$lb, sh$ub), min(sh$plb, sh$pub), max(sh$plb, sh$pub),
          r$tau2, r$I2, r$QE, r$QEp, "% reduction in time (1 - new/baseline)", "Knapp-Hartung t on log scale, back-transformed",
          paste0("Unweighted: no study reports a variance for a time ratio, so the REML weights are equal and the CI reflects between-unit spread only. ", note))
  invisible(r)
}
prim <- subset(cst, evidence_class == "A")
u_prim <- unit_level(prim)
cat("\nMA3 primary units:\n"); print(u_prim[, c("unit", "setting", "platform", "pct_red", "nrows")])
r3 <- ma3_row("MA3.time.primary", "Time reduction, digital vs manual: class A (author-measured/reported both rates), all settings", u_prim, "primary",
              "Units = paper x tool; E14 (Komyshev seed counting, no numeric manual baseline, 'no saving') excluded here and added back as ratio=1 in a sensitivity analysis.")
s3 <- summ(r3, function(x) 100 * (1 - exp(x)))
cat("MA3 primary: ", sprintf("%.1f%% [%.1f, %.1f] PI [%.1f, %.1f]\n", s3$est, min(s3$lb, s3$ub), max(s3$lb, s3$ub), min(s3$plb, s3$pub), max(s3$plb, s3$pub)))
# descriptive
desc_row("MA3.time.primary.median", "Time reduction: class A units, median (range)", "time savings", u_prim$pct_red, "% time reduction",
         "Median of unit-level percentage reductions; gives the robust central value because the log-mean is dominated by lab seed-counting ratios of 0.01-0.07.", length(unique(u_prim$paper)))
# paper level
pl <- do.call(rbind, lapply(split(u_prim, u_prim$paper), function(x) data.frame(unit = x$paper[1], paper = x$paper[1], setting = x$setting[1],
                                                                                 platform = x$platform[1], yi = mean(x$yi), nrows = sum(x$nrows),
                                                                                 pct_red = 100 * (1 - exp(mean(x$yi))), classes = "A")))
ma3_row("MA3.time.paperlevel", "Time reduction: class A, one estimate per paper", pl, note = "k = number of papers.")
ma3_row("MA3.time.field", "Time reduction: class A, field-plot data capture only", subset(u_prim, setting == "field plot"))
ma3_row("MA3.time.lab", "Time reduction: class A, lab/sample-level image analysis (seed counting, scoring images)", subset(u_prim, setting == "lab/sample"))
ma3_row("MA3.time.handheld", "Time reduction: class A, smartphone/handheld platforms (field + lab)", subset(u_prim, platform == "smartphone/handheld"))
ma3_row("MA3.time.vehicle", "Time reduction: class A, vehicle/UAV/robot platforms", subset(u_prim, platform == "vehicle/UAV/robot"))
ma3_row("MA3.time.fieldhandheld", "Time reduction: class A, field-plot AND smartphone/handheld (closest to the modeled tool)", subset(u_prim, setting == "field plot" & platform == "smartphone/handheld"),
        note = "Only Walter 2019 hand-held camera (class A). With the Lasdun 2024 row (class C) see MA3.time.fieldhandheld_withC.")
# sensitivity: add C/B classes and the null E14
allc <- subset(cst, evidence_class %in% c("A", "B", "C") & !is.na(unit) & effect_id != "E03")
allc <- subset(allc, !(effect_id %in% c("E04b")))      # keep E04 (capture only) as the unit's value; E04b analysed below
u_all <- unit_level(allc)
ma3_row("MA3.time.ABC", "Time reduction: classes A+B+C (adds Lasdun E01 derived, Tattaris assumption, Reynolds model)", u_all,
        note = "B/C rows rest on derived or assumed baselines.")
u_fh <- unit_level(subset(cst, effect_id %in% c("E04", "E01")))
ma3_row("MA3.time.fieldhandheld_withC", "Time reduction: field-plot smartphone/handheld incl. Lasdun E01 (class C)", u_fh,
        note = "Walter 2019 hand-held (capture only) and Lasdun 2024 (derived baseline). k = 2 -> descriptive.")
u_tot <- u_prim; u_tot$yi[u_tot$unit == "Walter2019_handheld"] <- log(180.3 / 293.1)
u_tot$pct_red[u_tot$unit == "Walter2019_handheld"] <- 100 * (1 - 180.3 / 293.1)
ma3_row("MA3.time.walter_total", "Time reduction: class A with Walter hand-held replaced by total time incl. image analysis (E04b)", u_tot,
        note = "Capture 55.6% falls to 38.5% once ~10 min operator image analysis per trial is counted.")
u_null <- rbind(u_prim, data.frame(unit = "Komyshev2017_count", paper = "Komyshev2017", setting = "lab/sample", platform = "smartphone/handheld",
                                   yi = 0, nrows = 1, pct_red = 0, classes = "A (ratio set to 1)"))
ma3_row("MA3.time.plusnull", "Time reduction: class A + Komyshev seed-counting null result (ratio set to 1)", u_null,
        note = "E14: authors state manual counting of ~50 grains 'may be a little less' time than the app; ratio assumed 1.0 (no saving).")
# leave-one-paper-out
for (pp in unique(u_prim$paper)) {
  x <- subset(u_prim, paper != pp); r <- ma3_fit(x); p <- summ(r, function(v) 100 * (1 - exp(v)))
  add_loo("MA3.time.primary", pp, p$est, min(p$lb, p$ub), max(p$lb, p$ub))
}
cat(sprintf("MA3: %d units from %d papers (k<10 papers: small-study tests not applicable; no variances available in any case)\n",
            nrow(u_prim), length(unique(u_prim$paper))))
# ---- cost (not time) outcomes: descriptive -------------------------------------------
cost_rows <- subset(cst, effect_id %in% c("E21", "E22", "E23") & ok)
desc_row("MA3.cost.descriptive", "Cost reduction (imaging cost per plot / total experiment cost), Reynolds 2019 cost model: E21 (-70%), E22 (0%), E23 (-7%)",
         "cost", num(cost_rows$pct_reduction), "% cost reduction",
         "All three rows come from ONE modelled study (Reynolds 2019, high-income infrastructure costs; UAV vs ground vehicle / hand-held), none is manual vs digital in an African breeding programme. k_papers = 1: NOT poolable.", 1)
desc_row("MA3.cost.share_nonhandling", "Share of total experiment cost that is NOT plot handling (E24, Reynolds 2019)", "cost",
         c(25.0, 30.0, 36.3), "% of total cost", "Reported range 25-36.3% (central ~30%), = upper bound on any saving from digital data capture alone in that cost structure.", 1)
desc_row("MA3.cost.plot_benchmark", "Baseline cost per plot: Juliana 2018 assumption (E26, USD 10) vs Reynolds 2019 plot handling (E25, USD 30-50 high-income)", "cost",
         c(10, 30, 40, 50), "USD per plot", "Not an effect size. USD 10 is an unreferenced wheat-nursery assumption; Das 2025 gives East African maize per-row minimum costs (USD 5.73 evaluation, 9.35 germplasm).", 2)
# implied total plot-cost reduction = (share of plot cost that is data capture/analysis) x (time reduction), bounding calculation
share <- c(0.25, 0.30, 0.363)                    # Reynolds 2019 (E24): non-plot-handling share of total experiment cost
tr_hh <- c(0.556, 0.792)                          # Walter 2019 handheld (A) and Lasdun 2024 (C, derived)
grid_cost <- outer(share, tr_hh)
add_row("MA3.cost.implied_handheld", "Implied reduction in TOTAL cost per plot = data-capture share x time reduction (smartphone/handheld, field)", "cost", "descriptive",
        "bounding arithmetic (not a pooled estimate)", length(grid_cost), NA, 0.30 * mean(tr_hh), min(grid_cost), max(grid_cost), NA, NA, NA, NA, NA, NA,
        "fraction of total plot cost", "min/max over share {0.25,0.30,0.363} x time reduction {0.556,0.792}",
        "Assumes plot handling (land, planting, harvest; 65-77% of cost in Reynolds 2019) is unchanged and the whole non-handling share falls in proportion to time. Reynolds' share is for high-income field infrastructure and an imaging baseline: East African bean plots are unmeasured. With the pooled all-platform time reduction (94.5%) the range is 24-34%.")
# ---- MA3 forest plot -------------------------------------------------------------------
u_plot <- u_prim[order(u_prim$setting, u_prim$platform, u_prim$pct_red), ]
ratio <- exp(u_plot$yi)
rows3 <- rbind(
  data.frame(label = "Field plot data capture", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  data.frame(label = gsub("_", " - ", u_plot$unit[u_plot$setting == "field plot"]), est = ratio[u_plot$setting == "field plot"], lb = NA, ub = NA, type = "study",
             annot = sprintf("%.0f%% less time", u_plot$pct_red[u_plot$setting == "field plot"]), size = 1.3, plb = NA, pub = NA),
  data.frame(label = "Lab / sample-level image analysis", est = NA, lb = NA, ub = NA, type = "header", annot = NA, size = NA, plb = NA, pub = NA),
  data.frame(label = gsub("_", " - ", u_plot$unit[u_plot$setting == "lab/sample"]), est = ratio[u_plot$setting == "lab/sample"], lb = NA, ub = NA, type = "study",
             annot = sprintf("%.0f%% less time", u_plot$pct_red[u_plot$setting == "lab/sample"]), size = 1.3, plb = NA, pub = NA),
  data.frame(label = "Pooled (unweighted RE, log ratio)", est = exp(r3$beta[1]), lb = exp(predict(r3)$ci.lb), ub = exp(predict(r3)$ci.ub), type = "pooled",
             annot = sprintf("%.0f%% [%.0f, %.0f]", s3$est, min(s3$lb, s3$ub), max(s3$lb, s3$ub)), size = NA,
             plb = exp(predict(r3)$pi.lb), pub = exp(predict(r3)$pi.ub)),
  data.frame(label = "", est = NA, lb = NA, ub = NA, type = "note", annot = sprintf("PI %.0f to %.0f%%; median %.0f%%", min(s3$plb, s3$pub), max(s3$plb, s3$pub), median(u_prim$pct_red)),
             size = NA, plb = NA, pub = NA))
forest_plot(file.path(FIG, "fig_ma3_forest_time_savings.png"), rows3, xlim = c(0.003, 2), logx = TRUE, ref = 1,
            xlab = "Time needed by digital method relative to manual (log scale; 1 = no saving)",
            title = sprintf("Time saving of digital vs manual phenotyping (%d units, %d papers)", nrow(u_prim), length(unique(u_prim$paper))),
            subtitle = c("Class A rows only (both rates reported by the authors); capture/processing time, not total programme cost.",
                         "No study reports a variance for a time ratio, so squares carry no CI and pooling is unweighted."),
            xticks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), xtick_labels = c("0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1"),
            secondary_axis = "ratio 0.1 = 90% less time; 0.5 = 50% less", label_cex = 0.82)

# =============================================================================
# write results
# =============================================================================
out <- do.call(rbind, RES)
# Q / I2 are meaningless when no sampling variances exist (nominal v = 1e-8): blank them
nom <- grepl("nominal v=1e-8", out$model)
out$Q[nom] <- NA; out$Q_p[nom] <- NA; out$I2_pct[nom] <- NA
num_cols <- c("estimate", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "I2_pct", "Q", "Q_p")
out[num_cols] <- lapply(out[num_cols], function(x) round(as.numeric(x), 4))
write.csv(out, file.path(LIT, "meta_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, LOO), file.path(LIT, "meta_leave_one_out.csv"), row.names = FALSE)
# harmonise columns (union) before writing
all_cols <- unique(unlist(lapply(USED, names)))
used <- do.call(rbind, lapply(USED, function(d) { for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA; d[, all_cols, drop = FALSE] }))
write.csv(used, file.path(LIT, "meta_inputs_used.csv"), row.names = FALSE)
cat("\nWrote", nrow(out), "result rows.\n")
print(out[, c("analysis_id", "role", "k_effects", "k_clusters", "estimate", "ci_lb", "ci_ub", "pi_lb", "pi_ub", "tau2", "I2_pct")], row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(LIT, "extract/synthesis/sessionInfo.txt"))
