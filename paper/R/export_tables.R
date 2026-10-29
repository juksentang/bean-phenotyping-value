# ============================================================================
# Export Table 2 (meta-analysis results used to set model inputs) -> paper/tables/table2_meta.tex
# ============================================================================
# Every value comes from paper/lit/meta_results.csv; the study lists in the table
# note come from meta_inputs_used.csv, extraction.csv and the DOI-to-BibTeX map.
# Run via `make tables` from the project root.

suppressMessages(library(data.table))

LIT <- "paper/lit"
mr  <- fread(file.path(LIT, "meta_results.csv"))
bk  <- fread(file.path(LIT, "extract/synthesis/bibkeys.tsv"))
bib <- readLines(file.path(LIT, "refs.bib"), warn = FALSE)
bib_keys <- sub("^@[A-Za-z]+\\{([^,]+),.*$", "\\1", grep("^@", bib, value = TRUE))
doi2key  <- setNames(bk$key, tolower(bk$doi))

minus <- function(x) sub("^-", "$-$", x)
f <- function(x, d = 2, scale = 1) minus(formatC(x * scale, format = "f", digits = d))
span_txt <- function(lo, hi, d, scale) paste(f(lo, d, scale), "to", f(hi, d, scale))

# id, label, digits, scale (1 or 100), role in the model, descriptive (median and range) flag
spec <- data.table(
  id = c("MA1.CB.primary", "MA1.CB.clustercollapsed", "MA1.CB.embrapaB08", "MA1.CB.exB02", "MA1.CB.imputedA", "MA1.OC.primary",
         "MA2.primary", "MA2.nopreprint", "MA2.trait.height", "MA2.trait.plant_count", "MA2.trait.disease_severity",
         "MA2.manualrel", "MA2.manualrel.validity", "MA2.ER.T1", "MA2.ER.T2", "MA2.ER.T2v", "MA2.ER.h2.median",
         "MA3.time.primary", "MA3.time.fieldhandheld_withC", "MA3.time.plusnull",
         "MA3.cost.share_nonhandling", "MA3.cost.implied_handheld"),
  label = c("Realized gain, common bean, three-level model (\\% of mean yield per year)",
            "One estimate per program",
            "Embrapa mixed-model estimate in place of the direct estimate",
            "Without IAC 1997 to 2007",
            "With one added estimate whose standard error is imputed (rule A)",
            "Realized gain, maize and rice in Africa (context)",
            "Tool-manual correlation, all traits (\\emph{r})",
            "Peer-reviewed studies only",
            "Plant height",
            "Plant count",
            "Disease severity",
            "Manual-scoring reliability as entered (\\emph{r})",
            "Manual-scoring validity after conversion (\\emph{r})",
            "Implied error reduction, reference error free, T1 (\\%)",
            "Implied error reduction, T2 with the reliability as entered (\\%)",
            "Implied error reduction, T2 with the converted validity (\\%)",
            "Implied error reduction, heritability pairs (\\%)",
            "Time reduction, digital versus manual (\\%)",
            "Time reduction, handheld tools in field plots, capture time (\\%)",
            "Time reduction with the seed-counting null result added (\\%)",
            "Share of plot cost that is not plot handling (\\%)",
            "Implied whole-plot cost reduction, handheld tools (\\%)"),
  d     = c(2, 2, 2, 2, 2, 2,  2, 2, 2, 2, 2,  2, 2, 0, 0, 0, 0,  0, 0, 0,  0, 0),
  scale = c(1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 100, 100, 100, 100,  1, 1, 1,  1, 100),
  role  = c("Calibration target", "Sensitivity", "Sensitivity", "Sensitivity", "Sensitivity", "Context only",
            "Reported", "Sensitivity", "Subgroup", "Subgroup", "Subgroup",
            "Used in T2 as entered", "Used in T2 converted", "Evidence route", "Evidence route", "Evidence route", "Central value (rounded)",
            "Reported", "Cost reduction", "Sensitivity",
            "Cost reduction", "Central cost reduction"),
  desc  = c(rep(FALSE, 13), TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE))
stopifnot(all(spec$id %in% mr$analysis_id))

fmt_row <- function(i) {
  s <- spec[i]; r <- mr[analysis_id == s$id]
  k <- if (is.na(r$k_clusters) || r$k_clusters == r$k_effects) as.character(r$k_effects)
       else sprintf("%d (%d)", r$k_effects, r$k_clusters)
  if (is.na(r$k_effects)) k <- "--"
  est <- f(r$estimate, s$d, s$scale)
  ci  <- if (is.na(r$ci_lb)) "--" else span_txt(r$ci_lb, r$ci_ub, s$d, s$scale)
  if (s$desc && !is.na(r$ci_lb)) ci <- paste0(ci, "\\textsuperscript{b}")
  pi  <- if (is.na(r$pi_lb)) "--" else span_txt(r$pi_lb, r$pi_ub, s$d, s$scale)
  i2  <- if (is.na(r$I2_pct)) "--" else formatC(r$I2_pct, format = "f", digits = 0)
  paste(s$label, k, est, ci, pi, i2, s$role, sep = " & ")
}
rows <- vapply(seq_len(nrow(spec)), fmt_row, character(1))
blocks <- list(c(1, 6), c(7, 17), c(18, 22))
head_rows <- c("\\multicolumn{7}{l}{\\textit{MA1: realized genetic gain}}\\\\",
               "\\multicolumn{7}{l}{\\textit{MA2: tool-manual agreement and implied error reduction}}\\\\",
               "\\multicolumn{7}{l}{\\textit{MA3: time and cost savings}}\\\\")

# --- study lists for the table note (included studies, by BibTeX key) ---
u  <- fread(file.path(LIT, "meta_inputs_used.csv"))
keyof <- function(dois) {
  k <- unname(doi2key[tolower(unique(dois))])
  stopifnot(!anyNA(k), all(k %in% bib_keys))
  paste(sort(k), collapse = ",")
}
ma1_cb <- u[analysis_family == "MA1" & included == TRUE & primary == "Y" & stratum3 == "common bean"]
ma1_oc <- u[analysis_family == "MA1" & included == TRUE & primary == "Y" & stratum3 %in% c("other crops", "other grain legumes")]
ma2    <- u[analysis_family == "MA2" & included == TRUE]
ex <- fread(file.path(LIT, "extract/meta_pheno_cost/extraction.csv"))
sa <- fread(file.path(LIT, "extract/meta_pheno_cost/studies_assessed.csv"))
ma3 <- ex[evidence_class == "A" & doi %in% sa[included == TRUE]$doi]

note <- paste0(
  "\\par\\smallskip\\noindent{\\footnotesize{Estimates are pooled by restricted maximum likelihood; the interval for the Monte Carlo and descriptive rows is the 2.5 to 97.5 percentile range or the minimum to maximum, and the time reductions are unweighted because no study reports a variance. ",
  "$k$: number of effects (number of programs or studies); CI: confidence interval; PI: prediction interval; $I^2$: heterogeneity. ",
  "\\textsuperscript{a} Included studies. MA1, common bean: \\cite{", keyof(ma1_cb$doi), "}; MA1, other crops: \\cite{", keyof(ma1_oc$doi),
  "}; MA2: \\cite{", keyof(ma2$doi), "}; MA3: \\cite{", keyof(ma3$doi), "}. ",
  "\\textsuperscript{b} Descriptive: median with minimum to maximum, not pooled.}}")

out <- c(
  "\\begin{table}[!htbp]",
  "\\caption{Meta-analysis results used to set model inputs.\\textsuperscript{a}\\label{tab:meta}}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabularx}{\\linewidth}{>{\\raggedright\\arraybackslash}X>{\\centering\\arraybackslash}p{1.2cm}>{\\centering\\arraybackslash}p{1.4cm}>{\\centering\\arraybackslash}p{2.4cm}>{\\centering\\arraybackslash}p{2.4cm}>{\\centering\\arraybackslash}p{1.0cm}>{\\raggedright\\arraybackslash}p{3.6cm}}",
  "\\toprule",
  "\\textbf{Quantity} & \\textbf{\\emph{k}} & \\textbf{Estimate} & \\textbf{95\\% CI} & \\textbf{95\\% PI} & \\textbf{\\emph{I}\\textsuperscript{2} (\\%)} & \\textbf{Use in model}\\\\",
  "\\midrule")
for (b in seq_along(blocks)) {
  out <- c(out, head_rows[b], paste0(rows[blocks[[b]][1]:blocks[[b]][2]], "\\\\"))
  if (b < length(blocks)) out <- c(out, "\\midrule")
}
out <- c(out, "\\bottomrule", "\\end{tabularx}", note, "\\end{table}")
dir.create("paper/tables", showWarnings = FALSE)
writeLines(out, "paper/tables/table2_meta.tex")
cat("Wrote paper/tables/table2_meta.tex\n")


# ============================================================================
# Table A1: appraisal of the MA2 effects (reference type, platform, unit of observation, preprint status)
# ============================================================================
ax <- fread(file.path(LIT, "extract/meta_pheno_accuracy/extraction.csv"))
ax <- ax[analysis_set == "primary" & verified_crossref == TRUE & fulltext_read == TRUE & !is.na(n)]
pre <- tolower(ax$doi) %in% tolower(cr_pre <- {
  cr <- fread(file.path(LIT, "extract/synthesis/crossref_verification.tsv"), quote = "")
  cr[type == "posted-content"]$doi })
ref_type <- function(x) {
  x <- tolower(x)
  ifelse(grepl("visual|graded", x), "visual score",
  ifelse(grepl("count", x), "manual count",
  ifelse(grepl("ruler|caliper|ground-measured|manual plant height|height", x), "instrument measurement", "other")))
}
metric_tx <- function(m) ifelse(grepl("R2", m), "$\\sqrt{R^2}$", ifelse(grepl("Spearman", m), "Spearman $\\rho$", "Pearson"))
tex_esc <- function(x) gsub("%", "\\%", gsub("_", "\\_", x, fixed = TRUE), fixed = TRUE)
tool_short <- function(x) {
  x <- tolower(x)
  ifelse(grepl("smartphone", x), "smartphone",
  ifelse(grepl("uav", x), "UAV camera", "camera or scanner"))
}
rowsA <- sprintf("\\cite{%s} & %s & %s & %s & %s & %s & %.2f (%s) & %d & %s\\\\",
                 vapply(ax$doi, function(d) unname(doi2key[tolower(d)]), ""),
                 ax$crop, ax$trait_class, tool_short(ax$tool_class), ref_type(ax$manual_reference), ax$unit_of_observation,
                 ax$r, metric_tx(ax$metric), as.integer(ax$n), ifelse(pre, "yes", "no"))
rowsA <- gsub("&amp;", "\\&", rowsA)
outA <- c(
  "\\begin{table}[!htbp]",
  "\\caption{Appraisal of the effects pooled in MA2: crop, trait, platform, type of manual reference, unit of observation (which defines $n$ in the Fisher $z$ variance $1/(n-3)$) and preprint status.\\label{tab:ma2appraisal}}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabularx}{\\linewidth}{>{\\raggedright\\arraybackslash}p{3.6cm}>{\\raggedright\\arraybackslash}p{2.0cm}>{\\raggedright\\arraybackslash}p{2.6cm}>{\\raggedright\\arraybackslash}p{2.6cm}>{\\raggedright\\arraybackslash}p{3.2cm}>{\\raggedright\\arraybackslash}p{1.8cm}>{\\raggedright\\arraybackslash}X>{\\raggedleft\\arraybackslash}p{1.0cm}>{\\raggedright\\arraybackslash}p{1.3cm}}",
  "\\toprule",
  "\\textbf{Study} & \\textbf{Crop} & \\textbf{Trait} & \\textbf{Platform} & \\textbf{Manual reference} & \\textbf{Unit} & \\textbf{Correlation entered} & \\textbf{\\emph{n}} & \\textbf{Preprint}\\\\",
  "\\midrule", rowsA, "\\bottomrule", "\\end{tabularx}",
  "\\par\\smallskip\\noindent{\\footnotesize{Manual-reference types were coded from the description of the reference in each source. The first column gives the study, and the unit column gives what $n$ counts (plots, plants, images, accessions, lines or cultivars), so $n$ is not the number of independent genotypes. Correlations are entered as reported, with square roots of $R^2$ and Spearman coefficients treated as correlations.}}",
  "\\end{table}")
writeLines(outA, "paper/tables/tableA1_ma2_appraisal.tex")
cat("Wrote paper/tables/tableA1_ma2_appraisal.tex\n")


# ============================================================================
# Table: image-based phenotyping tools, including deep learning, in the MA2 evidence
# -> paper/tables/table_ai_tools.tex and paper/tables/table_ai_tools.csv
# ============================================================================
# One row per group of primary MA2 effects (same inclusion rule as MA2.primary).
# Pooled r is read from meta_results.csv where a pooled subgroup row exists; the
# median and range of the effect-level r are computed here for every group.
ai <- fread(file.path(LIT, "extract/meta_pheno_accuracy/extraction.csv"))
ai <- ai[analysis_set == "primary" & verified_crossref == TRUE & fulltext_read == TRUE & !is.na(n) &
           !is.na(r) & abs(r) < 1]
ma2p <- mr[analysis_id == "MA2.primary"]
stopifnot(nrow(ai) == ma2p$k_effects, uniqueN(ai$study_id) == ma2p$k_clusters)

# Short names for the measured traits (keys are the trait strings of the extraction sheet)
trait_short <- c("plant stand count" = "stand count", "plant height" = "plant height",
                 "canopy height" = "canopy height", "pod number per plant" = "pods per plant",
                 "frogeye leaf spot severity (% leaf area)" = "frogeye leaf spot",
                 "seedling emergence count" = "seedling emergence count",
                 "ear density (ears per image)" = "ears per image", "ears per m2" = "ears per m$^2$",
                 "leaf area covered by lesions (STB severity)" = "Septoria tritici blotch",
                 "angular leaf spot severity" = "angular leaf spot", "seed area" = "seed area",
                 "canopy senescence" = "canopy senescence")
stopifnot(all(ai$trait %in% names(trait_short)))
ai[, trait_txt := unname(trait_short[trait])]
trait_lab <- c("plant count" = "Plant count", "height" = "Height", "disease severity" = "Disease severity",
               "seed size" = "Seed size", "seed/pod count" = "Seed or pod count", "senescence" = "Senescence")
stopifnot(all(ai$trait_class %in% names(trait_lab)))
ai[, class_txt := unname(sub("seed/pod count", "seed or pod count", trait_class, fixed = TRUE))]

# Crop names: winter and durum wheat are wheat; dry bean and common bean are Phaseolus vulgaris
ai[, crop_txt := ifelse(grepl("wheat", crop), "wheat", ifelse(grepl("dry bean|common bean", crop), "common bean", crop))]

# Platform from the tool class, using the same patterns as the MA2.tool.* subgroups
plat_of <- function(tc) {
  ifelse(grepl("smartphone", tc, ignore.case = TRUE), "smartphone",
  ifelse(grepl("UAV-RGB", tc, ignore.case = TRUE), "UAV RGB",
  ifelse(grepl("camera", tc, ignore.case = TRUE), "camera",
  ifelse(grepl("UAV-multispectral", tc, ignore.case = TRUE), "UAV multispectral",
  ifelse(grepl("scanner", tc, ignore.case = TRUE), "scanner", NA_character_)))))
}
ai[, platform := plat_of(tool_class)]
stopifnot(!anyNA(ai$platform),
          all(rowSums(sapply(c("smartphone", "UAV-RGB", "camera"), function(p) grepl(p, ai$tool_class, ignore.case = TRUE))) <= 1))
ai[, platform_grp := ifelse(platform %in% c("smartphone", "UAV RGB", "camera"), platform, "other platforms")]

# Deep learning: the tool column names a deep neural network
ai[, dl := grepl("R-CNN|YOLO|DeepLab|DeepCount|Siamese", tool) | grepl("DL|CNN", tool_class)]
ai[, method := ifelse(dl, "deep learning", "other")]
ai[, platform_dl := ifelse(dl, paste(platform, "+ DL"), platform)]

# Recording setting: controlled if the setting names a greenhouse, studio, bench, laboratory or chamber
# (images acquired away from the crop, even when the material was grown in the field); field if it names a
# field and none of those. Any other setting stops the script.
ai[, controlled := grepl("greenhouse|studio|bench|\\blab\\b|laborator|chamber", setting, ignore.case = TRUE)]
ai[, field := !controlled & grepl("field", setting, ignore.case = TRUE)]
stopifnot(all(ai$controlled | ai$field))

# Bibliography keys for the studies
ai[, key := vapply(doi, function(d) unname(doi2key[tolower(d)]), "")]
stopifnot(!anyNA(ai$key), all(ai$key %in% bib_keys))

# Pooled estimate lookup from meta_results.csv (pooled rows only; descriptive rows are not pooled)
pooled_row <- function(id) {
  if (is.na(id)) return(NULL)
  r <- mr[analysis_id == id]
  if (nrow(r) == 0 || r$role == "descriptive" || is.na(r$I2_pct)) return(NULL)
  r
}

make_row <- function(d, panel, group, id, col2, col4) {
  st <- uniqueN(d$study_id); ef <- nrow(d)
  pr <- pooled_row(id)
  if (!is.null(pr)) stopifnot(pr$k_effects == ef, pr$k_clusters == st)
  if (!is.na(id) && nrow(mr[analysis_id == id]) == 1 && mr[analysis_id == id]$role == "descriptive") {
    # a descriptive row of meta_results.csv must agree with the median and range computed here
    dr <- mr[analysis_id == id]
    stopifnot(abs(dr$estimate - median(d$r)) < 1e-3, abs(dr$ci_lb - min(d$r)) < 1e-3, abs(dr$ci_ub - max(d$r)) < 1e-3)
  }
  list(panel = panel, group = group,
       traits = col2(d),
       crops = paste(sort(unique(d$crop_txt)), collapse = ", "),
       detail = col4(d),
       n_studies = st, n_effects = ef,
       n_field_effects = sum(d$field), pct_field_effects = 100 * mean(d$field),
       n_field_studies = uniqueN(d[field == TRUE]$study_id),
       pooled_r = if (is.null(pr)) NA_real_ else pr$estimate,
       pooled_lo = if (is.null(pr)) NA_real_ else pr$ci_lb,
       pooled_hi = if (is.null(pr)) NA_real_ else pr$ci_ub,
       pooled_id = if (is.null(pr)) NA_character_ else id,
       median_r = median(d$r), min_r = min(d$r), max_r = max(d$r),
       keys = paste(sort(unique(d$key)), collapse = ","))
}
count_txt <- function(d, col) {
  tb <- unique(d[, c("study_id", col), with = FALSE])[, .N, by = col][order(-N)]
  paste0(tb[[col]], " (", tb$N, ")", collapse = ", ")
}
traits_txt  <- function(d) paste(unique(d$trait_txt), collapse = ", ")
classes_txt <- function(d) paste(unique(d$class_txt), collapse = ", ")
plat_txt    <- function(d) count_txt(d, "platform_dl")
method_txt  <- function(d) count_txt(d, "method_short")
cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
ai[, method_short := ifelse(dl, "DL", "other")]

res <- list()
for (tc in names(trait_lab)) {
  d <- ai[trait_class == tc]
  res[[length(res) + 1]] <- make_row(d, "A", trait_lab[[tc]], paste0("MA2.trait.", gsub(" ", "_", tc)), traits_txt, plat_txt)
}
plat_id <- c(smartphone = "MA2.tool.smartphone", "UAV RGB" = "MA2.tool.UAV-RGB", camera = "MA2.tool.camera",
             "other platforms" = NA_character_)
for (pg in names(plat_id)) {
  d <- ai[platform_grp == pg]
  lab <- if (pg == "other platforms") paste0("Other (", paste(sort(unique(d$platform)), collapse = ", "), ")") else
         if (pg == "UAV RGB") pg else cap1(pg)
  res[[length(res) + 1]] <- make_row(d, "B", lab, plat_id[[pg]], classes_txt, method_txt)
}
method_lab <- c("deep learning" = "Deep learning (DL)", "other" = "Other image analysis")
for (m in names(method_lab)) {
  res[[length(res) + 1]] <- make_row(ai[method == m], "C", method_lab[[m]], NA_character_, classes_txt, plat_txt)
}
res[[length(res) + 1]] <- make_row(ai, "All", "All primary effects", "MA2.primary",
                                   function(d) sprintf("%d trait classes", uniqueN(d$trait_class)),
                                   function(d) sprintf("%d platforms", uniqueN(d$platform)))
tab <- rbindlist(res)
tab[panel == "All", crops := sprintf("%d crops", uniqueN(ai$crop_txt))]
# every panel partitions the same effects
for (p in c("A", "B", "C")) stopifnot(sum(tab[panel == p]$n_effects) == nrow(ai), sum(tab[panel == p]$n_studies) == uniqueN(ai$study_id))

write.csv(transform(tab,
                    pct_field_effects = round(pct_field_effects, 1),
                    pooled_r = round(pooled_r, 4), pooled_lo = round(pooled_lo, 4), pooled_hi = round(pooled_hi, 4),
                    median_r = round(median_r, 4), min_r = round(min_r, 4), max_r = round(max_r, 4)),
          "paper/tables/table_ai_tools.csv", row.names = FALSE, na = "")
cat("Wrote paper/tables/table_ai_tools.csv\n")

# --- LaTeX ---
cell_pooled <- function(i) if (is.na(tab$pooled_r[i])) "--" else
  paste0(f(tab$pooled_r[i], 2), " (", span_txt(tab$pooled_lo[i], tab$pooled_hi[i], 2, 1), ")")
cell_median <- function(i) if (tab$n_effects[i] == 1) f(tab$median_r[i], 2) else
  paste0(f(tab$median_r[i], 2), " (", span_txt(tab$min_r[i], tab$max_r[i], 2, 1), ")")
# the studies behind each row are listed effect by effect in the MA2 appraisal table
cell_n <- function(i) paste0(tab$n_studies[i], " (", tab$n_effects[i], ")")
row_tex <- function(i) {
  paste0(paste(tab$group[i], tab$traits[i], tab$crops[i], tab$detail[i], cell_n(i),
               sprintf("%d/%d (%.0f)", tab$n_field_effects[i], tab$n_effects[i], tab$pct_field_effects[i]),
               cell_pooled(i), cell_median(i), sep = " & "), "\\\\")
}
panel_head <- c(A = "A. By trait class", B = "B. By platform", C = "C. By analysis method")
body <- character()
for (p in c("A", "B", "C", "All")) {
  body <- c(body, if (p != "A") "\\midrule",
            if (p != "All") sprintf("\\multicolumn{8}{l}{\\textit{%s}}\\\\", panel_head[[p]]),
            vapply(which(tab$panel == p), row_tex, ""))
}
i2_rng <- range(mr[analysis_id %in% tab$pooled_id[!is.na(tab$pooled_id)]]$I2_pct)
par_key <- unique(ai[study_id == "Parrella2013"]$key)
stopifnot(length(par_key) == 1)
par_dis <- ai[trait_class == "disease severity",
              .(n = .N, n_field = sum(field), n_par = sum(field & study_id == "Parrella2013"))]
stopifnot(par_dis$n_par >= 1)
ai_note <- paste0(
  "\\par\\smallskip\\noindent{\\footnotesize{Primary MA2 effects only (",
  ma2p$k_effects, " effects from ", ma2p$k_clusters, " studies); each panel splits the same effects. ",
  "Platforms (panels A and C) and analysis methods (panel B) are counted in studies; the studies are listed effect by effect in Table~\\ref{tab:ma2appraisal}. ",
  "Field: effects recorded in field conditions out of all effects in the group, classified from the setting recorded at extraction (rule in the caption); ",
  "material from field plots is classified as field when the source does not state where it was imaged; this applies to one disease-severity study of leaflets \\cite{", par_key, "} and raises the field share of disease severity from ",
  par_dis$n_field - par_dis$n_par, "/", par_dis$n, " to ", par_dis$n_field, "/", par_dis$n, ". ",
  "Pooled $r$: back-transformed random-effects mean of Fisher $z$ with its 95\\% confidence interval, from the MA2 subgroup models ($I^2$ from ",
  formatC(i2_rng[1], format = "f", digits = 1), " to ", formatC(i2_rng[2], format = "f", digits = 1),
  "\\%); it is shown only where such a subgroup estimate exists, and the median alone is given otherwise. ",
  "Median: median of the effect-level $r$ as entered (Pearson $r$, the square root of $R^2$ or Spearman $\\rho$), with the minimum to maximum in parentheses. ",
  "Plant count includes counts of plants, seedlings and ears. ",
  "DL: deep learning, that is, a deep neural network (R-CNN, YOLO, DeepLab, DeepCount or Siamese network); other: thresholding, image-analysis software, photogrammetry or regression on spectral data.}}")

ai_out <- c(
  "\\begin{table}[!htbp]",
  "\\caption{Image-based phenotyping tools, including deep-learning tools, in the MA2 evidence: traits, crops, platforms, share of effects recorded in the field, and agreement with manual scoring. The correlation $r$ is the agreement between the tool and the manual score of the same plants, plots or samples; it is not the accuracy of the tool against the true value, because manual scores contain error. A setting is field when the source records a field trial or field plots, and controlled when the images were acquired in a greenhouse, a studio, on a bench or in a laboratory, including material taken from field trials.\\label{tab:aitools}}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabularx}{\\linewidth}{>{\\raggedright\\arraybackslash}p{2.0cm}>{\\raggedright\\arraybackslash}X>{\\raggedright\\arraybackslash}p{2.4cm}>{\\raggedright\\arraybackslash}X>{\\raggedright\\arraybackslash}p{1.5cm}>{\\centering\\arraybackslash}p{1.5cm}>{\\raggedright\\arraybackslash}p{3.3cm}>{\\raggedright\\arraybackslash}p{3.0cm}}",
  "\\toprule",
  "\\textbf{Group} & \\textbf{Traits} & \\textbf{Crops} & \\textbf{Platforms or methods} & \\textbf{Studies (effects)} & \\textbf{Field, \\emph{n}/\\emph{N} (\\%)} & \\textbf{Pooled \\emph{r} (95\\% CI)} & \\textbf{Median \\emph{r} (range)}\\\\",
  "\\midrule",
  paste0(body, collapse = "\n"),
  "\\bottomrule", "\\end{tabularx}", ai_note, "\\end{table}")
writeLines(ai_out, "paper/tables/table_ai_tools.tex")
cat("Wrote paper/tables/table_ai_tools.tex\n")
