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
  "\\noindent{\\footnotesize{Estimates are pooled by restricted maximum likelihood; the interval for the Monte Carlo and descriptive rows is the 2.5 to 97.5 percentile range or the minimum to maximum, and the time reductions are unweighted because no study reports a variance. ",
  "$k$: number of effects (number of programs or studies); CI: confidence interval; PI: prediction interval; $I^2$: heterogeneity. ",
  "\\textsuperscript{a} Included studies. MA1, common bean: \\cite{", keyof(ma1_cb$doi), "}; MA1, other crops: \\cite{", keyof(ma1_oc$doi),
  "}; MA2: \\cite{", keyof(ma2$doi), "}; MA3: \\cite{", keyof(ma3$doi), "}. ",
  "\\textsuperscript{b} Descriptive: median with minimum to maximum, not pooled.}}")

out <- c(
  "\\begin{table}[H]",
  "\\caption{Meta-analysis results used to set model inputs.\\textsuperscript{a}\\label{tab:meta}}",
  "\\begin{adjustwidth}{-\\extralength}{0cm}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabularx}{\\fulllength}{>{\\raggedright\\arraybackslash}X>{\\centering\\arraybackslash}p{1.2cm}>{\\centering\\arraybackslash}p{1.3cm}>{\\centering\\arraybackslash}p{2.5cm}>{\\centering\\arraybackslash}p{2.5cm}>{\\centering\\arraybackslash}p{0.9cm}>{\\raggedright\\arraybackslash}p{2.7cm}}",
  "\\toprule",
  "\\textbf{Quantity} & \\textbf{\\emph{k}} & \\textbf{Estimate} & \\textbf{95\\% CI} & \\textbf{95\\% PI} & \\textbf{\\emph{I}\\textsuperscript{2} (\\%)} & \\textbf{Use in model}\\\\",
  "\\midrule")
for (b in seq_along(blocks)) {
  out <- c(out, head_rows[b], paste0(rows[blocks[[b]][1]:blocks[[b]][2]], "\\\\"))
  if (b < length(blocks)) out <- c(out, "\\midrule")
}
out <- c(out, "\\bottomrule", "\\end{tabularx}", "\\end{adjustwidth}", note, "\\end{table}")
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
  "\\begin{table}[H]",
  "\\caption{Appraisal of the effects pooled in MA2: crop, trait, platform, type of manual reference, unit of observation (which defines $n$ in the Fisher $z$ variance $1/(n-3)$) and preprint status.\\label{tab:ma2appraisal}}",
  "\\begin{adjustwidth}{-\\extralength}{0cm}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabularx}{\\fulllength}{>{\\raggedright\\arraybackslash}p{1.8cm}>{\\raggedright\\arraybackslash}p{1.8cm}>{\\raggedright\\arraybackslash}p{2.2cm}>{\\raggedright\\arraybackslash}p{2.2cm}>{\\raggedright\\arraybackslash}p{2.6cm}>{\\raggedright\\arraybackslash}p{1.6cm}>{\\raggedright\\arraybackslash}X>{\\raggedright\\arraybackslash}p{0.8cm}>{\\raggedright\\arraybackslash}p{1.2cm}}",
  "\\toprule",
  "\\textbf{Study} & \\textbf{Crop} & \\textbf{Trait} & \\textbf{Platform} & \\textbf{Manual reference} & \\textbf{Unit} & \\textbf{Correlation entered} & \\textbf{\\emph{n}} & \\textbf{Preprint}\\\\",
  "\\midrule", rowsA, "\\bottomrule", "\\end{tabularx}", "\\end{adjustwidth}",
  "\\noindent{\\footnotesize{Manual-reference types were coded from the description of the reference in each source. The first column gives the study, and the unit column gives what $n$ counts (plots, plants, images, accessions, lines or cultivars), so $n$ is not the number of independent genotypes. Correlations are entered as reported, with square roots of $R^2$ and Spearman coefficients treated as correlations.}}",
  "\\end{table}")
writeLines(outA, "paper/tables/tableA1_ma2_appraisal.tex")
cat("Wrote paper/tables/tableA1_ma2_appraisal.tex\n")
