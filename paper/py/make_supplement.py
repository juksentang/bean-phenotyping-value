"""Supplementary Figure S1 (PRISMA 2020-style flow diagrams) and Table S1 (search queries).

Run from the project root: python3 paper/py/make_supplement.py
Inputs: paper/lit/search/prisma/prisma_flow.csv and search_queries.csv (reconciled from the
on-disk search and screening records; see reconciliation_notes.md in the same folder).
"""
import os
import re
import sys
from pathlib import Path

import pandas as pd
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import ScaledTranslation

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import TOK, PCOL, FULL_W, setup, new_figure, claim, save

os.chdir(Path(__file__).resolve().parents[2])
setup()
PR = Path("paper/lit/search/prisma")
flow = pd.read_csv(PR / "prisma_flow.csv")
queries = pd.read_csv(PR / "search_queries.csv")
OUT = Path("paper/supplementary")
OUT.mkdir(exist_ok=True)

SOURCE_NAME = {
    "federated": "Federated search",
    "Crossref REST": "Crossref",
    "Europe PMC": "Europe PMC",
}


def source_label(s):
    for k, v in SOURCE_NAME.items():
        if k in str(s):
            return v
    return str(s)


def get(analysis, box, source=None):
    d = flow[(flow.analysis == analysis) & (flow.box == box)]
    if source is not None:
        d = d[d.source == source]
    if d.empty:
        return None
    v = d["count"].iloc[0]
    return None if pd.isna(v) or str(v) == "NA" else int(float(v))


def fmt(v):
    return "not recorded" if v is None else f"{v:,}"


# ---------------------------------------------------------------------------
# Figure S1: PRISMA 2020-style flow, one column per meta-analysis
# ---------------------------------------------------------------------------
# Main boxes follow the records down; grey side boxes give what left the flow at
# each step. Stage bands on the left follow PRISMA 2020 (identification,
# screening, included).
TITLES = {"MA1": "MA1 Realized genetic gain", "MA2": "MA2 Tool-manual agreement",
          "MA3": "MA3 Time and cost savings"}
COLORS = {"MA1": TOK["ink2"], "MA2": PCOL["error_reduction"], "MA3": PCOL["cost_reduction"]}
PRIMARY = {"MA1": "MA1.CB.primary", "MA2": "MA2.primary", "MA3": "MA3.time.primary"}

ROW_Y = [104, 87, 70, 53, 36, 17]          # box centers, top to bottom (axis units = mm)
BOX_H = 12
MAIN_W, SIDE_W, GAP = 28, 21, 4            # main box, side box, gap between them
COL_X = [18, 74, 130]                      # left edge of each meta-analysis column


def box(ax, x, yc, w, title, body, edge, fill, title_color=None):
    ax.add_patch(FancyBboxPatch((x, yc - BOX_H / 2), w, BOX_H, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=fill, ec=edge, lw=0.7))
    ax.text(x + w / 2, yc + (2.2 if body else 0), title, ha="center", va="center", fontsize=6.3,
            fontweight="bold", color=title_color or TOK["ink"])
    if body:
        ax.text(x + w / 2, yc - 2.6, body, ha="center", va="center", fontsize=5.6, color=TOK["ink2"],
                linespacing=1.15)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=TOK["muted"], lw=0.6, mutation_scale=6,
                                shrinkA=0, shrinkB=0))


def flow_column(ax, ma, x):
    col = COLORS[ma]
    db_total = get(ma, "Records identified from databases/registers", "all")
    other = get(ma, "Records identified via other methods", "other methods")
    dup = get(ma, "Duplicate records removed", "database search")
    after = get(ma, "Records after duplicates removed", "database search")
    screened = get(ma, "Records screened", "database search")
    excl = get(ma, "Records excluded", "database search")
    sought = get(ma, "Reports sought for retrieval", "all")
    notret = get(ma, "Reports not retrieved", "all")
    notass = get(ma, "Reports retrieved but not assessable", "all") or 0
    assessed = get(ma, "Reports assessed for eligibility (full text read)", "all")
    excl_ft = get(ma, "Reports excluded after full-text assessment", "all")
    incl = get(ma, "Studies included in review", "all")
    prim = PRIMARY[ma]
    k_st = get(ma, f"Studies contributing to primary pooled model ({prim})", "all")
    k_ef = get(ma, f"Effects contributing to primary pooled model ({prim})", "all")

    main = [(f"Records identified", f"databases {fmt(db_total)}\nother methods {fmt(other)}"),
            (f"After deduplication", f"n = {fmt(after)}"),
            (f"Records screened", f"n = {fmt(screened)}"),
            (f"Reports sought", f"n = {fmt(sought)}"),
            (f"Full texts assessed", f"n = {fmt(assessed)}"),
            (f"Studies included", f"n = {fmt(incl)}")]
    side = [None,
            ("Duplicates", f"n = {fmt(dup)}"),
            ("Excluded", f"n = {fmt(excl)}"),
            ("Not retrieved", f"n = {fmt(notret)}" + (f"\nnot assessable {notass}" if notass else "")),
            ("Excluded", f"n = {fmt(excl_ft)}"),
            ("Primary model", f"{fmt(k_ef)} effects,\n{fmt(k_st)} studies")]
    ax.text(x + (MAIN_W + GAP + SIDE_W) / 2, ROW_Y[0] + BOX_H / 2 + 5, TITLES[ma], ha="center", va="bottom",
            fontsize=6.8, fontweight="bold", color=col)
    for k, ((t, b), s) in enumerate(zip(main, side)):
        final = k == len(main) - 1
        box(ax, x, ROW_Y[k], MAIN_W, t, b, col, TOK["plane"] if final else TOK["surface"])
        if k:
            arrow(ax, x + MAIN_W / 2, ROW_Y[k - 1] - BOX_H / 2, x + MAIN_W / 2, ROW_Y[k] + BOX_H / 2)
        if s:
            sx = x + MAIN_W + GAP
            box(ax, sx, ROW_Y[k], SIDE_W, s[0], s[1], TOK["axis"], TOK["neutral"], title_color=TOK["ink2"])
            if not final:
                arrow(ax, x + MAIN_W, ROW_Y[k], sx, ROW_Y[k])
            else:
                ax.plot([x + MAIN_W, sx], [ROW_Y[k], ROW_Y[k]], color=TOK["muted"], lw=0.6, ls=(0, (2, 1.5)))


fig, axd = new_figure(FULL_W, 122, "a")
ax = axd["a"]
ax.set_xlim(0, 184); ax.set_ylim(8, 120); ax.axis("off")
for label, top, bottom in [("Identification", ROW_Y[0], ROW_Y[1]), ("Screening", ROW_Y[2], ROW_Y[4]),
                           ("Included", ROW_Y[5], ROW_Y[5])]:
    y0, y1 = bottom - BOX_H / 2, top + BOX_H / 2
    ax.add_patch(FancyBboxPatch((1, y0), 9, y1 - y0, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=TOK["neutral"], ec="none"))
    ax.text(5.5, (y0 + y1) / 2, label, rotation=90, ha="center", va="center", fontsize=6.3,
            fontweight="bold", color=TOK["ink2"])
for ma, x in zip(["MA1", "MA2", "MA3"], COL_X):
    flow_column(ax, ma, x)
issues = save(fig, "figS1_prisma", out=str(OUT))
print("figS1_prisma:", "ok" if not issues else issues)

# ---------------------------------------------------------------------------
# Table S2: search summary (the full query list is published as a CSV file)
# ---------------------------------------------------------------------------
SRC_SHORT = {"federated": "Federated search (OpenAlex, Semantic Scholar, Crossref)",
             "Crossref REST": "Crossref", "Europe PMC": "Europe PMC"}


def src_short(s):
    for k, v in SRC_SHORT.items():
        if k in str(s):
            return v
    return str(s)


q = queries.copy()
q["src"] = q.source.map(src_short)
q["failed"] = q.n_hits.astype(str).str.lower().eq("failed")
q["hits"] = pd.to_numeric(q.n_hits, errors="coerce")
rows = []
for ma, d in q.groupby("analysis", sort=True):
    first = True
    for src, e in d.groupby("src", sort=False):
        lead = TITLES.get(ma, ma) if first else ""
        first = False
        rows.append(f"{lead} & {src} & {len(e)} & {int(e.hits.sum()):,} & {int(e.failed.sum())}\\\\")
    rows.append(r"\addlinespace")
rows = rows[:-1]
date = pd.to_datetime(queries.date.dropna().unique()[0]).strftime("%-d %B %Y")
tab = [
    r"\begin{table}[!htbp]",
    rf"\caption{{Literature searches for the three meta-analyses, run on {date}. The federated search sent each "
    r"query to OpenAlex, Semantic Scholar and Crossref and returned at most 25 records per query and source. "
    r"Hits: records returned before deduplication; failed: queries with no usable response. The full list of "
    rf"{len(queries)} queries with their hit counts is archived with the code "
    r"(\texttt{paper/lit/search/prisma/search\_queries.csv}).}\label{tab:S1}",
    r"\small",
    r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}p{0.27\linewidth}X rrr@{}}",
    r"\toprule",
    r"Meta-analysis & Source & Queries & Hits & Failed\\",
    r"\midrule",
    *rows,
    r"\bottomrule",
    r"\end{tabularx}",
    r"\end{table}",
]
(OUT / "tableS1_queries.tex").write_text("\n".join(tab) + "\n")
print(f"tableS1_queries.tex: {len(queries)} queries summarized in {len([r for r in rows if '&' in r])} rows")
