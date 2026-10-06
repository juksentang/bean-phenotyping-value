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
# Figure S1: one flow diagram per meta-analysis
# ---------------------------------------------------------------------------
TITLES = {"MA1": "MA1 Realized genetic gain", "MA2": "MA2 Tool-manual agreement",
          "MA3": "MA3 Time and cost savings"}
COLORS = {"MA1": TOK["ink2"], "MA2": PCOL["error_reduction"], "MA3": PCOL["cost_reduction"]}
PRIMARY = {"MA1": "MA1.CB.primary", "MA2": "MA2.primary", "MA3": "MA3.time.primary"}


def box(ax, x, y, w, h, text, edge, fill=None, bold_first=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=fill or TOK["surface"], ec=edge, lw=0.7))
    lines = text.split("\n")
    top = ax.transData + ScaledTranslation(0, -3 / 72, ax.figure.dpi_scale_trans)
    body = ax.transData + ScaledTranslation(0, -12 / 72, ax.figure.dpi_scale_trans)
    ax.text(x + w / 2, y + h, lines[0], ha="center", va="top", fontsize=6.2, transform=top,
            fontweight="bold" if bold_first else "normal")
    if len(lines) > 1 and "\n".join(lines[1:]).strip():
        ax.text(x + w / 2, y + h, "\n".join(lines[1:]), ha="center", va="top", fontsize=5.7,
                color=TOK["ink2"], linespacing=1.2, transform=body)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=TOK["muted"], lw=0.6, mutation_scale=6))


def flow_panel(ax, ma):
    ax.set_xlim(0, 100); ax.set_ylim(0, 128); ax.axis("off")
    col = COLORS[ma]
    db_rows = flow[(flow.analysis == ma) & (flow.box == "Records identified from databases/registers")
                   & (flow.source != "all")]
    db_lines = [f"{source_label(s)}: {fmt(None if pd.isna(c) else int(c))}"
                for s, c in zip(db_rows.source, db_rows["count"]) if not pd.isna(c) and int(c) > 0]
    db_total = get(ma, "Records identified from databases/registers", "all")
    other = get(ma, "Records identified via other methods", "other methods")
    dup = get(ma, "Duplicate records removed", "database search")
    after = get(ma, "Records after duplicates removed", "database search")
    screened = get(ma, "Records screened", "database search")
    excl = get(ma, "Records excluded", "database search")
    sought_db = get(ma, "Reports sought for retrieval", "database search")
    sought_oth = get(ma, "Reports sought for retrieval", "other methods")
    notret = get(ma, "Reports not retrieved", "all")
    notass = get(ma, "Reports retrieved but not assessable", "all") or 0
    assessed = get(ma, "Reports assessed for eligibility (full text read)", "all")
    excl_ft = get(ma, "Reports excluded after full-text assessment", "all")
    incl = get(ma, "Studies included in review", "all")
    prim = PRIMARY[ma]
    k_st = get(ma, f"Studies contributing to primary pooled model ({prim})", "all")
    k_ef = get(ma, f"Effects contributing to primary pooled model ({prim})", "all")

    L, R, W, WR = 3, 60, 52, 37
    box(ax, L, 103, W, 24, f"Records from databases (n = {fmt(db_total)})\n" + "\n".join(db_lines), col)
    box(ax, R, 103, WR, 24, f"Other methods (n = {fmt(other)})\ncitation chasing, related\nproject searches", col)
    box(ax, L, 84, W, 14, f"After duplicates removed (n = {fmt(after)})\nremoved by DOI: {fmt(dup)}", col)
    box(ax, L, 65, W, 14, f"Records screened: {fmt(screened)}\nexcluded: {fmt(excl)}", col)
    box(ax, L, 46, W, 11, f"Reports sought (n = {fmt(sought_db)})\n", col)
    box(ax, R, 46, WR, 11, f"Reports sought (n = {fmt(sought_oth)})\n", col)
    box(ax, L, 26, W + 4 + WR, 14, f"Full texts assessed (n = {fmt(assessed)})\n"
        f"not retrieved: {fmt(notret)}; retrieved, not assessable: {fmt(notass)}; excluded: {fmt(excl_ft)}", col)
    box(ax, L, 6, W + 4 + WR, 14, f"Studies included (n = {fmt(incl)})\n"
        f"primary pooled model: {fmt(k_ef)} effects from {fmt(k_st)} studies or programs", col,
        fill=TOK["plane"])
    arrow(ax, L + W / 2, 103, L + W / 2, 98)
    arrow(ax, L + W / 2, 84, L + W / 2, 79)
    arrow(ax, L + W / 2, 65, L + W / 2, 57)
    arrow(ax, R + WR / 2, 103, R + WR / 2, 57)
    arrow(ax, L + W / 2, 46, L + W / 2, 40)
    arrow(ax, R + WR / 2, 46, R + WR / 2, 40)
    arrow(ax, 50, 26, 50, 20)


fig, axd = new_figure(FULL_W, 215, "a;b;c")
for key, ma in zip("abc", ["MA1", "MA2", "MA3"]):
    flow_panel(axd[key], ma)
    claim(axd[key], key, TITLES[ma])
issues = save(fig, "figS1_prisma", out=str(OUT))
print("figS1_prisma:", "ok" if not issues else issues)

# ---------------------------------------------------------------------------
# Table S1: search queries
# ---------------------------------------------------------------------------
def tex_escape(s):
    s = "" if pd.isna(s) else str(s)
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
                 ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\^{}")]:
        s = s.replace(a, b)
    return s


rows = []
for ma, d in queries.groupby("analysis", sort=True):
    rows.append(rf"\multicolumn{{4}}{{l}}{{\textbf{{{TITLES.get(ma, ma)}}}}}\\")
    for _, r in d.iterrows():
        q = "not recorded" if pd.isna(r.query_string) or str(r.query_string) == "NA" else r.query_string
        hits = r.n_hits if not pd.isna(r.n_hits) else "not recorded"
        rows.append(f"{tex_escape(r.query_id)} & {tex_escape(source_label(r.source))} & "
                    f"{tex_escape(q)} & {tex_escape(hits)}\\\\")
date = queries.date.dropna().unique()
tab = [
    r"\begin{longtable}{@{}>{\raggedright\arraybackslash}p{0.06\linewidth}>{\raggedright\arraybackslash}p{0.17\linewidth}>{\raggedright\arraybackslash}p{0.58\linewidth}r@{}}",
    rf"\caption{{Search queries for the three meta-analyses, run on {', '.join(pd.to_datetime(d).strftime('%-d %B %Y') for d in date)}. "
    r"Federated search: one query sent to OpenAlex, Semantic Scholar and Crossref, which "
    r"returned at most 25 records per query and source. Hits: records returned; failed: no usable "
    r"response.}\label{tab:S1}\\",
    r"\toprule ID & Source & Query & Hits\\ \midrule \endfirsthead",
    r"\toprule ID & Source & Query & Hits\\ \midrule \endhead",
    r"\bottomrule \endfoot",
    *rows,
    r"\end{longtable}",
]
(OUT / "tableS1_queries.tex").write_text("\n".join(tab) + "\n")
print(f"tableS1_queries.tex: {len(queries)} queries")
