"""Graphical abstract (MDPI: PNG, 2800 x 5500 px max, height x width).

Run from the project root: python3 paper/py/make_graphical_abstract.py
Numbers come from paper/numbers.csv and the robustness-block cell summaries.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from PIL import Image, ImageChops

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import TOK, PCOL, setup, check_layout

os.chdir(Path(__file__).resolve().parents[2])
setup()
mpl.rcParams["font.family"] = "Arial"

num = pd.read_csv("paper/numbers.csv", dtype=str).set_index("name")["value"]
rb = pd.read_csv("layer1_genetic_sim/outputs/sweep_v2_robust/cell_summaries.csv")
rb["uplift"] = 100 * rb.mean_dg_diff / rb.mean_dg_trad
CYCLE = 10
years = lambda u: CYCLE - CYCLE / (1 + np.asarray(u) / 100)
err_lo, err_hi = 0.10, rb.error_reduction.max()
h2_main = 0.30
m = rb[np.isclose(rb.h2_yield, h2_main)].groupby("error_reduction").uplift.mean()
yr_evid, yr_opt = float(years(m.loc[err_lo])), float(years(m.loc[err_hi]))

W_IN, H_IN, DPI = 11.0, 5.6, 500               # 5500 x 2800 px
fig = plt.figure(figsize=(W_IN, H_IN), facecolor=TOK["surface"])
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 110); ax.set_ylim(0, 56); ax.axis("off")

# --- left: the tool and its three levers -----------------------------------
ax.add_patch(FancyBboxPatch((6, 15), 12, 27, boxstyle="round,pad=0,rounding_size=2.2",
                            fc=TOK["plane"], ec=TOK["ink2"], lw=1.6))
ax.add_patch(FancyBboxPatch((7.4, 19), 9.2, 19.5, boxstyle="round,pad=0,rounding_size=0.8",
                            fc="#e7f0fb", ec="none"))
for k, h in enumerate([0.55, 0.85, 0.65, 0.95, 0.7]):          # a stylised crop row on the screen
    x = 8.6 + k * 1.75
    ax.plot([x, x], [21, 21 + 9 * h], color=PCOL["cost_reduction"], lw=2.4, solid_capstyle="round")
    ax.add_patch(Circle((x, 21 + 9 * h), 0.55, fc=PCOL["cost_reduction"], ec="none"))
ax.add_patch(Circle((12, 40.3), 0.6, fc=TOK["ink2"], ec="none"))
ax.text(12, 47.5, "Smartphone AI\nphenotyping tool", ha="center", va="center", fontsize=15,
        fontweight="bold", color=TOK["ink"], linespacing=1.1)
levers = [("error_reduction", "Accuracy"), ("cost_reduction", "Cost per plot"),
          ("cycle_compression", "Breeding-cycle time")]
for k, (p, lab) in enumerate(levers):
    y = 34 - k * 8
    ax.plot([18.4, 22.5], [y, y], color=TOK["axis"], lw=1.2)
    ax.add_patch(Circle((23.6, y), 1.1, fc=PCOL[p], ec=TOK["surface"], lw=1.5))
    ax.text(25.6, y, lab, va="center", fontsize=13.5, color=TOK["ink"],
            fontweight="bold" if p == "cycle_compression" else "normal")

# --- middle: the modelling chain --------------------------------------------
steps = [("Breeding simulation", f"{num['NcellsOK']} tool profiles"),
         ("Genetic gain", f"calibrated to {num['GainTargetPct']}%/yr"),
         ("Adoption and NPV", f"{num['CfgHorizon']}-year horizon")]
for k, (t, b) in enumerate(steps):
    y = 41 - k * 12
    ax.add_patch(FancyBboxPatch((46, y - 4.6), 20.5, 9.2, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=TOK["plane"], ec=TOK["grid"], lw=1.2))
    ax.text(56.25, y + 1.4, t, ha="center", va="center", fontsize=12.5, fontweight="bold", color=TOK["ink"])
    ax.text(56.25, y - 2.0, b, ha="center", va="center", fontsize=11, color=TOK["ink2"])
    if k < 2:
        ax.annotate("", xy=(56.25, y - 7.2), xytext=(56.25, y - 4.8),
                    arrowprops=dict(arrowstyle="-|>", color=TOK["muted"], lw=1.4, mutation_scale=14))
ax.annotate("", xy=(45.6, 29), xytext=(39.5, 29),
            arrowprops=dict(arrowstyle="-|>", color=TOK["muted"], lw=1.6, mutation_scale=16))
ax.annotate("", xy=(72, 29), xytext=(67, 29),
            arrowprops=dict(arrowstyle="-|>", color=TOK["muted"], lw=1.6, mutation_scale=16))

# --- right: the headline and the exchange rate -------------------------------
ax.text(73, 49.5, "Time beats accuracy", fontsize=23, fontweight="bold", color=TOK["ink"], va="center")
ax.text(73, 44.3, "Equivalent years of cycle compression", fontsize=11.5, color=TOK["ink2"], va="center")
x0, scale = 73, 22                                   # 1 year = 22 units of bar
bars = [(f"Accuracy, {int(err_lo * 100)}% error cut (evidence)", yr_evid, PCOL["error_reduction"], 0.45),
        (f"Accuracy, {int(err_hi * 100)}% error cut (optimistic)", yr_opt, PCOL["error_reduction"], 1.0),
        ("One year shorter breeding cycle", 1.0, PCOL["cycle_compression"], 1.0)]
for k, (lab, v, c, a) in enumerate(bars):
    y = 37 - k * 8.2
    ax.text(x0, y + 2.6, lab, fontsize=11.5, color=TOK["ink"], va="center")
    ax.add_patch(FancyBboxPatch((x0, y - 1.6), max(v * scale, 0.6), 3.2,
                                boxstyle="round,pad=0,rounding_size=0.8", fc=c, ec="none", alpha=a))
    ax.text(x0 + max(v * scale, 0.6) + 1.2, y, f"{12 * v:.1f} months" if v < 1 else "12 months",
            fontsize=12, fontweight="bold", color=TOK["ink"], va="center")
ax.text(73, 12.4, "Worth measuring first in a pilot (EVPPI)", fontsize=11.5, color=TOK["ink2"], va="center")
ax.text(73, 8.4, f"cycle time \\${num['EVPPIIndiffCyc']}M  vs  accuracy \\${num['EVPPIIndiffErrRed']}M",
        fontsize=12.5, fontweight="bold", color=TOK["ink"], va="center")
ax.text(6, 3.0, f"Simulated East African common-bean program; plot-level h² = {h2_main:.2f}",
        fontsize=9.5, color=TOK["muted"], va="center")

fig.canvas.draw()
issues = check_layout(fig, fig.canvas.get_renderer())
out = Path("paper/figures/graphical_abstract.png")
fig.savefig(out, dpi=DPI, facecolor=TOK["surface"])
plt.close(fig)
img = Image.open(out).convert("RGB")                 # trim surrounding blank space (MDPI requirement)
bg = Image.new("RGB", img.size, img.getpixel((0, 0)))
l, t, r, b = ImageChops.difference(img, bg).getbbox()
pad = 60
img.crop((max(l - pad, 0), max(t - pad, 0), min(r + pad, img.width), min(b + pad, img.height))).save(out)
print("graphical_abstract:", "ok" if not issues else issues)
