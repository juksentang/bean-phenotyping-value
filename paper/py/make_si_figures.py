"""Supplementary Figure S3: breeder's-equation prediction versus simulation.

Run from the project root: python3 paper/py/make_si_figures.py
Input: layer1_genetic_sim/outputs/<TOOL_ROBUST_SUBDIR>/be_vs_sim.csv
(layer1_genetic_sim/R/08_breeders_equation_check.R). Output: paper/supplementary/figS3_be_vs_sim.pdf/.png.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import TOK, PCOL, LW, MS, COL_W, setup, new_figure, claim, clean, point, save

os.chdir(Path(__file__).resolve().parents[2])
setup()
bv = pd.read_csv(Path("layer1_genetic_sim/outputs") / os.environ.get("TOOL_ROBUST_SUBDIR", "sweep_v3_gxe")
                 / "be_vs_sim.csv")
bv["kappa"] = bv.kappa.round(2)
ER, CR_ACC = bv.error_reduction.max(), bv.cost_reduction.min()      # accuracy panel: largest error cut, smallest cost cut
CR = bv.cost_reduction.max()                                        # cost panel: largest cost cut, no error cut

# Standard error of the simulated accuracy effect (difference of two independent cell means)
se0 = bv[bv.error_reduction == 0].set_index(["h2_yield", "kappa", "cost_reduction"]).sim_se
bv["sim_acc_se"] = np.sqrt(bv.sim_se ** 2 + se0.reindex(
    pd.MultiIndex.from_frame(bv[["h2_yield", "kappa", "cost_reduction"]])).values ** 2)


def panel(ax, d, be_col, sim_col, se_col, color):
    h2s, ks = np.sort(d.h2_yield.unique()), np.sort(d.kappa.unique())
    xt, xl = [], []
    for i, h in enumerate(h2s):
        for j, k in enumerate(ks):
            x = i * (len(ks) + 0.8) + j
            r = d[(d.h2_yield == h) & (d.kappa == k)].iloc[0]
            ax.plot([x, x], [r[sim_col] - 1.96 * r[se_col], r[sim_col] + 1.96 * r[se_col]], color=color, lw=LW)
            point(ax, x, r[sim_col], color)
            ax.plot(x, r[be_col], marker="o", ms=MS, mfc=TOK["surface"], mec=TOK["ink2"], mew=1.0, ls="none")
            xt.append(x); xl.append(f"{k:g}")
        ax.text(i * (len(ks) + 0.8) + (len(ks) - 1) / 2, -0.14, f"h² = {h:.2f}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5, fontweight="bold")
    ax.axhline(0, color=TOK["axis"], lw=0.5)
    ax.set_xticks(xt, xl, fontsize=6.5)
    ax.set_xlabel("GxE-to-genetic variance ratio κ, by plot-level heritability", labelpad=24)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:g}%")
    ax.plot([], [], marker="o", ls="none", mfc=TOK["surface"], mec=TOK["ink2"], ms=MS, label="breeder's equation")
    ax.plot([], [], marker="o", ls="none", mfc=color, mec="none", ms=MS, label="simulation (95% interval)")
    ax.legend(loc="upper right", ncols=2, handletextpad=0.2)
    clean(ax, "y")


fig, axd = new_figure(COL_W, 120, "a;b")
panel(axd["a"], bv[(bv.error_reduction == ER) & (bv.cost_reduction == CR_ACC)], "be_fixed", "sim_acc", "sim_acc_se",
      PCOL["error_reduction"])
axd["a"].set_ylabel(f"Per-cycle uplift from a\n{100 * ER:.0f}% error reduction")
claim(axd["a"], "a", "The breeder's equation overstates the accuracy effect")
panel(axd["b"], bv[(bv.error_reduction == 0) & (bv.cost_reduction == CR)], "be_opt", "sim_vsbase", "sim_se",
      PCOL["cost_reduction"])
axd["b"].set_ylabel(f"Per-cycle uplift from a\n{100 * CR:.0f}% cost reduction")
claim(axd["b"], "b", "It misses how the cost channel depends on GxE")
issues = save(fig, "figS3_be_vs_sim", out="paper/supplementary")
print("figS3_be_vs_sim:", "ok" if not issues else issues)
