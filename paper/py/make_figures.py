"""Manuscript figures 1-6. Run from the project root: python3 paper/py/make_figures.py

Reads the Layer-1 sweep summaries, the Layer-2 outputs and the meta-analysis
results; no value shown in a figure is entered by hand. Writes PDF + PNG to
paper/figures/ and the figure captions to paper/figures/captions.md. Figures are
drawn at their printed size (column or full width); a layout check reports any
text that leaves the canvas or overlaps other text, and the run fails if any
is found.

Storyline: (1) design and input colors, (2) time beats accuracy, (3) why:
per-cycle gain is near its ceiling, (4) value follows cycle time, (5) what a
pilot should measure, (6) evidence behind the inputs.
"""
import os
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import (TOK, ORD4, CMAP_BLUE, PCOL, TECH, ECON, LW, LW_HAIR, MS, COL_W, FULL_W, PARAM_LABEL,
                      KEY_LABEL, setup, clean, point, label_column, raincloud, dumbbell, new_figure, claim, save)

os.chdir(Path(__file__).resolve().parents[2])
setup()


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def cfg(name):
    with open(f"config/{name}.yaml") as fh:
        return yaml.safe_load(fh)


econ, breed, gen = cfg("economic_params"), cfg("breeding_params"), cfg("genetic_params")
OUT = Path(os.environ.get("TOOL_L2_OUT", "layer2_economics/outputs_v2"))
L1_SUBS = os.environ.get("TOOL_L1_SUBDIR", "sweep_v2,sweep_v2_extA,sweep_v2_extB").split(",")
LIT = Path("paper/lit")

l1 = pd.concat([pd.read_csv(f"layer1_genetic_sim/outputs/{s.strip()}/cell_summaries.csv")
                for s in L1_SUBS], ignore_index=True)
CYCLE = breed["cycle_years"]


def annual_uplift(per_cycle_pct, cyc):
    """Annual-gain uplift (%) from the per-cycle uplift and years of cycle compression."""
    return 100 * ((1 + np.asarray(per_cycle_pct) / 100) * CYCLE / (CYCLE - np.asarray(cyc)) - 1)


def years_equivalent(uplift_pct):
    """Years of cycle compression that raise annual gain by the same percentage."""
    return CYCLE - CYCLE / (1 + np.asarray(uplift_pct) / 100)


l1["uplift"] = 100 * l1.mean_dg_diff / l1.mean_dg_trad
l1["annual"] = annual_uplift(l1.uplift, l1.cycle_compression)
rb_path = Path("layer1_genetic_sim/outputs/sweep_v2_robust/cell_summaries.csv")
rb = pd.read_csv(rb_path) if rb_path.exists() else None
if rb is not None:
    rb["uplift"] = 100 * rb.mean_dg_diff / rb.mean_dg_trad

CYC_YEAR_PCT = float(annual_uplift(0, 1))
META_ERR = econ["tool_error_reduction_frac"]
META_COST = econ["tool_cost_reduction_frac"]
H2 = gen["h2_yield"]
is_meta = pd.Series(np.isclose(l1.error_reduction, META_ERR) & np.isclose(l1.cost_reduction, META_COST),
                    index=l1.index)

lk = pd.read_parquet(OUT / "lookup_table.parquet",
                     columns=["error_reduction", "cost_reduction", "cycle_compression", "bean_price",
                              "discount_rate", "adoption_ceiling", "total_area_ha", "npv", "eaa_per_ha_yr"])


def nearest(v, values):
    grid = np.sort(np.unique(values))
    return grid[np.argmin(np.abs(grid - v))]


BASE = dict(bean_price=nearest(econ["bean_price_usd_per_kg"], lk.bean_price),
            discount_rate=nearest(econ["discount_rate"], lk.discount_rate),
            adoption_ceiling=nearest(econ["adoption_ceiling"], lk.adoption_ceiling),
            total_area_ha=nearest(econ["total_area_ha"], lk.total_area_ha))
mr = pd.read_csv(LIT / "meta_results.csv").set_index("analysis_id")
mi = pd.read_csv(LIT / "meta_inputs_used.csv", low_memory=False)
CAPTIONS, ISSUES = {}, {}


def pct(v, d=0):
    return f"{100 * v:.{d}f}%"


def money(x):
    """USD in millions -> compact label ($0.7M, $17M, $10k)."""
    if abs(x) >= 10:
        return f"${x:,.0f}M"
    if abs(x) >= 0.1:
        return f"${x:,.1f}M"
    return f"${1000 * x:,.0f}k"


def finish(fig, name, caption):
    ISSUES[name] = save(fig, name)
    CAPTIONS[name] = caption


# ---------------------------------------------------------------------------
# Figure 1. Study design and color code for the model inputs
# ---------------------------------------------------------------------------
def figure1():
    fig, axd = new_figure(COL_W, 52, "a")
    ax = axd["a"]
    ax.set_xlim(0, 100); ax.set_ylim(-1.8, 40); ax.axis("off")
    n_econ = lk.groupby(["bean_price", "discount_rate", "adoption_ceiling", "total_area_ha"]).ngroups
    steps = [("Tool profile", "5 technical\ninputs"),
             ("Breeding\nsimulation", f"{len(l1):,} profiles\nx {int(l1.n_reps.iloc[0])} paired replicates"),
             ("Genetic\ngain", f"versus check varieties,\n{econ['gain_target_pct_per_yr']:.1f}% per year calibration"),
             ("Adoption\nand NPV", f"{n_econ} economic states,\n{econ['evaluation_horizon']}-yr horizon"),
             ("Attribution and\nvalue of information", "Sobol, Shapley,\nPAWN, EVPPI")]
    xs = np.linspace(9, 91, len(steps))
    ax.plot([xs[0], xs[-1]], [33, 33], color=TOK["axis"], lw=LW, zorder=1)
    for i, (x, (t, b)) in enumerate(zip(xs, steps)):
        ax.plot(x, 33, marker="o", ms=11, mfc=TOK["surface"], mec=TOK["ink2"], mew=0.9, zorder=2)
        ax.text(x, 33, str(i + 1), ha="center", va="center", fontsize=6.5, fontweight="bold", zorder=3)
        ax.text(x, 27.5, t, ha="center", va="top", fontsize=6.5, fontweight="bold", linespacing=1.05)
        ax.text(x, 19.5, b, ha="center", va="top", fontsize=6.5, color=TOK["ink2"], linespacing=1.15)
    for row, (label, params) in enumerate([("Technical", TECH[:3]), ("", TECH[3:]), ("Economic", ECON)]):
        y = 7.2 - row * 3.4
        ax.text(0.5, y, label, fontsize=6.5, fontweight="bold", va="center")
        for k, p in enumerate(params):
            x = 16 + k * 28
            ax.plot(x, y, marker="o", ms=4.2, mfc=PCOL[p], mec=TOK["surface"], mew=0.6)
            ax.text(x + 1.6, y, KEY_LABEL.get(p, PARAM_LABEL[p]), fontsize=6.5, va="center", color=TOK["ink2"])
    ax.plot([0.5, 99.5], [11.2, 11.2], color=TOK["grid"], lw=LW_HAIR)
    finish(fig, "fig1_design",
           "Study design. Each tool profile (five technical inputs) is simulated in a stochastic common-bean "
           "breeding program; its genetic gain over the check varieties is calibrated to realized progress, "
           "valued through variety adoption and discounting over the economic states, and attributed to its "
           "inputs. The colors identify the model inputs in all later figures.")


# ---------------------------------------------------------------------------
# Figure 2. Time beats accuracy
# ---------------------------------------------------------------------------
def figure2():
    fig, axd = new_figure(FULL_W, 78, "ab", width_ratios=[1.2, 1])

    ax = axd["b"]
    emax = None
    if rb is not None:
        t = rb.groupby(["h2_yield", "gxe_cor_onstation", "error_reduction"]).uplift.mean().reset_index()
        t["years"] = years_equivalent(t.uplift)
        h2s = np.sort(t.h2_yield.unique())
        emax = t.error_reduction.max()
        for k, h in enumerate(h2s):
            d = t[t.h2_yield == h]
            lo = d[np.isclose(d.error_reduction, META_ERR)].years
            hi = d[np.isclose(d.error_reduction, emax)].years
            dumbbell(ax, k, lo.mean(), hi.mean(), PCOL["error_reduction"])
            ax.scatter(np.r_[lo, hi], np.full(len(lo) + len(hi), k - 0.22), s=5, color=PCOL["error_reduction"],
                       alpha=0.5, lw=0)
            ax.text(hi.mean() + 0.035, k, f"{12 * hi.mean():.1f} mo", va="center", fontsize=6.5)
            ax.text(lo.mean(), k + 0.2, f"{12 * lo.mean():.1f} mo", va="bottom", ha="center", fontsize=6.5,
                    color=TOK["ink2"])
        ax.axvline(1, color=PCOL["cycle_compression"], lw=LW)
        ax.text(0.97, -0.5, "one year of\ncycle compression", fontsize=6.5, ha="right", va="bottom",
                color=PCOL["cycle_compression"])
        ax.set_yticks(range(len(h2s)), [f"h² = {h:.2f}" for h in h2s])
        ax.set_ylim(-0.55, len(h2s) + 0.3)
        ax.set_xlim(-0.08, 1.38)
        ax.set_xlabel("Equivalent years of cycle compression")
        ax.plot([], [], marker="o", ls="none", mfc=TOK["surface"], mec=PCOL["error_reduction"], ms=MS,
                label=f"error reduction {pct(META_ERR)} (evidence)")
        ax.plot([], [], marker="o", ls="none", mfc=PCOL["error_reduction"], mec="none", ms=MS,
                label=f"error reduction {pct(emax)} (optimistic)")
        ax.legend(loc="upper left", handletextpad=0.3)
        clean(ax, "x")
    claim(ax, "b", "Accuracy buys months, not years")

    ax = axd["a"]
    cycs = np.sort(l1.cycle_compression.unique())
    groups = [l1.loc[l1.cycle_compression == c, "annual"].values for c in cycs]
    hls = [is_meta[l1.cycle_compression == c].values for c in cycs]
    raincloud(ax, groups, list(range(len(cycs))), [PCOL["cycle_compression"]] * len(cycs), orient="h",
              highlight=hls, highlight_color=PCOL["error_reduction"])
    for k, gv in enumerate(groups):
        ax.text(np.median(gv), k + 0.47, f"{np.median(gv):.0f}%", fontsize=6.5, ha="center")
    ax.set_yticks(range(len(cycs)), [f"{c} yr" for c in cycs])
    ax.set_ylim(-0.5, len(cycs) - 0.15)
    ax.set_xlim(min(-3, min(g.min() for g in groups) - 1), max(g.max() for g in groups) + 2)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
    ax.set_xlabel("Annual genetic-gain uplift versus manual phenotyping\n(compression applied as a scenario, not simulated)")
    ax.set_ylabel("Cycle compression")
    ax.scatter([], [], s=9, color=PCOL["error_reduction"], label="evidence-based tool profile")
    ax.legend(loc="lower right", handletextpad=0.2)
    clean(ax, "x")
    claim(ax, "a", "Cycle time separates outcomes; accuracy and cost barely do")
    finish(fig, "fig2_time_beats_accuracy",
           "Time beats accuracy. (a) Annual genetic-gain uplift over manual phenotyping across all simulated tool "
           "profiles, by years of cycle compression, which is applied as a scenario and not simulated (cloud: "
           "density; rain: profiles; bar and ring: interquartile range and median). Orange: evidence-based "
           f"profile (error reduction {pct(META_ERR)}, cost reduction {pct(META_COST)}). (b) Exchange rate between "
           "accuracy and cycle time: years of breeding-cycle compression that raise annual genetic gain as much "
           "as the tool's per-cycle uplift, at the "
           f"evidence-based ({pct(META_ERR)}) and an optimistic ({pct(emax) if emax is not None else 'maximum'}) "
           "error reduction, by plot-level heritability (large markers: mean over on-station GxE parameter "
           "0.4-0.8; small dots: individual settings).")


# ---------------------------------------------------------------------------
# Figure 3. Why accuracy and cost move gain only a little
# ---------------------------------------------------------------------------
def figure3():
    fig, axd = new_figure(COL_W, 112, "ab;cc", width_ratios=[1, 1.25], height_ratios=[1, 0.95])
    ymin = min(l1.uplift.min(), -1) - 0.5
    ymax = max(l1.uplift.quantile(0.998), CYC_YEAR_PCT) + 1
    for key, param, lab, ttl in [("a", "error_reduction", "Error reduction (%)", "Accuracy lever"),
                                 ("b", "cost_reduction", "Cost reduction (%)", "Cost lever")]:
        ax = axd[key]
        levels = np.sort(l1[param].unique())
        groups = [l1.loc[np.isclose(l1[param], v), "uplift"].values for v in levels]
        meta_v = META_ERR if param == "error_reduction" else META_COST
        raincloud(ax, groups, list(range(len(levels))), [PCOL[param]] * len(levels), width=0.42,
                  point_size=1.0, max_points=200, jitter=0.06)
        j = int(np.argmin(np.abs(levels - meta_v)))
        ax.axvspan(j - 0.48, j + 0.48, color=TOK["neutral"], zorder=0, lw=0)
        ax.axhline(CYC_YEAR_PCT, color=PCOL["cycle_compression"], lw=LW)
        ax.axhline(0, color=TOK["axis"], lw=LW_HAIR)
        ax.set_xticks(range(len(levels)), [f"{100 * v:.0f}" for v in levels])
        ax.set_xlim(-0.6, len(levels) - 0.3)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(lab)
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
        if key == "a":
            ax.set_ylabel("Per-cycle gain uplift")
        else:
            ax.sharey(axd["a"])
            ax.tick_params(labelleft=False)
            ax.text(len(levels) - 0.35, CYC_YEAR_PCT + 0.25, "one year of cycle compression", fontsize=6.5,
                    ha="right", va="bottom", color=PCOL["cycle_compression"])
        clean(ax, "y")
        claim(ax, key, ttl)

    ax = axd["c"]
    if rb is not None:
        t = rb.groupby(["h2_yield", "gxe_cor_onstation", "error_reduction"]).uplift.mean().reset_index()
        e0, e1 = t.error_reduction.min(), t.error_reduction.max()
        h2s, gxes = np.sort(t.h2_yield.unique()), np.sort(t.gxe_cor_onstation.unique())
        xt, xl = [], []
        for i, h in enumerate(h2s):
            for j, gx in enumerate(gxes):
                x = i * (len(gxes) + 0.8) + j
                d = t[(t.h2_yield == h) & (t.gxe_cor_onstation == gx)]
                a = d[np.isclose(d.error_reduction, e0)].uplift.iloc[0]
                b = d[np.isclose(d.error_reduction, e1)].uplift.iloc[0]
                ax.plot([x, x], [a, b], color=PCOL["error_reduction"], lw=LW * 1.6, alpha=0.45,
                        solid_capstyle="round")
                ax.plot(x, a, marker="o", ms=MS, mfc=TOK["surface"], mec=PCOL["error_reduction"], mew=1.0)
                point(ax, x, b, PCOL["error_reduction"])
                xt.append(x); xl.append(f"{gx:.1f}")
            ax.text(i * (len(gxes) + 0.8) + (len(gxes) - 1) / 2, -0.17, f"h² = {h:.2f}",
                    transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5, fontweight="bold")
        ax.axhline(CYC_YEAR_PCT, color=PCOL["cycle_compression"], lw=LW)
        ax.axhline(0, color=TOK["axis"], lw=LW_HAIR)
        ax.set_xticks(xt, xl, fontsize=6.5)
        ax.set_xlabel("On-station GxE parameter, by plot-level heritability", labelpad=20)
        ax.set_ylabel(f"Per-cycle uplift, {pct(e0)} to {pct(e1)}")
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
        ax.plot([], [], marker="o", ls="none", mfc=TOK["surface"], mec=PCOL["error_reduction"], ms=MS,
                label=f"error reduction {pct(e0)}")
        ax.plot([], [], marker="o", ls="none", mfc=PCOL["error_reduction"], mec="none", ms=MS,
                label=f"error reduction {pct(e1)}")
        ax.legend(loc="upper right", ncols=2, handletextpad=0.2)
        clean(ax, "y")
    claim(ax, "c", "Only at low heritability does accuracy rival a year of compression")
    finish(fig, "fig3_gain_ceiling",
           "Accuracy and cost savings move per-cycle genetic gain by a few percent. (a, b) Per-cycle "
           "genetic-gain uplift of the tool over manual phenotyping across all simulated tool profiles, by "
           f"measurement-error reduction (a) and whole-plot cost reduction (b) (h2 = {H2:.2f}; cloud: density; "
           "rain: profiles; bar and ring: interquartile range and median); shaded columns: evidence-based "
           "values. (c) Robustness: per-cycle uplift with no and with the largest error reduction for each "
           "combination of plot-level heritability and on-station GxE parameter. Blue lines: the annual-gain "
           "increase from one year of cycle compression.")


# ---------------------------------------------------------------------------
# Figure 4. Value follows cycle time
# ---------------------------------------------------------------------------
def figure4():
    fig, axd = new_figure(FULL_W, 104, "ab;cd", width_ratios=[1.1, 1])

    ax = axd["a"]
    cycs = np.sort(lk.cycle_compression.unique())
    data = []
    for c in cycs:
        v = lk.loc[lk.cycle_compression == c, "npv"].values / 1e6
        v = v[v > 0]                                   # violin: positive values only (log scale)
        data.append(np.log10(v[v >= np.percentile(v, 1)]))
    parts = ax.violinplot(data, positions=range(len(cycs)), widths=0.8, showextrema=False)
    for body, col in zip(parts["bodies"], ORD4):
        body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(0.55); body.set_linewidth(0.6)
    for k, c in enumerate(cycs):
        allv = lk.loc[lk.cycle_compression == c, "npv"].values / 1e6      # bar, ring and label: all states
        q1, med, q3 = np.percentile(allv, [25, 50, 75])
        assert q1 > 0
        ax.plot([k, k], [np.log10(q1), np.log10(q3)], color=TOK["ink"], lw=2.2, solid_capstyle="butt")
        ax.plot(k, np.log10(med), marker="o", ms=3.2, mfc=TOK["surface"], mec=TOK["ink"], mew=0.8)
        ev = lk.loc[(lk.cycle_compression == c) & np.isclose(lk.error_reduction, META_ERR)
                    & np.isclose(lk.cost_reduction, META_COST), "npv"] / 1e6
        if len(ev):
            point(ax, k + 0.3, np.log10(ev.median()), PCOL["error_reduction"], marker="D", size=MS * 0.9)
        ax.text(k - 0.42, np.log10(med), money(med), ha="right", va="center", fontsize=6.5)
    lo_t, hi_t = int(np.floor(min(map(np.min, data)))), int(np.ceil(max(map(np.max, data))))
    ax.set_yticks(range(lo_t, hi_t + 1), [money(10.0 ** e) for e in range(lo_t, hi_t + 1)])
    ax.set_xlim(-1.0, len(cycs) - 0.4)
    ax.set_xticks(range(len(cycs)), [f"{c} yr" for c in cycs])
    ax.set_xlabel("Cycle compression")
    ax.set_ylabel("NPV over all states (log scale;\nviolins show positive values only)")
    ax.plot([], [], marker="D", ls="none", mfc=PCOL["error_reduction"], mec="none", ms=MS * 0.9,
            label="evidence-based profile")
    ax.legend(loc="lower right", handletextpad=0.2)
    clean(ax, "y")
    claim(ax, "a", "Each year of compression multiplies NPV")

    ax = axd["b"]
    sl = lk[(lk.discount_rate == BASE["discount_rate"]) & (lk.adoption_ceiling == BASE["adoption_ceiling"])
            & (lk.total_area_ha == BASE["total_area_ha"])]
    t = sl.groupby(["cycle_compression", "bean_price"]).eaa_per_ha_yr.median().unstack()
    norm = Normalize(0, t.values.max())
    ax.pcolormesh(np.arange(t.shape[1] + 1), np.arange(t.shape[0] + 1), t.values, cmap=CMAP_BLUE, norm=norm,
                  edgecolors=TOK["surface"], linewidth=1.2)
    for (i, j), v in np.ndenumerate(t.values):
        ax.text(j + 0.5, i + 0.5, f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                color=TOK["surface"] if norm(v) > 0.55 else TOK["ink"])
    ax.set_xticks(np.arange(t.shape[1]) + 0.5, [f"{p:.2f}" for p in t.columns])
    ax.set_yticks(np.arange(t.shape[0]) + 0.5, [f"{c} yr" for c in t.index])
    ax.set_xlabel("Bean farm-gate price (USD per kg)")
    ax.set_ylabel("Cycle compression")
    ax.grid(False); ax.spines[:].set_visible(False)
    claim(ax, "b", "Value per ha per year (USD), median over profiles")

    ax = axd["c"]
    dyn = pd.read_csv(OUT / "dynamic_sensitivity_ST_wide.csv")
    last = dyn.iloc[-1].drop("year").sort_values(ascending=False)
    shown = list(last.index[:4])
    for p in dyn.columns.drop("year"):
        ax.plot(dyn.year, dyn[p].clip(0, 1), color=PCOL[p], lw=LW if p in shown else 0.6,
                alpha=1 if p in shown else 0.5)
    label_column(ax, [(dyn.year.iloc[-1], dyn[p].clip(0, 1).iloc[-1], PARAM_LABEL[p]) for p in shown],
                 x_text=dyn.year.iloc[-1] + 0.8, min_gap=0.09, fontsize=6.5)
    ax.text(0.4, dyn.total_budget.iloc[0] - 0.07, PARAM_LABEL["total_budget"], fontsize=6.5, va="top",
            color=PCOL["total_budget"])
    ax.set_xlim(0, dyn.year.iloc[-1] + 9); ax.set_ylim(0, 1.02)
    ax.set_xticks(range(0, int(dyn.year.iloc[-1]) + 1, 5))
    ax.set_xlabel("Evaluation horizon (years)")
    ax.set_ylabel("Sobol total-order index")
    clean(ax, "y")
    claim(ax, "c", "Cycle time dominates early; prices and discounting catch up")

    ax = axd["d"]
    sc = pd.read_csv(OUT / "scenario_decomposition.csv").set_index("scen").loc[["pessimistic", "other", "optimistic"]]
    names = {"pessimistic": "Pessimistic", "other": "Mixed", "optimistic": "Optimistic"}
    cols = {"pessimistic": PCOL["discount_rate"], "other": TOK["ink2"], "optimistic": PCOL["cycle_compression"]}
    x_lo = sc.median_npv.min() / 1e6 / 2.5
    for k, (s_, r) in enumerate(sc.iterrows()):
        dumbbell(ax, k, r.median_npv / 1e6, r.mean_npv / 1e6, cols[s_])
        ax.text(x_lo * 1.15, k + 0.36, f"{100 * r.prob:.0f}% of states; P(loss) {100 * r.pr_neg:.1f}%",
                va="center", fontsize=6.5, color=TOK["ink2"])
    ax.set_xscale("log")
    ax.set_xlim(x_lo, sc.mean_npv.max() / 1e6 * 4)
    ax.xaxis.set_major_formatter(lambda v, _: money(v))
    ax.set_yticks(range(3), [names[s_] for s_ in sc.index]); ax.set_ylim(-0.6, 2.75)
    ax.set_xlabel("NPV, median (open) to mean (filled), log scale")
    clean(ax, "x")
    claim(ax, "d", "Negative NPV is rare; the upside is long")
    finish(fig, "fig4_value_follows_time",
           "Value follows cycle time. (a) NPV of the tool across all tool profiles and economic states by cycle "
           "compression (violins of positive values on a log scale; labels, bar and ring: median and interquartile "
           "range of all states; orange diamonds: evidence-based profile). (b) Median equivalent annual value per hectare by cycle compression and bean "
           f"price (discount rate {pct(BASE['discount_rate'])}, adoption ceiling {pct(BASE['adoption_ceiling'])}). "
           "(c) Sobol total-order indices by evaluation horizon; the four inputs with the largest index at the "
           "final year are labelled. (d) Median and mean NPV in the pessimistic and optimistic corners of the "
           "input space and in the remaining (mixed) states; P(loss): share of states with negative NPV.")


# ---------------------------------------------------------------------------
# Figure 5. What a pilot should measure
# ---------------------------------------------------------------------------
def figure5():
    fig, axd = new_figure(COL_W, 112, "ab;cc", width_ratios=[1, 1], height_ratios=[1, 0.9])
    sob = pd.read_csv(OUT / "sobol_indices_full.csv").set_index("param")
    order = sob.ST.sort_values().index

    ax = axd["a"]
    for k, p in enumerate(order):
        dumbbell(ax, k, max(sob.S1[p], 0), max(sob.ST[p], 0), PCOL[p])
    ax.set_yticks(range(len(order)), [KEY_LABEL.get(p, PARAM_LABEL[p]) for p in order], fontsize=6.5)
    ax.set_xlim(-0.02, sob.ST.max() * 1.12)
    ax.set_xlabel("Sobol index\n(open: first order; filled: total)")
    clean(ax, "x")
    claim(ax, "a", "Variance attribution")

    ax = axd["b"]
    ev = pd.read_csv(OUT / "evppi.csv")
    measures = {
        "Sobol": sob.ST,
        "Shap.\nind.": pd.read_csv(OUT / "shapley_indices_indep.csv").set_index("param").shapley_effect,
        "Shap.\ncop.": pd.read_csv(OUT / "shapley_indices_copula.csv").set_index("param").shapley_effect,
        "PAWN": pd.read_csv(OUT / "pawn_indices.csv").set_index("param").PAWN_median,
        "EVPPI": ev[ev.at_indifference].set_index("param").EVPPI,
    }
    minor = ["error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget"]
    major = [p for p in sob.index if p not in minor]
    for m, v in measures.items():            # the four minor inputs lie below every major input on every measure
        assert v[minor].max() < v[major].min(), m
    ranks = pd.DataFrame({m: v[major].rank(ascending=False, method="first") for m, v in measures.items()})
    xm = np.arange(len(measures))
    ax.axhspan(4.55, 8.5, color=TOK["neutral"], zorder=0, lw=0)
    for p, r in ranks.iterrows():
        strong = p == "cycle_compression"
        ax.plot(xm, r.values, color=PCOL[p], lw=LW * (1.5 if strong else 0.9), zorder=3)
        point(ax, xm, r.values, PCOL[p], size=MS)
    for k, p in enumerate(minor):
        point(ax, xm + (k - 1.5) * 0.11, np.full(len(xm), 6.5), PCOL[p], size=MS * 0.8)
    ax.set_xticks(xm, list(measures), fontsize=6.5)
    ax.set_yticks([1, 2, 3, 4, 6.5], ["1", "2", "3", "4", "5 to 8,\nnot resolved"], fontsize=6.5)
    ax.set_ylim(8.5, 0.5); ax.set_xlim(-0.3, xm[-1] + 0.3)
    ax.set_ylabel("Importance rank")
    clean(ax, "y", left=False)
    claim(ax, "b", "Top four ranks agree across measures")

    ax = axd["c"]
    if "hurdle_mult" not in ev:
        ev["hurdle_mult"] = ev.threshold / ev.loc[ev.at_indifference, "threshold"].iloc[0]
    ev["M"] = ev.EVPPI.clip(lower=0) / 1e6
    top = ev[ev.at_indifference].sort_values("EVPPI", ascending=False).param.iloc[:3].tolist()
    for p, d in ev.sort_values("hurdle_mult").groupby("param"):
        ax.plot(d.hurdle_mult, d.M, color=PCOL[p], lw=LW if p in top else 0.6, alpha=1 if p in top else 0.6)
    items = []
    for p in top:
        pk = ev[(ev.param == p) & ev.at_indifference].iloc[0]
        point(ax, pk.hurdle_mult, pk.M, PCOL[p])
        items.append((pk.hurdle_mult, pk.M, f"{PARAM_LABEL[p]}  {money(pk.M)}"))
    ymax = ev.M.max()
    label_column(ax, items, x_text=1 + 0.35 * (ev.hurdle_mult.max() - 1), min_gap=ymax * 0.12, fontsize=6.5)
    ax.axvline(1, color=TOK["axis"], lw=LW_HAIR)
    ax.text(1, ymax * 1.12, " hurdle = expected NPV", fontsize=6.5, color=TOK["muted"], va="top")
    ax.set_ylim(0, ymax * 1.15)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:g}x")
    ax.set_xlabel("Decision hurdle (multiple of expected NPV)")
    ax.set_ylabel("EVPPI (USD million)")
    clean(ax, "y")
    claim(ax, "c", "A pilot should measure cycle time first")
    finish(fig, "fig5_what_to_measure",
           "What a pilot should measure. (a) First-order (open) and total-order (filled) Sobol indices of NPV. "
           "(b) Importance ranks of the inputs under five measures: Sobol total-order index, Shapley effects with "
           "independent and copula-correlated inputs, PAWN median and EVPPI at the indifference hurdle; the four "
           "minor technical inputs (dots in the shaded band, colors as in panel a) lie below every other input on "
           "every measure, and their order is not resolved. (c) Expected value of perfect partial information "
           "of each input as the payoff of the alternative investment varies (multiple of expected NPV); labels "
           "give the three largest values at the indifference point.")


# ---------------------------------------------------------------------------
# Figure 6. Evidence behind the inputs
# ---------------------------------------------------------------------------
def pretty_study(s):
    pre = "_preprint" in s
    s = s.replace("_preprint", "")
    s = re.sub(r"([a-z]{3,})([A-Z])", r"\1-\2", s, count=1)
    s = re.sub(r"(\d{4})$", r" \1", s)
    return s + ("*" if pre else "")


def figure6():
    fig, axd = new_figure(COL_W, 120, "ab;cd", width_ratios=[1.15, 1], height_ratios=[1.1, 1])
    ex = pd.read_csv(LIT / "extract/meta_genetic_gain/extraction.csv")
    gain_col = TOK["ink2"]

    ax = axd["a"]
    m1 = mi[(mi.analysis_family == "MA1") & (mi.primary == "Y") & (mi.stratum3 == "common bean")
            & mi.sei_pct.notna()].merge(ex[["effect_id", "label"]], on="effect_id", how="left", suffixes=("", "_ex"))
    yr = np.where((m1.crossref_first_author == "Amongi") & (m1.crossref_year < 1990), 2023, m1.crossref_year)

    def short(a, y, l):
        l = str(l)
        if "PABRA" in l:
            m = re.search(r"(large|small)-seeded (bush|climbers?)", l)
            l = f"PABRA {m.group(1)} {m.group(2).rstrip('s')}" if m else "PABRA"
        else:
            l = re.sub(r"(\s*(beans?|cultivars|Brazil)\b)", "", l.split(" (")[0].split(",")[0]).strip()
        return f"{a} {int(y)}, {l}"
    m1["row"] = [short(a, y, l) for a, y, l in zip(m1.crossref_first_author, yr, m1.label_ex)]
    m1 = m1.sort_values("yi_pct").reset_index(drop=True)
    y = np.arange(len(m1)) + 1.4
    ax.hlines(y, m1.yi_pct - 1.96 * m1.sei_pct, m1.yi_pct + 1.96 * m1.sei_pct, color=gain_col, lw=0.6)
    point(ax, m1.yi_pct, y, gain_col, size=MS * 0.85)
    p1 = mr.loc["MA1.CB.primary"]
    ax.add_patch(Rectangle((p1.ci_lb, 0.15), p1.ci_ub - p1.ci_lb, 0.5, color=gain_col, alpha=0.18, lw=0))
    point(ax, p1.estimate, 0.4, TOK["ink"], marker="D", size=MS)
    ax.text(p1.ci_ub + 0.12, 0.4, f"{p1.estimate:.2f}%", va="center", fontsize=6.5)
    ax.axvline(0, color=TOK["axis"], lw=LW_HAIR)
    ax.set_yticks(np.r_[0.4, y], ["Pooled"] + list(m1.row), fontsize=6.5)
    ax.set_ylim(-0.3, len(m1) + 2)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:g}%")
    ax.set_xlabel("Realized yield gain per year")
    clean(ax, "x")
    claim(ax, "a", "Bean breeding gains ~0.5%/yr")

    ax = axd["b"]
    m2 = mi[(mi.analysis_family == "MA2") & (mi.analysis_set == "primary") & mi.r.notna()].copy()
    order = m2.groupby("trait_class").r.median().sort_values().index
    for k, tc in enumerate(order):
        gv = m2.loc[m2.trait_class == tc, "r"].values
        raincloud(ax, [gv], [k], [PCOL["error_reduction"]], orient="h", width=0.42, point_size=6, min_cloud=10)
    p2 = mr.loc["MA2.primary"]
    ax.axvspan(p2.ci_lb, p2.ci_ub, color=PCOL["error_reduction"], alpha=0.1, lw=0)
    ax.axvline(p2.estimate, color=PCOL["error_reduction"], lw=LW)
    ax.text(p2.ci_lb - 0.02, len(order) - 0.35, f"pooled\nr = {p2.estimate:.2f}", ha="right", va="center",
            fontsize=6.5)
    ax.set_yticks(range(len(order)), [t.capitalize() for t in order], fontsize=6.5)
    ax.set_ylim(-0.6, len(order)); ax.set_xlim(0, 1.02)
    ax.set_xlabel("Correlation with manual (r)")
    clean(ax, "x")
    claim(ax, "b", "Agreement varies by trait")

    ax = axd["c"]
    m3 = mi[(mi.analysis_family == "MA3") & mi.pct_reduction.notna() & mi.unit.notna()
            & mi.evidence_class.isin(["A", "C"]) & (mi.pct_reduction > -100)].drop_duplicates("unit")
    tot = m3[m3.unit.str.contains("_total")]          # handheld time when image-analysis time is counted
    m3 = m3[~m3.unit.str.contains("_total")]
    plats = ["Laboratory imaging", "Vehicle, UAV, robot", "Handheld, field"]
    platform = np.select([m3.unit.str.contains("handheld|phone", case=False),
                          m3.unit.str.contains("tractor|boom|robot|UAS|UAV", case=False)], [plats[2], plats[1]],
                         plats[0])
    for k, pl in enumerate(plats):
        gv = m3.pct_reduction.values[platform == pl]
        raincloud(ax, [gv], [k], [PCOL["cost_reduction"]], orient="h", width=0.42, point_size=6, min_cloud=10)
    if len(tot):
        ax.plot(tot.pct_reduction.values, np.full(len(tot), 2 - 0.16), ls="none", marker="D", ms=MS * 0.8,
                mfc=TOK["surface"], mec=PCOL["cost_reduction"], mew=1.0, zorder=5)
        ax.annotate("incl. image-\nanalysis time", xy=(tot.pct_reduction.iloc[0], 2 - 0.16),
                    xytext=(tot.pct_reduction.iloc[0] - 2, 2 + 0.34), fontsize=6.5, ha="right", va="center",
                    color=TOK["ink2"], arrowprops=dict(arrowstyle="-", color=TOK["axis"], lw=LW_HAIR))
    pt = mr.loc["MA3.time.primary"]
    ax.axvline(pt.estimate, color=PCOL["cost_reduction"], lw=LW)
    ax.text(pt.estimate - 1, len(plats) - 0.3, f"pooled\n{pt.estimate:.0f}%", ha="right", va="center",
            fontsize=6.5)
    ax.set_yticks(range(len(plats)), plats, fontsize=6.5)
    ax.set_ylim(-0.6, len(plats)); ax.set_xlim(0, 103)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
    ax.set_xlabel("Time saved vs manual")
    clean(ax, "x")
    claim(ax, "c", "Handheld tools save least time")

    ax = axd["d"]
    XLO, XHI = -1.15, 1.0
    r_er = [("T1: reference\nerror free", mr.loc["MA2.ER.T1"]), ("T2: reliability\nas entered", mr.loc["MA2.ER.T2"]),
            ("T2: validity\nconverted", mr.loc["MA2.ER.T2v"]), ("Heritability\npairs", mr.loc["MA2.ER.h2.median"])]
    ypos = {"model_err": 0.0, "model_cost": 6.2, "cost_ev": 7.2}
    ax.axvline(0, color=TOK["axis"], lw=LW_HAIR)
    ax.plot((l1.error_reduction.min(), l1.error_reduction.max()), [ypos["model_err"]] * 2, color=PCOL["error_reduction"],
            lw=5, alpha=0.25, solid_capstyle="butt")
    point(ax, META_ERR, ypos["model_err"], PCOL["error_reduction"])
    labels, ticks = ["Model range\n(error)"], [ypos["model_err"]]
    def interval(y, lo, hi, est, col, marker="o"):
        clo, chi = max(lo, XLO), min(hi, XHI)
        ax.plot([clo, chi], [y, y], color=col, lw=1.2)
        if lo < XLO:
            ax.plot(XLO, y, marker="<", ms=3.5, color=col, clip_on=False)
        for v in (lo, hi):
            if XLO <= v <= XHI:
                ax.plot([v, v], [y - 0.16, y + 0.16], color=col, lw=0.9)
        point(ax, est, y, col, marker=marker, size=MS * 0.9)
    for k, (nm, r) in enumerate(r_er):
        y = 1.2 + k
        interval(y, r.ci_lb, r.ci_ub, r.estimate, PCOL["error_reduction"])
        labels.append(nm); ticks.append(y)
    ax.plot((l1.cost_reduction.min(), l1.cost_reduction.max()), [ypos["model_cost"]] * 2, color=PCOL["cost_reduction"],
            lw=5, alpha=0.25, solid_capstyle="butt")
    point(ax, META_COST, ypos["model_cost"], PCOL["cost_reduction"])
    rc = mr.loc["MA3.cost.implied_handheld"]
    interval(ypos["cost_ev"], rc.ci_lb, rc.ci_ub, rc.estimate, PCOL["cost_reduction"])
    labels += ["Model range\n(cost)", "Time saving x\nnon-handling share"]; ticks += [ypos["model_cost"], ypos["cost_ev"]]
    ax.set_yticks(ticks, labels, fontsize=6.5)
    ax.set_ylim(8.0, -0.8); ax.set_xlim(XLO, XHI)
    ax.xaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
    ax.set_xlabel("Reduction vs manual")
    clean(ax, "x")
    claim(ax, "d", "Accuracy evidence spans zero")
    finish(fig, "fig6_evidence",
           "Evidence behind the model inputs. (a) Realized yield gain of common-bean breeding programs "
           f"(k = {int(p1.k_effects)} estimates, {int(p1.k_clusters)} programs, I2 = {p1.I2_pct:.0f}%; diamond "
           "and band: three-level REML estimate and 95% CI). (b) Agreement between image-based or AI tools and "
           f"manual reference measurements by trait (k = {int(p2.k_effects)} effects; line and band: pooled "
           "estimate and 95% CI). (c) Time saved by digital phenotyping relative to manual recording, by "
           "platform (line: pooled estimate; two descriptive units for handheld tools). (d) Ranges of "
           "measurement-error and whole-plot cost reduction explored by the model (bars; dots: central case) "
           "against the evidence (points: median or estimate; lines: 95% interval, or range for heritability "
           "pairs; arrows: interval continues beyond -115%).")


failed = []
for f in (figure1, figure2, figure3, figure4, figure5, figure6):
    try:
        f()
    except Exception as e:
        traceback.print_exc()
        failed.append(f"{f.__name__}: {e}")

with open("paper/figures/captions.md", "w") as fh:
    for k, v in CAPTIONS.items():
        fh.write(f"## {k}\n\n{v}\n\n")

n_issues = 0
for name, issues in ISSUES.items():
    print(f"{name}: {'ok' if not issues else f'{len(issues)} layout issue(s)'}")
    for i in issues:
        print("   ", i)
    n_issues += len(issues)
if failed or n_issues:
    print("\n".join(failed))
    sys.exit(1)
