"""Shared style for the manuscript figures: colors, fonts, sizes and helpers.

Colors are color-vision-deficiency safe on a white background. Figures are
drawn at their printed width (138.6 or 184.7 mm) with 7 pt text and bold lower-case panel letters.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

MM = 1 / 25.4
COL_W = 138.6 * MM     # MDPI text-column width (13.86 cm), the printed width of \linewidth
FULL_W = 184.7 * MM    # MDPI full width via adjustwidth (18.47 cm), the printed width

TOK = dict(
    surface="#ffffff",
    plane="#f9f9f7",     # page plane: schematic boxes, bands
    ink="#0b0b0b",
    ink2="#52514e",
    muted="#898781",
    grid="#e1e0d9",
    axis="#c3c2b7",      # baseline / axis rule; also de-emphasis marks
    neutral="#f0efec",   # diverging midpoint / neutral band
)
CAT = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100")
ACCENT = CAT["blue"]
DEEMPH = TOK["axis"]
SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
ORD3 = ["#86b6ef", "#3987e5", "#1c5cab"]
ORD4 = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
CMAP_BLUE = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)

LW = 1.0          # data lines (pt) ~ 2 px at print scale
LW_HAIR = 0.4     # grid / rules
MS = 4.2          # marker size (pt)

# One fixed color per model input, used in every figure. Order follows the
# categorical palette so neighbouring inputs stay distinguishable.
PCOL = {
    "cycle_compression": "#2a78d6",   # blue
    "error_reduction":   "#eb6834",   # orange
    "cost_reduction":    "#1baf7a",   # aqua
    "adoption_ceiling":  "#eda100",   # yellow
    "bean_price":        "#e87ba4",   # magenta
    "tool_fixed_cost":    "#008300",   # green
    "discount_rate":     "#4a3aa7",   # violet
    "total_budget":      "#e34948",   # red
}
TECH = ["cycle_compression", "error_reduction", "cost_reduction", "tool_fixed_cost", "total_budget"]
ECON = ["bean_price", "discount_rate", "adoption_ceiling"]

PARAM_LABEL = {
    "cycle_compression": "Cycle compression",
    "discount_rate": "Discount rate",
    "bean_price": "Bean price",
    "adoption_ceiling": "Adoption ceiling",
    "error_reduction": "Measurement-error reduction",
    "cost_reduction": "Whole-plot cost reduction",
    "tool_fixed_cost": "Tool fixed cost",
    "total_budget": "Program budget",
}


def setup():
    mpl.rcParams.update({
        "font.family": "Lato",
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "text.color": TOK["ink"],
        "axes.labelcolor": TOK["ink2"],
        "xtick.color": TOK["ink2"],
        "ytick.color": TOK["ink2"],
        "axes.edgecolor": TOK["axis"],
        "axes.linewidth": LW_HAIR,
        "axes.facecolor": TOK["surface"],
        "figure.facecolor": TOK["surface"],
        "axes.grid": True,
        "grid.color": TOK["grid"],
        "grid.linewidth": LW_HAIR,
        "grid.linestyle": "-",
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "xtick.major.pad": 3,
        "ytick.major.pad": 3,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.2,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
    })


def panel_label(ax, letter, dx=-0.02, dy=1.02):
    """Bold lower-case panel letter at the top-left of the axes' bounding box."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="right", va="bottom", color=TOK["ink"])


def clean(ax, grid="y", left=True):
    ax.grid(False)
    if grid in ("x", "both"):
        ax.xaxis.grid(True)
    if grid in ("y", "both"):
        ax.yaxis.grid(True)
    ax.spines["left"].set_visible(left)


def point(ax, x, y, color, size=MS, zorder=4, marker="o", **kw):
    """Filled marker with a white surface ring."""
    return ax.plot(x, y, marker=marker, ls="none", ms=size, mfc=color, mec=TOK["surface"],
                   mew=0.8, zorder=zorder, **kw)


def label_column(ax, items, x_text, min_gap, x_anchor=None, fontsize=6.5, ha="left"):
    """Place labels in a column at x_text, spread vertically by at least
    min_gap (data units), with hairline leaders back to their anchors.
    items: list of (anchor_x, anchor_y, text)."""
    items = sorted(items, key=lambda t: t[1])
    ys = [it[1] for it in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + min_gap)
    for (ax_x, ax_y, text), y in zip(items, ys):
        ax.annotate(text, xy=(ax_x, ax_y), xytext=(x_text, y), fontsize=fontsize, ha=ha,
                    va="center", color=TOK["ink"],
                    arrowprops=dict(arrowstyle="-", color=TOK["axis"], lw=LW_HAIR,
                                    shrinkA=1, shrinkB=2))


def raincloud(ax, groups, positions, colors, width=0.38, jitter=0.07, point_size=2.2, highlight=None,
              highlight_color=None, orient="v", seed=0, max_points=400, min_cloud=6):
    """Half violin ('cloud') + jittered points ('rain') + slim box (median, IQR).
    groups: list of 1-D arrays; positions: category positions; colors: one per group.
    highlight: optional list of boolean masks (same length as each group) drawn on top."""
    from scipy.stats import gaussian_kde
    import numpy as np
    rng = np.random.default_rng(seed)
    for k, (v, pos, col) in enumerate(zip(groups, positions, colors)):
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v)]
        if len(v) < 3:                      # too few for a box: show the points only
            if orient == "v":
                point(ax, np.full(len(v), pos), v, col, size=MS * 0.85)
            else:
                point(ax, v, np.full(len(v), pos), col, size=MS * 0.85)
            continue
        lo, hi = np.percentile(v, [0.5, 99.5])
        pad = 0.05 * (hi - lo + 1e-9)
        grid = np.linspace(lo - pad, hi + pad, 200)
        dens = gaussian_kde(v)(grid) if len(v) >= min_cloud else np.zeros_like(grid)
        dens = width * dens / max(dens.max(), 1e-12)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        idx = rng.choice(len(v), size=min(len(v), max_points), replace=False)
        jit = rng.uniform(-jitter, jitter, len(idx))
        hl = None if highlight is None else np.asarray(highlight[k])[np.isfinite(np.asarray(groups[k], float))]
        if orient == "v":
            ax.fill_betweenx(grid, pos, pos + dens, color=col, alpha=0.28, lw=0)
            ax.plot(pos + dens, grid, color=col, lw=0.7)
            ax.scatter(pos - 0.16 + jit, v[idx], s=point_size, color=col, alpha=0.35, lw=0, zorder=2)
            if hl is not None and hl.any():
                ax.scatter(pos - 0.16 + rng.uniform(-jitter, jitter, hl.sum()), v[hl], s=point_size * 2.2,
                           color=highlight_color, lw=0.3, edgecolor=TOK["surface"], zorder=3)
            ax.plot([pos, pos], [q1, q3], color=TOK["ink"], lw=2.2, solid_capstyle="butt", zorder=4)
            ax.plot(pos, med, marker="o", ms=3.2, mfc=TOK["surface"], mec=TOK["ink"], mew=0.8, zorder=5)
        else:
            ax.fill_between(grid, pos, pos + dens, color=col, alpha=0.28, lw=0)
            ax.plot(grid, pos + dens, color=col, lw=0.7)
            ax.scatter(v[idx], pos - 0.16 + jit, s=point_size, color=col, alpha=0.35, lw=0, zorder=2)
            if hl is not None and hl.any():
                ax.scatter(v[hl], pos - 0.16 + rng.uniform(-jitter, jitter, hl.sum()), s=point_size * 2.2,
                           color=highlight_color, lw=0.3, edgecolor=TOK["surface"], zorder=3)
            ax.plot([q1, q3], [pos, pos], color=TOK["ink"], lw=2.2, solid_capstyle="butt", zorder=4)
            ax.plot(med, pos, marker="o", ms=3.2, mfc=TOK["surface"], mec=TOK["ink"], mew=0.8, zorder=5)


def dumbbell(ax, y, x0, x1, color, size=MS, open_start=True, lw=None):
    """Dumbbell from x0 (open marker) to x1 (filled marker) at height y."""
    import numpy as np
    for yi, a, b in zip(np.atleast_1d(y), np.atleast_1d(x0), np.atleast_1d(x1)):
        ax.plot([a, b], [yi, yi], color=color, lw=lw or LW * 1.6, alpha=0.45, solid_capstyle="round", zorder=2)
        if open_start:
            ax.plot(a, yi, marker="o", ms=size, mfc=TOK["surface"], mec=color, mew=1.0, zorder=3)
        else:
            point(ax, a, yi, color, size=size, zorder=3)
        point(ax, b, yi, color, size=size, zorder=4)


KEY_LABEL = {"error_reduction": "Error reduction", "cost_reduction": "Cost reduction"}


def param_key(ax, params, y, x0, fontsize=6.0, char_w=0.0092, gap=0.035):
    """Horizontal color key for model inputs (dot + name), spaced by label length (axes coords)."""
    xs = x0
    for p in params:
        ax.plot(xs, y, marker="o", ms=4.2, mfc=PCOL[p], mec=TOK["surface"], mew=0.6, transform=ax.transAxes,
                clip_on=False)
        lab = KEY_LABEL.get(p, PARAM_LABEL[p])
        ax.text(xs + 0.01, y, lab, fontsize=fontsize, va="center", color=TOK["ink2"], transform=ax.transAxes)
        xs += 0.01 + char_w * len(lab) + gap


def new_figure(width, height_mm, mosaic, **kw):
    """Figure with constrained layout and a named-panel mosaic (e.g. "ab;cc")."""
    fig = plt.figure(figsize=(width, height_mm * MM), layout="constrained")
    fig.get_layout_engine().set(w_pad=1.5 / 72, h_pad=1.5 / 72, wspace=0.03, hspace=0.04)
    axd = fig.subplot_mosaic(mosaic, **kw)
    fig._panel_titles = {}
    return fig, axd


def claim(ax, letter, text):
    """Panel title: the claim the panel supports. The bold letter is placed after layout."""
    ax.set_title(text, loc="left", fontsize=7, fontweight="bold", color=TOK["ink"], pad=5)
    ax.figure._panel_titles[ax] = letter


def _place_letters(fig, renderer):
    """Bold panel letter at the outer-left edge of each panel, level with its title. Where the
    panel has no y-axis labels (shared axis) there is no room, so the letter leads the title."""
    inv = fig.transFigure.inverted()
    for ax, letter in fig._panel_titles.items():
        bb = ax.get_tightbbox(renderer)
        t = ax._left_title.get_window_extent(renderer)
        x0 = inv.transform((bb.x0, 0))[0]
        y = inv.transform((0, t.y0 + t.height / 2))[1]
        lt = fig.text(x0, y, letter, fontsize=9, fontweight="bold", ha="left", va="center", color=TOK["ink"])
        if lt.get_window_extent(renderer).x1 > t.x0 - 3:
            lt.remove()
            ax.set_title(f"{letter}   {ax.get_title(loc='left')}", loc="left", fontsize=7, fontweight="bold",
                         color=TOK["ink"], pad=5)


def _texts(fig):
    """Text artists that are actually drawn (tick labels outside the view limits are skipped)."""
    out = list(fig.texts)
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        out += [ax.title, ax._left_title, ax._right_title] + list(ax.texts)
        if not ax.axison:
            continue
        out += [ax.xaxis.label, ax.yaxis.label]
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis._update_ticks():
                if tick.label1.get_visible():
                    out.append(tick.label1)
        leg = ax.get_legend()
        if leg is not None:
            out += list(leg.get_texts())
    return [t for t in out if t.get_visible() and t.get_text().strip()]


def check_layout(fig, renderer, tol=0.5):
    """Report text that leaves the canvas or overlaps other text."""
    W, H = fig.bbox.width, fig.bbox.height
    boxes = []
    for t in _texts(fig):
        try:
            # an Annotation's extent includes its leader line; only the text itself can collide with text
            bb = (mpl.text.Text.get_window_extent(t, renderer) if isinstance(t, mpl.text.Annotation)
                  else t.get_window_extent(renderer))
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        boxes.append((t, bb))
    issues = []
    for t, bb in boxes:
        if bb.x0 < -tol or bb.y0 < -tol or bb.x1 > W + tol or bb.y1 > H + tol:
            issues.append(f"outside canvas: {t.get_text()!r}")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > 1.0 and oy > 1.0:
                issues.append(f"overlap: {boxes[i][0].get_text()!r} / {boxes[j][0].get_text()!r}")
    return issues


def save(fig, name, out="paper/figures"):
    Path(out).mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    if getattr(fig, "_panel_titles", None):
        _place_letters(fig, renderer)
        fig.canvas.draw()
    issues = check_layout(fig, renderer)
    fig.savefig(Path(out) / f"{name}.pdf")
    fig.savefig(Path(out) / f"{name}.png", dpi=600)
    plt.close(fig)
    return issues
