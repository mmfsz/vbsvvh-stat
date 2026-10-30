"""Expected-limit plots in CMS style (matplotlib + mplhep).

Two figures from ``limits.json``:

``limits_<signal>``  one row per fit (each run, then the combination);
``limits_<fit>``     one row per signal hypothesis.
"""

import re

from .config import RUNS

#: CMS band colours (ROOT kGreen+1, kOrange), band half-height, figure size and header font size.
BAND_68, BAND_95 = "#00cc00", "#ffcc00"
HALF_HEIGHT = 0.35
WIDTH, HEIGHT_PER_ROW, HEIGHT_EXTRA = 12, 1.4, 3.5
HEADER_SIZE = 22


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    return plt, hep


def plot_rows(points, labels, outstem, lumi_title, divider_before=None, title=None, cms_label="Preliminary"):
    """``points``: ordered ``{key: quantiles}``, drawn top to bottom; ``labels``: key -> axis label.

    ``title`` is written inside the frame, above the rows (the signal hypothesis).

    ``divider_before`` is the index of the first row that combines the ones
    above it, marked by a dashed line since it is not another measurement.
    """
    plt, hep = _pyplot()
    keys = list(points)
    n = len(keys)
    fig, ax = plt.subplots(figsize=(WIDTH, HEIGHT_EXTRA + HEIGHT_PER_ROW * n))

    for i, key in enumerate(keys):
        q, y = points[key], n - 1 - i
        ax.fill_betweenx([y - HALF_HEIGHT, y + HALF_HEIGHT], q["exp-2"], q["exp+2"], color=BAND_95, lw=0)
        ax.fill_betweenx([y - HALF_HEIGHT, y + HALF_HEIGHT], q["exp-1"], q["exp+1"], color=BAND_68, lw=0)
        ax.plot([q["exp0"]] * 2, [y - HALF_HEIGHT, y + HALF_HEIGHT], color="black", ls="--", lw=2)
    if divider_before:
        ax.axhline(n - divider_before - 0.5, color="gray", ls="--", lw=1.5)

    lo = min(q["exp-2"] for q in points.values())
    hi = max(q["exp+2"] for q in points.values())
    ax.set_xscale("log")
    ax.set_xlim(lo / 3.0, hi * 3.0)
    ax.set_ylim(-0.6, n - 0.4 + (1.3 if title else 0.8))
    if title:
        ax.text(0.03, n - 0.2, title, transform=ax.get_yaxis_transform(), va="center")
    ax.set_yticks(range(n), [labels.get(k, k) for k in reversed(keys)])
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_xlabel(r"95% CL upper limit on $\mu$")

    handles = [
        plt.Line2D([], [], color="black", ls="--", lw=2, label="Median expected"),
        plt.Rectangle((0, 0), 1, 1, color=BAND_68, label="68% expected"),
        plt.Rectangle((0, 0), 1, 1, color=BAND_95, label="95% expected"),
    ]
    ax.legend(handles=handles, loc="upper center", ncol=3, fontsize="small", handlelength=1.5, columnspacing=1.0)
    hep.cms.label(cms_label, data=True, rlabel=lumi_title, loc=0, fontsize=HEADER_SIZE, ax=ax)
    outstem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{outstem}.{ext}", bbox_inches="tight")
    plt.close(fig)


ENERGY = {"Run2": "13 TeV", "Run3": "13.6 TeV"}


def lumi_title(eras, analysis):
    """``138 fb^-1 (13 TeV) + 221 fb^-1 (13.6 TeV)`` for a set of eras."""
    parts = []
    for run, energy in ENERGY.items():
        total = sum(analysis.luminosity.get(e, 0.0) for e in eras if e in RUNS[run])
        if total:
            parts.append(f"{total:.0f} fb$^{{-1}}$ ({energy})")
    return " + ".join(parts)


def run_label(run):
    """``Run2`` -> ``Run 2``; ``Run3_2024`` -> ``Run 3 (2024)``."""
    m = re.fullmatch(r"Run(\d)(?:_(.+))?", run)
    if not m:
        return run
    return f"Run {m[1]}" + (f" ({m[2]})" if m[2] else "")


def fit_label(runs):
    """A fit's runs joined: ``Run 2 + Run 3``."""
    return " + ".join(run_label(r) for r in runs)


def plot_limits(limits, fits, fit_eras, analysis, out_dir):
    """Every figure of a result from its ``limits.json`` content."""
    combined = [f for f in fits if len(fits[f]) > 1]
    for signal in sorted({s for per in limits.values() for s in per}):
        points = {f: limits[f][signal] for f in fits if signal in limits.get(f, {})}
        first_combined = next((i for i, f in enumerate(points) if f in combined), None)
        all_eras = {e for f in points for e in fit_eras[f]}
        labels = {f: fit_label(fits[f]) for f in points}
        plot_rows(points, labels, out_dir / f"limits_{signal}", lumi_title(all_eras, analysis), first_combined,
                  title=analysis.signals.get(signal, signal))
    for fit, per in limits.items():
        points = {s: per[s] for s in analysis.signals if s in per}
        plot_rows(points, analysis.signals, out_dir / f"limits_{fit}", lumi_title(fit_eras[fit], analysis))
