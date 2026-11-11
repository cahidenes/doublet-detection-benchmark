"""Shared names, colors and data loading for the analysis scripts."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results/figures"
STATS = ROOT / "results/stats"

NAMES = {"libsize": "Library size", "ngenes": "# genes", "scrublet": "Scrublet",
         "doubletdetection": "DoubletDetection", "solo": "Solo", "vaeda": "vaeda",
         "scdblfinder": "scDblFinder", "dbldensity": "DoubletDensity", "cxds": "cxds",
         "bcds": "bcds", "hybrid": "hybrid", "doubletfinder": "DoubletFinder", "chord": "Chord"}
BASELINES = ["libsize", "ngenes"]
# labeled doublet technology (README, Data)
TECH = {"pbmc-ch": "Cell hashing", "cline-ch": "Cell hashing", "mkidney-ch": "Cell hashing",
        "hm-12k": "Species mix", "hm-6k": "Species mix",
        **{d: "demuxlet" for d in ["pbmc-1A-dm", "pbmc-1B-dm", "pbmc-1C-dm", "pbmc-2ctrl-dm",
                                   "pbmc-2stim-dm", "J293t-dm"]},
        **{d: "MULTI-seq" for d in ["pdx-MULTI", "HMEC-orig-MULTI", "HMEC-rep-MULTI",
                                    "HEK-HMEC-MULTI", "nuc-MULTI"]}}

# reference palette (dataviz skill, validated): categorical slots 1-4, blue ramp, blue<->red
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list(
    "blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("blue_red", ["#184f95", "#6da7ec", "#f0efec",
                                                     "#ec8f8e", "#e34948"])


def setup():
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "figure.facecolor": SURF, "axes.facecolor": SURF,
                         "savefig.facecolor": SURF})


def clean(ax, grid="y"):
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def main_metrics():
    """Main run: every tool at the labeled rate (1x), seed 0."""
    d = pd.read_csv(ROOT / "results/metrics.tsv", sep="\t")
    return d[(d.rate_mult == 1.0) & (d.status == "ok")].copy()


def text_on(v, vmin, vmax):
    """Ink that stays readable on a sequential cell."""
    return "white" if (v - vmin) / (vmax - vmin + 1e-12) > 0.55 else INK
