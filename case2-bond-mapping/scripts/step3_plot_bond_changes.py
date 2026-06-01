#!/usr/bin/env python3
"""
Step 3 - Plot the intramolecular bond-order comparison.

Grouped bars show each intramolecular bond's DDEC6 bond order in the gas phase
and on the two surfaces, so weakening (activation) or strengthening on Cu(111)
vs Ni(111) is read directly against the gas-phase reference.

Input : results/intramolecular_comparison.csv
Output: results/bond_order_comparison.png
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FONT_FAMILY = "Times New Roman"
LABEL_FONTSIZE = 26
TICKLABEL_FONTSIZE = 18
LEGEND_FONTSIZE = 20
FIGURE_DPI = 600
COLORS = {"gas": "#A9A9A9", "Cu": "#FF8C00", "Ni": "#0096FF"}


def main():
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = [FONT_FAMILY, "DejaVu Serif", "serif"]

    df = pd.read_csv(ROOT / "results" / "intramolecular_comparison.csv")
    labels = df["bond"].tolist()
    x = np.arange(len(labels))
    w = 0.27

    fig, ax = plt.subplots(figsize=(max(8, 1.6 * len(labels)), 7), dpi=FIGURE_DPI)
    ax.bar(x - w, df["bo_gas"], w, label="gas", color=COLORS["gas"], edgecolor="black")
    ax.bar(x,     df["bo_Cu"],  w, label="Cu(111)", color=COLORS["Cu"], edgecolor="black")
    ax.bar(x + w, df["bo_Ni"],  w, label="Ni(111)", color=COLORS["Ni"], edgecolor="black")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=TICKLABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel("DDEC6 bond order", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.tick_params(axis="y", labelsize=TICKLABEL_FONTSIZE)
    for lab in ax.get_yticklabels():
        lab.set_fontweight("bold")
    ax.legend(fontsize=LEGEND_FONTSIZE, framealpha=0.0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", ls=":", alpha=0.4)

    out = ROOT / "results" / "bond_order_comparison.png"
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
