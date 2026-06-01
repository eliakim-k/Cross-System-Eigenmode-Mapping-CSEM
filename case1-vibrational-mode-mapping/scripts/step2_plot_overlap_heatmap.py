#!/usr/bin/env python3
"""
Step 2 - Plot the mode-overlap heatmap.

Reads the overlap matrix from step 1 and draws it as a heatmap: rows are the
reference (Cu) modes, columns the moved (Ni) modes, color is the absolute
eigenvector overlap in [0, 1]. Each cell is annotated with its overlap value,
and the greedy one-to-one matches are outlined, so the matched pairing reads as
a near-diagonal track of boxed cells.

Input : results/overlap_matrix.csv, results/matches.csv
Output: results/overlap_heatmap.png
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = Path(__file__).resolve().parent.parent
REF, MOV = "Cu", "Ni"

FONT_FAMILY = "Times New Roman"
LABEL_FONTSIZE = 26
TICKLABEL_FONTSIZE = 16
ANNOT_FONTSIZE = 12
TITLE_FONTSIZE = 26
FIGURE_DPI = 600
CMAP = "viridis"


def main():
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = [FONT_FAMILY, "DejaVu Serif", "serif"]

    S = pd.read_csv(ROOT / "results" / "overlap_matrix.csv", index_col=0)
    M = pd.read_csv(ROOT / "results" / "matches.csv")
    A = S.to_numpy(float)
    nref, nmov = A.shape

    fig, ax = plt.subplots(figsize=(10, 9), dpi=FIGURE_DPI)
    im = ax.imshow(A, cmap=CMAP, vmin=0.0, vmax=1.0, aspect="equal", origin="upper")

    # Annotate every cell with its overlap value (black on bright, white on dark).
    for i in range(nref):
        for j in range(nmov):
            v = A[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=ANNOT_FONTSIZE, color="black" if v > 0.5 else "white")

    # Outline each matched (reference -> moved) cell.
    for _, r in M.iterrows():
        i, j = int(r[f"{REF}_mode"]), int(r[f"{MOV}_mode"])
        ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                        fill=False, edgecolor="red", linewidth=2.0))

    ax.set_xticks(range(nmov)); ax.set_yticks(range(nref))
    ax.set_xticklabels(range(nmov), fontsize=TICKLABEL_FONTSIZE)
    ax.set_yticklabels(range(nref), fontsize=TICKLABEL_FONTSIZE)
    ax.set_xlabel(f"{MOV}(111) mode index", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel(f"{REF}(111) mode index", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.set_title("Eigenmode overlap  |v$_{Cu}$ . v$_{Ni}$|",
                 fontsize=TITLE_FONTSIZE, fontweight="bold")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("overlap", fontsize=LABEL_FONTSIZE - 4, fontweight="bold")
    cb.ax.tick_params(labelsize=TICKLABEL_FONTSIZE)

    out = ROOT / "results" / "overlap_heatmap.png"
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
