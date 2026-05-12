"""
Aim 1 plotting script (matplotlib only; no seaborn)

Inputs (in current working dir or provide full paths):
  - Aim1_summary_metrics.csv
  - Aim1_pocket_residue_deviation.csv

Outputs (saved into ./analysis_out/ by default):
  1) Aim1_Fig1_Global_vs_Pocket_RMSD.png (+ .pdf)
  2) Aim1_Fig2_PocketResidue_dCA_Heatmap.png (+ .pdf)
  3) Aim1_Fig2_PocketResidue_dCA_Heatmap_TOPN.png (+ .pdf)  [optional]
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# User settings
# -------------------------
SUMMARY_CSV = "Aim1_summary_metrics.csv"
RESDEV_CSV = "Aim1_pocket_residue_deviation.csv"
OUT_DIR = "analysis_out"

# Optional: make a second heatmap showing only the TOP-N most variable pocket residues (by mean dCA)
TOPN_RESIDUES = 25  # set to None to disable

# Structure ordering (optional). If None, uses alphabetical order from file.
STRUCTURE_ORDER = None  # e.g. ["2DN3", "6BB5", "AncA", "AncB", "AncMH"]

# -------------------------
os.makedirs(OUT_DIR, exist_ok=True)


def _read_inputs():
    df_sum = pd.read_csv(SUMMARY_CSV)
    df_res = pd.read_csv(RESDEV_CSV)
    return df_sum, df_res


def fig1_barplot_global_vs_pocket(df_sum: pd.DataFrame):
    # Ensure numeric
    for col in ["global_RMSD_CA", "pocket_RMSD_CA"]:
        df_sum[col] = pd.to_numeric(df_sum[col], errors="coerce")

    # Order
    if STRUCTURE_ORDER is not None:
        df_sum["structure"] = pd.Categorical(df_sum["structure"], categories=STRUCTURE_ORDER, ordered=True)
        df_sum = df_sum.sort_values("structure")
    else:
        df_sum = df_sum.sort_values("structure")

    x = np.arange(len(df_sum))
    w = 0.38

    fig = plt.figure(figsize=(10, 4.6))
    ax = plt.gca()

    ax.bar(x - w/2, df_sum["global_RMSD_CA"].values, width=w, label="Global RMSD (Cα)")
    ax.bar(x + w/2, df_sum["pocket_RMSD_CA"].values, width=w, label="Pocket RMSD (Cα)")

    ax.set_xticks(x)
    ax.set_xticklabels(df_sum["structure"].tolist(), rotation=35, ha="right")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Global vs Heme-pocket backbone deviation (β-chain)")

    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.25)

    # annotate values (optional but useful)
    for i, (g, p) in enumerate(zip(df_sum["global_RMSD_CA"].values, df_sum["pocket_RMSD_CA"].values)):
        if np.isfinite(g):
            ax.text(i - w/2, g + 0.02, f"{g:.2f}", ha="center", va="bottom", fontsize=8)
        if np.isfinite(p):
            ax.text(i + w/2, p + 0.02, f"{p:.2f}", ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT_DIR, "Aim1_Fig1_Global_vs_Pocket_RMSD.png")
    out_pdf = os.path.join(OUT_DIR, "Aim1_Fig1_Global_vs_Pocket_RMSD.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)
    print("Saved:", out_png)
    print("Saved:", out_pdf)


def _make_residue_labels(df_res: pd.DataFrame) -> pd.Series:
    # Use reference index + reference residue name to make stable labels
    # Example: "LEU_28" (index is 0-based in the script output; that is fine as long as consistent)
    # If you prefer 1-based, add +1 below.
    ref_idx = pd.to_numeric(df_res["ref_res_index"], errors="coerce")
    ref_name = df_res["ref_resname"].astype(str)
    labels = ref_name + "_" + (ref_idx.astype("Int64") + 1).astype(str)  # +1 => human-friendly indexing
    return labels


def fig2_heatmap_residue_dca(df_res: pd.DataFrame, topn: int | None = None):
    df = df_res.copy()
    df["dCA_A"] = pd.to_numeric(df["dCA_A"], errors="coerce")

    # Create residue label (ref-based) and pivot
    df["res_label"] = _make_residue_labels(df)

    # Optional: pick top-n residues by mean deviation across structures
    if topn is not None:
        means = df.groupby("res_label")["dCA_A"].mean().sort_values(ascending=False)
        keep = set(means.head(topn).index)
        df = df[df["res_label"].isin(keep)]

    # Order structures
    if STRUCTURE_ORDER is not None:
        df["structure"] = pd.Categorical(df["structure"], categories=STRUCTURE_ORDER, ordered=True)

    # Pivot: rows=structures, cols=residues, values=dCA
    mat = df.pivot_table(index="structure", columns="res_label", values="dCA_A", aggfunc="mean")

    # Sort columns by reference index (extracted from label suffix)
    def _col_key(label: str) -> int:
        # label like "LEU_28"
        try:
            return int(label.split("_")[-1])
        except Exception:
            return 10**9

    mat = mat.reindex(sorted(mat.columns, key=_col_key), axis=1)

    # Sort rows if not categorical
    if STRUCTURE_ORDER is None:
        mat = mat.sort_index()

    data = mat.to_numpy(dtype=float)

    fig_w = max(10, 0.28 * data.shape[1])
    fig_h = max(4.5, 0.42 * data.shape[0])
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = plt.gca()

    im = ax.imshow(data, aspect="auto", interpolation="nearest")

    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index.tolist())

    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns.tolist(), rotation=90, fontsize=7)

    ax.set_xlabel("Pocket residues (reference identity)")
    ax.set_ylabel("Structure")
    title = "Pocket residue backbone deviation (ΔCα, Å)"
    if topn is not None:
        title += f" (TOP {topn} most variable residues)"
    ax.set_title(title)

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("ΔCα (Å)")

    fig.tight_layout()

    suffix = f"_TOP{topn}" if topn is not None else ""
    out_png = os.path.join(OUT_DIR, f"Aim1_Fig2_PocketResidue_dCA_Heatmap{suffix}.png")
    out_pdf = os.path.join(OUT_DIR, f"Aim1_Fig2_PocketResidue_dCA_Heatmap{suffix}.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)
    print("Saved:", out_png)
    print("Saved:", out_pdf)


def main():
    df_sum, df_res = _read_inputs()

    # Figure 1
    fig1_barplot_global_vs_pocket(df_sum)

    # Figure 2 (full pocket heatmap)
    fig2_heatmap_residue_dca(df_res, topn=None)

    # Optional: top-N heatmap (often more readable)
    if TOPN_RESIDUES is not None and isinstance(TOPN_RESIDUES, int) and TOPN_RESIDUES > 0:
        fig2_heatmap_residue_dca(df_res, topn=TOPN_RESIDUES)


if __name__ == "__main__":
    main()