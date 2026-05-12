import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv("master_docking_summary.csv")

# Ensure consistent ordering
protein_order = [
    "1A3N", "2DN3", "6BB5",
    "mut_D99N", "mut_K82N", "mut_N102T", "mut_Y35F",
    "AncA", "AncB", "AncMH_nat"
]

ligand_order = ["heme_Co", "heme_Fe", "heme_O2", "heme", "2,3-BPG", "IHP"]

df["protein"] = pd.Categorical(df["protein"], protein_order)
df["ligand"] = pd.Categorical(df["ligand"], ligand_order)

# Pivot tables
energy_map = df.pivot(index="protein", columns="ligand", values="mean_energy")
cluster_map = df.pivot(index="protein", columns="ligand", values="members_fraction")

# -------- Heatmap A: Mean binding energy --------
font1 = {'family': 'sans-serif',
        'color': 'black',
        'size': 16,
        'weight': 'semibold',
        }

plt.figure(figsize=(10, 6))
sns.heatmap(
    energy_map,
    annot=True,
    fmt=".2f",
    cmap="viridis",
    cbar_kws={"label": "Mean ΔG (kcal/mol)"}
)
plt.title("Mean Binding Energy of Best Cluster", fontdict=font1)
plt.ylabel("Protein", fontdict=font1, labelpad=20)
plt.xlabel("Ligand", fontdict=font1, labelpad=20)
plt.tight_layout()
plt.show()

# -------- Heatmap B: Cluster reproducibility --------
plt.figure(figsize=(10, 6))
sns.heatmap(
    cluster_map,
    annot=True,
    fmt=".2f",
    cmap="magma",
    cbar_kws={"label": "Cluster Fraction (members / 100)"}
)
plt.title("Reproducibility of Best Binding Mode", fontdict=font1)
plt.ylabel("Protein", fontdict=font1, labelpad=20)
plt.xlabel("Ligand", fontdict=font1, labelpad=20)
plt.tight_layout()
plt.show()
