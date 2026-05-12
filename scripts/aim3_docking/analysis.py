import re
import csv
from pathlib import Path
from collections import defaultdict

# --------- CONFIG ---------
ROOT = Path(".")  # set to your docking results root folder
TOTAL_RUNS = 100  # GA runs

PROTEINS = [
    "1A3N", "2DN3", "6BB5",
    "mut_D99N", "mut_K82N", "mut_N102T", "mut_Y35F",
    "AncA", "AncB", "AncMH_nat"
]

LIGANDS = ["2,3-BPG", "heme_Co", "heme_Fe", "heme_O2", "heme", "IHP"]

# For alternate spellings, add aliases here:
ALIASES = {
    "2,3-BPG": ["2,3-BPG", "BPG", "23BPG", "2_3-BPG", "2,3BPG", "2_3_BPG"],
    "heme_Co": ["heme_Co", "heme_CO", "hemeCO"],
    "heme_Fe": ["heme_Fe", "heme_FE", "hemeFe"],
    "heme_O2": ["heme_O2", "heme_O2", "hemeO2"],
}

# --------- HELPERS ---------
def infer_label(path_str, targets):
    """Return first matching target found in path_str, else None."""
    for t in targets:
        if t in path_str:
            return t
    return None

def infer_ligand(path_str):
    # Try exact ligand tokens first
    lig = infer_label(path_str, LIGANDS)
    if lig:
        return lig
    # Try aliases
    for canon, variants in ALIASES.items():
        for v in variants:
            if v in path_str:
                return canon
    return None

def parse_clusters_from_histogram(text):
    """
    Parses ASCII CLUSTERING HISTOGRAM rows:
      1 |  -8.97 | 43 | -8.29 | 24 | ####
    Returns list of dicts.
    """
    start = text.find("CLUSTERING HISTOGRAM")
    if start == -1:
        return []

    segment = text[start:start+12000]

    row_re = re.compile(
        r"^\s*(\d+)\s*\|\s*([+\-]?\d+(?:\.\d+)?)\s*\|\s*(\d+)\s*\|\s*([+\-]?\d+(?:\.\d+)?)\s*\|\s*(\d+)",
        re.MULTILINE
    )
    rows = row_re.findall(segment)
    clusters = []
    for r in rows:
        clusters.append({
            "cluster_rank": int(r[0]),
            "lowest_energy": float(r[1]),
            "best_run": int(r[2]),
            "mean_energy": float(r[3]),
            "members": int(r[4]),
        })
    return clusters

def choose_best_cluster(clusters, delta_kcal=0.5):
    """
    Selection logic aligned with our manual rule:
    - Use multi-member clusters if available.
    - Find best (most negative) mean energy.
    - Consider all clusters within delta 0.5 kcal of the best mean.
    - Among those, pick the largest cluster (members).
    """
    if not clusters:
        return None, False

    multi = [c for c in clusters if c["members"] >= 2]
    use_multi = True
    candidates = multi if multi else clusters
    if not multi:
        use_multi = False

    best_mean = min(c["mean_energy"] for c in candidates)
    near_best = [c for c in candidates if c["mean_energy"] <= best_mean + delta_kcal]

    near_best.sort(key=lambda c: (-c["members"], c["mean_energy"], c["cluster_rank"]))
    return near_best[0], use_multi

# --------- MAIN ---------
results = []
coverage = defaultdict(lambda: {"found": False, "path": ""})

dlg_files = list(ROOT.rglob("*.dlg"))

for dlg in dlg_files:
    pstr = str(dlg)
    protein = infer_label(pstr, PROTEINS)
    ligand = infer_ligand(pstr)

    # If we cannot infer metadata, skip but keep traceable output
    if protein is None or ligand is None:
        results.append({
            "protein": protein or "UNKNOWN",
            "ligand": ligand or "UNKNOWN",
            "dlg_path": str(dlg),
            "error": "Could not infer protein/ligand from path/filename",
        })
        continue

    text = dlg.read_text(errors="ignore")
    clusters = parse_clusters_from_histogram(text)

    best, used_multi = choose_best_cluster(clusters)
    if best is None:
        results.append({
            "protein": protein, "ligand": ligand, "dlg_path": str(dlg),
            "error": "No CLUSTERING HISTOGRAM found or no parsable cluster rows",
        })
        continue

    members_frac = best["members"] / TOTAL_RUNS
    low_conf = best["members"] < 5  # adjustable

    results.append({
        "protein": protein,
        "ligand": ligand,
        "dlg_path": str(dlg),
        "best_cluster_rank": best["cluster_rank"],
        "best_run": best["best_run"],
        "members": best["members"],
        "members_fraction": round(members_frac, 3),
        "mean_energy": best["mean_energy"],
        "lowest_energy": best["lowest_energy"],
        "multi_member_clusters_present": used_multi,
        "low_confidence": low_conf,
        "error": "",
    })

    coverage[(protein, ligand)] = {"found": True, "path": str(dlg)}

# Write master summary
out_master = Path("master_docking_summary.csv")
fields = sorted({k for r in results for k in r.keys()})
with out_master.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(results)

# Write coverage matrix (expected combinations)
out_cov = Path("coverage_matrix.csv")
with out_cov.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["protein"] + LIGANDS)
    for prot in PROTEINS:
        row = [prot]
        for lig in LIGANDS:
            row.append("YES" if coverage[(prot, lig)]["found"] else "NO")
        w.writerow(row)

print(f"Wrote: {out_master.resolve()}")
print(f"Wrote: {out_cov.resolve()}")
print(f"DLGs scanned: {len(dlg_files)}")
