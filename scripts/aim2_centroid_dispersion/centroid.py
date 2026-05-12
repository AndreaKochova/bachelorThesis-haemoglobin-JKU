import re
import csv
from pathlib import Path
from statistics import mean, median

# ---------- Patterns (edit here if your DLG differs) ----------
RE_MODEL_START = re.compile(r'^\s*DOCKED:\s*MODEL\s+(\d+)', re.IGNORECASE)
RE_MODEL_END   = re.compile(r'^\s*DOCKED:\s*ENDMDL', re.IGNORECASE)

# Coordinates lines are often "DOCKED: ATOM ..." or "DOCKED: HETATM ..."
# PDB-format columns: atom serial, atom name, res name, chain, res seq, x, y, z ...
# We'll extract the last 3 floats on the line (x y z).
RE_COORD_LINE  = re.compile(r'^\s*DOCKED:\s*(ATOM|HETATM)\s+', re.IGNORECASE)

# Energy line often includes "Estimated Free Energy of Binding"
RE_ENERGY_LINE = re.compile(r'Estimated Free Energy of Binding\s*=\s*([-0-9.]+)', re.IGNORECASE)

def extract_last_three_floats(line: str):
    # robust float capture
    floats = re.findall(r'[-+]?\d*\.\d+|[-+]?\d+', line)
    if len(floats) < 3:
        return None
    x, y, z = map(float, floats[-3:])
    return x, y, z

def parse_dlg_models(dlg_path: Path):
    """
    Returns a list of dicts, one per MODEL block:
    {
        "model_id": int or None,
        "energy": float or None,
        "coords": [(x,y,z), ...]
    }
    """
    models = []
    current = None

    with dlg_path.open("r", errors="ignore") as f:
        for line in f:
            m_start = RE_MODEL_START.match(line)
            if m_start:
                # start new model
                current = {"model_id": int(m_start.group(1)), "energy": None, "coords": []}
                continue

            if current is not None:
                # energy
                m_e = RE_ENERGY_LINE.search(line)
                if m_e:
                    current["energy"] = float(m_e.group(1))

                # coordinates
                if RE_COORD_LINE.match(line):
                    xyz = extract_last_three_floats(line)
                    if xyz:
                        current["coords"].append(xyz)

                # end
                if RE_MODEL_END.match(line):
                    # only keep models that actually have coords
                    if current["coords"]:
                        models.append(current)
                    current = None

    return models

def centroid(coords):
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return (mean(xs), mean(ys), mean(zs))

def euclid(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5

def summarize_centroids(centroids_list):
    # centroids_list: [(cx,cy,cz), ...]
    cx_bar = mean([c[0] for c in centroids_list])
    cy_bar = mean([c[1] for c in centroids_list])
    cz_bar = mean([c[2] for c in centroids_list])
    c_bar = (cx_bar, cy_bar, cz_bar)

    rs = [euclid(c, c_bar) for c in centroids_list]
    rs_sorted = sorted(rs)
    # IQR
    q1 = rs_sorted[len(rs_sorted)//4]
    q3 = rs_sorted[(3*len(rs_sorted))//4]
    iqr = q3 - q1

    return {
        "N": len(centroids_list),
        "mean_r": mean(rs),
        "median_r": median(rs),
        "iqr_r": iqr
    }

def infer_protein_ligand_from_path(dlg_path: Path):
    """
    Assumes folder structure like: Parent/Protein/Ligand/dock.dlg
    Adjust if yours differs.
    """
    parts = dlg_path.parts
    if len(parts) >= 3:
        ligand = dlg_path.parent.name
        protein = dlg_path.parent.parent.name
        return protein, ligand
    return "UNKNOWN_PROTEIN", "UNKNOWN_LIGAND"

def main(parent_dir: str, out_centroids="centroids.csv", out_summary="summary.csv"):
    parent = Path(parent_dir)

    centroid_rows = []
    summary_rows = []

    dlg_files = list(parent.rglob("dock.dlg"))
    print(f"Found {len(dlg_files)} dock.dlg files under {parent}")

    for dlg in dlg_files:
        protein, ligand = infer_protein_ligand_from_path(dlg)
        models = parse_dlg_models(dlg)

        # sanity: how many MODEL blocks?
        print(f"{protein}/{ligand}: models parsed = {len(models)}  (expected ~100 if GA runs=100)")

        centroids_list = []
        for m in models:
            c = centroid(m["coords"])
            centroids_list.append(c)
            centroid_rows.append({
                "protein": protein,
                "ligand": ligand,
                "model_id": m["model_id"],
                "energy": m["energy"],
                "cx": c[0],
                "cy": c[1],
                "cz": c[2],
            })

        if centroids_list:
            summ = summarize_centroids(centroids_list)
            summary_rows.append({
                "protein": protein,
                "ligand": ligand,
                **summ
            })

    # Write outputs
    with open(out_centroids, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["protein","ligand","model_id","energy","cx","cy","cz"])
        w.writeheader()
        w.writerows(centroid_rows)

    with open(out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["protein","ligand","N","mean_r","median_r","iqr_r"])
        w.writeheader()
        w.writerows(summary_rows)

    print(f"Wrote {out_centroids} and {out_summary}")

# Example usage:
# main("DockingResults")

if __name__ == "__main__":
    main("DockingResults")

