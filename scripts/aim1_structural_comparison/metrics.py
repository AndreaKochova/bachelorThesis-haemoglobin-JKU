import os
import glob
import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd

from Bio import pairwise2
from Bio.SeqUtils import seq1

# --------------------------
# User settings
# --------------------------
ALN_DIR = "chimerax_aligned"
OUT_DIR = "analysis_out"

REF_FILE = os.path.join(ALN_DIR, "ref_aln.pdb")  # adjust name if different

SHELL1 = 5.0
SHELL2 = 7.0

HEME_NAMES = {"HEM"}  # reference has HEM
# --------------------------

os.makedirs(OUT_DIR, exist_ok=True)


def get_chain_residues_and_sequence(u: mda.Universe) -> tuple[list[mda.core.groups.Residue], str]:
    """Return protein residues list (in order) and 1-letter sequence."""
    prot = u.select_atoms("protein")
    residues = list(prot.residues)
    seq = ""
    for r in residues:
        # MDAnalysis gives 3-letter res names; convert to 1-letter
        try:
            seq += seq1(r.resname)
        except Exception:
            seq += "X"
    return residues, seq


def find_reference_heme(u_ref: mda.Universe):
    heme = u_ref.select_atoms("resname " + " ".join(HEME_NAMES) + " and not name H*")
    if len(heme) == 0:
        raise RuntimeError("Reference PDB contains no HEM heavy atoms. Reference must retain HEM.")
    return heme


def residues_within_shell(u: mda.Universe, heme_atoms, rmin: float, rmax: float | None = None) -> set[int]:
    """
    Return a set of residue indices (0-based residue index in the protein residue list)
    that have any heavy atom within distance to heme atoms.
    NOTE: returns indices into protein residue list, not PDB res ids.
    """
    prot = u.select_atoms("protein and not name H*")
    prot_res = list(prot.residues)
    prot_res_atoms = [r.atoms.select_atoms("not name H*") for r in prot_res]

    heme_pos = heme_atoms.positions

    idxs = set()
    for i, atoms in enumerate(prot_res_atoms):
        if len(atoms) == 0:
            continue
        # minimum heavy-atom distance residue ↔ heme
        dmin = np.min(np.linalg.norm(atoms.positions[:, None, :] - heme_pos[None, :, :], axis=2))
        if rmax is None:
            if dmin <= rmin:
                idxs.add(i)
        else:
            if (dmin > rmin) and (dmin <= rmax):
                idxs.add(i)
    return idxs


def sequence_map_ref_to_target(ref_seq: str, tgt_seq: str) -> dict[int, int]:
    """
    Map reference residue index -> target residue index via global alignment.
    Indices are 0-based positions in their respective residue lists.
    """
    # Global alignment with affine gaps; conservative defaults
    aln = pairwise2.align.globalms(ref_seq, tgt_seq, 2, -1, -5, -0.5, one_alignment_only=True)[0]
    ref_aln, tgt_aln = aln.seqA, aln.seqB

    mapping = {}
    ref_i = -1
    tgt_i = -1
    for r_char, t_char in zip(ref_aln, tgt_aln):
        if r_char != "-":
            ref_i += 1
        if t_char != "-":
            tgt_i += 1
        if (r_char != "-") and (t_char != "-"):
            mapping[ref_i] = tgt_i
    return mapping


def ca_positions_for_residue_indices(u: mda.Universe, residue_indices: list[int]) -> np.ndarray:
    """Return CA positions for a list of protein residue indices (0-based)."""
    prot_res, _ = get_chain_residues_and_sequence(u)
    coords = []
    for idx in residue_indices:
        if idx < 0 or idx >= len(prot_res):
            coords.append([np.nan, np.nan, np.nan])
            continue
        ca = prot_res[idx].atoms.select_atoms("name CA")
        if len(ca) != 1:
            coords.append([np.nan, np.nan, np.nan])
        else:
            coords.append(ca.positions[0])
    return np.array(coords, dtype=float)


# --------------------------
# Load reference
# --------------------------
u_ref = mda.Universe(REF_FILE)
ref_residues, ref_seq = get_chain_residues_and_sequence(u_ref)

ref_ca = u_ref.select_atoms("protein and name CA")
ref_heme = find_reference_heme(u_ref)

# Define PocketSet on the reference (indices into reference residue list)
shell1_ref = residues_within_shell(u_ref, ref_heme, SHELL1, None)
shell2_ref = residues_within_shell(u_ref, ref_heme, SHELL1, SHELL2)
pocket_ref = sorted(shell1_ref.union(shell2_ref))

# Precompute reference pocket CA coords
ref_pocket_ca = ca_positions_for_residue_indices(u_ref, pocket_ref)
ref_pocket_ca = ref_pocket_ca[~np.isnan(ref_pocket_ca).any(axis=1)]

rows = []
per_res_rows = []

# --------------------------
# Process all aligned PDBs (including modern apo, mutants, ancestors)
# --------------------------
for f in sorted(glob.glob(os.path.join(ALN_DIR, "*.pdb"))):

    tag = os.path.basename(f).replace(".pdb", "")
    u = mda.Universe(f)

    # --- build target sequence ---
    tgt_residues, tgt_seq = get_chain_residues_and_sequence(u)

    # --- sequence mapping reference → target ---
    mapping = sequence_map_ref_to_target(ref_seq, tgt_seq)

    # --- Global RMSD using mapped residues ---
    ref_match = []
    tgt_match = []

    for ref_idx, tgt_idx in mapping.items():

        ref_ca_atom = ref_residues[ref_idx].atoms.select_atoms("name CA")
        tgt_ca_atom = tgt_residues[tgt_idx].atoms.select_atoms("name CA")

        if len(ref_ca_atom) == 1 and len(tgt_ca_atom) == 1:
            ref_match.append(ref_ca_atom.positions[0])
            tgt_match.append(tgt_ca_atom.positions[0])

    ref_match = np.array(ref_match)
    tgt_match = np.array(tgt_match)

    if len(ref_match) > 10 and ref_match.shape == tgt_match.shape:
        global_r = rmsd(tgt_match, ref_match, center=True, superposition=True)
    else:
        global_r = np.nan
    

    # Sequence mapping ref->target
    tgt_residues, tgt_seq = get_chain_residues_and_sequence(u)
    mapping = sequence_map_ref_to_target(ref_seq, tgt_seq)

    # Transfer pocket residues by index mapping
    pocket_tgt = [mapping[i] for i in pocket_ref if i in mapping]
    pocket_tgt = sorted(set(pocket_tgt))

    # Pocket RMSD (CA) using mapped indices
    mob_pocket_ca = ca_positions_for_residue_indices(u, pocket_tgt)
    mob_ok = mob_pocket_ca[~np.isnan(mob_pocket_ca).any(axis=1)]

    # For RMSD we need the corresponding reference positions for those mapped residues.
    # Build the reference indices that successfully mapped:
    ref_mapped = [i for i in pocket_ref if i in mapping]
    ref_mapped = sorted(ref_mapped)
    ref_mapped_ca = ca_positions_for_residue_indices(u_ref, ref_mapped)
    ref_ok = ref_mapped_ca[~np.isnan(ref_mapped_ca).any(axis=1)]

    if len(mob_ok) == len(ref_ok) and len(mob_ok) > 0:
        pocket_r = rmsd(mob_ok, ref_ok, center=True, superposition=True)
    else:
        pocket_r = np.nan

    # Per-residue CA displacement (in aligned coordinates; no extra superposition here)
    # We report displacements for mapped pocket residues only.
    for ref_idx in ref_mapped:
        tgt_idx = mapping[ref_idx]
        ref_ca_atom = ref_residues[ref_idx].atoms.select_atoms("name CA")
        tgt_ca_atom = tgt_residues[tgt_idx].atoms.select_atoms("name CA")
        if len(ref_ca_atom) == 1 and len(tgt_ca_atom) == 1:
            disp = float(np.linalg.norm(tgt_ca_atom.positions[0] - ref_ca_atom.positions[0]))
            per_res_rows.append({
                "structure": tag,
                "ref_res_index": ref_idx,
                "ref_resname": ref_residues[ref_idx].resname,
                "tgt_res_index": tgt_idx,
                "tgt_resname": tgt_residues[tgt_idx].resname,
                "dCA_A": disp
            })

    rows.append({
        "structure": tag,
        "n_ref_res": len(ref_residues),
        "n_tgt_res": len(tgt_residues),
        "global_RMSD_CA": global_r,
        "pocket_RMSD_CA": pocket_r,
        "pocket_size_ref": len(pocket_ref),
        "pocket_size_mapped": len(pocket_tgt)
    })

df = pd.DataFrame(rows).sort_values("structure")
df.to_csv(os.path.join(OUT_DIR, "Aim1_summary_metrics.csv"), index=False)

df_res = pd.DataFrame(per_res_rows)
df_res.to_csv(os.path.join(OUT_DIR, "Aim1_pocket_residue_deviation.csv"), index=False)

print("Wrote:")
print(" -", os.path.join(OUT_DIR, "Aim1_summary_metrics.csv"))
print(" -", os.path.join(OUT_DIR, "Aim1_pocket_residue_deviation.csv"))