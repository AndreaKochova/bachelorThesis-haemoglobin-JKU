import csv

modern = "modern_centr.csv"
ancestral = "ancestral_centr.csv"
mutant = "mutant_centr.csv"

INPUT = mutant
OUTPUT = "centroids_K82N_2_3_BPG.pdb"

with open(INPUT) as f, open(OUTPUT, "w") as out:
    reader = csv.DictReader(f)
    i = 1
    for row in reader:
        x = float(row["cx"])
        y = float(row["cy"])
        z = float(row["cz"])
        out.write(
            f"HETATM{i:5d}  CEN CEN A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
        i += 1
