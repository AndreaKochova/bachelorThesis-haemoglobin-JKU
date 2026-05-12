# Aim 1: Structural Comparison of Haemoglobin Variants

## Purpose

The purpose of this analysis was to compare the structural similarity and variability of different haemoglobin variants. The aim was to determine whether modern, mutant, and ancestral haemoglobin structures remain globally conserved or show localized structural differences that could influence ligand-binding behaviour.

This aim provides the structural baseline for the later docking analyses.

## Protein Structures

The analysed haemoglobin variants included modern, mutant, and ancestral protein structures.

### Modern haemoglobins

- 1A3N
- 2DN3
- 6BB5

### Mutant haemoglobins

- mut_D99N
- mut_K82N
- mut_N102T
- mut_Y35F

### Ancestral haemoglobins

- AncA
- AncB
- AncMH_nat

## Reference Structure

The modern haemoglobin structure `1A3N` was used as the main structural reference. All other structures were compared against this reference to quantify structural deviation.

## Method Summary

The structural comparison consisted of the following steps:

1. Protein structures were prepared and cleaned.
2. Non-standard residues, ligands, and solvent molecules were removed where necessary.
3. Structures were aligned to the reference haemoglobin structure.
4. RMSD values were calculated to quantify structural deviation.
5. Residue-level variability was examined to identify regions of increased structural difference.

## RMSD

Root-mean-square deviation, or RMSD, measures the average spatial deviation between equivalent atoms after structural alignment.

A low RMSD indicates that two structures are globally similar, whereas a higher RMSD indicates greater structural deviation.

In this project, RMSD was used to assess whether the haemoglobin variants differed strongly at the global structural level or whether their differences were more localized.

## Interpretation

The structural comparison showed that haemoglobin variants were generally conserved at the global fold level. However, localized differences were still relevant because small structural changes may affect ligand accessibility, docking pose distribution, and central cavity organization.

This distinction is important because ligand-binding behaviour may change even when the overall protein fold remains similar.

## Main Outputs

Expected output files from this analysis include:

```text
results/aim1/tables/
results/aim1/figures/