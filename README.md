# Comparative Analysis of Structure and Ligand Bindings in Haemoglobin Variants
Fully computational analysis of ligand-binding behaviour across modern, mutant, and ancestral haemoglobin variants.

## Overview
The repository contains the computational workflow, analysis scripts, processed results, and selected figures from the bachelor thesis on ligand docking behaviour across haemoglobin variants.

The project investigated how different haemoglobin structures interact with biologically relevant ligands using molecular docking, structural comparison, and Python-based post-processing. The analysis focused on docking energy patterns, cluster reproducibility, centroid dispersion of docked poses, and structural variability between haemoglobin variants.

## Research Aims

### Aim 1 - Structural Comparison
Quantify structural conservation and variability across modern, mutant, and ancestral haemoglobin variants using RMSD analysis.

### Aim 2 - Ligand Spatial Behaviour
Analyze docking pose convergence and ligand dispersion inside the haemoglobin cavity using centroid-based metrics.

### Aim 3 - Binding Energy vs. Reproducibility
Compare thermodynamic favourability and docking reproducibility across different ligand families.

## Materials and Methods

The workflow combined protein structure preparation, ligand preparation, blind docking, and downstream computational analysis.

Main tools used:
- UCSF ChimeraX for protein cleaning, structure inspection, and visualization.
- Avogadro for ligand preparation and geometry optimization.
- MGLTools / AutoDockTools for PDBQT preparation, charge assignment, grid setup, and docking configuration.
- AutoGrid and AutoDock 4 for blind docking.
- Python for parsing docking outputs, calculating summary metrics, and generating figures.
- MDAnalysis / Python-based workflows for structural comparison and RMSD analysis.

## Key Findings
- Strong binding affinity did not necessarily correspond to high docking reproducibility.
- Heme-derived ligands frequently showed highly favourable energies but moderate geometric convergence.
- 2,3-BPG displayed comparatively stable spatial clustering despite weaker average binding energies.
- IHP exhibited highly constrained docking distributions consistent with strong electrostatic guidance.
- Structural variability between haemoglobin classes remained globally low but localized differences affected ligand behaviour.


## Repository Structure

- manuscript/    Final thesis document and thesis presentation 
- docs/          Methodological notes and workflow documentation   
- data/          Processed result tables used for analysis 
- scripts/       Python scripts for parsing and analysis   
- results/       Resulting tables and figures 