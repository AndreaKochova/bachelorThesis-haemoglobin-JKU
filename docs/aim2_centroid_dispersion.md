# Aim 2: Centroid Dispersion Analysis

## Purpose

The purpose of this analysis was to evaluate how consistently ligands occupied spatial regions within haemoglobin docking results.

Instead of relying only on docking energy, this aim focused on the geometric distribution of docking poses. The central question was whether ligands repeatedly converged toward similar regions or whether they were broadly dispersed across the protein structure.

This aim was especially relevant for evaluating ligand behaviour in the haemoglobin central cavity.

## Core Idea

Each AutoDock run produces ligand poses with 3D atomic coordinates.

For each docking pose, a ligand centroid was calculated. The centroid represents the geometric center of the ligand pose.

In simplified form:

    centroid = average position of all ligand atoms

For a ligand with multiple atoms, the centroid is calculated as:

    cx = mean of all x coordinates
    cy = mean of all y coordinates
    cz = mean of all z coordinates

Each docking pose therefore becomes one point in 3D space.

This makes it possible to compare pose localization across many docking runs and across many protein-ligand systems.


## Why Centroids Were Used

Docking produces many ligand poses, and each pose contains many atoms. Direct atom-by-atom comparison across all poses would be unnecessarily complex for this analysis.

Centroids provide a compact spatial representation of each docking solution.

They allow the analysis to ask:

- Do poses cluster in a similar region?
- Are poses widely spread?
- Does one protein variant show more ligand dispersion than another?
- Are some ligands more spatially constrained?
- Does spatial convergence match binding energy or not?

This makes centroid dispersion a useful complement to docking energy.


## Method Summary

The input data for this aim came from AutoDock docking output files. For each protein-ligand pair, docking runs produced multiple ligand poses. These poses were parsed and converted into centroid coordinates.

The centroid analysis followed this workflow:

1. AutoDock docking output files were parsed.
2. Ligand pose coordinates were extracted.
3. One centroid was calculated for each docking pose.
4. Centroid coordinates were stored in a table.
5. Radial distances and dispersion metrics were calculated.
6. Results were summarized per protein-ligand pair.
7. Centroid distributions were visualized using plots and 3D figures.

## Centroid Table

The full centroid table stores one row per docking pose.

Table file:

    data/processed/aim2_centroids/centroids.csv

Columns:

    protein
    ligand
    pose_id
    energy
    cx
    cy
    cz
    r

### Column Meaning

| Column | Meaning |
|---|---|
| `protein` | Haemoglobin variant |
| `ligand` | Docked ligand |
| `pose_id` | Docking pose identifier |
| `energy` | Predicted docking energy for the pose |
| `cx` | Centroid x-coordinate |
| `cy` | Centroid y-coordinate |
| `cz` | Centroid z-coordinate |
| `r` | Radial distance used for dispersion analysis |

---

## Summary Table

The summary table stores one row per protein-ligand pair.

Summary file:

    data/processed/aim2_centroids/summary.csv

 With columns:

    protein
    ligand
    N
    mean_r
    sd_r
    median_r
    iqr_r

### Summary Metrics

| Metric | Meaning |
|---|---|
| `N` | Number of docking poses included |
| `mean_r` | Mean radial distance of centroids |
| `sd_r` | Standard deviation of radial distances |
| `median_r` | Median radial distance |
| `iqr_r` | Interquartile range of radial distances |



## Interpretation of Dispersion

Centroid dispersion was interpreted as a measure of spatial consistency.

    low dispersion  = poses converge into a similar region
    high dispersion = poses are broadly distributed

A low dispersion value suggests that docking repeatedly placed the ligand in a similar spatial region.

A high dispersion value suggests that the ligand had several possible docking regions or less constrained docking behaviour.

Importantly, centroid dispersion is not the same as binding energy. A ligand may show favourable predicted binding energy but still be spatially variable.



## Relation to Binding Energy

Binding energy describes predicted interaction favourability.

Centroid dispersion describes spatial consistency.

These two quantities answer different questions:

    Binding energy       → how favourable is the predicted interaction?
    Centroid dispersion  → how consistently does the ligand localize?

For this reason, Aim 2 provides information that binding energy alone cannot show.



## Main Outputs

Output locations:

    results/aim2/tables/
    results/aim2/figures/

