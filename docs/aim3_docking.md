# Aim 3: Docking Energy and Cluster Reproducibility Patterns

## Purpose

The purpose of this analysis was to compare docking behaviour across all protein-ligand combinations using two complementary docking metrics:

1. Predicted binding energy
2. Cluster reproducibility

This aim was designed to evaluate whether energetically favourable docking results also showed reproducible geometric convergence.


## Core Idea

Molecular docking should not be interpreted only by the most negative binding energy value.
A docking result can appear energetically favourable while still producing poses that are spread across different regions or clusters. In that case, the predicted interaction may be less geometrically stable or less specific.

Therefore, this aim compared binding energy with cluster reproducibility.

    Binding energy          = predicted interaction favourability
    Cluster reproducibility = consistency of pose convergence across docking runs

Together, these metrics give a broader view of docking behaviour.



## Input Data

The input data for this aim came from AutoDock docking log files.

Each docking log contained information about:

- predicted docking energies
- ranked poses
- docking clusters
- cluster sizes
- representative poses

Typical raw input location:

    data/raw/aim3_docking_summary_inputs/

Typical processed output location:

    data/processed/aim3_heatmap_tables/



## Binding Energy

Binding energy was extracted from AutoDock output files. In AutoDock-style interpretation, more negative binding energy values indicate more favourable predicted ligand-protein interaction. However, binding energy alone does not show whether the ligand repeatedly converged to the same region. This is why cluster reproducibility was analysed separately.



## Cluster Reproducibility

Cluster reproducibility describes how consistently docking runs produced similar ligand poses. A highly reproducible docking result means that many independent docking runs converged into the same or similar cluster.A less reproducible result means that docking poses were distributed across several clusters. Cluster reproducibility is useful because it provides geometric information that binding energy alone cannot capture.



## Why Energy and Reproducibility Were Compared

Binding energy and cluster reproducibility answer different questions.

| Metric | Main Question |
|---|---|
| Binding energy | Is the predicted interaction energetically favourable? |
| Cluster reproducibility | Do docking runs converge consistently? |

These two measures do not always agree.

Possible cases:

| Pattern | Interpretation |
|---|---|
| Strong energy + high reproducibility | Favourable and geometrically consistent docking |
| Strong energy + low reproducibility | Favourable predicted energies, but multiple possible poses |
| Moderate energy + high reproducibility | Less strong energy, but stable pose convergence |
| Weak energy + low reproducibility | Weak and inconsistent docking behaviour |

This comparison was central to Aim 3.



## Method Summary

The analysis followed this workflow:

1. AutoDock docking log files were parsed.
2. Binding energy values were extracted.
3. Cluster information was extracted.
4. Summary statistics were calculated per protein-ligand pair.
5. Binding energy tables were generated.
6. Cluster reproducibility tables were generated.
7. Heatmaps were produced for visual comparison.
8. Energy and reproducibility patterns were interpreted together.



### Processed Tables

 THe processed tables include:

    data/processed/aim3_heatmap_tables/binding_energy_summary.csv
    data/processed/aim3_heatmap_tables/cluster_reproducibility_summary.csv




## Heatmaps

Heatmaps were used to summarize large docking patterns across the full dataset. They allowed broad trends to be interpreted visually across modern, mutant, and ancestral haemoglobin groups. See results for figures.



## Main Observations

The analysis showed that predicted binding energy and cluster reproducibility were not always directly coupled.

General observations:

- Heme-derived ligands often showed strongly favourable predicted binding energies.
- Strong binding energy did not always correspond to high cluster reproducibility.
- 2,3-BPG often showed more interpretable spatial convergence despite more moderate predicted energies.
- IHP showed narrow spatial dispersion, consistent with strong electrostatic constraints.
- Docking behaviour differed not only between proteins, but also between ligand families.

These observations supported the interpretation that docking results should not be judged by energy values alone.



## Main Outputs

 osutput locations:

    results/aim3/tables/
    results/aim3/figures/




### Related Scripts

The scripts for this aim are stored in:

    scripts/aims3_docking/

