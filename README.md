# Modeling Ontogenetic Modulation of Drug–Drug Interactions in Pediatric Populations

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)

## Description

This repository contains the scripts used for the validation of a Generalized Additive Model (GAM) for the detection of Drug-Drug Interaction (DDI) signals with ontogenetic dynamics within spontaneous adverse event reporting systems in the pediatric population

## Key Concepts

- **Triplet**: combination of `drugA`, `drugB`, and `meddra`.
- **Ontogenic modulation**: change in the disproportionality pattern across NICHD stages.
- **Dynamics**: shape of the injected pattern across the 7 stages (`uniform`, `increase`, `decrease`, `plateau`, `inverse_plateau`).
- **IOR**: interaction measure on the multiplicative scale.
- **RERI**: interaction measure on the additive scale.
- **Null distribution**: empirical distribution built through stage-stratified permutation.
- **Power-calibrated subset**: subset of triplets retained to ensure fairer comparisons across methods.

## Repository Structure

```text
00_functions.R               # Core pipeline functions, global configuration, and utilities
01_theme.R                   # Shared plotting theme
10_augmentation.R            # Semi-synthetic generation, positives, negatives, and sensitivity
20_null.R                    # Empirical null distribution construction
30_metrics.R                 # Metric calculation, power, and method comparisons
40_network.R                 # Network analysis, biological support, and comparison with TWOSIDES
41_graphs.R                  # Final faceted figure generation
dataset_drugbank.R           # DrugBank input preparation

drug.csv                     # Drug information
drug_gene.csv                # Drug-gene/protein relationships
ade_raw.csv                  # Curated spontaneous reporting dataset
twosides/                    # Reference data for external comparison
vocabulary/                  # OMOP/MedDRA/ATC/RxNorm vocabulary
```

## Input Data

### Main File

`ade_raw.csv` must contain, at a minimum:

- `safetyreportid`: unique report identifier.
- `atc_concept_id`: ATC drug identifier.
- `meddra_concept_id`: MedDRA event identifier.
- `nichd`: NICHD pediatric development stage.

### Other Required Inputs

- `drug.csv`: drug names and metadata.
- `drug_gene.csv`: drug-gene/protein pairs for biological support.
- `vocabulary/concept.csv`: OMOP, MedDRA, ATC, and RxNorm mappings.
- `vocabulary/concept_relationship.csv`: vocabulary relationships.
- `twosides/TWOSIDES.csv.gz`: external comparison dataset.

## Global Configuration

The central configuration lives in `00_functions.R`.

Active parameters in the current pipeline version:

```r
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
include_nichd <- FALSE
nichd_spline <- FALSE
k_spline <- 7
bs_type <- "cs"
select <- FALSE
method <- "fREML"
percentil <- "p95"
```

With this combination, the default results suffix is:

```r
suffix <- "sics"
```

Therefore, outputs from the current configuration are written to `results/sics/`.

## Methodological Summary

### 1. Semi-Synthetic Generation (`10_augmentation.R`)

- harmonizes drug IDs that correspond to the same active compound;
- builds candidate triplets from reports with at least two drugs and one event;
- filters positives with at least `2` reports and presence in at least `2` stages;
- selects `500` base positive triplets;
- injects ontogenic dynamics using five profiles (`uniform` included as comparator);
- assigns effect sizes using a negative exponential distribution (`lambda_fc = 0.75`);
- fits the GAM and computes `log(IOR)` and `RERI` for each triplet;
- generates a negative set of up to `10000` triplets using drugs and events from the same universe;
- performs iterative dataset reduction (`10%` to `90%`) for sensitivity analysis;
- selects `100000` reports to build the pool used for the null distribution.

### 2. Null Distribution (`20_null.R`)

- takes the pool saved by `10_augmentation.R`;
- permutes drugs and events within each developmental stage while preserving stage-specific structure;
- generates permuted triplets and reinserts those permutations into the original dataset;
- fits the GAM for triplets in the null universe;
- computes, by stage, empirical percentiles (`p90`, `p95`, `p99`) of the lower bound of the 90% CI for `log(IOR)` and `RERI`.

### 3. Metrics and Power Calibration (`30_metrics.R`)

- expands positive and negative results by stage;
- compares observed distributions against the null distribution;
- calculates sensitivity, specificity, `PPV`, `NPV`, `F1`, and `AUC`;
- estimates intervals through nonparametric bootstrap (`n_boot = 2000`);
- defines subsets calibrated to a target power of `80%`;
- evaluates results globally, by dynamics, and by stage;
- analyzes three scenarios: `original`, `filtered`, and `intersection`.

### 4. Networks and External Validation (`40_network.R`)

- fits candidates from the original dataset with at least `5` reports per triplet;
- detects positive signals using the GAM-IOR criterion and stage-specific null thresholds;
- integrates drug-gene/protein relationships derived from DrugBank;
- evaluates biological support through shared genes and the Jaccard index;
- contrasts the observed network against null networks generated by partial rewiring;
- compares positive pairs and triplets against `TWOSIDES`;
- separates concordant, potentially novel, and missed signals.

### 5. Final Figures (`41_graphs.R`)

- generates figures for simulated dynamics;
- produces faceted figures by stage, metric, and dataset version;
- exports main-result plots in `png` and `svg` format.

## Detection Criteria

### GAM

A triplet is classified as positive if there is at least one stage where:

- the lower bound of the 90% CI for `log(IOR)` or `RERI` is greater than `0`; and
- it also exceeds the empirical threshold from the null distribution for that stage.

By default, the percentile used is `p95`.

### Stratified Methods

Classic methods use a nominal per-stage criterion:

- `log(IOR)_lower90 > 0`
- `RERI_lower90 > 0`


## Main Outputs

Under the current configuration, results are organized in `results/sics/`:

```text
results/sics/
  augmentation_results/
  null_distribution_results/
  metrics_results/
  network/
```

Relevant files generated by the pipeline:

- `augmentation_results/positive_triplets_metadata.csv`
- `augmentation_results/positive_triplets_results*.csv`
- `augmentation_results/negative_triplets_results*.csv`
- `augmentation_results/null_pool_reports_metadata.csv`
- `null_distribution_results/null_distribution.csv`
- `null_distribution_results/null_thresholds.csv`
- `null_distribution_results/null_thresholds_reri.csv`
- `metrics_results/metrics_global_original.csv`
- `metrics_results/metrics_dynamic_original.csv`
- `metrics_results/metrics_stage_original.csv`

## Dependencies

The pipeline loads libraries through `pacman::p_load()`. Main dependencies include:

- `data.table`
- `tidyverse`
- `mgcv`
- `parallel`
- `doParallel`
- `foreach`
- `igraph`
- `ggraph`

## Practical Considerations

- `ade_raw.csv` is large, so working with sufficient RAM is recommended (16GB used).
- Parallelization uses `75%` of available CPU cores by default.
- The scripts are intended to be executed from the project root.
- The pipeline stores intermediate checkpoints to avoid losing progress in long runs.


