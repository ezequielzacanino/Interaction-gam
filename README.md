# Modeling Ontogenic Dynamics of Pharmacological Interactions in Pediatrics

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)


## Description

This repository contains the scripts used for the validation of a Generalized Additive Model (GAM) for the detection of Drug-Drug Interaction (DDI) signals with ontogenetic dynamics within spontaneous adverse event reporting systems in the pediatric population

## Concepts

- **Triplet**: DrugA-DrugB-Event combination.
- **Ontogenic modulation**: Variations in the disproportionality signal between developmental stages
- **Dynamic**: Temporal pattern of the ontogenic modulation (uniform, increase, decrease, plateau, inverse-plateau)
- **IOR (Interaction Odds Ratio)**: Metric quantifying the synergistic or antagonistic interaction between drugs.
- **RERI (Relative Excess Risk Due to Interaction)**: Metric quantifying the additive interaction between drugs.
- **Stratified methods**: Refers to the arithmetic calculation of IOR/RERI interaction measures by stratifying the dataset by stages.

## Workflow

```
00_functions.R               # Pipeline functions. 
01_theme.R                   # Plot theme.
02_descriptive.R             # Descriptive analysis of the dataset.
10_augmentation.R            # Data generation and results.
20_null.R                    # Generation of the null distribution.
30_metrics.R                 # Results analysis.
31_metrics_graphs.R          # Generation of plots.
40_network.R                 # Network analysis
```

## Methodology 

### 1. Semi-synthetic data generation (`10_augmentation.R`)

**Positive control set:**
- Assembly of triplets that meet minimum requirements (final works requires at least 2 A-B reports with A-B reports in at least 2 different ontogenic stages)
- Random selection of 500 candidate triplets
- Signal injection with simulated ontogenic dynamics (increase, decrease, plateau, inverse-plateau, uniform)
- Assignment of fold-changes sampled from a negative exponential distribution (λ = 0.75)
- Model fitting of the positive pool in an independent copy of the original dataset 
- 5 sets of 500 positive triplets (uniform excluded from final analysis) 

**Negative control set:**
- Random selection of 500 candidate triplets (mutually exclusivy from the positive set), with similar reporting characteristics as the positive set (same drugs and events)
- Model fitting of the negative set

**Injection Formula:**
```
P(event|stage) = fold_change × e_j + dynamic(stage)

Where e_j = P(A) + P(B) - P(A)×P(B)  
```

### 2. Null Distribution (`20_null.R`)

- Random selection of 100,000 raws 
- Stratified permutation of drugs and events by developmental stage (breaking drug-event associations and keeping stage reporting characteristics)
- Generation of 10000 triplets (maximum 25 per permutation) 
- Percentile calculation (P90, P95, P99) of lower bound 90IC by stage

### 3. Detection Criteria

A signal is considered **positive** if at least one of the stages meets the following criteria:
1. **Nominal criterion**: lower bound 90IC of the log-IOR/RERI > 0 (only criterion used in stratified methods)
2. **Null model criterion**: lower bound 90IC of the log-IOR/RERI > 95 Percentile of the null distribution 

### 4. GAM formula

**Base formula:**
```r
event ~ s(nichd_num, by=drugA, bs="cs", k=7) +
        s(nichd_num, by=drugB, bs="cs", k=7) +
        s(nichd_num, by=drugA_drugB, bs="cs", k=7) 
         
```

**Interaction Metrics:**
```
log(IOR) = lp₁₁ - lp₁₀ - lp₀₁ + lp₀₀           # predictions obtained with predict()
RERI = p₁₁ - p₁₀ - p₀₁ + p₀₀                   # probabilities obtained with plogis()
```

Where:
- ₁₁: DrugA + DrugB
- ₁₀: DrugA alone
- ₀₁: DrugB alone
- ₀₀: None


### Formula adjustment parameters

Configuration used in final work:
```r
spline_individuales <- TRUE    # Splines individual components
include_sex <- FALSE           # Includes sex as covariable
include_stage_sex <- FALSE     # Spline for sex
nichd_spline <- FALSE          # Spline for stage
bs_type <- "cs"                # Type of spline: "cs", "tp", "cr"
k_spline <- 7                  # Number of knots (= stages)
select <- FALSE                # Penalization allowing terms to shrink to zero
```

## Database structure

Columns required:
- `safetyreportid`: Unique report ID
- `atc_concept_id`: ATC drug code
- `meddra_concept_id`: MedDRA adverse event code
- `nichd`: stages according to the NICHD classification(factor ordenado)
- `nichd_num`: numeric NICHD stage (1-7)
- `sex`: Sex (opctional, if `include_sex = TRUE`)


