################################################################################
# Functions script
# Script: 00_functions.R
# To use: source("00_functions.R", local = TRUE)
################################################################################

################################################################################
# General pipeline configuration
################################################################################

# Project setup
setwd("D:/Bioestadística/Interaction-gam")
set.seed(7113)

# Load all pipeline libraries using pacman
library(pacman)
pacman::p_load(tidyverse, data.table, pbapply, parallel, doParallel, foreach,
  mgcv, MASS, akima, doRNG, pROC, svglite, DHARMa,
  igraph, ggraph, tidygraph, scales, RColorBrewer, patchwork, graphlayouts, ggrepel, networkD3, htmlwidgets, ggalluvial)

# NICHD stage level ordering used throughout all scripts
niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

# Number of cores for parallelization (0.80 crashes with 16GB RAM)
n_cores <- max(1, floor(detectCores() * 0.50)) 

# General file paths
ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_gene <- "./drug_gene.csv"
ruta_drug_info <- "./drug.csv" 
ruta_concept <- "./vocabulary/concept.csv"
ruta_rel <- "./vocabulary/concept_relationship.csv"
ruta_twosides <- "./twosides/TWOSIDES.csv.gz"

# GAM formula parameters
spline_individuales <- TRUE 
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7    
include_nichd <- FALSE
nichd_spline <- FALSE 
bs_type <- "cs"
select <- FALSE
method <- "fREML" 

# Formula encoding string for output file naming
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

Z90 <- qnorm(0.95)  # 90th percentile quantile for confidence intervals

# Chosen percentile from the null distribution
percentil <- "p95"

source("01_theme.R", local = TRUE)
theme_set(theme_base())               

################################################################################
# Function to build triplets
################################################################################

# Generates triplets (drugA, drugB, event) for a single report
#
# Parameters:
# drug: drug IDs in the report (atc_concept_id)
# event: event IDs in the report (meddra_concept_id)
# report_id: report ID (safetyreportid)
# nichd_stage: numeric NICHD stage of the report
# 
# Return:
# data.table with columns: safetyreportid, drugA, drugB, meddra, nichd_num
# 
# Implementation:
# Generates all drug pair combinations (if >= 2 drugs)
# Crosses each pair with every event in the report
# Enforces ordering: drugA <= drugB (to avoid duplicates from reversed order)
# Returns NULL if the report does NOT have >= 2 drugs or at least 1 event

make_triplets <- function(drug, event, report_id, nichd_stage) {
  
  if (length(drug) < 2 || length(event) < 1) return(NULL)
  
  drug <- unique(drug)
  event <- unique(event)
  
  # All drug pair combinations
  if (length(drug) == 2) {
    combination <- matrix(c(min(drug), max(drug)), nrow = 1)
  } else {
    combination <- t(combn(drug, 2))
    # Ordering: min first
    combination <- t(apply(combination, 1, function(x) c(min(x), max(x))))
  }
  
  n_combination <- nrow(combination)
  n_events <- length(event)
  
  data.table(
    safetyreportid = report_id,
    drugA = rep(combination[,1], times = n_events),
    drugB = rep(combination[,2], times = n_events),
    meddra = rep(event, each = n_combination),
    nichd_num = nichd_stage
  )
}

################################################################################
# Effect size function (fold-change)
################################################################################

# Samples fold-changes from a negative exponential distribution
#
# Parameters:
# n: Number of fold-changes to generate
# lambda: Rate parameter for the exponential distribution (0.75)
# 
# Return: numeric vector of fold-changes >= 1
# 
# Implementation:
# - FC ~ 1 + exp(λ = 0.75)
# - Typical range: [1, 10] with right skew

fold_change <- function(n, lambda = 0.75) {
  1 + rexp(n, rate = lambda)
}

################################################################################
# Helper function: compute co-administration counts by stage for a triplet
################################################################################

# Computes the number of A+B co-administration reports per stage for each triplet
# Used to compute the superset for the classical method
# 
# Return:
# n_coadmin_stage: number of A+B reports for the event in each stage
#
# Implementation:
# Identifies A+B reports
# Counts reports per stage
# Fills stages with no reports with 0

coadmin_by_stage <- function(drugA, drugB, meddra, ade_data) {
  
  # Reports containing each drug
  reports_A <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  reports_B <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  
  # Reports with both drugs (co-administration)
  reports_AB <- intersect(reports_A, reports_B)
  
  if (length(reports_AB) == 0) {
    # If no co-administration exists, return empty structure
    return(data.table(
      nichd_num = 1:7,
      nichd = niveles_nichd,
      n_coadmin_stage = 0L
    ))
  }
  
  # NICHD stage of each co-administered report
  stage_counts <- unique(ade_data[
    safetyreportid %in% reports_AB,
    .(safetyreportid, nichd, nichd_num)
  ])[, .N, by = .(nichd, nichd_num)]
  
  # Fill stages with no reports with zeros
  full_stages <- data.table(
    nichd_num = 1:7,
    nichd = niveles_nichd
  )
  
  stage_counts <- merge(
    full_stages,
    stage_counts,
    by = c("nichd_num", "nichd"),
    all.x = TRUE
  )
  
  stage_counts[is.na(N), N := 0L]
  setnames(stage_counts, "N", "n_coadmin_stage")
  
  return(stage_counts[order(nichd_num)])
}

################################################################################
# Batch co-administration counts by NICHD stage
################################################################################

# Computes A+B co-administration counts by stage for a batch of triplets
# using a single set of joins instead of iterating pair-by-pair.
#
# Parameters:
# pairs_dt: data.table with triplet_id, drugA, drugB, and meddra
# ade_dt: pre-computed unique report-drug-stage table
#
# Return:
# data.table with one row per triplet_id x nichd_num
# and the corresponding n_coadmin_stage count

compute_coadmin_batch <- function(pairs_dt, ade_dt) {
  
  pair_meta <- unique(pairs_dt[, .(triplet_id, drugA, drugB, meddra)])
  
  if (nrow(pair_meta) == 0) {
    return(data.table(
      triplet_id = integer(),
      drugA = character(),
      drugB = character(),
      meddra = character(),
      nichd_num = integer(),
      nichd = character(),
      n_coadmin_stage = integer()
    ))
  }
  
  # Join each unique drug only once to avoid exploding rows when many
  # triplets share the same drugA or drugB.
  unique_pairs <- unique(pair_meta[, .(drugA, drugB)])
  unique_a <- unique(pair_meta[, .(drugA)])
  unique_b <- unique(pair_meta[, .(drugB)])
  
  reports_a <- ade_dt[
    unique_a,
    on = .(atc_concept_id = drugA),
    nomatch = 0L,
    .(drugA = i.drugA, safetyreportid, nichd_num)
  ]
  
  reports_b <- ade_dt[
    unique_b,
    on = .(atc_concept_id = drugB),
    nomatch = 0L,
    .(drugB = i.drugB, safetyreportid)
  ]
  
  # Build report-level co-administrations at the drug-pair level first.
  coadmin_drug <- reports_a[
    reports_b,
    on = .(safetyreportid),
    nomatch = 0L,
    allow.cartesian = TRUE,
    .(drugA, drugB = i.drugB, safetyreportid, nichd_num)
  ]
  
  # Keep only the drug pairs that were requested in pairs_dt.
  coadmin_drug <- coadmin_drug[
    unique_pairs,
    on = .(drugA, drugB),
    nomatch = 0L
  ]
  
  stage_counts <- coadmin_drug[
    , .(n_coadmin_stage = .N),
    by = .(drugA, drugB, nichd_num)
  ]
  
  full_grid <- CJ(
    row_id = seq_len(nrow(unique_pairs)),
    nichd_num = seq_along(niveles_nichd),
    unique = TRUE
  )
  full_grid[, `:=`(
    drugA = unique_pairs$drugA[row_id],
    drugB = unique_pairs$drugB[row_id]
  )]
  full_grid[, row_id := NULL]
  
  result_drug <- merge(
    full_grid,
    stage_counts,
    by = c("drugA", "drugB", "nichd_num"),
    all.x = TRUE
  )
  
  result_drug[is.na(n_coadmin_stage), n_coadmin_stage := 0L]
  result_drug[, nichd := niveles_nichd[nichd_num]]
  
  result <- merge(
    pair_meta,
    result_drug,
    by = c("drugA", "drugB"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  setcolorder(
    result,
    c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "nichd", "n_coadmin_stage")
  )
  
  return(result[order(triplet_id, nichd_num)])
}

################################################################################
# Dynamic pattern generation function
################################################################################

# Generates normalized dynamic patterns for signal injection
#
# Parameters:
# type: Dynamic type: "uniform", "increase", "decrease", "plateau", "inverse_plateau"
# N: Number of stages (7 for NICHD)
# 
# Return: numeric vector normalized to [-1, 1] representing the temporal pattern
# 
# Dynamic shapes:
# uniform: constant signal (0 across all stages)
# increase: monotonic increase from -1 to +1
# decrease: monotonic decrease from +1 to -1
# plateau: peak at central stages (bell shape)
# inverse_plateau: trough at central stages (U shape)

generate_dynamic <- function(type, N = 7) {
  type <- as.character(type)
  if (type == "uniform") {
    return(rep(0, N))
  }
  if (type == "increase") {
    return(tanh(seq(-pi, pi, length.out = N)))
  }
  if (type == "decrease") {
    return(-tanh(seq(-pi, pi, length.out = N)))
  }
  if (type == "plateau") {
    return(c(
      tanh(seq(-pi, pi, length.out = floor(N/2))),
      tanh(seq(pi, -pi, length.out = ceiling(N/2)))
    ))
  }
  if (type == "inverse_plateau") {
    return(c(
      tanh(seq(pi, -pi, length.out = floor(N/2))),
      tanh(seq(-pi, pi, length.out = ceiling(N/2)))
    ))
  }
}

################################################################################
# Signal injection function
################################################################################

# Injects a drug-drug interaction signal
#
# Parameters:
# drugA_id: Drug A ID (ATC concept_id)
# drugB_id: Drug B ID (ATC concept_id)
# event_id: Adverse event ID (MedDRA concept_id)
# dynamic_type: Dynamic type (calls generate_dynamic)
# fold_change: Effect magnitude (calls fold_change)
# ade_raw_dt: data.table with the original dataset
# 
# Return: List with: success, n_injected, n_coadmin, ade_aug, diagnostics
# 
# Implementation:
# 
# 1 Base rate (e_j):
# Baseline value to which the fold-change is applied
# Must be calculated accounting for individual drug components
# The baseline choice affects injection consistency across scenarios
# Script "simulacion_inyeccion" illustrates consistency under different scenarios
# Methods assuming event "independence" are used only as a proxy;
# they do not assume true independence
#
# 2 Fold_change (FC):
# Multiplicative effect size factor
# 
# 3 Dynamic f(j):
# Stage-normalized function in [0, 1]
# f(j) = calls generate_dynamic(type, N=7)  <--- N always 7 for NICHD stages
# 
# 4 Dynamic probability: (final reporting probability formula)
# p_dynamic(j) = e_j × FC + f(j)
# 
# 5 Simulation:
# For each report containing drugA + drugB
# Y_new ~ Bernoulli(p_dynamic(stage_j))
# Injects ONLY into events that do not already exist (if the event is present, it is kept;
# this attempts to model the dynamic on top of pre-existing events)

inject_signal <- function(drugA_id, drugB_id, event_id, 
                          dynamic_type, fold_change, 
                          ade_raw_dt) {
  
  if (drugA_id > drugB_id) {   # enforce ordering (already done in other functions, but added as a safety measure)
    temp <- drugA_id
    drugA_id <- drugB_id
    drugB_id <- temp
  }
  
  # Independent copy (so injected data does not corrupt the original dataset)
  ade_aug <- copy(ade_raw_dt)
  ade_aug[, simulated_event := FALSE]
  
  # 1- Reports containing both drugs (co-administration)
  reports_A <- unique(ade_aug[atc_concept_id == drugA_id, safetyreportid])
  reports_B <- unique(ade_aug[atc_concept_id == drugB_id, safetyreportid])
  reports_AB <- intersect(reports_A, reports_B)
  
  if (length(reports_AB) <= 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      ade_aug = NULL,
      message = sprintf(
        "coadministración insuficiente: %d reportes",
        length(reports_AB),
      ),
      diagnostics = list(
        reason = "insufficient_coadmin",
        n_coadmin = length(reports_AB),
        drugA = drugA_id,
        drugB = drugB_id,
        event = event_id
      )
    ))
  }
  
  # 2- Target reports (co-administration reports only)
  target_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_AB, 
    .(safetyreportid, nichd, nichd_num)
  ])
  
  # Identifies whether the event ALREADY exists in each report
  event_in_report <- unique(ade_raw_dt[
    meddra_concept_id == event_id, 
    safetyreportid
  ])
  
  target_reports[, e_old := as.integer(safetyreportid %in% event_in_report)]
  
  # 3- Compute base rate (e_j)

  # Compute individual base rates
  reports_A_clean <- setdiff(reports_A, reports_AB)
  reports_B_clean <- setdiff(reports_B, reports_AB)
  p_baseA <- mean(reports_A_clean %in% event_in_report)
  p_baseB <- mean(reports_B_clean %in% event_in_report)
  p_base0 <- length(event_in_report) / length(unique(ade_raw_dt$safetyreportid))

  # Note: these are not probabilities but odds; the extrapolation is valid when p < 0.1
 
  # Several methods exist to compute the base rate;
  # each behaves differently depending on the scenario
 
  # - Additive method:
  # Consistent; produces high IOR when individual risks are low and baseline risk is high
  # e_j = P(event | A ∪ B) assuming independence
  # Formula: P(A ∪ B) = P(A) + P(B) - P(A) × P(B)
  # p_baseA + p_baseB - (p_baseA * p_baseB)
 
  # - Multiplicative method:
  # Consistent; produces high IOR only when baseline risk is low and individual risks are high
  # Accounts for all components (individual and global risks)
  # e_j = (p_A × p_B) / p_0
  # Formula: P(A) × P(B) / P(0)
  # (p_baseA * p_baseB) / p_base0
  e_j <- p_baseA + p_baseB - (p_baseA * p_baseB)

  # t_ij = fold_change * e_j (effect size)
  t_ij <- fold_change * e_j

  # 4- Stage-specific probabilities
  # bprobs = rep(tij, N)
  # dy = tanh(...) * tij  (the dynamic is scaled by tij)
  # rprobs = bprobs + dy
  
  N <- 7
  bprobs <- rep(t_ij, N)
  
  if (dynamic_type == "uniform") {
    dy <- rep(0, N)
  } else if (dynamic_type == "increase") {
    dy <- tanh(seq(-pi, pi, length.out = N)) * t_ij
  } else if (dynamic_type == "decrease") {
    dy <- -tanh(seq(-pi, pi, length.out = N)) * t_ij
  } else if (dynamic_type == "plateau") {
    dy <- c(
      tanh(seq(-pi, pi, length.out = floor(N/2))) * t_ij,
      tanh(seq(pi, -pi, length.out = ceiling(N/2))) * t_ij
    )
  } else if (dynamic_type == "inverse_plateau") {
    dy <- c(
      tanh(seq(pi, -pi, length.out = floor(N/2))) * t_ij,
      tanh(seq(-pi, pi, length.out = ceiling(N/2))) * t_ij
    )
  }
  
  rprobs <- bprobs + dy

  # Clip probabilities to valid range
  rprobs <- pmax(pmin(rprobs, 0.999), 0.001)

  # 5- Stage probability table
  stage_probs <- data.table(
    nichd_num = 1:N,
    bprobs = bprobs,
    dy = dy,
    p_dynamic = rprobs
  )
  
  # 6- Merge and generate new events
  target_reports <- merge(
    target_reports, 
    stage_probs[, .(nichd_num, p_dynamic)], 
    by = "nichd_num", 
    all.x = TRUE
  )
  
  # Bernoulli simulation per report
  target_reports[, e_new := rbinom(.N, 1, p_dynamic)]
  
  # Combine with pre-existing events
  target_reports[, e_final := pmax(e_old, e_new)]
  
  # 7- Flag reports to inject
  # Identify reports that SHOULD have the event (simulated)
  reports_to_mark <- target_reports[e_old == 0 & e_final == 1, safetyreportid]
  
  # Validation: at least 1 event must have been injected
  if (length(reports_to_mark) == 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      ade_aug = NULL,
      message = sprintf(
        "Inyección fallida: 0 eventos generados (prob. media = %.4f, max = %.4f)",
        mean(target_reports$p_dynamic),
        max(target_reports$p_dynamic)
      ),
      diagnostics = list(
        reason = "zero_events_injected",
        low_probability_injection = TRUE,
        e_j = e_j,
        t_ij = t_ij,
        fold_change = fold_change,
        dynamic_type = dynamic_type,
        mean_p_dynamic = mean(target_reports$p_dynamic),
        max_p_dynamic = max(target_reports$p_dynamic),
        min_p_dynamic = min(target_reports$p_dynamic),
        n_eligible = nrow(target_reports[e_old == 0]),
        n_already_with_event = sum(target_reports$e_old),
        stage_probs = stage_probs
      )
    ))
  }
  
  # 8a- Change: (experimental) Validate that at least 1 injected event falls in high reporting stages
  # for the intended dynamic.
  # Reason: a nominally successful injection (n_injected >= 1) but with all
  # events in low reporting stages generate a structural false negative
  # the GAM cannot detect signal in the stages that the design intends to test.
  # The failure here is different from a method failure; It is a failure of the injection itself.
  high_stages_by_dynamic <- list(
    "uniform" = 1:7,      # all relevant stages
    "increase" = c(6L, 7L),
    "decrease" = c(1L, 2L),
    "plateau" = c(3L, 4L, 5L),
    "inverse_plateau"  = c(1L, 7L)
  )
  
  high_stages <- high_stages_by_dynamic[[dynamic_type]]
  
  injection_by_stage_temp <- target_reports[
    safetyreportid %in% reports_to_mark, .N, by = nichd_num]
  
  n_injected_high <- injection_by_stage_temp[nichd_num %in% high_stages, sum(N, na.rm = TRUE)]
  
  #sum() on empty table returns NA, not 0
  if (length(n_injected_high) == 0 || is.na(n_injected_high)) n_injected_high <- 0L
  
  if (n_injected_high == 0L) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = length(reports_to_mark),
      n_coadmin = length(reports_AB),
      ade_aug = NULL,
      message = sprintf(
        "Inyección sin señal en etapas clave: 0 eventos en etapas %s (dinámica: %s, total inyectado: %d)",
        paste(high_stages, collapse = ","),
        dynamic_type,
        length(reports_to_mark)
      ),
      diagnostics = list(            
      reason = "zero_events_in_high_stages",
      n_injected_high = n_injected_high,
      high_stages = high_stages,
      n_injected_total = length(reports_to_mark),
      e_j = e_j,
      t_ij = t_ij,
      fold_change = fold_change,
      dynamic_type = dynamic_type,
      stage_probs = stage_probs
    )
    ))
  }

  # 8b- Flag reports that now carry the simulated event (injected)
  ade_aug[
    safetyreportid %in% reports_to_mark,
    `:=`(
      simulated_event = TRUE,
      simulated_meddra = event_id,
      simulated_drugA = drugA_id,
      simulated_drugB = drugB_id
    )
  ]

  injection_rate <- length(reports_to_mark) / nrow(target_reports[e_old == 0])
  
  # 9- Diagnostics
  diagnostics <- list(
    e_j = e_j,
    t_ij = t_ij,
    fold_change = fold_change,
    dynamic_type = dynamic_type,
    stage_probs = stage_probs,
    mean_p_dynamic = mean(target_reports$p_dynamic),
    max_p_dynamic = max(target_reports$p_dynamic),
    min_p_dynamic = min(target_reports$p_dynamic),
    n_eligible = nrow(target_reports),
    n_already_with_event = sum(target_reports$e_old),
    n_without_event = nrow(target_reports[e_old == 0]),
    n_new_events = length(reports_to_mark),
    injection_rate = length(reports_to_mark) / nrow(target_reports[e_old == 0]),
    injection_by_stage = injection_by_stage_temp
  )
  
  diagnostics$n_injected_high <- n_injected_high
  diagnostics$high_stages <- high_stages
  diagnostics$n_injected_total <- length(reports_to_mark)

  return(list(
    success = TRUE,
    injection_success = TRUE,
    n_injected = length(reports_to_mark),
    n_coadmin = length(reports_AB),
    ade_aug = ade_aug,
    message = sprintf(
      "Inyección exitosa: %d eventos en %d reportes (tasa: %.2f%%)",
      length(reports_to_mark),
      length(reports_AB),
      injection_rate * 100
    ),
    diagnostics = diagnostics
  ))
}

################################################################################
# Helper function: compute basic counts
################################################################################

# Computes basic counts (events, co-administrations) for a drug-event pair
#
# Parameters:
#   ade_data: augmented data.table with reports (columns: atc_concept_id, meddra_concept_id, safetyreportid).
#   drugA: Drug A ID (atc_concept_id).
#   drugB: Drug B ID (atc_concept_id).
#   meddra: Event ID (meddra_concept_id).
#
# Return:
#   List with: n_events, n_events_coadmin, and n_coadmin

calc_basic_counts <- function(ade_data, drugA, drugB, meddra) {
  r_a <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  r_b <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  r_coadmin <- intersect(r_a, r_b)
  r_ea <- unique(ade_data[meddra_concept_id == meddra, safetyreportid])
  if("simulated_event" %in% names(ade_data)) {
    r_ea_sim <- unique(ade_data[simulated_event == TRUE & simulated_meddra == meddra, safetyreportid])
    r_ea <- union(r_ea, r_ea_sim)
  }
  list(
    n_events = length(r_ea),
    n_events_coadmin = length(intersect(r_coadmin, r_ea)),
    n_coadmin = length(r_coadmin)
  )
}

################################################################################
# GAM fitting function
################################################################################

# Fits a GAM model for drug-drug interaction
#
# Parameterized to facilitate iterative testing
#
# Parameters:
# drugA_id: Drug A ID
# drugB_id: Drug B ID
# event_id: Adverse event ID
# ade_data: data.table with the original dataset
# 
# include_nichd: If TRUE, adds the NICHD stage as a baseline covariate
# nichd_spline: If TRUE, fits the NICHD effect as a spline.
#               If FALSE, uses a linear coefficient (default: TRUE)
# spline_individuales: If TRUE, uses splines for individual drug effects (smooths baseline risks for drugA and drugB separately)
# bs_type: Spline basis type: "cs", "tp", or "cr"
# select: If TRUE, allows shrinkage to zero (useful for assessing whether a term contributes)
# include_sex: If TRUE, includes sex as a covariate
# include_stage_sex: If TRUE, includes a stage-by-sex interaction
# k_spline: Number of knots for splines (should always be 7 for NICHD stages)
# method: GAM fitting method (keep as "fREML")
# 
# Return:
# List with: success, n_events, n_coadmin, log_ior, reri, etc.

fit_gam <- function(drugA_id, drugB_id, event_id, ade_data,
                                 nichd_spline = TRUE,
                                 include_nichd = TRUE,
                                 spline_individuales = FALSE,
                                 bs_type = "cs",
                                 select = FALSE,
                                 include_sex = FALSE,
                                 include_stage_sex = FALSE,
                                 k_spline = 7,
                                 method = "fREML") {
  ###########
  # 1- Identify reports
  ###########
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  reportes_coadmin <- intersect(reportes_droga_a, reportes_droga_b)
  
  # Include simulated reports flagged during injection
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  # Pre-existing real events
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Combine real + simulated events
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  n_events_total <- length(reportes_ea)   # total events (with or without drug)
  n_coadmin <- length(reportes_coadmin)   # total A+B reports (with or without event)
  n_events_coadmin <- length(intersect(reportes_coadmin, reportes_ea))  # events in A+B reports

  ###########
  # 2- Build the modeling dataset
  ###########
  
  # Required base columns
  cols_necesarias <- c("safetyreportid", "nichd", "nichd_num")
  
  # Sex covariate
  if (include_sex) {
    cols_necesarias <- c(cols_necesarias, "sex")
  }
  
  datos_modelo <- unique(ade_data[, ..cols_necesarias])
  
  # Exposure variables
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # If sex is included
  if (include_sex) {
    # Standardize sex values
    datos_modelo[, sex := toupper(trimws(sex))]
    datos_modelo[sex == "M", sex := "MALE"]
    datos_modelo[sex == "F", sex := "FEMALE"]
    # Convert to factor with standard levels
    datos_modelo[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }
  
  ###########
  # 4- Build formula from parameters
  ###########
  
  # Response variable
  formula_parts <- "ea_ocurrio ~ "
  
  # Option A: linear individual drug effects
  if (!spline_individuales) {
    formula_parts <- paste0(formula_parts, "droga_a + droga_b + ")
  } else {
    # Option B: individual drug effects with splines
    formula_parts <- paste0(
      formula_parts,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", 
              k_spline, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", 
              k_spline, bs_type)
    )
  }

  # NICHD effect - spline or linear (if include_nichd = TRUE)
  if (include_nichd) {
    if (nichd_spline) {
    # Spline basis for NICHD effect
      formula_parts <- paste0(
        formula_parts,
        sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline, bs_type)
      )
    } else {
    # Linear NICHD effect
      formula_parts <- paste0(formula_parts, "nichd_num + ")
    }
  }

  # Interaction spline (do not modify this term)
  formula_parts <- paste0(
    formula_parts,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline, bs_type)
  )
  
  # If sex is included
  if (include_sex) {
    if (include_stage_sex) {
      # Sex-by-stage spline interaction
      formula_parts <- paste0(
        formula_parts,
        sprintf(" + s(nichd_num, k = %d, bs = '%s', by = sex)", 
                k_spline, bs_type)
      )
    } else {
      # Linear sex effect only
      formula_parts <- paste0(formula_parts, " + sex")
    }
  }
  
  # Final model formula
  formula_final <- as.formula(formula_parts)
  
  ###########
  # 5- Model fitting
  ###########
  
  tryCatch({
    
    modelo <- bam(
      formula = formula_final,
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = method,
      select = select,    
      discrete = TRUE,
      nthreads = 1  # Setting nthreads = 1 prevents conflicts with external parallelization in 10_augmentation
    )
    
    ###########
    # 6- Compute log-IOR per NICHD stage
    ###########
    
    # Prediction grid (all exposure combinations)
    grid_dif <- CJ(
      nichd_num = 1:7, 
      droga_a = c(0, 1), 
      droga_b = c(0, 1)
    )
    grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
    
    # If sex is in the formula, set reference level
    if (include_sex) {
      # Male as reference level
      grid_dif[, sex := factor("MALE", levels = c("MALE", "FEMALE"))]
    }
    
    # Predictions on the link scale
    pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
    grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
    
    # Pivot for contrast computation
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                  value.var = c("lp", "se"))
    
    # Final computation: log(IOR) = log(OR₁₁) - log(OR₁₀) - log(OR₀₁) + log(OR₀₀)
    # This is equivalent to log( OR₁₁ / (OR₁₀ · OR₀₁) )
    log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
    
    ###########
    # 7- Compute log-IOR standard error using the covariance matrix
    ###########
    
    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE) # unconditional = TRUE applies correction to include smoothing uncertainty
    cov_link <- Xp %*% Vb %*% t(Xp)
    
    log_ior_se <- numeric(7)
    for (stage in 1:7) {
      idx_00 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 0)
      idx_01 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 1)
      idx_10 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 0)
      idx_11 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 1)
      
      # Contrast vector
      cvec <- rep(0, nrow(grid_dif))
      cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)
      
      # SE = sqrt(c' Σ c)
      log_ior_se[stage] <- sqrt(max(
        as.numeric(t(cvec) %*% cov_link %*% cvec), 
        0
      ))
    }
    
    ###########
    # 8- Compute confidence intervals and metrics
    ###########
    
    # z90 for confidence intervals
    if (!exists("Z90")) {
      Z90 <- qnorm(0.95)
    }
    
    log_ior_lower90 <- log_ior - Z90 * log_ior_se
    log_ior_upper90 <- log_ior + Z90 * log_ior_se
    ior_values <- exp(log_ior)
    
    n_stages_significant <- sum(log_ior_lower90 > 0)
    max_ior <- max(ior_values)
    mean_ior <- mean(ior_values)

    ###########
    # 9- Compute RERI per stage (relative risk scale) with 90% CI
    ###########
    
    # Stages
    stages <- sort(unique(datos_modelo$nichd_num))
    
    # Prediction data with the 4 exposure combinations per stage
    nd_reri <- rbindlist(lapply(stages, function(s) {
      data.table(
        nichd_num = s,
        droga_a   = c(0, 1, 0, 1),
        droga_b   = c(0, 0, 1, 1),
        droga_ab  = c(0, 0, 0, 1)
      )
    }), use.names = TRUE)
    
    # Additional covariates (if present)
    if (include_sex) {
      nd_reri[, sex := factor(levels(datos_modelo$sex)[1], 
                             levels = levels(datos_modelo$sex))]
    }
    if (include_nichd && !nichd_spline) {
      nd_reri[, nichd := factor(niveles_nichd[nichd_num],
                                levels = niveles_nichd,
                                ordered = TRUE)]
    }
    
    # Prediction on the risk scale
    pred_reri <- predict(modelo, newdata = nd_reri, type = "link", se.fit = TRUE)
    nd_reri[, `:=`(
      eta = pred_reri$fit,
      se  = pred_reri$se.fit
    )]
    
    # Parametric bootstrap --> unlike the delta method, this captures non-linearity
    # Design matrix and coefficients
    X_reri <- predict(modelo, newdata = nd_reri, type = "lpmatrix")
    beta_hat <- coef(modelo)
    V_beta <- vcov(modelo, unconditional = TRUE)
    
    # Number of bootstrap simulations
    B <- 2000
    
    # Simulate coefficients from their joint distribution
    # beta_sim ~ MVN(beta_hat, V_beta)
    beta_sims <- mvrnorm(n = B, mu = beta_hat, Sigma = V_beta)
    
    # Compute predictions for each set of simulated coefficients
    p_sims <- matrix(NA, nrow = nrow(nd_reri), ncol = B)

    p_sims <- plogis(X_reri %*% t(beta_sims))

    # Helper function for RERI calculation
    calc_reri <- function(p) {
      # p is a vector of 4 probabilities per stage: [p00, p10, p01, p11]
      p11 <- p[4]; p10 <- p[2]; p01 <- p[3]; p00 <- p[1]
      if (p00 <= 0) return(NA_real_)       # against division by zero
      p11/p00 - p10/p00 - p01/p00 + 1
    }
    
    # Per-stage computation
    reri_dt <- nd_reri[, {
      idx <- .I
      p_mat <- p_sims[idx, , drop = FALSE]
      reri_sim <- apply(p_mat, 2, calc_reri)
      
      data.table(
        RERI = mean(reri_sim),
        RERI_lower90 = quantile(reri_sim, 0.05),
        RERI_upper90 = quantile(reri_sim, 0.95)
      )
    }, by = nichd_num]
    
    # Return vectors
    reri_values <- reri_dt$RERI
    reri_lower90 <- reri_dt$RERI_lower90
    reri_upper90 <- reri_dt$RERI_upper90
    
    ###########
    # 10- Results
    ###########

    return(list(
      success = TRUE,
      n_events = n_events_total,
      n_coadmin = n_coadmin,
      n_events_coadmin = n_events_coadmin,
      log_ior = log_ior,
      log_ior_lower90 = log_ior_lower90,
      log_ior_upper90 = log_ior_upper90,
      log_ior_se = log_ior_se,
      ior_values = ior_values,
      n_stages_significant = n_stages_significant,
      max_ior = max_ior,
      mean_ior = mean_ior,
      reri_values = reri_values,
      reri_lower90 = reri_lower90,
      reri_upper90 = reri_upper90,
      n_stages_reri_significant = sum(reri_lower90 > 0),
      model_aic = AIC(modelo),
      model_deviance = deviance(modelo),
      formula_used = formula_parts,  # Store the formula that was used
      nichd_spline = nichd_spline,
      include_nichd = include_nichd,
      spline_individuales = spline_individuales,
      bs_type = bs_type,              
      select = select,                
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE, 
      n_events = n_events_total, 
      n_coadmin = n_coadmin,
      error_msg = e$message,
      formula_attempted = formula_parts
    ))
  })
}

################################################################################
# Classical IOR calculation
################################################################################

# Computes classical IOR using 2x2 contingency tables per stage
#
# Parameters:
# drugA_id: Drug A ID (ATC concept_id)
# drugB_id: Drug B ID (ATC concept_id)
# event_id: Adverse event ID (MedDRA concept_id)
# ade_data: data.table with the dataset
# 
# Return: List with stage, ior_classic, ior_classic_lower90, ior_classic_upper90
# 
# Implementation:
# For each stage j:
# 2x2 contingency table:
# a: event + co-administration
# b: no event + co-administration
# c: event + no co-administration
# d: no event + no co-administration
# OR_11 = (a/b) / (c/d) = ad/bc
# OR_10 = reports with drug A only
# OR_01 = reports with drug B only
# OR_00 = reports with neither A nor B
# IOR = (OR_11 × OR_00) / (OR_10 × OR_01)
# 90% CI computed using the Woolf method (log scale)

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Report identification
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reports with the event (real/observed)
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Reports with simulated event
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  
  # Real + simulated events combined
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # One row per report with exposure indicators
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Per-stage calculation
  stage_results <- datos_unicos[, {
    
    # 2x2 table for each exposure combination
    # Group 11: A + B (co-administration)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_no_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 0)
    
    # Group 10: drug A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_no_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 0)
    
    # Group 01: drug B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_no_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 0)
    
    # Group 00: neither A nor B
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_no_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 0)
    
    # Compute OR for each exposure group
    or_11 <- (n_11_evento / n_11_no_evento) / (n_00_evento / n_00_no_evento)
    or_10 <- (n_10_evento / n_10_no_evento) / (n_00_evento / n_00_no_evento)
    or_01 <- (n_01_evento / n_01_no_evento) / (n_00_evento / n_00_no_evento)
    or_00 <- 1 # by definition
    
    # Since OR₀₀ = 1 by definition, this simplifies to IOR = OR₁₁ / (OR₁₀ · OR₀₁)
    # IOR computation
    ior_val <- (or_11 * or_00) / (or_10 * or_01)
    log_ior <- log(ior_val)
    
    # Variance on the log scale (Woolf method)
    var_log_ior <- (1/n_11_evento + 1/n_11_no_evento +
                1/n_10_evento + 1/n_10_no_evento +
                1/n_01_evento + 1/n_01_no_evento +
                1/n_00_evento + 1/n_00_no_evento)
    se_log_ior <- sqrt(var_log_ior)
    
    # 90% CI on the log scale
    z90 <- qnorm(0.95)
    log_ior_lower90 <- log_ior - z90 * se_log_ior
    log_ior_upper90 <- log_ior + z90 * se_log_ior
    
    # IOR back-transformed to the original scale
    ior_lower90 <- exp(log_ior_lower90)
    ior_upper90 <- exp(log_ior_upper90)
    
    data.table(
      stage = nichd_num[1],
      ior_classic = ior_val,
      log_ior_classic = log_ior,
      ior_classic_lower90 = ior_lower90,
      ior_classic_upper90 = ior_upper90,
      log_ior_classic_lower90 = log_ior_lower90,
      log_ior_classic_upper90 = log_ior_upper90,
      se_log_ior_classic = se_log_ior,
      # Diagnostics
      n_11_evento = n_11_evento,
      n_11_total = n_11_evento + n_11_no_evento
    )
    
  }, by = nichd_num]
  
  # Order by stage
  setorder(stage_results, nichd_num)
  
  return(list(
    success = TRUE,
    results_by_stage = stage_results
  ))
}

################################################################################
# Stratified RERI calculation
################################################################################

# Computes stratified RERI using 2x2 contingency tables per stage
#
# Parameters:
# drugA_id: Drug A ID (ATC concept_id)
# drugB_id: Drug B ID (ATC concept_id)
# event_id: Adverse event ID (MedDRA concept_id)
# ade_data: data.table with the dataset (may be augmented)
# 
# Return: List with success, results_by_stage (stage, RERI, RERI_lower90, RERI_upper90)
# 
# Implementation:
# For each stage j:
# Builds a 2x2 table:
# R11: risk with A+B
# R10: risk with A only
# R01: risk with B only
# R00: risk with neither A nor B (reference)
# RERI = R11 - R10 - R01 + 1
# 90% CI using the delta method (variance propagation)
# Variance of each risk: Var(R) = p(1-p)/n
# RERI variance: sum of individual variances (assuming independence)

calculate_classic_reri <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Report identification
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reports with the event (real/observed)
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Reports with simulated event (if present)
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  
  # Combine real + simulated events
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # One row per report with exposure indicators
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Per-stage calculation
  stage_results <- datos_unicos[, {
    
    # Counts for each exposure combination
    # Group 11: A + B (co-administration)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_total <- sum(droga_a == 1 & droga_b == 1)
    
    # Group 10: drug A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_total <- sum(droga_a == 1 & droga_b == 0)
    
    # Group 01: drug B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_total <- sum(droga_a == 0 & droga_b == 1)
    
    # Group 00: neither A nor B (for reference)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_total <- sum(droga_a == 0 & droga_b == 0)
    
    # Validation: all groups must have observations
    if (n_11_total == 0 || n_10_total == 0 || n_01_total == 0 || n_00_total == 0) {
      data.table(      
        stage = nichd_num[1],
        RERI_classic = NA_real_,
        RERI_classic_lower90 = NA_real_,
        RERI_classic_upper90 = NA_real_,
        RERI_classic_se = NA_real_,
        n_11_evento = n_11_evento, 
        n_11_total = n_11_total,
        n_10_evento = n_10_evento, 
        n_10_total = n_10_total,
        n_01_evento = n_01_evento, 
        n_01_total = n_01_total,
        n_00_evento = n_00_evento, 
        n_00_total = n_00_total,
        insufficient_data = TRUE
      )
  }
    else {  # Compute risks as proportions
      R11 <- n_11_evento / n_11_total
      R10 <- n_10_evento / n_10_total
      R01 <- n_01_evento / n_01_total
      R00 <- n_00_evento / n_00_total
    
      if (R00 == 0) {
        data.table(
          stage = nichd_num[1],
          RERI_classic = NA_real_,
          RERI_classic_lower90 = NA_real_,
          RERI_classic_upper90 = NA_real_,
          RERI_classic_se = NA_real_,
          n_11_evento = n_11_evento, n_11_total = n_11_total,
          n_10_evento = n_10_evento, n_10_total = n_10_total,
          n_01_evento = n_01_evento, n_01_total = n_01_total,
          n_00_evento = n_00_evento, n_00_total = n_00_total,
          insufficient_data = TRUE
        )
      } else {
        reri_val <- R11/R00 - R10/R00 - R01/R00 + 1

        # Binomial variances for each risk
        var_r <- function(r, n) ifelse(r > 0 & r < 1, r*(1-r)/n, 0.25/n)
        var_R11 <- var_r(R11, n_11_total)
        var_R10 <- var_r(R10, n_10_total)
        var_R01 <- var_r(R01, n_01_total)    
        var_R00 <- var_r(R00, n_00_total)

        # RERI variance by delta method:
        # RERI = (R11 - R10 - R01) / R00 + 1
        # dRERI/dR11 =  1/R00
        # dRERI/dR10 = -1/R00
        # dRERI/dR01 = -1/R00
        # dRERI/dR00 = -(R11 - R10 - R01) / R00^2
        numerador <- R11 - R10 - R01
        var_reri <- (var_R11 + var_R10 + var_R01) / R00^2 +
          numerador^2 * var_R00 / R00^4
        se_reri <- sqrt(var_reri)

        # 90% CI
        z90 <- qnorm(0.95)
        reri_lower90 <- reri_val - z90 * se_reri
        reri_upper90 <- reri_val + z90 * se_reri

        data.table(
          stage = nichd_num[1],
          RERI_classic = reri_val,
          RERI_classic_lower90 = reri_lower90,
          RERI_classic_upper90 = reri_upper90,
          RERI_classic_se = se_reri,
          # Individual risks for diagnostics
          R11 = R11, R10 = R10, R01 = R01, R00 = R00,
          # Diagnostic counts
          n_11_evento = n_11_evento, n_11_total = n_11_total,
          n_10_evento = n_10_evento, n_10_total = n_10_total,
          n_01_evento = n_01_evento, n_01_total = n_01_total,
          n_00_evento = n_00_evento, n_00_total = n_00_total,
          insufficient_data = FALSE
        ) 
      }
    }
  }, by = nichd_num]
  # Order by stage
  setorder(stage_results, nichd_num)
  
  # Check whether there is sufficient data in at least one stage
  if (all(is.na(stage_results$RERI_classic))) {
    return(list(
      success = FALSE,
      message = "Datos insuficientes en todas las etapas",
      results_by_stage = stage_results
    ))
  }
  
  return(list(
    success = TRUE,
    results_by_stage = stage_results
  ))
}

################################################################################
# Bootstrap function by dynamic type and stage
################################################################################

# Bootstraps the difference in log-IOR between a dynamic type and the baseline (uniform)
#
# Parameters:
# data: data.table with columns `dynamic`, `stage`, and `log_ior`
# dynamic_type: name of the dynamic type
# stage_num: NICHD stage number
# n_boot: number of bootstrap replicates (default 100)
#
# Return:
# data.table with bootstrap statistics for the difference (mean, sd, CI)

bootstrap_dynamic_diff <- function(data, dynamic_type, stage_num, n_boot = 100) {
  
  # Data for the target dynamic
  target_data <- data[dynamic == dynamic_type & stage == stage_num, log_ior]
  
  # Data for the uniform baseline
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, log_ior]
  
  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(
      mean_diff = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  # Bootstrap
  boot_diffs <- replicate(n_boot, {
    target_sample <- sample(target_data, replace = TRUE)
    uniform_sample <- sample(uniform_data, replace = TRUE)
    mean(target_sample) - mean(uniform_sample)
  })
  data.table(
    mean_diff = mean(boot_diffs, na.rm = TRUE),
    ci_lower = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}

################################################################################
# Bootstrap function by dynamic type and stage (RERI)
################################################################################

# Bootstraps the difference in RERI between a dynamic type and the baseline (uniform)
#
# Parameters:
# data: data.table with columns `dynamic`, `stage`, and `reri`
# dynamic_type: name of the dynamic type
# stage_num: NICHD stage number
# n_boot: number of bootstrap replicates (default 100)
#
# Return:
# data.table with bootstrap statistics for the difference (mean, sd, CI)

bootstrap_dynamic_diff_reri <- function(data, dynamic_type, stage_num, n_boot = 100) {
  
  target_data <- data[dynamic == dynamic_type & stage == stage_num, reri]
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, reri]
  
  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(mean_diff = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  
  boot_diffs <- replicate(n_boot, {
    median(sample(target_data,  replace = TRUE)) -   # median
    median(sample(uniform_data, replace = TRUE))
  })
  
  data.table(
    mean_diff = median(boot_diffs, na.rm = TRUE),    #  median
    ci_lower  = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper  = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}

################################################################################
# Function to process a single positive triplet with sensitivity analysis
################################################################################

# Processes a single positive triplet across all reduction levels
# Designed for batch parallel execution (memory-intensive)
#
# Parameters:
# idx: triplet index in pos_meta
# pos_meta: data.table with positive triplet metadata
# ade_raw_dt: data.table with the original dataset
# reduction_levels: numeric vector of reduction percentages to apply
# GAM formula configuration parameters
# base_seed: base seed for reproducibility
#
# Return:
#  Combined data.table with results for all reduction levels
#
# Implementation:
# inject_signal() to create an independent augmented dataset
# For each reduction level --> reduce dataset, fit GAM, compute classical IOR/RERI
# fit_reduced_model() as a wrapper for fitting

process_single_positive <- function(idx, pos_meta, ade_raw_dt, reduction_levels, 
                                    spline_individuales, include_sex, include_stage_sex,
                                    k_spline, bs_type, select, nichd_spline, z90, base_seed = 9427) {
  # Unique seed per triplet (same scheme as 10_augmentation)
  set.seed(9427 + idx)
  
  rowt <- pos_meta[idx]
  rowt$type <- "positive"
  
  # Create INDEPENDENT augmented dataset
  inj_result <- tryCatch({
    inject_signal(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      dynamic_type = rowt$dynamic,
      fold_change = rowt$fold_change,
      ade_raw_dt = ade_raw_dt
    )
  }, error = function(e) {
    list(
      success = FALSE, 
      injection_success = FALSE,
      n_injected = 0, 
      n_coadmin = 0,
      ade_aug = NULL,
      message = paste("Error en inyección:", e$message),
      diagnostics = list(reason = "exception", error = e$message)
    )
  })

  inj_success <- inj_result$success
  n_injected_val <- inj_result$n_injected
  n_coadmin_val <- inj_result$n_coadmin
  diag_data <- list(inj_result$diagnostics)
  inj_message <- if(!is.null(inj_result$message)) inj_result$message else NA_character_
  
  t_ij_val <- if(inj_success && !is.null(inj_result$diagnostics$t_ij)) {
    inj_result$diagnostics$t_ij
  } else {
    NA_real_
  }

  rowt$t_ij <- t_ij_val

  if (!inj_success) {
    # Failure result with no sensitivity analysis
    base_result <- data.table(
      triplet_id = idx,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "positive",
      reduction_pct = 0,
      N = rowt$N,
      dynamic = rowt$dynamic,
      fold_change = rowt$fold_change,
      t_ij = t_ij_val,
      model_success = FALSE,
      injection_success = FALSE,
      n_injected = n_injected_val,
      n_coadmin = n_coadmin_val,
      n_events = NA_integer_,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = FALSE,
      log_ior_classic = list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = list(rep(NA_real_, 7)),
      ior_classic = list(rep(NA_real_, 7)),
      reri_classic_success = FALSE,
      reri_values = list(rep(NA_real_, 7)),
      reri_lower90 = list(rep(NA_real_, 7)),
      reri_upper90 = list(rep(NA_real_, 7)),
      n_stages_reri_significant = NA_integer_,
      RERI_classic = list(rep(NA_real_, 7)),
      RERI_classic_lower90 = list(rep(NA_real_, 7)),
      RERI_classic_upper90 = list(rep(NA_real_, 7)),
      RERI_classic_se = list(rep(NA_real_, 7)),
      diagnostics = diag_data,
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = NA_character_,
      error_msg = inj_message
    )
    rm(inj_result); gc(verbose = FALSE)
    return(base_result)
  }

  # Sensitivity loop
  all_results <- list()
  
  # Base result (0% reduction)
  base_result <- fit_reduced_model(inj_result$ade_aug, rowt, 0)
  base_result$n_injected <- n_injected_val
  base_result$injection_success <- TRUE
  base_result$diagnostics <- diag_data
  
  all_results[[1]] <- base_result
  
  # Loop over reduction levels
  for (red_pct in reduction_levels) {
    # Reduce dataset
    ade_reduced <- reduce_dataset_by_stage(
      inj_result$ade_aug, 
      red_pct,
      seed = base_seed + idx)   # Seed propagation for greater reproducibility
    
    # Fit model on reduced dataset
    reduced_result <- fit_reduced_model(ade_reduced, rowt, red_pct)
    reduced_result$n_injected <- n_injected_val
    reduced_result$injection_success <- TRUE
    reduced_result$diagnostics <- diag_data
    
    all_results[[length(all_results) + 1]] <- reduced_result
    
    rm(ade_reduced); gc(verbose = FALSE)
  }
  
  rm(inj_result); gc(verbose = FALSE)
  
  # Combine all sensitivity results
  combined_results <- rbindlist(all_results, fill = TRUE)
  
  return(combined_results)
}

################################################################################
# Helper function: Reduce dataset by stage
################################################################################

# Reduces an augmented dataset by randomly removing a percentage of rows per stage
# 
# Parameters:
# ade_aug: augmented data.table
# reduction_pct: percentage to remove (e.g. 10 for 10%)
# nichd_col: column containing the NICHD stage
# seed: provided for reproducibility
#
# Return:
# Reduced data.table

reduce_dataset_by_stage <- function(ade_aug, reduction_pct, nichd_col = "nichd", seed = NULL) {
  
  # For each stage, remove the specified percentage of rows
  stages <- unique(ade_aug[[nichd_col]])
  
  reduced_rows <- lapply(stages, function(stage) {
    stage_data <- ade_aug[get(nichd_col) == stage]
    n_rows <- nrow(stage_data)
    n_keep <- ceiling(n_rows * (1 - reduction_pct / 100))
    
    if (n_keep < 1 && n_rows > 0) {
      n_keep <- 1  # Keep at least one row if data is present
    }
    
    if (n_rows > 0) {
      # Deterministic seed based on stage + reduction_pct + dimensions
      if (!is.null(seed)) {
        local_seed <- seed + as.integer(stage) * 100 + as.integer(reduction_pct) * 10
        set.seed(local_seed)
      }
      stage_data[sample(.N, min(n_keep, n_rows))]
    } else {
      stage_data
    }
  })
  
  rbindlist(reduced_rows, use.names = TRUE)
}

################################################################################
# Model fitting function on a reduced dataset
################################################################################

# Fits GAM and classical IOR/RERI models on a reduced dataset
# Wrapper function for model fitting and dataset reduction
#
# Parameters:
# ade_reduced: reduced dataset
# rowt: triplet metadata (drugA, drugB, meddra, etc.)
# reduction_pct: reduction percentage applied
# 
# Return:
# data.table with fitting results

fit_reduced_model <- function(ade_reduced, rowt, reduction_pct) {
  
  counts_reduced <- calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra)
  
  # GAM fitting
  model_res <- tryCatch({
    fit_gam(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced,
      spline_individuales = spline_individuales,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      nichd_spline = nichd_spline
    )
  }, error = function(e) {
    list(
      success = FALSE, 
      n_events_total = counts_reduced$n_events_total, 
      n_coadmin = counts_reduced$n_coadmin,
      error_msg = paste("Error en modelo reducido:", e$message)
    )
  })
  
  # Classical IOR computation
  classic_res <- tryCatch({
    calculate_classic_ior(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
  }, error = function(e) {
    list(success = FALSE)
  })
  
  # Classical RERI computation
  classic_reri <- tryCatch({
    calculate_classic_reri(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
  }, error = function(e) {
    list(success = FALSE)
  })
  
  # Prepare result
  if (!model_res$success) {
    result <- data.table(
      triplet_id = rowt$triplet_id,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = rowt$type,
      reduction_pct = reduction_pct,
      N = counts_reduced$n_events_coadmin,
      dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
      fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
      t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
      model_success = FALSE,
      injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
      n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
      n_coadmin = counts_reduced$n_coadmin,
      n_events = counts_reduced$n_events,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = classic_res$success,
      log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
      ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
      reri_classic_success = classic_reri$success,
      reri_values = list(rep(NA_real_, 7)),
      reri_lower90 = list(rep(NA_real_, 7)),
      reri_upper90 = list(rep(NA_real_, 7)),
      n_stages_reri_significant = NA_integer_,
      RERI_classic = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic) else list(rep(NA_real_, 7)),
      RERI_classic_lower90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_lower90) else list(rep(NA_real_, 7)),
      RERI_classic_upper90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_upper90) else list(rep(NA_real_, 7)),
      RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7)),
      diagnostics = list(list(error = model_res$error_msg)),
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = if(!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
      error_msg = if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_
    )
    return(result)
  }
  
  # Successful result
  result <- data.table(
    triplet_id = rowt$triplet_id,
    drugA = rowt$drugA,
    drugB = rowt$drugB,
    meddra = rowt$meddra,
    type = rowt$type,
    reduction_pct = reduction_pct,
    N = model_res$n_events_coadmin,
    dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
    fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
    t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
    model_success = TRUE,
    injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
    n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
    n_coadmin = model_res$n_coadmin,
    n_events = model_res$n_events,
    n_stages_significant = model_res$n_stages_significant,
    max_ior = model_res$max_ior,
    mean_ior = model_res$mean_ior,
    model_aic = model_res$model_aic,
    stage = list(1:7),
    log_ior = list(model_res$log_ior),
    log_ior_lower90 = list(model_res$log_ior_lower90),
    ior_values = list(model_res$ior_values),
    classic_success = classic_res$success,
    log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
    log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
    ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
    reri_classic_success = classic_reri$success,
    reri_values = list(model_res$reri_values),
    reri_lower90 = list(model_res$reri_lower90),
    reri_upper90 = list(model_res$reri_upper90),
    n_stages_reri_significant = model_res$n_stages_reri_significant,
    RERI_classic = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic) else list(rep(NA_real_, 7)),
    RERI_classic_lower90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_lower90) else list(rep(NA_real_, 7)),
    RERI_classic_upper90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_upper90) else list(rep(NA_real_, 7)),
    RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7)),
    diagnostics = list(list()),
    spline_individuales = spline_individuales,
    nichd_spline = nichd_spline,
    include_sex = include_sex,
    include_stage_sex = include_stage_sex,
    k_spline = k_spline,
    bs_type = bs_type,
    select = select,
    formula_used = model_res$formula_used
  )
  
  return(result)
}

################################################################################
# Permutation function for the null distribution
################################################################################

# Permutes events across reports
# 
# Parameters:
# perm_events: if TRUE, permutes events (meddra_concept_id)
# perm_drugs: if TRUE, also permutes drugs across reports
#
# Uses the pool of reports selected in 10_augmentation (pool_meta)
# Breaks the drug-event association to generate a negative ground truth

permute_pool <- function(pool_meta, niveles_nichd, 
                         perm_events = TRUE, 
                         perm_drugs = FALSE, 
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  pool_copy <- copy(pool_meta)
  pool_copy[, drugs_perm := vector("list", .N)]
  pool_copy[, events_perm := vector("list", .N)]
  
  for (stage in niveles_nichd) {
    idx <- which(pool_copy$nichd == stage)
    
    if (length(idx) <= 1) {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
      next
    }
    
    if (perm_drugs) {
      perm_idx_drugs <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx[perm_idx_drugs]]
    } else {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
    }
    
    if (perm_events) {
      perm_idx_events <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$events_perm[idx] <- pool_copy$events[idx[perm_idx_events]]
    } else {
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
    }
  }
  
  pool_copy[, .(safetyreportid, nichd, nichd_num, drugs_perm, events_perm)]
}

# Function to reintroduce permuted reports into the original dataset before model fitting

reintroduce_permuted_reports <- function(ade_original, permuted_pool) {
  
  pool_report_ids <- unique(permuted_pool$safetyreportid)
  ade_without_pool <- ade_original[!safetyreportid %in% pool_report_ids]
  
  permuted_rows <- permuted_pool[, {
    drugs_vec <- drugs_perm[[1]]
    events_vec <- events_perm[[1]]
    
    if (length(drugs_vec) > 0 && length(events_vec) > 0) {
      CJ(safetyreportid = safetyreportid,
         atc_concept_id = drugs_vec,
         meddra_concept_id = events_vec,
         nichd = nichd,
         nichd_num = nichd_num)
    } else {
      data.table()
    }
  }, by = safetyreportid]
  
  rbindlist(list(ade_without_pool, permuted_rows), use.names = TRUE, fill = TRUE)
}

################################################################################
# Expansion function
################################################################################

# Expands triplet results to long format with all metrics
#
# Parameters:
# dt: data.table with triplet results (with list-type columns)
# label_val: classification label (1 = positive, 0 = negative)
# null_thresholds_dt: data.table with null distribution thresholds per stage
# use_threshold_ior: if TRUE, applies the IOR null distribution threshold
# use_threshold_reri: if TRUE, applies the RERI null distribution threshold
# 
# Return:
# Expanded data.table with one row per triplet-stage
# 
# Implementation:
# Unpacks list-type columns (stage, log_ior, etc.) into individual rows
# Merges with null distribution thresholds
# Retains triplet-level metadata (dynamic, t_ij, n_coadmin)

expand_clean_all_metrics <- function(dt, label_val, null_thresholds_dt,
                                     use_threshold_ior = TRUE, use_threshold_reri = TRUE) {

  has_dynamic <- "dynamic" %in% names(dt)
  by_cols <- "triplet_id"
  if ("dynamic" %in% names(dt)) {
    by_cols <- c(by_cols, "dynamic", "t_ij")
  }
  if ("n_coadmin" %in% names(dt)) {  
    by_cols <- c(by_cols, "n_coadmin")
  }
  if ("N" %in% names(dt)) {       
    by_cols <- c(by_cols, "N")
  }

  expanded <- dt[, {
    stages <- unlist(stage)
    
    # GAM
    gam_log_ior <- unlist(log_ior)
    gam_log_ior_lower90 <- unlist(log_ior_lower90)
    gam_reri <- unlist(reri_values)
    gam_reri_lower90 <- unlist(reri_lower90)
    
    # Classical stratified metrics
    cls_log_ior <- unlist(log_ior_classic)
    cls_log_ior_lower90 <- unlist(log_ior_classic_lower90)
    cls_reri <- unlist(RERI_classic)
    cls_reri_lower90 <- unlist(RERI_classic_lower90)
    
    n <- min(length(stages), length(gam_log_ior), length(gam_log_ior_lower90),
             length(gam_reri), length(gam_reri_lower90),
             length(cls_log_ior), length(cls_log_ior_lower90),
             length(cls_reri), length(cls_reri_lower90))
    
    if (n > 0) {
      data.table(
        stage_num = stages[1:n],
        # GAM
        gam_log_ior = gam_log_ior[1:n],
        gam_log_ior_lower90 = gam_log_ior_lower90[1:n],
        gam_reri = gam_reri[1:n],
        gam_reri_lower90 = gam_reri_lower90[1:n],
        # Stratified
        classic_log_ior = cls_log_ior[1:n],
        classic_log_ior_lower90 = cls_log_ior_lower90[1:n],
        classic_reri = cls_reri[1:n],
        classic_reri_lower90 = cls_reri_lower90[1:n]
      )
    }
  }, by = by_cols]
  
  if (!has_dynamic) {
    expanded[, `:=`(dynamic = "control", t_ij = 0)]
  }
  
  expanded[, nichd := niveles_nichd[stage_num]]
  expanded[, label := label_val]
  
  # Merge with null distribution thresholds
  expanded <- merge(expanded, null_thresholds_dt, 
                   by.x = "stage_num", by.y = "stage", all.x = TRUE)
  
  # Store threshold usage parameters
  expanded[, `:=`(
    use_threshold_ior = use_threshold_ior,
    use_threshold_reri = use_threshold_reri
  )]
  return(expanded)
}

################################################################################
# Statistical power calculation function — GAM method
################################################################################

# Computes statistical power by filtering on effect size and co-administration report count
#
# Parameters:
# data_pos: data.table with expanded positive triplets
# data_neg: data.table with expanded negative triplets (not used)
# target_power: target power level (0.80)
# metric_n: column name to use as count filter ("n_coadmin" or "n_events")
# grid_resolution: number of steps for the search grid (default 30x30)
# use_threshold_ior: if TRUE, applies the IOR null distribution threshold
# use_threshold_reri: if TRUE, applies the RERI null distribution threshold
# detection: detection mode ("log-ior", "reri", or "double")
#
# Return:
# power_surface: data.table with t_threshold, n_threshold, power, len
# t_star: optimal t_ij threshold (effect size in the injection formula)
# n_star: optimal n_coadmin threshold
# superset_pos: positive triplets passing the optimal filters
# achieved_power: power achieved at the optimal point
#
# Implementation:
# Generates a 2D grid (t_ij × n_coadmin) using quantiles
# For each grid point, computes statistical power (TP / (TP + FN))
# Identifies the optimal configuration achieving target_power with maximum retention
# Aggregates to triplet level: detected if ANY stage detects

calculate_power_gam <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,  
  use_threshold_ior = TRUE,   
  use_threshold_reri = TRUE,   
  detection = "double"  ) {    # "ior", "reri", or "double"
  
    library(data.table)

    ###########  
    # 1. Prepare positive data
    ###########
    pos_clean <- copy(data_pos)
  
    # Basic cleaning
    pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
    # Merge thresholds if provided and specific columns are missing
    if (!is.null(null_thresholds)) {
      # Check whether the threshold columns we need are missing
      need_merge <- FALSE
      if (use_threshold_ior && !"threshold_ior" %in% names(pos_clean)) need_merge <- TRUE
      if (use_threshold_reri && !"threshold_reri" %in% names(pos_clean)) need_merge <- TRUE
      if (need_merge) {
        # Ensure null_thresholds has a 'stage' column
        if (!"stage" %in% names(null_thresholds)) {
          stop("null_thresholds debe tener una columna llamada 'stage'")
        }
        pos_clean <- merge(
          pos_clean,
          null_thresholds,
          by.x = "stage_num",
          by.y = "stage",
          all.x = TRUE
        )
      }
    }
  
  ###########
  # 2. 2D grid (t_ij x n_coadmin)
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Power calculation
  ###########
  
  # Pre-compute detection flags according to the "detection" parameter
  
  # IOR criterion
  if (detection %in% c("ior", "double")) {
    if (use_threshold_ior) {
      if (!"threshold_ior" %in% names(pos_clean)) {
        stop("columna threshold_ior no existe en los datos ni en null_thresholds")
      }
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & 
          gam_log_ior_lower90 > 0 & 
          gam_log_ior_lower90 > threshold_ior
      )]
    } else {
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & gam_log_ior_lower90 > 0
      )]
    }
  }
  
  # RERI criterion
  if (detection %in% c("reri", "double")) {
    if (use_threshold_reri) {
      if (!"threshold_reri" %in% names(pos_clean)) {
        stop("columna threshold_reri no existe en los datos ni en null_thresholds")
      }
      pos_clean[, reri_detected := (
        !is.na(gam_reri_lower90) & 
          gam_reri_lower90 > 0 & 
          gam_reri_lower90 > threshold_reri
      )]
    } else {
      pos_clean[, reri_detected := (
        !is.na(gam_reri_lower90) & gam_reri_lower90 > 0
      )]
    }
  }
  
  # Final detection flag based on the chosen mode
  if (detection == "ior") {
    pos_clean[, is_detected := ior_detected]
  } else if (detection == "reri") {
    pos_clean[, is_detected := reri_detected]
  } else {  # "double"
    pos_clean[, is_detected := ior_detected | reri_detected]
  }

  # Aggregate to triplet level: detected if ANY stage detects
  triplet_detection_gam <- pos_clean[, .(
    triplet_detected = any(is_detected, na.rm = TRUE),
    t_ij_triplet = unique(t_ij)[1],
    n_coadmin_triplet = unique(get(metric_n))[1]
  ), by = triplet_id]
  
  # Iterate over the search grid
  power_surface <- search_grid[, {
    
    # Triplet-level filter: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      triplet_detection_gam$t_ij_triplet >= t_threshold & 
      triplet_detection_gam$n_coadmin_triplet >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(triplet_detection_gam$triplet_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identify optimal point
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("No se alcanzó poder objetivo (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # Informational message based on detection mode
  detection_label <- switch(detection,
    "ior" = "solo IOR",
    "reri" = "solo RERI", 
    "double" = "IOR O RERI"
  )
  
  message(sprintf("\nUMBRALES ÓPTIMOS (GAM - %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Poder alcanzado: %.1f%%", achieved_power * 100))
  message(sprintf("  Tripletes retenidos: %d / %d (%.1f%%)",
                  n_retained,
                  uniqueN(triplet_detection_gam$triplet_id),
                  100 * n_retained / uniqueN(triplet_detection_gam$triplet_id)))
  
  ###########
  # 5. Build supersets
  ###########
  
  # Positive superset: ALL observations from triplets passing the filters
  triplets_passed <- triplet_detection_gam[
    t_ij_triplet >= t_star & n_coadmin_triplet >= n_star,
    triplet_id
  ]
  
  superset_pos <- pos_clean[triplet_id %in% triplets_passed]
  
  power_surface[, method := paste0("GAM-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "gam",
    detection_mode = detection,
    metric_n_used = metric_n
  ))
}

################################################################################
# Statistical power calculation function — Classical stratified method
################################################################################

# Computes statistical power by filtering on effect size and co-administration report count
#
# Classical methods are more prone to producing NAs or infinite CIs, so NA handling is parameterized
#
# Parameters:
# data_pos: data.table with expanded positive triplets
# data_neg: data.table with expanded negative triplets
# target_power: target power level (0.80)
# metric_n: column name to use as count filter ("n_coadmin" or "n_events")
# grid_resolution: number of steps for the search grid (default 30x30)
# detection: criterion for the superset — IOR, RERI, or both
# na_remove: controls NA handling.
#  TRUE: superset includes ALL triplets meeting minimum detection criteria
#  FALSE: superset excludes NAs even if they meet the threshold criteria
#
# Return:
# power_surface: data.table with t_threshold, n_threshold, power, len
# t_star: optimal t_ij threshold
# n_star: optimal n_coadmin threshold
# superset_pos: positive triplets passing the optimal filters
# achieved_power: power achieved at the optimal point
#
# Implementation:
# See calculate_power_gam
# Key difference: power is computed at the stage level

calculate_power_classic <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,
  detection = "double",  # "ior", "reri", or "double"
  na_remove = "TRUE") {
    
  ###########
  # 1. Prepare positive data
  ###########

  detection <- match.arg(detection, choices = c("ior", "reri", "double"))

  pos_clean <- copy(data_pos)
  
  n_triplets_original_total <- uniqueN(pos_clean$triplet_id)
  n_obs_original_total <- nrow(pos_clean)

  # Remove NAs according to the detection criterion  
  if (na_remove) {
    n_obs_before_na <- nrow(pos_clean)
    n_triplets_before_na <- uniqueN(pos_clean$triplet_id)
    if (detection == "ior") {
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90)]
    } else if (detection == "reri") {
      pos_clean <- pos_clean[!is.na(classic_reri_lower90)]
    } else {  # double
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90) & !is.na(classic_reri_lower90)]
    }
  }
  # Log how many were removed
  n_obs_after_na <- nrow(pos_clean)
  n_triplets_after_na <- uniqueN(pos_clean$triplet_id)
  n_triplets_lost_na <- n_triplets_before_na - n_triplets_after_na
 
  # Basic cleaning
  pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
  if (!is.null(null_thresholds) && !"threshold" %in% names(pos_clean)) {
    pos_clean <- merge(
      pos_clean,
      null_thresholds,
      by.x = "stage_num",
      by.y = "stage",
      all.x = TRUE
    )
  }
  
  ###########
  # 2. 2D search grid
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Power calculation
  ###########
  
  # Each row (triplet_id + stage) is an independent observation
  # because the method is stratified by stage
  
  # Pre-compute detection flags according to the 'detection' parameter
  
  if (detection == "ior") {
    # IOR only: 90% CI lower bound > 0
    pos_clean[, is_detected := (
      !is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0
    )]
    
  } else if (detection == "reri") {
    # RERI only: 90% CI lower bound > 0
    pos_clean[, is_detected := (
      !is.na(classic_reri_lower90) & classic_reri_lower90 > 0
    )]
    
  } else {  # detection == "double"
    # Double criterion: IOR OR RERI, 90% CI lower bound > 0
    pos_clean[, is_detected := (
      (!is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0) |
      (!is.na(classic_reri_lower90) & classic_reri_lower90 > 0)
    )]
  }
    
  # Save a copy with detection flags
  pos_all_with_detection <- copy(pos_clean)

  # Computed at the observation level (row): t_ij, n_coadmin, detection
  # Not aggregated to triplet level
  
  # Iterate over the search grid
  power_surface <- search_grid[, {
    
    # Observation-level filter: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      pos_clean$t_ij >= t_threshold & 
        pos_clean[[metric_n]] >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(pos_clean$is_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identify optimal point
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("No se alcanzó poder objetivo (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # Detection mode for the message
  detection_label <- switch(detection,
    "ior" = "solo IOR",
    "reri" = "solo RERI", 
    "double" = "IOR O RERI"
  )

  ###########
  # 5. Build supersets
  ###########
  
  # Positive superset: observations passing the filters
  superset_pos <- pos_all_with_detection[t_ij >= t_star & get(metric_n) >= n_star]
  
  # Superset metrics
  n_retained_total <- nrow(superset_pos)
  n_detected_total <- sum(superset_pos$is_detected, na.rm = TRUE)
  power_total <- ifelse(n_retained_total > 0, n_detected_total / n_retained_total, 0)

  # Count of unique triplets meeting the retention criterion
  n_triplets_retained <- uniqueN(superset_pos$triplet_id)
  n_triplets_total <- uniqueN(pos_all_with_detection$triplet_id)

  n_triplets_retained_vs_original <- n_triplets_retained
  pct_retained_vs_original <- 100 * n_triplets_retained_vs_original / n_triplets_original_total

  message(sprintf("\nUmbrales (Método Estratificado, %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Poder alcanzado: %.1f%%", achieved_power * 100))
  message(sprintf("  Tripletes retenidos: %d / %d (%.1f%%)",
                  n_triplets_retained, n_triplets_total,
                  100 * n_triplets_retained / n_triplets_total))
  message(sprintf("  Poder en muestra completa: %.1f%%", power_total * 100))

  # Lines to assess true retention, accounting for NAs excluded at the start of the function
  message(sprintf("  Tripletes retenidos (vs total original): %d / %d (%.1f%%) [INCLUYE NA]",
                  n_triplets_retained_vs_original, n_triplets_original_total,
                  pct_retained_vs_original))
  
  power_surface[, method := paste0("Estratificado-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "classic",
    metric_n_used = metric_n,
    detection_mode = detection,
    na_remove = na_remove,
    achieved_power_total = power_total,  # power in full sample
    n_retained_grid = n_retained,  # retained in grid
    n_retained_total = n_retained_total,  # total retained
    n_detected_total = n_detected_total,  # total detected
    n_triplets_retained = n_triplets_retained,
    n_triplets_total = n_triplets_total
  ))
}

################################################################################
# Function to visualize the power surface over the full grid
################################################################################

# Generates a power surface heatmap with interpolation to a regular grid
# Faceted by method
#
# Parameters:
# power_result: list returned by calculate_power_gam() or calculate_power_classic()
# target_power: target power level (for visual reference)
# detection: detection type ("IOR", "RERI", or "double")
# t_range: X-axis range (effect size)
# n_range: Y-axis range (co-administration report count)
# grid_size: interpolated grid size (default: 50x50 for smoothness)
# 
# Return:
# ggplot object with the full interpolated surface
#
# Implementation:
# Removes NAs and Inf values
# Interpolates to a regular grid using akima::interp (for missing data)
# Generates heatmap using geom_raster
# Returns optimal point statistics

plot_power_surface <- function(
  power_results_list,
  facet_by = "method", 
  target_power = 0.80, 
  detection = detection,
  t_range = c(0, 0.5),      
  n_range = c(0, 300),
  grid_size = 30) 
  {
  
  library(ggplot2)
  library(data.table)
  library(scales)
  library(akima)
  
  ###########
  # 1- Process multiple surfaces
  ###########
  
  all_surfaces <- rbindlist(lapply(names(power_results_list), function(met_name) {
    surface <- as.data.table(power_results_list[[met_name]]$power_surface)
    surface[, method_label := met_name]
    return(surface)
  }), fill = TRUE)
  
  # Extract optimal parameters for the subtitle
  opt_params <- lapply(names(power_results_list), function(met_name) {
    res <- power_results_list[[met_name]]
    data.table(
      method_label = met_name,
      t_star = res$t_star,
      n_star = res$n_star,
      achieved_power = res$achieved_power
    )
  })
  opt_params_dt <- rbindlist(opt_params)
  
  ###########
  # 2- Clean original data
  ###########
  
  # Remove NAs and non-finite values
  surface_clean <- all_surfaces[is.finite(t_threshold) & is.finite(n_threshold) & is.finite(power)]
  
  ###########
  # 3- Interpolate to a full regular grid
  ###########
  
  # Regular sequences with the same range for X and Y (to keep a square aspect)
  # Use the wider range across both axes to maintain proportion
  common_range <- c(0, max(t_range[2], n_range[2]))

  # Create regular grid sequences
  t_seq <- seq(t_range[1], t_range[2], length.out = grid_size)
  n_seq <- seq(n_range[1], n_range[2], length.out = grid_size)
  
  surfaces_interp <- lapply(unique(surface_clean$method_label), function(met) {
    surface_met <- surface_clean[method_label == met]
    
    interp_result <- tryCatch({
      interp(
        x = surface_met$t_threshold,
        y = surface_met$n_threshold,
        z = surface_met$power,
        xo = t_seq,
        yo = n_seq,
        linear = TRUE,
        extrap = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(interp_result)) return(NULL)
    
    dt <- data.table(expand.grid(
      t_threshold = interp_result$x,
      n_threshold = interp_result$y
    ))
    dt[, power := as.vector(interp_result$z)]
    dt[, method_label := met]
    return(dt)
  })
  
  surface_plot <- rbindlist(surfaces_interp, fill = TRUE)
  surface_plot <- surface_plot[!is.na(power)]
  surface_plot[, power := pmin(pmax(power, 0), 1)]
  
  ###########
  # 5- Method labels
  ###########

  # Subtitle with per-method info
  subtitulo <- paste(opt_params_dt[, sprintf("%s: t≥%.4f, N≥%.0f (%.1f%%)", method_label, t_star, n_star, achieved_power*100)], collapse = " | ")

  method_label_clean <- switch(detection,
    "IOR" = "Detección IOR",
    "RERI" = "Detección RERI",
    detection
  )
  
  ###########
  # 6- Build the plot
  ###########
  
  p <- ggplot(surface_plot, aes(x = t_threshold, y = n_threshold, fill = power)) +
    
    # geom_raster() for regular grids
    geom_raster(interpolate = FALSE) +  # interpolate=TRUE smooths visually

    # Faceting
    facet_wrap(~ method_label, ncol = 2, scales = "fixed") +  

    # Color scale
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1), 
      name = "Poder",
      option = "plasma",
      na.value = "white",
     guide = guide_colorbar(theme = theme( legend.text = element_text(angle = 45, hjust = 1, vjust = 1)))
    ) +
    
    # Axis scales
    scale_x_continuous(
      limits = t_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 5)
    ) +
    
    scale_y_continuous(
      limits = n_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 6)
    ) +

    # Square coordinate system
    coord_fixed(ratio = diff(t_range)/diff(n_range)) +

    # Labels
    labs(
      title = sprintf("Superficie de Poder - %s", method_label_clean),
      subtitle = sprintf("Objetivo: %.0f%% | %s", target_power * 100, subtitulo),
      x = expression("Tamaño de Efecto " * (t[italic(ij)])),
      y = "Número de reportes A-B",
      caption = sprintf("Grilla: %dx%d celdas", grid_size, grid_size)
    ) 
  return(p)
}

################################################################################
# Bootstrap-based metrics calculation function
################################################################################

# Computes classification metrics with 95% bootstrap confidence intervals
# 
# Parameters:
# dt: data.table with columns: triplet_id, detected, label
# n_boot: number of bootstrap replications
# aggregate_triplet: if TRUE, aggregates to triplet level using any(); if FALSE, keeps current granularity (stage or dynamic)
# 
# Return:
# data.table with metrics and 95% CIs

calculate_metrics <- function(dt, n_boot = 2000, aggregate_triplet = TRUE, score_type, score_type_auc) {
  
  # Aggregate to triplet level if requested
  if (aggregate_triplet) {
    dt_eval <- dt[, .(
      detected = any(detected, na.rm = TRUE),
      label = unique(label),
      score = max(get(score_type), na.rm = TRUE),
      score_auc = max(get(score_type_auc), na.rm = TRUE)
    ), by = triplet_id]
  } else {
    dt_eval <- dt[, .(triplet_id, detected, label, score = get(score_type), score_auc = get(score_type_auc))]
  }
  
  n_pos <- sum(dt_eval$label == 1)
  n_neg <- sum(dt_eval$label == 0)
  n_total <- nrow(dt_eval)
  
  auc_result <- tryCatch({
    # Filter NA and infinite values before computing ROC
    dt_roc <- dt_eval[is.finite(score_auc) & !is.na(score_auc) & !is.na(label)]
    
    if (nrow(dt_roc) < 2 || length(unique(dt_roc$label)) < 2) {
    stop("Datos insuficientes")
    }
    
    roc_data <- pROC::roc(
      response = dt_roc$label,
      predictor = dt_roc$score_auc,
      quiet = TRUE,
      direction = "<"
    )
    auc <- as.numeric(pROC::auc(roc_data))
    
    # Bootstrap for AUC
    b_auc <- replicate(n_boot, {
      b_idx <- sample(nrow(dt_roc), replace = TRUE)
      b_dt <- dt_roc[b_idx]
      
      # Check that both classes are present in the bootstrap sample
      if (length(unique(b_dt$label)) < 2) {
        return(NA_real_)
      }
      
      roc_boot <- pROC::roc(
        response = b_dt$label,
        predictor = b_dt$score_auc,
        quiet = TRUE,
        direction = "<"
      )
      as.numeric(pROC::auc(roc_boot))
    })
    
    list(
      auc = auc,
      auc_lower = unname(quantile(b_auc, 0.025, na.rm = TRUE)),
      auc_upper = unname(quantile(b_auc, 0.975, na.rm = TRUE))
    )
  }, error = function(e) {
    warning(sprintf("Error en cálculo de AUC: %s", e$message))
    list(auc = NA_real_, auc_lower = NA_real_, auc_upper = NA_real_)
  })
  
  # Bootstrap
  if (aggregate_triplet) {
    # Identify unique IDs for triplet-level bootstrap
    unique_pos_ids <- unique(dt_eval[label == 1, triplet_id])
    unique_neg_ids <- unique(dt_eval[label == 0, triplet_id])
  
    boot_stats <- replicate(n_boot, {
      boot_pos_ids <- sample(unique_pos_ids, length(unique_pos_ids), replace = TRUE)
      boot_neg_ids <- sample(unique_neg_ids, length(unique_neg_ids), replace = TRUE)
      ids_dt <- data.table(triplet_id = c(boot_pos_ids, boot_neg_ids))
      boot_dt <- dt_eval[ids_dt, on = "triplet_id", nomatch = NULL]
    
      b_tp <- sum(boot_dt$detected & boot_dt$label == 1)
      b_fn <- sum(!boot_dt$detected & boot_dt$label == 1)
      b_fp <- sum(boot_dt$detected & boot_dt$label == 0)
      b_tn <- sum(!boot_dt$detected & boot_dt$label == 0)
    
      b_sens <- ifelse((b_tp + b_fn) > 0, b_tp / (b_tp + b_fn), 0)
      b_spec <- ifelse((b_tn + b_fp) > 0, b_tn / (b_tn + b_fp), 0)
      b_ppv <- ifelse((b_tp + b_fp) > 0, b_tp / (b_tp + b_fp), 0)
      b_npv <- ifelse((b_tn + b_fn) > 0, b_tn / (b_tn + b_fn), 0)
      b_acc <- (b_tp + b_tn) / nrow(boot_dt)
      b_f1 <- ifelse((b_tp + b_fp) > 0 && (b_tp + b_fn) > 0, 2 * b_tp / (2 * b_tp + b_fp + b_fn), 0)
    
      c(b_sens, b_spec, b_ppv, b_npv, b_acc, b_f1, b_tp, b_fn, b_fp, b_tn)
    })

  } else {
    # Row-level bootstrap for stage-level analysis
    pos_idx <- which(dt_eval$label == 1)
    neg_idx <- which(dt_eval$label == 0)

    boot_stats <- replicate(n_boot, {
      boot_pos_idx <- sample(pos_idx, n_pos, replace = TRUE)
      boot_neg_idx <- sample(neg_idx, n_neg, replace = TRUE)
      boot_idx <- c(boot_pos_idx, boot_neg_idx)
      boot_dt <- dt_eval[boot_idx]
    
      b_tp <- sum(boot_dt$detected & boot_dt$label == 1)
      b_fn <- sum(!boot_dt$detected & boot_dt$label == 1)
      b_fp <- sum(boot_dt$detected & boot_dt$label == 0)
      b_tn <- sum(!boot_dt$detected & boot_dt$label == 0)
    
      b_sens <- ifelse((b_tp + b_fn) > 0, b_tp / (b_tp + b_fn), 0)
      b_spec <- ifelse((b_tn + b_fp) > 0, b_tn / (b_tn + b_fp), 0)
      b_ppv <- ifelse((b_tp + b_fp) > 0, b_tp / (b_tp + b_fp), 0)
      b_npv <- ifelse((b_tn + b_fn) > 0, b_tn / (b_tn + b_fn), 0)
      b_acc <- (b_tp + b_tn) / nrow(boot_dt)
      b_f1 <- ifelse((b_tp + b_fp) > 0 && (b_tp + b_fn) > 0, 2 * b_tp / (2 * b_tp + b_fp + b_fn), 0)
    
      c(b_sens, b_spec, b_ppv, b_npv, b_acc, b_f1, b_tp, b_fn, b_fp, b_tn)
    })
  }
  
  # Point estimates as bootstrap means (to avoid misalignment with CIs)
  sens_boot <- boot_stats[1, ]
  spec_boot <- boot_stats[2, ]
  ppv_boot <- boot_stats[3, ]
  npv_boot <- boot_stats[4, ]
  acc_boot <- boot_stats[5, ]
  f1_boot <- boot_stats[6, ]
  
  sens <- mean(sens_boot, na.rm = TRUE)
  spec <- mean(spec_boot, na.rm = TRUE)
  ppv <- mean(ppv_boot, na.rm = TRUE)
  npv <- mean(npv_boot, na.rm = TRUE)
  acc <- mean(acc_boot, na.rm = TRUE)
  f1 <- mean(f1_boot, na.rm = TRUE)
  
  # CI calculation
  calc_ci <- function(valores_boot) {
    c(
      point = mean(valores_boot, na.rm = TRUE),
      lower = unname(quantile(valores_boot, 0.025, na.rm = TRUE)),
      upper = unname(quantile(valores_boot, 0.975, na.rm = TRUE))
    )
  }
  
  sens_ci <- calc_ci(sens_boot)
  spec_ci <- calc_ci(spec_boot)
  ppv_ci <- calc_ci(ppv_boot)
  npv_ci <- calc_ci(npv_boot)
  acc_ci <- calc_ci(acc_boot)
  f1_ci <- calc_ci(f1_boot)
  
  # Counts from the original (non-bootstrap) dataset for reference
  tp_orig <- sum(dt_eval$detected & dt_eval$label == 1)
  fn_orig <- sum(!dt_eval$detected & dt_eval$label == 1)
  fp_orig <- sum(dt_eval$detected & dt_eval$label == 0)
  tn_orig <- sum(!dt_eval$detected & dt_eval$label == 0)
  
  data.table(
    sensitivity = sens_ci["point"], sensitivity_lower = sens_ci["lower"], sensitivity_upper = sens_ci["upper"],
    specificity = spec_ci["point"], specificity_lower = spec_ci["lower"], specificity_upper = spec_ci["upper"],
    PPV = ppv_ci["point"], PPV_lower = ppv_ci["lower"], PPV_upper = ppv_ci["upper"],
    NPV = npv_ci["point"], NPV_lower = npv_ci["lower"], NPV_upper = npv_ci["upper"],
    Accuracy = acc_ci["point"],
    F1 = f1_ci["point"], F1_lower = f1_ci["lower"], F1_upper = f1_ci["upper"],
    TP = tp_orig, FN = fn_orig, FP = fp_orig, TN = tn_orig,
    n_positives = n_pos, n_negatives = n_neg, n_total = n_total,
    AUC = auc_result$auc,
    AUC_lower = auc_result$auc_lower,
    AUC_upper = auc_result$auc_upper
  )
}

################################################################################
# Data expansion function with reduction
################################################################################

# Expands data to long format
#
# Implementation:
# Loads data and creates objects according to the reduction level
# Filters out failed injections
# Expands to long format
# Filters to stages with high reporting according to the dynamic

expand <- function(red_pct) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)  # suffix based on reduction level
  
  ruta_pos <- paste0(ruta_base_sensitivity, "positive_triplets_results", suffix_file, ".rds")
  ruta_neg <- paste0(ruta_base_sensitivity, "negative_triplets_results", suffix_file, ".rds")
  
  pos_raw <- readRDS(ruta_pos)
  neg_raw <- readRDS(ruta_neg)
  
  # Successful injections only
  pos_valid <- pos_raw[injection_success == TRUE]
  
  # Expand to long format
  pos_exp <- expand_clean_all_metrics(pos_valid, 1, null_thresholds, use_threshold_ior, use_threshold_reri)
  neg_exp <- expand_clean_all_metrics(neg_raw, 0, null_thresholds, use_threshold_ior, use_threshold_reri)
  
  # Merge with stage classification
  pos_exp <- merge(pos_exp, stage_class, by = c("nichd", "dynamic"), all.x = TRUE)
  neg_exp[, class := 0]
  
  # Merge with co-administration data
  pos_exp <- merge(pos_exp, coadmin_stage_pos[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  neg_exp <- merge(neg_exp, coadmin_stage_neg[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  
  # High-reporting dataset (excluding uniform dynamic)
  pos_high <- pos_exp[class == 1]
  neg_high <- neg_exp[class == 0]
  
  list(
    pos_all = pos_exp,
    pos_high = pos_high,
    neg_high = neg_high,
    reduction_pct = red_pct
  )
}

################################################################################
# Signal detection function
################################################################################

# Detects signals according to the specified method
#
# Parameters:
# data: expanded triplet results (long format)
# thresholds: per-stage thresholds
# method: detection method ("gam" or "classic")
# criterion: detection criterion ("ior", "reri", or "double")
# use_null: if TRUE, applies null distribution thresholds
#

detect_signal <- function(dt, method_name, detection_type, use_null) {
  
  is_gam <- grepl("GAM", method_name) # checks for "GAM" in the method string
  
  # Determine columns based on the method
  if (is_gam) {
    ior_col <- "gam_log_ior_lower90"
    reri_col <- "gam_reri_lower90"
    thresh_ior_col <- "threshold_ior"  # null distribution thresholds
    thresh_reri_col <- "threshold_reri"
  } else {
    ior_col <- "classic_log_ior_lower90"
    reri_col <- "classic_reri_lower90"
    # Classic null thresholds. populated when 20_null runs classic methods
    thresh_ior_col  <- "threshold_classic_ior"
    thresh_reri_col <- "threshold_classic_reri"
  }
  
  # Compute detection flags by criterion type
  # use_null applies to both GAM (threshold_ior/reri) and classic (threshold_classic_ior/reri)
  if (detection_type == "IOR") {
    if (use_null && thresh_ior_col %in% names(dt)) {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)]
    } else {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0]}
  } else if (detection_type == "RERI") {
    if (use_null && thresh_reri_col %in% names(dt)) {
      dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0 & get(reri_col) > get(thresh_reri_col)]
    } else {
      dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0]}
  } else {
    # Double criterion: IOR OR RERI
    ior_det <- if (use_null && thresh_ior_col %in% names(dt)) {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0 & dt[[ior_col]] > dt[[thresh_ior_col]]
    } else {
      !is.na(dt[[ior_col]]) & dt[[ior_col]] > 0}
    reri_det <- if (use_null && thresh_reri_col %in% names(dt)) {
      !is.na(dt[[reri_col]]) & dt[[reri_col]] > 0 & dt[[reri_col]] > dt[[thresh_reri_col]]
    } else {
      !is.na(dt[[reri_col]]) & dt[[reri_col]] > 0}
    dt[, detected := ior_det | reri_det]
  }
  
  # Replace NA with FALSE (i.e., not detected)
  dt[is.na(detected), detected := FALSE]
  
  return(dt)
}

