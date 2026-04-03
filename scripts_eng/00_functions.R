################################################################################
# Script de funciones
# Script: 00_functions.R
# para usar, source("00_functions.R", local = TRUE)
################################################################################

################################################################################
# General pipeline configuration
################################################################################

# project
setwd("D:/Bioestadística/gam-farmacovigilancia")
set.seed(7113)

# loading with libraries for the entire pipeline with pacman
library(pacman)
pacman::p_load(tidyverse, data.table, pbapply, parallel, doParallel, foreach,
  mgcv, MASS, akima, doRNG, pROC, svglte, DHARMa,
  igraph, ggraph, tidygraph, scales, RColorBrewer, patchwork, graphlayouts, ggrepel, networkD33, htmlwidgets, ggalluvial)

# level order for the rest of the scripts
niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

# cores for parallelization
n_cores <- max(1, floor(detectCores() * 0.75)) 

# general routes
ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_gene <- "./drug_gene.csv"
ruta_drug_info <- "./drug.csv" 
ruta_concept <- "./vocabulary/concept.csv"
ruta_rel <- "./vocabulary/concept_relationship.csv"
ruta_twosides <- "./twosides/TWOSIDES.csv.gz"

# Formula parameters for GAM
spline_individuales <- TRUE 
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7    
include_nichd <- FALSE
nichd_spline <- FALSE 
bs_type <- "cs"
select <- FALSE
method <- "fREML" 

# formula coding for saving
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

Z90 <- qnorm(0.95)  # 90% quantile for confidence intervals

# chosen percentile of the null distribution
percentil <- "p95"

source("01_theme.R", local = TRUE)
theme_set(theme_base())               

################################################################################
# Function to construct triplets
################################################################################

# Generate triplets (drugA, drugB, event) for an individual report
#
# Parameters:
# drug: drug ids in the report (atc_concept_id)
# event: ids of events in the report (meddra_concept_id)
# report_id: report id (safetyreportid)
# nichd_stage: numerical NICHD stage of the report
#
# Return:
# data.table with columns: safetyreportid, drugA, drugB, meddra, nichd_num
#
# Implementation:
# Generate all combinations of drug pairs (if >= 2 drugs)
# Cross each pair with each event in the report
# Apply order: drugA <= drugB (to avoid duplicates in reverse order)
# Returns NULL if the report does NOT have >= 2 drugs or events

make_triplets <- function(drug, event, report_id, nichd_stage) {
  
  if (length(drug) < 2 || length(event) < 1) return(NULL)
  
  drug <- unique(drug)
  event <- unique(event)
  
  # drug pair combinations
  if (length(drug) == 2) {
    combination <- matrix(c(min(drug), max(drug)), nrow = 1)
  } else {
    combination <- t(combn(drug, 2))
    # Order: min first
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

# Sample fold-changes following negative exponential distribution
#
#Parameters:
# n: Number of fold-changes to generate
# lambda: rate parameter for exponential distribution (0.75)
#
# return: numeric vector with fold-changes >= 1
#
# Implementation:
# -FC ~ 1 + exp(λ = 0.75)
# -typical range: [1, 10] right-skewed

fold_change <- function(n, lambda = 0.75) {
  1 + rexp(n, rate = lambda)
}

################################################################################
# Helper function: calculate co-management per stage for a triplet
################################################################################

# Calculates the number of A+B reports per stage for each triplet
# Allows calculating superset for classic method
#
# return:
# n_coadmin_stage: number of A+B reports for event in stage
#
# Implementation:
# Identifies A+B reports
# Count reports by stage
# Complete stages without reports with 0

coadmin_by_stage <- function(drugA, drugB, meddra, ade_data) {
  
  # reports with each drug
  reports_A <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  reports_B <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  
  # Reports with both drugs (coadministration)
  reports_AB <- intersect(reports_A, reports_B)
  
  if (length(reports_AB) == 0) {
    # If there is no co-administration, return empty structure
    return(data.table(
      nichd_num = 1:7,
      nichd = niveles_nichd,
      n_coadmin_stage = 0L
    ))
  }
  
  # NICHD stage of each report with coadministration
  stage_counts <- unique(ade_data[
    safetyreportid %in% reports_AB,
    .(safetyreportid, nichd, nichd_num)
  ])[, .N, by = .(nichd, nichd_num)]
  
  # zeros for stages without reports
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
# Dynamics generation function
################################################################################

# Generates normalized dynamic patterns for signal injection
#
#Parameters:
# type: Dynamic type: "uniform", "increase", "decrease", "plateau", "inverse_plateau"
# N: Number of stages (7 for NICHD)
#
# return: normalized numeric vector to [-1, 1] with the temporal pattern
#
# forms of dynamics
# uniform: constant signal (0 in all stages)
# increase: monotonic growth from -1 to +1
# decrease: monotonic decrease from +1 to -1
# plateau: peak in central stages (bell-shaped)
# inverse_plateau: valley in central stages (inverted U shape)

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

# Injects drug-drug interaction signal
#
#Parameters:
# drugA_id: id drug A (ATC concept_id)
# drugB_id: drug B id (ATC concept_id)
# event_id: id of the adverse event (MedDRA concept_id)
# dynamic_type: dynamic type (call generate_dynamic)
# fold_change: effect magnitude (called fold_change)
# ade_raw_dt: data.table of original dataset
#
# return List with: success, n_injected, n_coadmin, ade_aug, diagnostics
#
# Implementation:
#
#1 Base rate(e_j):
# Base value to apply fold-change
# The base value must be calculated taking into account individual components
# The base value alters the consistency of the injection in different scenarios
# Script "simulation_injeccion" shows consistency in different scenarios
# Methods that assume "independence" from events are only used as a proxy
# do not imply real independence
#
#2 Fold_change (FC):
# Multiplicative factor of the effect
#
#3 Dynamics f(j):
# Function normalized to [0, 1] per stage
# f(j) = call generate_dynamic(type, N=7) <---N always 7 by nichd stages
#
#4 Dynamic probability: (final report probability formula)
# p_dynamic(j) = e_j × FC + f(j)
#
#5 Simulation:
# For each report with drugA + drugB
# Y_new ~ Bernoulli(p_dynamic(stage_j))
# Inject ONLY into events that do not exist (if there is already an event, it is left, it tries to model dynamics on existing ones)

inject_signal <- function(drugA_id, drugB_id, event_id, 
                          dynamic_type, fold_change, 
                          ade_raw_dt) {
  
  if (drugA_id > drugB_id) {   # add sorting (it is already done in other functions but just in case)
    temp <- drugA_id
    drugA_id <- drugB_id
    drugB_id <- temp
  }
  
  # independent copy (so that injected data does not corrupt the original dataset)
  ade_aug <- copy(ade_raw_dt)
  ade_aug[, simulated_event := FALSE]
  
  #1-reports with both drugs (coadministration)
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
  
  #2-objective reporting (co-management only)
  target_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_AB, 
    .(safetyreportid, nichd, nichd_num)
  ])
  
  # Identifies if the event ALREADY exists in each report
  event_in_report <- unique(ade_raw_dt[
    meddra_concept_id == event_id, 
    safetyreportid
  ])
  
  target_reports[, e_old := as.integer(safetyreportid %in% event_in_report)]
  
  # 3-Calculate basic fee (e_j)

  # I calculate individual base rates
  reports_A_clean <- setdiff(reports_A, reports_AB)
  reports_B_clean <- setdiff(reports_B, reports_AB)
  p_baseA <- mean(reports_A_clean %in% event_in_report)
  p_baseB <- mean(reports_B_clean %in% event_in_report)
  p_base0 <- length(event_in_report) / length(unique(ade_raw_dt$safetyreportid))

  # keep in mind that they are not probabilities, they are odds, but the extrapolation is valid if p <0.1

  # There are various methods to consider base rate
  # Everyone behaves differently depending on the scenario

  # -Additive method:
  # consistent, generates high IORs when low individual risks + high baseline risk
  # e_j = P(event | A ∪ B) assuming independence
  # formula: P(A ∪ B) = P(A) + P(B) -P(A) × P(B)
  # p_baseA + p_baseB -(p_baseA *p_baseB)

  # -Multiplicative method
  # consistent, only generates high IORs when low basal risk + high individual risks
  # takes into account all components (individual and global risks)
  # e_j = (p_A × p_B) /p_0
  # formula: P(A) × P(B) /P(0)
  # (p_baseA *p_baseB) /p_base0
  e_j <- p_baseA + p_baseB - (p_baseA * p_baseB)

  # t_ij = fold_change *e_j (effect size)
  t_ij <- fold_change * e_j

  #4-Odds by stage
  # bprobs = rep(his, N)
  # dy = tanh(...) *age (dynamics scales by age)
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

  # probability clipping
  rprobs <- pmax(pmin(rprobs, 0.999), 0.001)

  #5-probability table by stage
  stage_probs <- data.table(
    nichd_num = 1:N,
    bprobs = bprobs,
    dy = dy,
    p_dynamic = rprobs
  )
  
  #6-Merge and generation of new events
  target_reports <- merge(
    target_reports, 
    stage_probs[, .(nichd_num, p_dynamic)], 
    by = "nichd_num", 
    all.x = TRUE
  )
  
  # Bernoulli simulation per report
  target_reports[, e_new := rbinom(.N, 1, p_dynamic)]
  
  # I combine with existing events
  target_reports[, e_final := pmax(e_old, e_new)]
  
  #7-frame reports to inject
  # I identify reports that SHOULD have the event (simulated)
  reports_to_mark <- target_reports[e_old == 0 & e_final == 1, safetyreportid]
  
  # validation of at least 1 injected event
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
  
  #8-frame reports that now have the event simulated (injected)
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
  
  #9-Diagnostics
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
    injection_by_stage = target_reports[
      safetyreportid %in% reports_to_mark, 
      .N, 
      by = nichd_num
    ]
  )
  
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
# Helper function: calculate basic counts
################################################################################

# Calculates basic counts (events, co-administrations) for a drug-event pair
#
#Parameters:
# ade_data: data.table augmented with reports (columns atc_concept_id, meddra_concept_id, safetyreportid).
# drugA: id of drug A (atc_concept_id).
# drugB: id of drug B (atc_concept_id).
# meddra: event id (meddra_concept_id).
#
#Return:
# list with: n_events, n_events_coadmin and n_coadmin

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
# GAM adjustment function
################################################################################

# Adjust GAM model for drug-drug interaction
#
# parameterized to run tests
#
#Parameters:
# drugA_id: id Drug A
# drugB_id: id Drug B
# event_id: id of the adverse event
# ade_data: data.table with original dataset
#
#include_nichd: If TRUE, add stage base effect as a covariate
# nichd_spline: If TRUE, use spline for stage base effect.
# If FALSE, use linear coefficient (default: TRUE)
# individual_spline: If TRUE, use splines for individual effects (smooth baseline risks drugA and B separately)
# bs_type: choice of spline type: "cs", "tp", "cr"
# select: If TRUE, allows a penalty up to zero (To see if any coefficient does not contribute)
#include_sex: If TRUE, include sex as a covariate
#include_stage_sex: If TRUE, include stage-sex interaction
# k_spline: number of knots for splines (should always be 7 per nichd)
#method; GAM adjustment method (leave "fREML")
#
#Return:
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
  #1-I identify reports
  ###########
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  reportes_coadmin <- intersect(reportes_droga_a, reportes_droga_b)
  
  # add "flagged" mock reports based on injection
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  # Previous events
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Combine real + simulated
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  n_events_total <- length(reportes_ea)  # total events (with or without drug)
  n_coadmin <- length(reportes_coadmin)  # total A+B reports (with or without event)
  n_events_coadmin <- length(intersect(reportes_coadmin, reportes_ea))  # eventos A+B

  ###########
  #2-Construction of dataset to fit
  ###########
  
  # Required base columns
  cols_necesarias <- c("safetyreportid", "nichd", "nichd_num")
  
  # sex
  if (include_sex) {
    cols_necesarias <- c(cols_necesarias, "sex")
  }
  
  datos_modelo <- unique(ade_data[, ..cols_necesarias])
  
  # exposure variables
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # if sex is present
  if (include_sex) {
    # I standardize values
    datos_modelo[, sex := toupper(trimws(sex))]
    datos_modelo[sex == "M", sex := "MALE"]
    datos_modelo[sex == "F", sex := "FEMALE"]
    # convert to factor with standard levels
    datos_modelo[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }
  
  ###########
  #4-Formula Construction with Parameters
  ###########
  
  # answer
  formula_parts <- "ea_ocurrio ~ "
  
  # option A: linear individual effects
  if (!spline_individuales) {
    formula_parts <- paste0(formula_parts, "droga_a + droga_b + ")
  } else {
    # option B: individual effects with splines
    formula_parts <- paste0(
      formula_parts,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", 
              k_spline, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", 
              k_spline, bs_type)
    )
  }

  # Effect of nichd -spline or linear (if include_nichd = TRUE)
  if (include_nichd) {
    if (nichd_spline) {
    # base spline
      formula_parts <- paste0(
        formula_parts,
        sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline, bs_type)
      )
    } else {
    # nichd linear effect
      formula_parts <- paste0(formula_parts, "nichd_num + ")
    }
  }

  # interaction spline (do not modify this)
  formula_parts <- paste0(
    formula_parts,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline, bs_type)
  )
  
  # yes sex TRUE
  if (include_sex) {
    if (include_stage_sex) {
      # Yes sex with spline by nichd
      formula_parts <- paste0(
        formula_parts,
        sprintf(" + s(nichd_num, k = %d, bs = '%s', by = sex)", 
                k_spline, bs_type)
      )
    } else {
      # Yes only sex lineal
      formula_parts <- paste0(formula_parts, " + sex")
    }
  }
  
  # Formula final
  formula_final <- as.formula(formula_parts)
  
  ###########
  #5-Template Tuning
  ###########
  
  tryCatch({
    
    modelo <- bam(
      formula = formula_final,
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = method,
      select = select,    
      discrete = TRUE,
      nthreads = 1      # if I don't put this, it gives me problems in 10_augmentation with the use of parallelization
    )
    
    ###########
    #6-Log-IOR calculation by nichd
    ###########
    
    # Prediction Grid (all combinations)
    grid_dif <- CJ(
      nichd_num = 1:7, 
      droga_a = c(0, 1), 
      droga_b = c(0, 1)
    )
    grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
    
    # if sex in the formula, choose reference level
    if (include_sex) {
      # MALE reference
      grid_dif[, sex := factor("MALE", levels = c("MALE", "FEMALE"))]
    }
    
    # predictions
    pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
    grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
    
    # Pivot for contrast
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                  value.var = c("lp", "se"))
    
    # Calculo final: log(IOR) = log(OR₁₁) -log(OR₁₀) -log(OR₀₁) + log(OR₀₀)
    # This is the same as doing log( OR₁₁ /OR₁₀ . OR₀₁)
    log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
    
    ###########
    #7-Calculation of log-IOR standard error with CoVAR matrix
    ###########
    
    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE)  # unconditional = TRUE applies correction to include smoothing uncertainty
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
      
      # contrast vector
      cvec <- rep(0, nrow(grid_dif))
      cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)
      
      # SE = sqrt(c' Σ c)
      log_ior_se[stage] <- sqrt(max(
        as.numeric(t(cvec) %*% cov_link %*% cvec), 
        0
      ))
    }
    
    ###########
    #8-IC calculation and metrics
    ###########
    
    # z90 para IC
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
    #9-Calculation of RERI per stage (absolute risk) with IC90
    ###########
    
    # stages
    stages <- sort(unique(datos_modelo$nichd_num))
    
    # newdata with the 4 A/B combinations per stage
    nd_reri <- rbindlist(lapply(stages, function(s) {
      data.table(
        nichd_num = s,
        droga_a   = c(0, 1, 0, 1),
        droga_b   = c(0, 0, 1, 1),
        droga_ab  = c(0, 0, 0, 1)
      )
    }), use.names = TRUE)
    
    # Additional covariates (if any)
    if (include_sex) {
      nd_reri[, sex := factor(levels(datos_modelo$sex)[1], 
                             levels = levels(datos_modelo$sex))]
    }
    
    if (include_nichd && !nichd_spline) {
      nd_reri[, nichd := factor(niveles_nichd[nichd_num],
                                levels = niveles_nichd,
                                ordered = TRUE)]
    }
    
    # Risk scale prediction
    pred_reri <- predict(modelo, newdata = nd_reri, type = "link", se.fit = TRUE)
    nd_reri[, `:=`(
      eta = pred_reri$fit,
      se  = pred_reri$se.fit
    )]
    
    # Parametric Bootstrap --> unlike the delta method, it allows capturing nonlinearity of the method
    # design matrix and coefficients
    X_reri <- predict(modelo, newdata = nd_reri, type = "lpmatrix")
    beta_hat <- coef(modelo)
    V_beta <- vcov(modelo, unconditional = TRUE)
    
    # n simulations
    B <- 2000
    
    # simulation of coefficients from their joint distribution
    # beta_sim ~ MVN(beta_hat, V_beta)
    beta_sims <- mvrnorm(n = B, mu = beta_hat, Sigma = V_beta)
    
    # I calculate predictions for each set of the simulated coefficients
    p_sims <- matrix(NA, nrow = nrow(nd_reri), ncol = B)

    p_sims <- plogis(X_reri %*% t(beta_sims))

    # Helper Para Reri
    calc_reri <- function(p) {
      #p is a vector of 4 probabilities per stage: [p00, p10, p01, p11]
      p11 <- p[4]; p10 <- p[2]; p01 <- p[3]; p00 <- p[1]
      if (p00 <= 0) return(NA_real_)          # save against division by zero
      p11/p00 - p10/p00 - p01/p00 + 1
    }
    
    # Calculation per stage
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
    
    # return vectors
    reri_values <- reri_dt$RERI
    reri_lower90 <- reri_dt$RERI_lower90
    reri_upper90 <- reri_dt$RERI_upper90
    
    ###########
    #10-Results
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
      formula_used = formula_parts,  # I save used formula
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
# Classic IOR calculation
################################################################################

# Calculate classical IOR using 2x2 tables per stage
#
# Parameters:
# drugA_id: drug id A (ATC concept_id)
#drugB_id: Drug B id (ATC concept_id)
# event_id: adverse event id (MedDRA concept_id)
# ade_data : data . table with dataset
#
# return : List with stage , ior_classic , ior_classic_lower90 , ior_classic_upper90
#
# Implementation:
# For each stage j:
# 2x2 table:
# a : event + coadmin
# b : unevent + coadmin
# c : event + no coadmin
# d: no event + no coadmin
# OR_11 = (a/b) /(c/d) = ad/bc
# OR_10 = reports with A only
# OR_01 = reports with B only
# OR_00 = reports without A or B
#IOR = (OR_11 × OR_00) /(OR_10 × OR_01)
# calculate IC90 by Woolf method (log scale)

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Identification of reports
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reports with the event (real)
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
  
  # real + simulated
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # Unique dataset per report with exposures
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Calculation per stage
  resultados_por_etapa <- datos_unicos[, {
    
    # 2x2 table for each exposure combination
    #Group 11: A + B (co-management)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_no_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 0)
    
    #Group 10: A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_no_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 0)
    
    #Group 01: B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_no_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 0)
    
    # Group 00: neither A nor B
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_no_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 0)
    
    # OR calculation for each group
    or_11 <- (n_11_evento / n_11_no_evento) / (n_00_evento / n_00_no_evento)
    or_10 <- (n_10_evento / n_10_no_evento) / (n_00_evento / n_00_no_evento)
    or_01 <- (n_01_evento / n_01_no_evento) / (n_00_evento / n_00_no_evento)
    or_00 <- 1  # by definition
    
    # as by definition OR₀₀ = 1, it is the same as doing IOR = OR₁₁ /OR₁₀ . OR₀₁
    # calculation of IOR
    ior_val <- (or_11 * or_00) / (or_10 * or_01)
    log_ior <- log(ior_val)
    
    # Variance on log scale (Woolf method)
    var_log_ior <- (1/n_11_evento + 1/n_11_no_evento +
                1/n_10_evento + 1/n_10_no_evento +
                1/n_01_evento + 1/n_01_no_evento +
                1/n_00_evento + 1/n_00_no_evento)
    se_log_ior <- sqrt(var_log_ior)
    
    # IC90 in log scale
    z90 <- qnorm(0.95)
    log_ior_lower90 <- log_ior - z90 * se_log_ior
    log_ior_upper90 <- log_ior + z90 * se_log_ior
    
    # IOR in original scale
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
      # diagnostics
      n_11_evento = n_11_evento,
      n_11_total = n_11_evento + n_11_no_evento
    )
    
  }, by = nichd_num]
  
  # sort by stage
  setorder(resultados_por_etapa, nichd_num)
  
  return(list(
    success = TRUE,
    results_by_stage = resultados_por_etapa
  ))
}

################################################################################
# Classic RERI calculation
################################################################################

# Calculate classic RERI using 2x2 tables per stage
#
#Parameters:
# drugA_id: id drug A (ATC concept_id)
# drugB_id: drug B id (ATC concept_id)
# event_id: id of the adverse event (MedDRA concept_id)
# ade_data: data.table with dataset (can be augmented)
#
# Return: List with success, results_by_stage (stage, RERI, RERI_lower90, RERI_upper90)
#
# Implementation:
# For each stage j:
# Build 2x2 table:
# R11: risk with A+B
#R10: risk with only A
#R01: risk with only B
# R00: risk without A or B (reference)
# RERI = R11 -R10 -R01 + 1
#90%CI using delta method (variance propagation)
# Variance of each risk: Var(R) = p(1-p)/n
#5RERI variance: sum of variances (assuming independence)

calculate_classic_reri <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Identification of reports
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reports with the event (real)
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Reports with simulated event (if they exist)
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  
  # Combine real + simulated
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # Unique dataset per report with exposures
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Calculation per stage
  resultados_por_etapa <- datos_unicos[, {
    
    # Counts for each exposure combination
    #Group 11: A + B (co-management)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_total <- sum(droga_a == 1 & droga_b == 1)
    
    #Group 10: A only
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_total <- sum(droga_a == 1 & droga_b == 0)
    
    #Group 01: B only
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_total <- sum(droga_a == 0 & droga_b == 1)
    
    # Group 00: neither A nor B (reference)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_total <- sum(droga_a == 0 & droga_b == 0)
    
    # Validation: all groups must have data
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
    else { # risk calculation, proportions
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

        # Binomial variances of each risk
        var_r <- function(r, n) ifelse(r > 0 & r < 1, r*(1-r)/n, 0.25/n)
        var_R11 <- var_r(R11, n_11_total)
        var_R10 <- var_r(R10, n_10_total)
        var_R01 <- var_r(R01, n_01_total)    
        var_R00 <- var_r(R00, n_00_total)

        # Variance of RERI by delta method:
        # RERI = (R11 -R10 -R01) /R00 + 1
        # dRERI/dR11 = 1/R00
        # dRERI/dR10 = -1/R00
        # dRERI/dR01 = -1/R00
        # dRERI/dR00 = -(R11 -R10 -R01) /R00^2
        numerador <- R11 - R10 - R01
        var_reri <- (var_R11 + var_R10 + var_R01) / R00^2 +
          numerador^2 * var_R00 / R00^4
        se_reri <- sqrt(var_reri)

        # IC90%
        z90 <- qnorm(0.95)
        reri_lower90 <- reri_val - z90 * se_reri
        reri_upper90 <- reri_val + z90 * se_reri

        data.table(
          stage = nichd_num[1],
          RERI_classic = reri_val,
          RERI_classic_lower90 = reri_lower90,
          RERI_classic_upper90 = reri_upper90,
          RERI_classic_se = se_reri,
          # Individual risks for diagnosis
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
  # Sort by stage
  setorder(resultados_por_etapa, nichd_num)
  
  # Check if there is enough data in at least one stage
  if (all(is.na(resultados_por_etapa$RERI_classic))) {
    return(list(
      success = FALSE,
      message = "Datos insuficientes en todas las etapas",
      results_by_stage = resultados_por_etapa
    ))
  }
  
  return(list(
    success = TRUE,
    results_by_stage = resultados_por_etapa
  ))
}

################################################################################
# Normalizing function
################################################################################

# Normalize triplet result with empty lists by default
#
#Return:
# Normalized list with default vectors instead of NULL
#
# Implementation:
# Prevent rbindlist errors when some list fields are NULL
# Apparently necessary when using parallelization in case convergence of some triplets fails

normalize_triplet_result <- function(result) {
  
  # field type list of columns that must exist
  list_fields <- c("stage", "log_ior", "log_ior_lower90", "ior_values",
                  "log_ior_classic", "log_ior_classic_lower90", "ior_classic")
  
  for (field in list_fields) {
    if (is.null(result[[field]])) {
      result[[field]] <- list(numeric(0))
    }
  }
  
  # If lists are empty, fill with default values
  if (length(result$stage[[1]]) == 0) {
    result$stage <- list(1:7)
  }
  
  if (length(result$log_ior[[1]]) == 0) {
    result$log_ior <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_lower90[[1]]) == 0) {
    result$log_ior_lower90 <- list(rep(NA_real_, 7))
  }
  
  if (length(result$ior_values[[1]]) == 0) {
    result$ior_values <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_classic[[1]]) == 0) {
    result$log_ior_classic <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_classic_lower90[[1]]) == 0) {
    result$log_ior_classic_lower90 <- list(rep(NA_real_, 7))
  }
  
  if (length(result$ior_classic[[1]]) == 0) {
    result$ior_classic <- list(rep(NA_real_, 7))
  }
  
  return(result)
}

################################################################################
# Bootstrap function by dynamics and stage
################################################################################

# Bootstrap log-IOR difference between dynamic and base (uniform)
#
#Parameters:
# data: data.table with `dynamic`, `stage` and `log_ior` columns
# dynamic_type: name of the dynamic
# stage_num: NICHD stage number.
# n_boot: number of bootstrap replicas (default 100)
#
#Return:
# data.table with bootstrap statistics for difference (mean, sd, CI)

bootstrap_dynamic_diff <- function(data, dynamic_type, stage_num, n_boot = 100) {
  
  # Data for target dynamics
  target_data <- data[dynamic == dynamic_type & stage == stage_num, log_ior]
  
  # Data for uniform (baseline)
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
# Bootstrap function by dynamics and stage (RERI)
################################################################################

# Bootstrap RERI difference between dynamic and base (uniform)
#
#Parameters:
# data: data.table with `dynamic`, `stage` and `reri` columns
# dynamic_type: name of the dynamic
# stage_num: NICHD stage number
# n_boot: number of bootstrap replicas (default 100)
#
#Return:
# data.table with bootstrap statistics for difference (mean, sd, CI)

bootstrap_dynamic_diff_reri <- function(data, dynamic_type, stage_num, n_boot = 100) {
  
  target_data <- data[dynamic == dynamic_type & stage == stage_num, reri]
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, reri]
  
  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(mean_diff = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  
  boot_diffs <- replicate(n_boot, {
    median(sample(target_data,  replace = TRUE)) -   # <--median
    median(sample(uniform_data, replace = TRUE))
  })
  
  data.table(
    mean_diff = median(boot_diffs, na.rm = TRUE),    # <--median
    ci_lower  = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper  = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}

################################################################################
# Function to process a single positive triplet with sensitivity
################################################################################

# Processes a single triplet with all reduction levels
# Designed for batch parallel use (memory intensive)
#
#Parameters:
# idx: id of the triplet in pos_meta
#pos_meta: data.table with positive triplet metadata
# ade_raw_dt: data.table with original dataset
# reduction_levels: numerical vector with reduction percentages to apply
# Configuration parameters for the GAM formula
#base_seed: base seed for reproducibility
#
#Return:
# data.table combined with results for all reduction levels
#
# Implementation:
# inject_signal() to create independent augmented dataset
# For each reduction level --> reduce dataset, adjust GAM, calculate classic IOR/RERI
# fit_reduced_model() as a fit wrapper

process_single_positive <- function(idx, pos_meta, ade_raw_dt, reduction_levels, 
                                    spline_individuales, include_sex, include_stage_sex,
                                    k_spline, bs_type, select, nichd_spline, z90, base_seed = 9427) {
  # Unique seed for the triplet (same scheme as 10_augmentation)
  set.seed(9427 + idx)
  
  rowt <- pos_meta[idx]
  rowt$type <- "positive"
  
  # INDEPENDENT augmented dataset creation
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
    # Failure result without sensitivity
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
      seed = base_seed + idx)  # Seed propagation for greater reproducibility
    
    # Fit model in reduced dataset
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
# Helper function: Reduce dataset in stages
################################################################################

# Reduce an augmented dataset by removing a random percentage of rows per stage
#
#Parameters:
# ade_aug: data.table increased
# reduction_pct: percentage to remove (ex: 10 for 10%)
# nichd_col: column with NICHD stage
# seed: implemented to make reproducible
#
#Return:
# reduced data.table

reduce_dataset_by_stage <- function(ade_aug, reduction_pct, nichd_col = "nichd", seed = NULL) {
  
  # For each stage, remove the specified percentage
  stages <- unique(ade_aug[[nichd_col]])
  
  reduced_rows <- lapply(stages, function(stage) {
    stage_data <- ade_aug[get(nichd_col) == stage]
    n_rows <- nrow(stage_data)
    n_keep <- ceiling(n_rows * (1 - reduction_pct / 100))
    
    if (n_keep < 1 && n_rows > 0) {
      n_keep <- 1  # Keep at least one row if there is data
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
# Model fitting function on reduced dataset
################################################################################

# Fit classic GAM and IOR/RERI model on a reduced dataset
# Function wrapper for adjustment and reduction of the dataset
#
#Parameters:
# ade_reduced: reduced dataset
# rowt: metadata of the triplet (drugA, drugB, meddra, etc.)
# reduction_pct: reduction percentage applied
#
#Return:
# data.table with adjustment results

fit_reduced_model <- function(ade_reduced, rowt, reduction_pct) {
  
  counts_reduced <- calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra)
  
  # Adjust GAM
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
  
  # Classic IOR calculation
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
  
  # Classic RERI calculation
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
      n_events = counts_reduced$n_events_total,
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
    n_events = model_res$n_events_total,
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
# Permutation formula for null distribution
################################################################################

# Swap events between reports
#
#Parameters:
# perm_events: if TRUE, swap events (meddra_concept_id)
# perm_drugs: if TRUE, also swap drugs between reports
#
# Use pool of reports selected in 10_augmentation (pool_meta)
# Break drug-event association to have a negative ground truth

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

# Function to put the swaps back into the original dataset before adjusting

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
#Parameters:
# dt: data.table with triplet results (with list-like columns)
# label_val: label for classification (1 = positive, 0 = negative)
#null_thresholds_dt: data.table with null distribution thresholds per stage
# use_threshold_ior: if TRUE, apply null distribution IOR threshold
# use_threshold_reri: if TRUE, apply RERI threshold of the null distribution
#
#Return:
# expanded data.table with one row per stage-triple
#
# Implementation:
# Breaks list-like columns (stage, log_ior, etc.) into individual rows
# Joins with thresholds of the null distribution
# Maintains triplet metadata (dynamic, t_ij, n_coadmin)

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
    
    # Layered
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
        # Layered
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
  
  # Merge con thresholds
  expanded <- merge(expanded, null_thresholds_dt, 
                   by.x = "stage_num", by.y = "stage", all.x = TRUE)
  
  # Saving Usage Parameters
  expanded[, `:=`(
    use_threshold_ior = use_threshold_ior,
    use_threshold_reri = use_threshold_reri
  )]
  
  return(expanded)
}

################################################################################
# GAM statistical power calculation function
################################################################################

# Calculate statistical power by filtering by effect size and number of coadmin reports
#
#Parameters:
# data_pos: data.table with expanded positive triplets
# data_neg: data.table with expanded negative triplets (not used)
#target_power: target power level (0.80)
# metric_n: name of the column to use as quantity filter ("n_coadmin" or "n_events")
# grid_resolution: number of steps for the search grid (30x30 by default)
# use_threshold_ior: if TRUE, use threshold of the null distribution
# use_threshold_reri: if TRUE, use threshold of the null distribution
# detection: detection mode ("log-ior", "reri", or "double")
#
# return:
# power_surface: data.table with t_threshold, n_threshold, power, len
# t_star: optimal threshold of t_ij (power size in injection formula)
# n_star: n_coadmin optimal threshold
# superset_pos: positive triplets that satisfy the optimal filters
# achieved_power: power achieved at the optimal point
#
# Implementation:
# Generate a 2D grid (according to t_ij × n_coadmin) using quantiles
# For each grid point, calculate the statistical power (TP /(TP + FN))
# Identifies optimal configuration that achieves target_power with maximum retention
# Add to triplet level: detection if ANY stage detects

calculate_power_gam <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,  
  use_threshold_ior = TRUE,    # Explicit parameter
  use_threshold_reri = TRUE,   # Explicit parameter
  detection = "double"  ) {    # "ior", "reri", a "double"
  
    library(data.table)

    ###########
    #1. Preparation of positive data
    ###########
  
    pos_clean <- copy(data_pos)
  
    # Check existence of columns
    if (!metric_n %in% names(pos_clean)) {
      stop(sprintf("La columna métrica '%s' no existe en los datos.", metric_n))
    }
  
    # Basic cleaning
    pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
    # Merge thresholds if specific columns are passed and missing
    if (!is.null(null_thresholds)) {
      # Check if the threshold columns that we are going to use are missing
      need_merge <- FALSE
      if (use_threshold_ior && !"threshold_ior" %in% names(pos_clean)) need_merge <- TRUE
      if (use_threshold_reri && !"threshold_reri" %in% names(pos_clean)) need_merge <- TRUE
      if (need_merge) {
        # Ensure null_thresholds has the 'stage' column
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
  # 2. 2D grid (day x n_coadmin)
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  #3. Power calculation
  ###########
  
  # Detection pre-calculation according to "detection" parameter
  
  # IOR criteria
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
  
  # Final detection according to the chosen mode
  if (detection == "ior") {
    pos_clean[, is_detected := ior_detected]
  } else if (detection == "reri") {
    pos_clean[, is_detected := reri_detected]
  } else {  # "double"
    pos_clean[, is_detected := ior_detected | reri_detected]
  }

  # Aggregation at triplet level: IF any stage detects, triplet detected
  triplet_detection_gam <- pos_clean[, .(
    triplet_detected = any(is_detected, na.rm = TRUE),
    t_ij_triplet = unique(t_ij)[1],
    n_coadmin_triplet = unique(get(metric_n))[1]
  ), by = triplet_id]
  
  # Iterate over the grid
  power_surface <- search_grid[, {
    
    # Filtro a nivel triplete: t_ij >= t_thresh AND n >= n_thresh
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
  #4. Sweet Spot Identification
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("ADVERTENCIA (GAM): No se alcanzó poder objetivo (%.0f%%). Max: %.1f%%",
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
  
  # Informative message according to detection mode
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
  #5. Superset construction
  ###########
  
  # Positive superset: ALL triplet observations that satisfy filters
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
# Statistical power calculation function for Stratified method
################################################################################

# Calculate statistical power by filtering by effect size and number of coadmin reports
#
# As classical methods tend to return NAs or ICs with Inf, I parameterize handling of NAs
#
#Parameters:
# data_pos: data.table with expanded positive triplets
# data_neg: data.table with expanded negative triplets
#target_power: target power level (0.80)
# metric_n: name of the column to use as quantity filter ("n_coadmin" or "n_events")
# grid_resolution: number of steps for the search grid (30x30 by default)
# detection: to superset IOR, RERI or both
# na_remove: choose NA handling.
# TRUE: superset includes ALL triplets with minimum detection criteria
# FALSE: final superset excludes NA even if they meet minimum threshold criteria
#
# return:
# power_surface: data.table with t_threshold, n_threshold, power, len
# t_star: optimal threshold of t_ij
# n_star: n_coadmin optimal threshold
# superset_pos: positive triplets that satisfy the optimal filters
# achieved_power: power achieved at the optimal point
#
# Implementation:
# See calculate_power_gam
# changes: calculate power at stage level

calculate_power_classic <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,
  detection = "double",  # "ior", "reri", a "double"
  na_remove = "TRUE") {
    
  ###########
  #1. Preparation of positive data
  ###########

  detection <- match.arg(detection, choices = c("ior", "reri", "double"))

  pos_clean <- copy(data_pos)
  
  n_triplets_original_total <- uniqueN(pos_clean$triplet_id)
  n_obs_original_total <- nrow(pos_clean)

  # NA cleaning according to detection criteria
  
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
  # record those that were deleted
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
  #2. 2D Grid
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  #3. Power calculation
  ###########
  
  # each row (triplet_id + stage) is a separate observation
  # because the method is stratified by stage
  
  # Detection pre-calculation according to 'detection' parameter
  
  if (detection == "ior") {
    # Solo IOR: IC90 > 0
    pos_clean[, is_detected := (
      !is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0
    )]
    
  } else if (detection == "reri") {
    # RELATED ONLY: IC90 > 0
    pos_clean[, is_detected := (
      !is.na(classic_reri_lower90) & classic_reri_lower90 > 0
    )]
    
  } else {  # detection == "double"
    # Double criterion: IOR O RERI, IC90 > 0
    pos_clean[, is_detected := (
      (!is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0) |
      (!is.na(classic_reri_lower90) & classic_reri_lower90 > 0)
    )]
  }
    
  # I keep a copy with detections
  pos_all_with_detection <- copy(pos_clean)

  # is calculated at the observation level for each row: t_ij, n_coadmin, detection
  # Not added at triplet level
  
  # Iterate over the grid
  power_surface <- search_grid[, {
    
    # Observation level filter: t_ij >= t_thresh AND n >= n_thresh
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
  #4. Sweet Spot Identification
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
  
  # detection mode for message
  detection_label <- switch(detection,
    "ior" = "solo IOR",
    "reri" = "solo RERI", 
    "double" = "IOR O RERI"
  )

  ###########
  #5. Superset construction
  ###########
  
  # Positive superset: observations that satisfy filters
  superset_pos <- pos_all_with_detection[t_ij >= t_star & get(metric_n) >= n_star]
  
  # Superset metrics
  n_retained_total <- nrow(superset_pos)
  n_detected_total <- sum(superset_pos$is_detected, na.rm = TRUE)
  power_total <- ifelse(n_retained_total > 0, n_detected_total / n_retained_total, 0)

  # count of unique triplets that meet retention criteria
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

  # Lines to see true retention, counting NAs that were excluded at the beginning of the function
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
    n_retained_total = n_retained_total,  # total withheld
    n_detected_total = n_detected_total,  # total detected
    n_triplets_retained = n_triplets_retained,
    n_triplets_total = n_triplets_total
  ))
}

################################################################################
# Function to display power surface with complete grid
################################################################################

# Generate power surface heatmap with regular grid interpolation
# Faceted by method
#
#Parameters:
# power_result: list returned by calculate_power_gam() or calculate_power_classic()
#target_power: target power level (for visual reference)
# detection: detection type ("IOR", "RERI", "double")
# t_range: X axis range (effect size)
# n_range: Y axis range (co-management reports)
# grid_size: interpolated grid size (default: 50x50 for smoothness)
#
#Return:
# ggplot object with full interpolated surface
#
# Implementation:
# Remove NAs and Inf
# Interpolate to regular grid using akima::interp (for missing data)
# Generate heatmap with geom_raster
# return optimal point statistics

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
  #1-Process multiple surfaces
  ###########
  
  all_surfaces <- rbindlist(lapply(names(power_results_list), function(met_name) {
    surface <- as.data.table(power_results_list[[met_name]]$power_surface)
    surface[, method_label := met_name]
    return(surface)
  }), fill = TRUE)
  
  # Extract optimal Parameters for subtitle
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
  #2-Cleaning original data
  ###########
  
  # Remove NA and non-finite values
  surface_clean <- all_surfaces[is.finite(t_threshold) & is.finite(n_threshold) & is.finite(power)]
  
  ###########
  #3-Full regular grid interpolation
  ###########
  
  # regular sequences with same range for X and Y (to make it square)
  # Use the widest range of both axes to maintain proportion
  common_range <- c(0, max(t_range[2], n_range[2]))

  # Create regular sequences for the grid
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
  #5-Method Tags
  ###########

  # subtitle with information about each method
  subtitulo <- paste(opt_params_dt[, sprintf("%s: t≥%.4f, N≥%.0f (%.1f%%)", method_label, t_star, n_star, achieved_power*100)], collapse = " | ")

  method_label_clean <- switch(detection,
    "IOR" = "Detección IOR",
    "RERI" = "Detección RERI",
    detection
  )
  
  ###########
  #6-Graph Construction
  ###########
  
  p <- ggplot(surface_plot, aes(x = t_threshold, y = n_threshold, fill = power)) +
    
    # geom_raster() for regular grids
    geom_raster(interpolate = FALSE) +  # interpolate=TRUE visually smoothes

    # faceted
    facet_wrap(~ method_label, ncol = 2, scales = "fixed") +  

    # Color scale
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1), 
      name = "Poder\n(TPR)",
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

    # square coordinates
    coord_fixed(ratio = diff(t_range)/diff(n_range)) +

    # Tags
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
# Metric calculation function by bootstrap
################################################################################

# Calculate classification metrics with 95%CI by bootstrap
#
#Parameters:
# dt data.table with columns: triplet_id, detected, label
# n_boot number of bootstrap replications
# aggregate_triplet if TRUE, aggregate at the triplet level with any(). if FALSE, maintains current granularity (stage or dynamic)
#
#Return:
# data.table with metrics and 95%CI

calculate_metrics <- function(dt, n_boot = 2000, aggregate_triplet = TRUE, score_type, score_type_auc) {
  
  # Aggregation at triplet level if requested
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
    # filter NA and infinite values ​​before calculating ROC
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
    
    # Bootstrap para AUC
    b_auc <- replicate(n_boot, {
      b_idx <- sample(nrow(dt_roc), replace = TRUE)
      b_dt <- dt_roc[b_idx]
      
      # check that there are both classes in the bootstrap sample
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
    # Identify unique IDs for bootstrap by triplet
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
    # Bootstrap by ROW for analysis by stage
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
  
  # Calculation of specific metrics as the average of the bootstrap (this way I avoid misaligning with IC)
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
  
  # IC calculation
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

# Function that expands data to long format
#
# Implementation:
# Load data and create objects according to reduction level
# Filter failed injections
# Expand to long format
# Filter by stages with high reporting according to dynamics

expand <- function(red_pct) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)  # object named according to reduction level
  
  ruta_pos <- paste0(ruta_base_sensitivity, "positive_triplets_results", suffix_file, ".rds")
  ruta_neg <- paste0(ruta_base_sensitivity, "negative_triplets_results", suffix_file, ".rds")
  
  pos_raw <- readRDS(ruta_pos)
  neg_raw <- readRDS(ruta_neg)
  
  # only successful injections
  pos_valid <- pos_raw[injection_success == TRUE]
  
  # expand to long format
  pos_exp <- expand_clean_all_metrics(pos_valid, 1, null_thresholds, use_threshold_ior, use_threshold_reri)
  neg_exp <- expand_clean_all_metrics(neg_raw, 0, null_thresholds, use_threshold_ior, use_threshold_reri)
  
  # Merge with stage sorting
  pos_exp <- merge(pos_exp, stage_class, by = c("nichd", "dynamic"), all.x = TRUE)
  neg_exp[, class := 0]
  
  # Merge with co-managed data
  pos_exp <- merge(pos_exp, coadmin_stage_pos[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  neg_exp <- merge(neg_exp, coadmin_stage_neg[, .(triplet_id, stage_num, n_coadmin_stage)], by = c("triplet_id", "stage_num"), all.x = TRUE)
  
  # High reporting dataset (exclude uniform)
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

# Detect signals according to method
#
#Parameters:
# data:expanded triplet results (long format)
# thresholds: thresholds per stage
# method: detection method ("gam" or "classic")
# criterion: detection criterion ("ior", "reri", or "double")
# use_null: if TRUE, apply thresholds from the null distribution
#

detect_signal <- function(dt, method_name, detection_type, use_null) {
  
  is_gam <- grepl("GAM", method_name) # look for GAM in the method string
  
  # Determine columns according to method
  if (is_gam) {
    ior_col <- "gam_log_ior_lower90"
    reri_col <- "gam_reri_lower90"
    thresh_ior_col <- "threshold_ior"   # null distribution thresholds
    thresh_reri_col <- "threshold_reri"
  } else {
    ior_col <- "classic_log_ior_lower90"
    reri_col <- "classic_reri_lower90"
    thresh_ior_col <- NULL
    thresh_reri_col <- NULL
  }
  
  # Calculate detection according to type
  if (detection_type == "IOR") {
    if (is_gam && use_null) {   # detection for gam with null distribution threshold
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)]
    } else { dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0]}
  } else if (detection_type == "RERI") {
    if (is_gam && use_null) { dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0 & get(reri_col) > get(thresh_reri_col)]
    } else {dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0]}
  } else { # Doble
    if (is_gam && use_null) {
      dt[, detected := (!is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)) | # both criteria
                       (!is.na(get(reri_col)) & get(reri_col) > 0 & get(reri_col) > get(thresh_reri_col))]
    } else { dt[, detected := (!is.na(get(ior_col)) & get(ior_col) > 0) | (!is.na(get(reri_col)) & get(reri_col) > 0)]}
  }
  
  # Replace NA with FALSE (that is, not detected)
  dt[is.na(detected), detected := FALSE]
  
  return(dt)
}

