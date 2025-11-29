################################################################################
# Script de generación de distribución nula
# Script 02_nulldistribution
################################################################################

library(data.table)
library(mgcv)

setwd("D:/Bioestadística/gam-farmacovigilancia")

################################################################################
# Configuración
################################################################################

ruta_ade_raw <- "./ade_raw.csv"
ruta_null_pool_meta <- "./augmentation_results/null_pool_reports_metadata.csv"
output_dir <- "./null_distribution_results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Parámetros para permutación
perm_events <- TRUE   
perm_drugs <- TRUE   

# Parámetros para muestreo
max_triplets_per_permutation <- 25  
min_reports_triplet <- 5            
target_total_triplets <- 2500        # Objetivo de tripletes
max_permutation_attempts <- 1000    # intentos máximos

seed_base <- 2025

# Parámetros de fórmula GAM
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML"


################################################################################
# Carga de funciones
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Carga de datos
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)

# columnas requeridas
cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")

# si se usa sex
if (include_sex) {
  cols_req <- c(cols_req, "sex")
}

# preprocesado
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

# Procesar sex si está presente
if (include_sex) {

  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  
  # factor con niveles estándar
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  
  # distribucion
  sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
  message("\n Distribución:")
  print(sex_summary)
}

null_pool_meta <- fread(ruta_null_pool_meta)

# Preparo listas de drogas y eventos
null_pool_data <- ade_raw_dt[safetyreportid %in% null_pool_meta$safetyreportid]

drugs_by_report <- unique(null_pool_data[!is.na(atc_concept_id), 
                                         .(safetyreportid, atc_concept_id)])
events_by_report <- unique(null_pool_data[!is.na(meddra_concept_id), 
                                          .(safetyreportid, meddra_concept_id)])

drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), 
                               by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), 
                                 by = safetyreportid]

pool_reports_meta <- merge(null_pool_meta, drugs_list, by = "safetyreportid", all.x = TRUE)
pool_reports_meta <- merge(pool_reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
pool_reports_meta[is.na(drugs), drugs := list(integer(0))]
pool_reports_meta[is.na(events), events := list(integer(0))]

message(sprintf("Pool nulo: %s reportes", format(nrow(pool_reports_meta), big.mark = ",")))

################################################################################
# LOOP PRINCIPAL SECUENCIAL
################################################################################

# no estoy logrando hacer que esto funcione con paralelización sin que crashee pc
# posiblemente muy ineficiente
# tiempo aprox para 5000 tripletes: 1 día y medio

all_results <- list()
triplets_collected <- 0
permutation_attempt <- 0
failed_attempts <- 0

while (triplets_collected < target_total_triplets && 
       permutation_attempt < max_permutation_attempts) {
  
  permutation_attempt <- permutation_attempt + 1
  
  if (permutation_attempt %% 5 == 1) {
    message(sprintf("\n--- Permutación %d | Tripletes: %d/%d (%.1f%%) ---",
                    permutation_attempt, 
                    triplets_collected, 
                    target_total_triplets,
                    100 * triplets_collected / target_total_triplets))
  }
  
  set.seed(seed_base + permutation_attempt)
  
  ###########
  # Permutación del pool
  ###########
  
  permuted_pool <- tryCatch({
    permute_pool(
      pool_reports_meta, niveles_nichd,
      perm_events = perm_events, 
      perm_drugs = perm_drugs,
      seed = permutation_attempt
    )
  }, error = function(e) {
    message(sprintf(" error %s", e$message))
    NULL
  })
  
  if (is.null(permuted_pool) || nrow(permuted_pool) == 0) {
    failed_attempts <- failed_attempts + 1
    next
  }
  
  ###########
  # Construcción de tripletes
  ###########
  
  triplets_perm <- tryCatch({
    permuted_pool[, {
      drugs_vec <- drugs_perm[[1]]
      events_vec <- events_perm[[1]]
      
      if (length(drugs_vec) >= 2 && length(events_vec) >= 1) {
        make_triplets_per_report(drugs_vec, events_vec, safetyreportid, nichd_num)
      } else {
        data.table()
      }
    }, by = safetyreportid]
  }, error = function(e) NULL)
  
  if (is.null(triplets_perm) || nrow(triplets_perm) == 0) {
    rm(permuted_pool)
    gc(verbose = FALSE)
    failed_attempts <- failed_attempts + 1
    next
  }
  
  ###########
  # Filtrado y selección aleatoria de los tripletes a ajustar
  ###########
  
  trip_counts <- unique(triplets_perm[, .(drugA, drugB, meddra, safetyreportid)])[
    , .N, by = .(drugA, drugB, meddra)
  ]
  candidate_triplets <- trip_counts[N >= min_reports_triplet]
  
  if (nrow(candidate_triplets) == 0) {
    rm(permuted_pool, triplets_perm, trip_counts)
    gc(verbose = FALSE)
    failed_attempts <- failed_attempts + 1
    next
  }
  
  n_to_sample <- min(nrow(candidate_triplets), max_triplets_per_permutation)
  selected_triplets <- candidate_triplets[sample(.N, n_to_sample)]
  
  ###########
  # Reintroducción de los tripletes a copia del dataset para ajustar
  ###########
  
  ade_modified <- tryCatch({
    reintroduce_permuted_reports(ade_raw_dt, permuted_pool)
  }, error = function(e) {
    message(sprintf(" error %s", e$message))
    NULL
  })
  
  if (is.null(ade_modified) || nrow(ade_modified) == 0) {
    rm(permuted_pool, triplets_perm, trip_counts, candidate_triplets, selected_triplets)
    gc(verbose = FALSE)
    failed_attempts <- failed_attempts + 1
    next
  }
  
  ###########
  # Validación de tripletes seleccionados (para que tenga requisitos mínimos para que modelo converja)
  ###########
  
  reports_by_drug <- ade_modified[!is.na(atc_concept_id), 
                                  .(reports = list(unique(safetyreportid))), 
                                  by = atc_concept_id]
  setkey(reports_by_drug, atc_concept_id)
  
  reports_by_event <- ade_modified[!is.na(meddra_concept_id), 
                                   .(reports = list(unique(safetyreportid))), 
                                   by = meddra_concept_id]
  setkey(reports_by_event, meddra_concept_id)
  
  validation_results <- selected_triplets[, {
    rA <- reports_by_drug[.(drugA), reports][[1]]
    rB <- reports_by_drug[.(drugB), reports][[1]]
    rE <- reports_by_event[.(meddra), reports][[1]]
    
    n_coadmin <- if (!is.null(rA) && !is.null(rB)) {
      length(intersect(rA, rB))
    } else 0
    
    n_events <- if (!is.null(rE)) length(rE) else 0
    
    data.table(
      drugA, drugB, meddra,
      valid_gam = (n_coadmin >= 3 && n_events >= 5)
    )
  }, by = .I]
  
  valid_triplets <- validation_results[valid_gam == TRUE]

  }
  
  ###########
  # Ajuste secuencial de modelo
  ###########
  
  n_valid <- nrow(valid_triplets)
  n_successful_gams <- 0
  
  for (ti in seq_len(n_valid)) {
    
    rowt <- valid_triplets[ti]
    
    model_res <- tryCatch({
      fit_differential_gam(
    drugA_id = rowt$drugA,
    drugB_id = rowt$drugB,
    event_id = rowt$meddra,
    ade_data = ade_modified,
    spline_individuales = spline_individuales,
    include_sex = include_sex,
    nichd_spline = nichd_spline,
    include_stage_sex = include_stage_sex,
    bs_type = bs_type,
    select = select,
    k_spline = k_spline)
    }, error = function(e) list(success = FALSE))
    
    if (!model_res$success) next
    
    if (any(is.na(model_res$log_ior)) || 
        any(is.infinite(model_res$log_ior)) ||
        any(abs(model_res$log_ior) > 20)) next
    
    # guardo resultados
    result_dt <- data.table(
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      stage = 1:7,
      log_lower90 = model_res$log_ior_lower90,
      log_ior = model_res$log_ior,
      permutation = permutation_attempt,
      spline_individuales = spline_individuales,
      include_sex = include_sex,
      nichd_spline = nichd_spline,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = if(model_res$success) model_res$formula_used else NA_character_
    )
    
    all_results[[length(all_results) + 1]] <- result_dt
    n_successful_gams <- n_successful_gams + 1
  }
  
  if (n_successful_gams > 0) {
    triplets_collected <- triplets_collected + n_successful_gams
    
    message(sprintf("     %d/%d ajustes exitosos | total: %d/%d (%.1f%%)",
                    n_successful_gams, n_valid,
                    triplets_collected, target_total_triplets,
                    100 * triplets_collected / target_total_triplets))
  } else {
    message("  Error ningun modelo esta convergiendo")
    failed_attempts <- failed_attempts + 1
  }
  
  ###########
  # Limpieza de memoria
  ###########
  
  rm(permuted_pool, triplets_perm, trip_counts, candidate_triplets,
     selected_triplets, ade_modified, reports_by_drug, reports_by_event,
     validation_results, valid_triplets)
  gc(verbose = FALSE)
  
  # Libero memoria cada 10 permutaciones
  if (permutation_attempt %% 10 == 0) {
    gc(full = TRUE)
    Sys.sleep(1)  # Pausa para estabilidad
  }
  
  # checkpoint cada 50 permutaciones
  if (permutation_attempt %% 50 == 0 && triplets_collected > 0) {
    message(sprintf("\n  Guardado de resultados parciales"))
    checkpoint_data <- rbindlist(all_results, fill = TRUE)
    fwrite(checkpoint_data, 
           paste0(output_dir, "null_distribution_checkpoint.csv"))
    message(sprintf("  guardados %d tripletes", nrow(checkpoint_data) / 7))
  }
}

################################################################################
# Resultados
################################################################################


message(sprintf("total de permutaciones: %d", permutation_attempt))
message(sprintf("permutaciones exitosas: %d", permutation_attempt - failed_attempts))
message(sprintf("tasa de éxito: %.1f%%",
                100 * (permutation_attempt - failed_attempts) / permutation_attempt))
message(sprintf("\nTripletes recolectados: %d", triplets_collected))
message(sprintf("Objetivo: %d (%.1f%% completado)",
                target_total_triplets,
                100 * triplets_collected / target_total_triplets))

null_all <- rbindlist(all_results, fill = TRUE)
message(sprintf("Total de observaciones: %d", nrow(null_all)))

rm(all_results)
gc(full = TRUE)

################################################################################
# Limpieza de outliers
################################################################################

# no se si hace falta esto, creo que es relevante solo en modelos con muchos parámetros que dan muchos valores inestables

null_all[, `:=`(
  q005 = quantile(log_lower90, 0.005, na.rm = TRUE),
  q995 = quantile(log_lower90, 0.995, na.rm = TRUE)
), by = stage]

null_cleaned <- null_all[
  log_lower90 >= q005 & log_lower90 <= q995
][, .(stage, log_lower90, log_ior, permutation)]

pct_removed <- 100 * (1 - nrow(null_cleaned) / nrow(null_all))
message(sprintf(" outliers removidos: %.1f%%", pct_removed))

################################################################################
# Calculo de umbrales
################################################################################


null_thresholds <- null_cleaned[, .(
  threshold_p90 = quantile(log_lower90, 0.90, na.rm = TRUE),
  threshold_p95 = quantile(log_lower90, 0.95, na.rm = TRUE),
  threshold_p99 = quantile(log_lower90, 0.99, na.rm = TRUE),
  threshold_p999 = quantile(log_lower90, 0.999, na.rm = TRUE),
  n_samples = .N,
  mean_null = mean(log_lower90, na.rm = TRUE),
  sd_null = sd(log_lower90, na.rm = TRUE)
), by = stage]

null_thresholds[, stage_name := niveles_nichd[stage]]

cat("umbrales por etapa:\n")
print(null_thresholds[, .(stage, stage_name, threshold_p99, n_samples)])

################################################################################
# Guardado de resultados 
################################################################################

fwrite(null_cleaned, paste0(output_dir, "null_distribution.csv"))
fwrite(null_thresholds, paste0(output_dir, "null_thresholds.csv"))

execution_summary <- data.table(
  parameter = c("perm_events", "perm_drugs", "max_triplets_per_permutation",
                "target_total_triplets", "total_permutations", 
                "triplets_collected", "success_rate"),
  value = c(perm_events, perm_drugs, max_triplets_per_permutation,
            target_total_triplets, permutation_attempt,
            triplets_collected,
            100 * (permutation_attempt - failed_attempts) / permutation_attempt)
)
fwrite(execution_summary, paste0(output_dir, "execution_summary.csv"))

