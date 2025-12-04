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

# codifico cómo se guardan los resultados de las distintas formulas
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

ruta_null_pool_meta <- paste0("./results/", suffix, "/augmentation_results/null_pool_reports_metadata.csv")

output_dir <- paste0("./results/", suffix, "/null_distribution_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


n_cores <- max(1, floor(detectCores() * 0.5))
batch_size <- 5  # para que no colapse por llenado de memoria
checkpoint_frequency <- 2  # guardo resultados cada 5 lotes

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
# LOOP PRINCIPAL PARALELIZADO
################################################################################

# Configuración de cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# objetos necesarios
clusterExport(cl, c(
  "ade_raw_dt", "pool_reports_meta", "niveles_nichd",
  "permute_pool", "reintroduce_permuted_reports", "make_triplets_per_report",
  "fit_differential_gam",
  "perm_events", "perm_drugs", "min_reports_triplet", 
  "max_triplets_per_permutation", "seed_base",
  "spline_individuales", "include_sex", "include_stage_sex",
  "k_spline", "nichd_spline", "bs_type", "select", "method"
), envir = environment())

# librerías
clusterEvalQ(cl, {
  library(data.table)
  library(mgcv)
})
  
# Variables de control
n_batches <- ceiling(max_permutation_attempts / batch_size)
triplets_collected <- 0
permutation_attempt <- 0
failed_attempts <- 0
batch_files <- character()  # para rastreo de archivos temporales

for (batch in 1:n_batches) {
  
  # Verificar si ya alcanzamos el objetivo
  if (triplets_collected >= target_total_triplets) {
    message("\n¡Objetivo alcanzado!")
    break
  }
  
  start_perm <- (batch - 1) * batch_size + 1
  end_perm <- min(batch * batch_size, max_permutation_attempts)
  batch_perms <- start_perm:end_perm
  
  message(sprintf("\n=== Lote %d/%d | Permutaciones %d-%d ===", 
                  batch, n_batches, start_perm, end_perm))
  
  # Procesamiento paralelo del lote
  batch_results <- foreach(
    perm_id = batch_perms,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE
  ) %dopar% {
    
    set.seed(seed_base + perm_id)
    
    ###########
    # Permutación del pool
    ###########
    
    permuted_pool <- tryCatch({
      permute_pool(
        pool_reports_meta, niveles_nichd,
        perm_events = perm_events, 
        perm_drugs = perm_drugs,
        seed = perm_id
      )
    }, error = function(e) NULL)
    
    if (is.null(permuted_pool) || nrow(permuted_pool) == 0) {
      return(list(success = FALSE, reason = "permutation_failed"))
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
      return(list(success = FALSE, reason = "no_triplets"))
    }
    
    ###########
    # Filtrado y selección de tripletes
    ###########
    
    trip_counts <- unique(triplets_perm[, .(drugA, drugB, meddra, safetyreportid)])[
      , .N, by = .(drugA, drugB, meddra)
    ]
    candidate_triplets <- trip_counts[N >= min_reports_triplet]
    
    if (nrow(candidate_triplets) == 0) {
      return(list(success = FALSE, reason = "no_candidates"))
    }
    
    n_to_sample <- min(nrow(candidate_triplets), max_triplets_per_permutation)
    selected_triplets <- candidate_triplets[sample(.N, n_to_sample)]
    
    ###########
    # Reintroducción de reportes permutados
    ###########
    
    ade_modified <- tryCatch({
      reintroduce_permuted_reports(ade_raw_dt, permuted_pool)
    }, error = function(e) NULL)
    
    if (is.null(ade_modified) || nrow(ade_modified) == 0) {
      return(list(success = FALSE, reason = "reintroduction_failed"))
    }
    
    ###########
    # Validación de tripletes
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
    
    if (nrow(valid_triplets) == 0) {
      return(list(success = FALSE, reason = "no_valid_triplets"))
    }
    
    ###########
    # Ajuste de modelos
    ###########
    
    triplet_results <- list()
    
    for (ti in seq_len(nrow(valid_triplets))) {
      
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
          k_spline = k_spline
        )
      }, error = function(e) list(success = FALSE))
      
      if (!model_res$success) next
      
      if (any(is.na(model_res$log_ior)) || 
          any(is.infinite(model_res$log_ior)) ||
          any(abs(model_res$log_ior) > 20)) next
      
      # guardado de resultados
      result_dt <- data.table(
        drugA = rowt$drugA,
        drugB = rowt$drugB,
        meddra = rowt$meddra,
        stage = 1:7,
        log_lower90 = model_res$log_ior_lower90,
        log_ior = model_res$log_ior,
        permutation = perm_id,
        spline_individuales = spline_individuales,
        include_sex = include_sex,
        nichd_spline = nichd_spline,
        include_stage_sex = include_stage_sex,
        k_spline = k_spline,
        bs_type = bs_type,
        select = select,
        formula_used = if(model_res$success) model_res$formula_used else NA_character_
      )
      
      triplet_results[[length(triplet_results) + 1]] <- result_dt
    }
    
    # Limpieza de memoria en worker
    rm(permuted_pool, triplets_perm, trip_counts, candidate_triplets,
       selected_triplets, ade_modified, reports_by_drug, reports_by_event,
       validation_results, valid_triplets)
    gc(verbose = FALSE, full = TRUE)
    
    # resultados del worker
    if (length(triplet_results) > 0) {
      return(list(
        success = TRUE, 
        n_triplets = length(triplet_results),
        results = rbindlist(triplet_results, fill = TRUE)
      ))
    } else {
      return(list(success = FALSE, reason = "no_convergence"))
    }
  }
  
  ###########
  # Procesamiento de resultados del lote
  ###########
  
  # filtrado de errores
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)
  
  # conteo de éxitos y fallos
  batch_successes <- sum(sapply(batch_results_clean, function(x) x$success))
  batch_failures <- length(batch_results_clean) - batch_successes
  
  permutation_attempt <- end_perm
  failed_attempts <- failed_attempts + batch_failures
  
  # extraigo y guardo resultados exitosos
  successful_results <- Filter(function(x) x$success, batch_results_clean)
  
  if (length(successful_results) > 0) {
    
    batch_triplets <- sum(sapply(successful_results, function(x) x$n_triplets))
    
    # combino resultados del lote
    batch_data <- rbindlist(lapply(successful_results, function(x) x$results), fill = TRUE)
    
    # guardo archivo temporal de lote para no colapsar memoria
    batch_file <- paste0(output_dir, "batch_", sprintf("%04d", batch), ".csv")
    fwrite(batch_data, batch_file)
    batch_files <- c(batch_files, batch_file)
    
    triplets_collected <- triplets_collected + batch_triplets
    
    message(sprintf("  Lote completado: %d/%d exitosos | Total: %d/%d tripletes (%.1f%%)",
                    batch_successes, length(batch_perms),
                    triplets_collected, target_total_triplets,
                    100 * triplets_collected / target_total_triplets))
    message(sprintf("  Guardado: %s (%d filas)", basename(batch_file), nrow(batch_data)))
    
  } else {
    message(sprintf("  Lote completado: 0/%d exitosos", length(batch_perms)))
  }
  
  # Limpieza de memoria después de cada lote
  rm(batch_results, batch_results_clean, successful_results)
  if (exists("batch_data")) rm(batch_data)
  gc(full = TRUE)
  
  # pausa para estabilidad del sistema
  Sys.sleep(2)
}

# Cerrar cluster
stopCluster(cl)

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

null_all <- rbindlist(lapply(batch_files, fread), fill = TRUE)
message(sprintf("Total de observaciones: %d", nrow(null_all)))


# Leer en chunks para evitar saturar memoria
chunk_size <- 10  # Leer 10 archivos a la vez
n_chunks <- ceiling(length(batch_files) / chunk_size)

null_all_chunks <- list()

for (chunk in 1:n_chunks) {
  
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, length(batch_files))
  chunk_files <- batch_files[start_idx:end_idx]
  
  chunk_data <- rbindlist(lapply(chunk_files, fread), fill = TRUE)
  null_all_chunks[[chunk]] <- chunk_data
  
  rm(chunk_data)
  gc()
}

# combino todos los chunks
null_all <- rbindlist(null_all_chunks, fill = TRUE)
rm(null_all_chunks)
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


