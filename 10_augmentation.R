################################################################################
# Script de generación de datos semisintéticos con análisis de sensibilidad
# Script 10_augmentation
################################################################################

library(data.table)
library(tidyverse)
library(mgcv)
library(parallel)
library(pbapply)
library(doParallel)
set.seed(9427)

################################################################################
# Configuración
################################################################################

setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)

ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_gene <- "./drug_gene.csv"
ruta_drug_info <- "./drug.csv" 

# Parámetros para inyección
n_pos <- 500
n_neg <- 5000
lambda_fc <- 0.75
dinamicas <- c("uniform","increase","decrease","plateau","inverse_plateau")
n_cores <- max(1, floor(detectCores() * 0.50)) 
z90 <- qnorm(0.95)

# Parámetros de filtrado para construcción de tripletes
min_reports_triplet <- 2         
min_nichd_with_rep <- 2          
all_nichd_rep <- FALSE           
max_events_per_pair <- 10000

# Parámetros de fórmula para GAM
spline_individuales <- TRUE 
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7    
include_nichd <- FALSE
nichd_spline <- FALSE 
bs_type <- "cs"
select <- FALSE
method <- "fREML" 

# Parámetros de sensibilidad
reduction_levels <- seq(10, 90, by = 10)  # 10%, 20% [...] 90%

# Parámetros de procesamiento por lotes
batch_size_pos <- 25  # Tamaño de lote para positivos 
save_interval <- 1    # Guardar cada X lotes completados

suffix <- paste0(
  "sens_",
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

output_dir <- paste0("./results/", suffix, "/augmentation_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
source("01_theme.R")
theme_set(theme_base())

################################################################################
# Carga de datos
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)

cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")

# por si incluyo sexo en la formula gam
if (include_sex) {
  cols_req <- c(cols_req, "sex")
}
if (include_sex) {
  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  
  sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
  message("\nDistribución de sexo:")
  print(sex_summary)
}

message(sprintf("Dataset %s filas", format(nrow(ade_raw_dt), big.mark = ",")))

################################################################################
# Preprocesamiento 
################################################################################

drug_info_original <- fread(ruta_drug_info)
drug_info_original[, atc_concept_id := as.character(atc_concept_id)]
drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]

canonical_map <- drug_info_original[, .(
  canonical_id = min(atc_concept_id)
), by = base_name]

translation_table <- merge(
  drug_info_original[, .(atc_concept_id, base_name)],
  canonical_map,
  by = "base_name"
)

cat(sprintf("IDs originales: %d\n", uniqueN(translation_table$atc_concept_id)))
cat(sprintf("IDs únicos: %d\n", uniqueN(translation_table$canonical_id)))

ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

# unificación de drogas con mismo compuesto activo pero distinto vehiculo
ade_raw_dt <- merge(
  ade_raw_dt, 
  translation_table[, .(atc_concept_id, canonical_id)], 
  by = "atc_concept_id", 
  all.x = TRUE
)

ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]

nrow_before <- nrow(ade_raw_dt)
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

################################################################################
# Construcción de tripletes candidatos 
################################################################################

# obtengo drogas por reporte y evento
drugs_by_report <- unique(ade_raw_dt[!is.na(atc_concept_id), .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[!is.na(meddra_concept_id), .(safetyreportid, meddra_concept_id)])

reports <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), by = safetyreportid]

reports_meta <- merge(reports, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
reports_meta[is.na(drugs), drugs := list(integer(0))]
reports_meta[is.na(events), events := list(integer(0))]

report_combo <- copy(reports_meta)  

triplets_list <- pblapply(seq_len(nrow(report_combo)), function(i) {
  rowi <- report_combo[i]
  make_triplets(
    drug = rowi$drugs[[1]], 
    event = rowi$events[[1]], 
    report_id = rowi$safetyreportid, 
    nichd_stage = rowi$nichd
  )
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

trip_counts_by_stage <- unique(triplets_dt[, .(drugA, drugB, meddra, nichd_num, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra, nichd_num)
]

trip_summary <- trip_counts_by_stage[, .(
  N = sum(N),
  n_stages = uniqueN(nichd_num),
  stages_with_data = list(nichd_num)
), by = .(drugA, drugB, meddra)]

if (all_nichd_rep) {
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet & 
    n_stages == 7
  ]
  message(sprintf("Filtro estricto: todas las etapas deben tener reportes"))
} else {
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet & 
    n_stages >= min_nichd_with_rep
  ]
  message(sprintf("Filtro flexible: mínimo %d etapas con reportes", min_nichd_with_rep)) # opción de poner filtro con todas las etapas
}

# datos de resultado de filtrado 
message("\nEstadísticas de filtrado:")
message(sprintf("  Total candidatos iniciales: %s", 
                format(nrow(trip_summary), big.mark = ",")))
message(sprintf("  Candidatos que cumplen filtro de reportes (N >= %d): %s",
                min_reports_triplet,
                format(sum(trip_summary$N >= min_reports_triplet), big.mark = ",")))
message(sprintf("  Candidatos que cumplen filtro de etapas: %s",
                format(nrow(candidatos_pos), big.mark = ",")))

stage_distribution <- candidatos_pos[, .N, by = n_stages][order(n_stages)]
message("\nDistribución de candidatos por número de etapas:")
print(stage_distribution)

candidatos_pos <- candidatos_pos[, 
  .SD[sample(.N, min(.N, max_events_per_pair))], 
  by = .(drugA, drugB)
]
message(sprintf("Tripletes post diversificación: %d", nrow(candidatos_pos)))

n_pos_final <- min(n_pos, nrow(candidatos_pos))
positivos_sel <- candidatos_pos[sample(.N, n_pos_final)]

# características de candidatos
message(sprintf("\nTripletes positivos seleccionados: %d", nrow(positivos_sel)))
message(sprintf("  Reportes promedio: %.1f", mean(positivos_sel$N)))
message(sprintf("  Etapas promedio: %.1f", mean(positivos_sel$n_stages)))

################################################################################
# Selección de positivos 
################################################################################

pos_meta_base <- positivos_sel[, .(drugA, drugB, meddra, N)]
pos_meta_base[, base_triplet_id := 1:.N]

pos_meta <- pos_meta_base[, {
  data.table(
    drugA = drugA,
    drugB = drugB,
    meddra = meddra,
    N = N,
    base_triplet_id = base_triplet_id,
    dynamic = dinamicas
  )
}, by = base_triplet_id]

pos_meta[, fold_change := fold_change(.N, lambda = lambda_fc)]
pos_meta[, triplet_id := 1:.N]

top30_ids <- pos_meta[order(-N)][1:min(.N, 30), triplet_id]
pos_meta[, is_top30 := triplet_id %in% top30_ids]

fwrite(pos_meta, paste0(output_dir, "positive_triplets_metadata.csv"))

################################################################################
# Cálculo de coadministración por etapa para positivos 
################################################################################

coadmin_stage_pos_list <- lapply(1:nrow(pos_meta), function(i) {
  row <- pos_meta[i]
  result <- coadmin_by_stage(
    drugA = row$drugA, 
    drugB = row$drugB, 
    meddra = row$meddra, 
    ade_data = ade_raw_dt
  )
  result[, `:=`(
    triplet_id = row$triplet_id,
    drugA = row$drugA,
    drugB = row$drugB,
    meddra = row$meddra
  )]
  return(result)
})

coadmin_stage_pos <- rbindlist(coadmin_stage_pos_list)
setcolorder(coadmin_stage_pos, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "nichd", "n_coadmin_stage"))

fwrite(coadmin_stage_pos, paste0(output_dir, "positive_coadmin_by_stage.csv"))

rm(coadmin_stage_pos_list, coadmin_stage_pos)
gc()

################################################################################
# Procesado en lotes de tripletes positivos con análisis de sensibilidad
################################################################################

n_batches <- ceiling(nrow(pos_meta) / batch_size_pos)
positives_scores_list <- list()

for (batch in 1:n_batches) {
  
  start_idx <- (batch - 1) * batch_size_pos + 1
  end_idx <- min(batch * batch_size_pos, nrow(pos_meta))
  batch_indices <- start_idx:end_idx

  message(sprintf("Lote %d / %d (tripletes %d-%d)", batch, n_batches, start_idx, end_idx))  

  # preparo cluster para lote
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # exporto funciones y datos necesarios
  clusterExport(cl, c("process_single_positive", "fit_reduced_model", 
  "reduce_dataset_by_stage", "inject_signal", "fit_gam",
  "calculate_classic_ior", "calculate_classic_reri", "calc_basic_counts",
  "pos_meta", "ade_raw_dt", "reduction_levels", "niveles_nichd",
  "spline_individuales", "include_sex", "include_stage_sex",
  "k_spline", "bs_type", "select", "nichd_spline", "z90",
  "generate_dynamic", "fold_change"), 
  envir = environment())
  
  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library(MASS)
    library(doRNG)
  })
  
  # Procesar lote en paralelo
  batch_results <- foreach(
    idx = batch_indices,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE,
    .options.RNG = 9427
  ) %dorng% {   # cambio de %dopar% a %dorng% para reproductibilidad del loop
    process_single_positive(
      idx, pos_meta, ade_raw_dt, reduction_levels,
      spline_individuales, include_sex, include_stage_sex,
      k_spline, bs_type, select, nichd_spline, z90,
      base_seed = 9427)
  }
  
  stopCluster(cl)
  # limpieza y normalización de resultados del lote
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)
  
  if (length(batch_results_clean) > 0) {
    batch_dt <- rbindlist(batch_results_clean, fill = TRUE)
    positives_scores_list[[batch]] <- batch_dt
    # info de lotes que se van procesando
    message(sprintf("Lote %d completado: %d tripletes exitosos (base)", batch, sum(batch_dt$model_success & batch_dt$reduction_pct == 0, na.rm = TRUE)))
    
    # guardo checkpoint cada save_interval lotes o al final
    if (batch %% save_interval == 0 || batch == n_batches) {
      checkpoint_file <- paste0(output_dir, "checkpoint_positives_batch_", batch, ".rds")
      saveRDS(batch_dt, checkpoint_file)
      message(sprintf("Checkpoint guardado: %s", checkpoint_file))
    }
  }
  # libero memoria
  rm(batch_results, batch_results_clean, batch_dt)
  gc(verbose = FALSE)
}

# Combinar todos los lotes
positives_scores <- rbindlist(positives_scores_list, fill = TRUE)
rm(positives_scores_list)
gc()

# datos de injección
message(sprintf("Total tripletes procesados: %d", nrow(positives_scores[reduction_pct == 0])))
message(sprintf("Exitosos (base): %d", sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE)))

# Guardar resultados por nivel de reducción
for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)
  
  subset_data <- positives_scores[reduction_pct == red_pct]
  
  if (nrow(subset_data) > 0) {
    # guardo RDS
    saveRDS(subset_data, paste0(output_dir, "positive_triplets_results", suffix_file, ".rds"))
    
    # versión CSV
    subset_csv <- copy(subset_data)
    
    if ("diagnostics" %in% names(subset_csv)) {
      subset_csv[, diagnostics := NULL]
    }
    
    # Aplanar listas
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]
    
    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }
    fwrite(subset_csv, paste0(output_dir, "positive_triplets_results", suffix_file, ".csv"))
  }
}

gc()

pos_success <- positives_scores[reduction_pct == 0 & injection_success == TRUE]

# datos de injecciones exitosas
message("\nPositivos exitosos (base): ", sum(positives_scores$model_success & positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE))

################################################################################
# Selección de datos negativos
################################################################################

drugs_from_pos <- unique(c(pos_meta$drugA, pos_meta$drugB))
events_from_pos <- unique(pos_meta$meddra)

# características iniciales
message(sprintf("Pool de drogas: %d", length(drugs_from_pos)))
message(sprintf("Pool de eventos: %d", length(events_from_pos)))

pos_triplet_ids <- paste(
  pmin(pos_meta$drugA, pos_meta$drugB),
  pmax(pos_meta$drugA, pos_meta$drugB),
  pos_meta$meddra,
  sep = "_"
)
pos_triplet_set <- unique(pos_triplet_ids)

chunk_size <- 50
n_drugs <- length(drugs_from_pos)
n_chunks <- ceiling(n_drugs / chunk_size)

candidatos_neg_list <- list()
pb <- txtProgressBar(max = n_chunks, style = 3)

for (chunk_idx in 1:n_chunks) {
  start_idx <- (chunk_idx - 1) * chunk_size + 1
  end_idx <- min(chunk_idx * chunk_size, n_drugs)
  
  drugs_chunk <- drugs_from_pos[start_idx:end_idx]
  
  chunk_combinations <- CJ(
    drugA = drugs_chunk,
    drugB = drugs_from_pos,
    meddra = events_from_pos
  )
  
  chunk_combinations[, `:=`(
    drugA_ord = pmin(drugA, drugB),
    drugB_ord = pmax(drugA, drugB)
  )]
  chunk_combinations <- chunk_combinations[
    drugA == drugA_ord
  ]
  chunk_combinations[, `:=`(
    drugA = drugA_ord,
    drugB = drugB_ord,
    drugA_ord = NULL,
    drugB_ord = NULL
  )]
  chunk_combinations <- chunk_combinations[drugA != drugB]
  
  chunk_combinations[, triplet_id := paste(drugA, drugB, meddra, sep = "_")]
  chunk_combinations <- chunk_combinations[!triplet_id %in% pos_triplet_set]

  chunk_candidates <- merge(
    chunk_combinations[, .(drugA, drugB, meddra, triplet_id)],
    trip_summary,
    by = c("drugA", "drugB", "meddra"),
    all.x = FALSE
  )

  chunk_candidates <- chunk_candidates[N >= min_reports_triplet]
  
  if (nrow(chunk_candidates) > 0) {
    candidatos_neg_list[[chunk_idx]] <- chunk_candidates
  }
  
  rm(chunk_combinations, chunk_candidates)
  gc(verbose = FALSE)
  
  setTxtProgressBar(pb, chunk_idx)
}
close(pb)

candidatos_neg <- rbindlist(candidatos_neg_list, use.names = TRUE)
rm(candidatos_neg_list)
gc()

candidatos_neg <- unique(candidatos_neg, by = "triplet_id")

neg_ids <- paste(candidatos_neg$drugA, candidatos_neg$drugB, candidatos_neg$meddra, sep = "_")

neg_counts_by_stage <- unique(
  triplets_dt[paste(drugA, drugB, meddra, sep = "_") %in% neg_ids, 
              .(drugA, drugB, meddra, nichd_num, safetyreportid)]
)[, .(n_reports = .N), by = .(drugA, drugB, meddra, nichd_num)]

# renombro
setnames(neg_counts_by_stage, "nichd_num", "nichd")

# orden de niveles_nich (objeto se crea en 00_functions)
nichd_to_num <- setNames(1:7, niveles_nichd)
neg_counts_by_stage[, nichd_num := nichd_to_num[nichd]]

# 4. Calcular n_stages_with_data directamente (sin dcast)
stage_counts_summary <- neg_counts_by_stage[, .(
  n_stages_with_data = uniqueN(nichd_num[n_reports > 0])
), by = .(drugA, drugB, meddra)]

# merge con candidatos_neg
setkey(candidatos_neg, drugA, drugB, meddra)
setkey(stage_counts_summary, drugA, drugB, meddra)

candidatos_neg_full <- merge(
  candidatos_neg,
  stage_counts_summary,
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)

candidatos_neg_filtered <- candidatos_neg_full[n_stages_with_data >= min_nichd_with_rep]

# características de candidatos
message(sprintf("\nCandidatos negativos filtrados: %d (de %d)", nrow(candidatos_neg_filtered), nrow(candidatos_neg)))

n_neg_final <- min(n_neg, nrow(candidatos_neg_filtered))
selected_negatives <- candidatos_neg_filtered[sample(.N, n_neg_final)]

selected_negatives[, triplet_id := 1:.N]

# características de pool negativo seleccionado
message(sprintf("\nNegativos seleccionados: %d", nrow(selected_negatives)))
message(sprintf("  Reportes promedio: %.1f", mean(selected_negatives$N)))
message(sprintf("  Etapas promedio: %.1f", mean(selected_negatives$n_stages_with_data)))

coadmin_stage_neg_list <- lapply(1:nrow(selected_negatives), function(i) {
  row <- selected_negatives[i]
  result <- coadmin_by_stage(
    drugA = row$drugA, 
    drugB = row$drugB, 
    meddra = row$meddra, 
    ade_data = ade_raw_dt
  )
  result[, `:=`(
    triplet_id = row$triplet_id,
    drugA = row$drugA,
    drugB = row$drugB,
    meddra = row$meddra
  )]
  return(result)
})

coadmin_stage_neg <- rbindlist(coadmin_stage_neg_list)
setcolorder(coadmin_stage_neg, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "n_coadmin_stage"))

fwrite(coadmin_stage_neg, paste0(output_dir, "negative_coadmin_by_stage.csv"))

rm(coadmin_stage_neg_list, coadmin_stage_neg)
gc()

################################################################################
# Procesamiento de negativos por lotes con sensibilidad
################################################################################

batch_size_neg <- 50
n_batches <- ceiling(nrow(selected_negatives) / batch_size_neg)
negatives_scores_list <- list()

for (batch in 1:n_batches) {
  
  start_idx <- (batch - 1) * batch_size_neg + 1
  end_idx <- min(batch * batch_size_neg, nrow(selected_negatives))
  batch_indices <- start_idx:end_idx
  
  message("\nLote ", batch, " / ", n_batches, 
          " (tripletes ", start_idx, "-", end_idx, ")")
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  clusterExport(cl, c("fit_gam", "ade_raw_dt", "z90",
                      "selected_negatives", "normalize_triplet_result", 
                      "spline_individuales", "include_sex",
                      "include_stage_sex", "k_spline", "nichd_spline", "include_nichd",
                      "bs_type", "select", "method",
                      "calculate_classic_ior", "calculate_classic_reri",
                      "reduce_dataset_by_stage",
                      "reduction_levels", "calc_basic_counts"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library(MASS)
    library(doRNG)
  })
  
  batch_results <- foreach(
    idx = batch_indices,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE,
    .options.RNG = 9427
  ) %dorng% {

    set.seed(9427 + idx)
    
    rowt <- selected_negatives[idx]
    rowt$type <- "negative"
    
    counts_neg <- tryCatch({
      calc_basic_counts(ade_raw_dt, rowt$drugA, rowt$drugB, rowt$meddra)
    }, error = function(e) list(n_events_coadmin = NA, n_events_total = NA, n_coadmin = NA))
    
    # Loop de sensibilidad para negativos
    all_results <- list()
    
    # Resultado base (0% reducción) - con dataset original completo
    base_result <- tryCatch({
      
      # Ajuste en dataset original
      model_res <- fit_gam(
        drugA_id = rowt$drugA,
        drugB_id = rowt$drugB,
        event_id = rowt$meddra,
        ade_data = ade_raw_dt,
        spline_individuales = spline_individuales,
        include_sex = include_sex,
        include_stage_sex = include_stage_sex,
        k_spline = k_spline,
        bs_type = bs_type,
        select = select,
        nichd_spline = nichd_spline
      )
      
      classic_res <- calculate_classic_ior(
        drugA_id = rowt$drugA,
        drugB_id = rowt$drugB,
        event_id = rowt$meddra,
        ade_data = ade_raw_dt
      )
      
      classic_reri <- calculate_classic_reri(
        drugA_id = rowt$drugA,
        drugB_id = rowt$drugB,
        event_id = rowt$meddra,
        ade_data = ade_raw_dt  
      )
      
      if (!model_res$success) {
        data.table(
          triplet_id = rowt$triplet_id,
          drugA = rowt$drugA,
          drugB = rowt$drugB,
          meddra = rowt$meddra,
          type = "negative",
          reduction_pct = 0,
          N = counts_neg$n_events_coadmin,
          model_success = FALSE,
          n_events = counts_neg$n_events_total,
          n_coadmin = counts_neg$n_coadmin,
          n_stages_significant = NA_integer_,
          max_ior = NA_real_,
          mean_ior = NA_real_,
          model_aic = NA_real_,
          stage = list(1:7),
          log_ior = list(rep(NA_real_, 7)),
          log_ior_lower90 = list(rep(NA_real_, 7)),
          ior_values = list(rep(NA_real_, 7)),
          formula_used = if(!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
          message = if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_,
          classic_success = classic_res$success,
          log_ior_classic = list(rep(NA_real_, 7)),
          log_ior_classic_lower90 = list(rep(NA_real_, 7)),
          ior_classic = list(rep(NA_real_, 7)),
          reri_classic_success = FALSE,
          RERI_classic = list(rep(NA_real_, 7)),
          RERI_classic_lower90 = list(rep(NA_real_, 7)),
          RERI_classic_upper90 = list(rep(NA_real_, 7)),
          RERI_classic_se = list(rep(NA_real_, 7))
        )
      } else {
        data.table(
          triplet_id = rowt$triplet_id,
          drugA = rowt$drugA,
          drugB = rowt$drugB,
          meddra = rowt$meddra,
          type = "negative",
          reduction_pct = 0,
          N = model_res$n_events_coadmin,
          model_success = TRUE,
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
          RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7))
        )
      }
    }, error = function(e) {
      data.table(
        triplet_id = rowt$triplet_id,
        drugA = rowt$drugA,
        drugB = rowt$drugB,
        meddra = rowt$meddra,
        type = "negative",
        reduction_pct = 0,
        N = counts_neg$n_events_coadmin,
        model_success = FALSE,
        n_events = counts_neg$n_events_total,
        n_coadmin = counts_neg$n_coadmin,
        error_msg = paste("Error:", e$message)
      )
    })
    
    all_results[[1]] <- base_result
    
    # Loop sobre niveles de reducción
    for (red_pct in reduction_levels) {
      
      # recalcular conteos con el dataset reducido
      reduced_result <- tryCatch({
        # Reducir dataset original
        ade_reduced <- reduce_dataset_by_stage(ade_raw_dt, red_pct, seed = 9427 + idx)
        counts_red <- tryCatch({
             calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra)
        }, error = function(e) list(n_events_coadmin=NA, n_events_total=NA, n_coadmin=NA))
        
        model_res <- fit_gam(
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
        
        classic_res <- calculate_classic_ior(
          drugA_id = rowt$drugA,
          drugB_id = rowt$drugB,
          event_id = rowt$meddra,
          ade_data = ade_reduced
        )
        
        classic_reri <- calculate_classic_reri(
          drugA_id = rowt$drugA,
          drugB_id = rowt$drugB,
          event_id = rowt$meddra,
          ade_data = ade_reduced
        )
        
        if (!model_res$success) {
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = counts_red$n_events_coadmin,
            model_success = FALSE,
            n_events = counts_red$n_events_total,
            n_coadmin = counts_red$n_coadmin,
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
            formula_used = if(!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
            error_msg = if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_
          )
        } else {
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = model_res$n_events_coadmin,
            model_success = TRUE,
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
            RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7))
          )
        }
      }, error = function(e) {
        counts_red <- tryCatch(
          calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra),
          error = function(e2) list(n_events_coadmin = NA, n_events_total = NA, n_coadmin = NA)
        )
        data.table(
          triplet_id = rowt$triplet_id,
          drugA = rowt$drugA,
          drugB = rowt$drugB,
          meddra = rowt$meddra,
          type = "negative",
          reduction_pct = red_pct,
          N = counts_red$n_events_coadmin,
          model_success = FALSE,
          n_events = counts_red$n_events_total,
          n_coadmin = counts_red$n_coadmin,
          error_msg = paste("Error en reducción", red_pct, ":", e$message)
        )
      })
      
      all_results[[length(all_results) + 1]] <- reduced_result
      
      rm(ade_reduced); gc(verbose = FALSE)
    }
    
    combined_results <- rbindlist(all_results, fill = TRUE)
    return(combined_results)
  }
  
  stopCluster(cl)
  
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)
  batch_results_normalized <- lapply(batch_results_clean, function(x) {
    if (is.data.table(x) && "reduction_pct" %in% names(x)) {
      x
    } else {
      normalize_triplet_result(x)
    }
  })
  batch_dt <- rbindlist(batch_results_normalized, fill = TRUE)
  
  negatives_scores_list[[batch]] <- batch_dt
  
  rm(batch_results, batch_results_clean, batch_results_normalized, batch_dt)
  gc()
  
  message("Lote completado: ", sum(negatives_scores_list[[batch]]$model_success & 
                                     negatives_scores_list[[batch]]$reduction_pct == 0, na.rm = TRUE), " exitosos (base)")
}

negatives_scores <- rbindlist(negatives_scores_list, fill = TRUE)
rm(negatives_scores_list)
gc()

# guardado de resultados negativos por nivel de reducción
for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)
  
  subset_data <- negatives_scores[reduction_pct == red_pct]
  
  if (nrow(subset_data) > 0) {
    # Preparar CSV
    subset_csv <- copy(subset_data)
    
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]
    
    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }
    
    fwrite(subset_csv, paste0(output_dir, "negative_triplets_results", suffix_file, ".csv"))
    saveRDS(subset_data, paste0(output_dir, "negative_triplets_results", suffix_file, ".rds"))
  }
}

gc()

message("\nNegativos exitosos (base): ", sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0, na.rm = TRUE))

################################################################################
# Comparación 
################################################################################

pos_success <- positives_scores[model_success == TRUE & injection_success == TRUE & reduction_pct == 0]
neg_success <- negatives_scores[model_success == TRUE & reduction_pct == 0]

message("\nTriplets exitosos (base):")
message(" Positivos: ", nrow(pos_success))
message(" Negativos: ", nrow(neg_success))

if (nrow(pos_success) > 0 && nrow(neg_success) > 0) {
  
  pos_expanded <- pos_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type, dynamic)]
  
  neg_expanded <- neg_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type)]
  
  message("\nEstadísticas de log-IOR (base):")
  message("  Positivos. Media: ", round(mean(pos_expanded$log_ior, na.rm = TRUE), 4))
  message("  Negativos.  Media: ", round(mean(neg_expanded$log_ior, na.rm = TRUE), 4))
  
  if (nrow(pos_expanded) > 30 && nrow(neg_expanded) > 30) {
    test_result <- wilcox.test(
      pos_expanded$log_ior,
      neg_expanded$log_ior,
      alternative = "greater"
    )
    message("\nWilcoxon Test (Positivos > Negativos):")
    message("  p-value: ", format.pval(test_result$p.value, digits = 3))
  }
  
  all_expanded <- rbind(
    pos_expanded[, .(triplet_id, type, log_ior, stage)],
    neg_expanded[, .(triplet_id, type, log_ior, stage)]
  )
  
  p1 <- ggplot(all_expanded, aes(x = type, y = log_ior, fill = type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Distribución de Log-IOR por tipo de Triplete (Base)",
      x = "Tipo",
      y = "Log-IOR"
    )
  
  ggsave(paste0(output_dir, "comparison_log_ior_distribution.png"),
         p1, width = 8, height = 6, dpi = 300)
  
}

################################################################################
# Resumen final
################################################################################

summary_stats <- data.table(
  metric = c(
    "config_min_reports_triplet",
    "config_min_nichd_with_rep",
    "config_all_nichd_rep",
    "reduction_levels_tested",
    "n_positive_total",
    "n_positive_injected",
    "n_positive_modeled_base",
    "n_negative_total",
    "n_negative_modeled_base",
    "mean_ior_positive_base",
    "mean_ior_negative_base",
    "mean_stages_sig_positive_base",
    "mean_stages_sig_negative_base"
  ),
  value = c(
    min_reports_triplet,
    min_nichd_with_rep,
    as.numeric(all_nichd_rep),
    paste(reduction_levels, collapse = ","),
    nrow(pos_meta),
    sum(positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    nrow(selected_negatives),
    sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE)
  )
)

print(summary_stats)
fwrite(summary_stats, paste0(output_dir, "summary_statistics.csv"))

################################################################################
# Análisis de sensibilidad agregado 
################################################################################

sensitivity_summary <- rbind(
  positives_scores[, .(
    type = "positive",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct],
  negatives_scores[, .(
    type = "negative",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct]
)

print(sensitivity_summary)
fwrite(sensitivity_summary, paste0(output_dir, "sensitivity_analysis_summary.csv"))

################################################################################
# Verificación de reproductibilidad
################################################################################

# esto es para ver reproductibilidad con 11_augmentation_base que es versión light sin reducción
# Hash de verificación solo para casos base (reducción 0%)
verification_data <- positives_scores[
  reduction_pct == 0 & injection_success == TRUE,
  .(
    triplet_id,
    drugA,
    drugB, 
    meddra,
    dynamic,
    fold_change,
    n_injected,
    n_coadmin,
    # Promedio de log-IOR para comparación
    mean_log_ior = sapply(log_ior, function(x) mean(unlist(x), na.rm = TRUE))
  )
][order(triplet_id)]

# Guardar
fwrite(verification_data, paste0(output_dir, "injection_verification_hash.csv"))

################################################################################
# Detección de dinámicas
################################################################################
# Para ver si el modelo efectivamente detecta la forma de la dinámica inyectada

# Para no correr todo el análisis de 0 
#ruta_pos_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
#positives_scores <- readRDS(ruta_pos_results)

# positivos exitosos
pos_for_dynamics <- positives_scores[
  model_success == TRUE & 
  injection_success == TRUE & 
  !is.na(dynamic)
]

# Expando datos por etapa
pos_dynamics_expanded <- pos_for_dynamics[, {
  stages <- unlist(stage)
  log_iors <- unlist(log_ior)
  
  # valido longitudes
  n <- min(length(stages), length(log_iors))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      log_ior = log_iors[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

# agrego stage_name después de expandir
pos_dynamics_expanded[, stage_name := niveles_nichd[stage]]

# calculo el promedio de log-IOR, clasificado por dinámica y etapa
dynamics_summary <- pos_dynamics_expanded[, .(
  mean_log_ior = mean(log_ior, na.rm = TRUE),
  sd_log_ior = sd(log_ior, na.rm = TRUE),
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

# agregar stage_name a summary
dynamics_summary[, stage_name := niveles_nichd[stage]]

# Tomo uniform como baseline para calcular la diferencia
uniform_baseline <- dynamics_summary[dynamic == "uniform", .(
  stage, 
  baseline_log_ior = mean_log_ior
)]

# diferencias vs uniform
dynamics_diff <- merge(
  dynamics_summary[dynamic != "uniform"],
  uniform_baseline,
  by = "stage",
  all.x = TRUE
)

dynamics_diff[, log_ior_diff := mean_log_ior - baseline_log_ior]

print(dynamics_diff[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_log_ior, 
  baseline = baseline_log_ior,
  difference = log_ior_diff
)])

################################################################################
# Bootstrap para intervalos de confianza
################################################################################

n_boot <- 100

# aplico bootstrap a todas las combinaciones
dynamics_nonuniform <- unique(pos_dynamics_expanded[dynamic != "uniform", dynamic])
stages <- 1:7

bootstrap_results <- rbindlist(pblapply(dynamics_nonuniform, function(dyn) {
  rbindlist(lapply(stages, function(s) {
    boot_res <- bootstrap_dynamic_diff(pos_dynamics_expanded, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

# merge con datos principales
dynamics_with_ci <- merge(
  dynamics_diff,
  bootstrap_results,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci, paste0(output_dir, "dynamics_recovery_analysis.csv"))

recovery_stats <- dynamics_with_ci[, .(
  mean_difference = mean(log_ior_diff, na.rm = TRUE),
  max_difference = max(abs(log_ior_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats)

################################################################################
# Gráfico
################################################################################

# colores para dinámicas
dynamic_colors <- c(
  "increase" = "#E41A1C",
  "decrease" = "#377EB8", 
  "plateau" = "#4DAF4A",
  "inverse_plateau" = "#984EA3"
)

nichd_labels <- c(
  "term_neonatal" = "Neonato a term.",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temp.",
  "late_adolescence" = "Adolescencia tardía"
)

# Diferencias vs uniform con intervalo de confianza
p_dynamics_diff <- ggplot(
  dynamics_with_ci[!is.na(mean_diff)],
  aes(x = stage_name, y = log_ior_diff, color = dynamic, group = dynamic)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  labs(
    x = "Etapa",
    y = "Δ Log(IOR) (vs uniforme)",
    color = "Dinámica inyectada",
    fill = "Dinámica inyectada"
  ) 

ggsave(
  paste0(output_dir, "fig_dynamics_recovery.png"),
  p_dynamics_diff,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_dynamics_diff)

################################################################################
# Detección de dinámicas con RERI
################################################################################

# Para ver si el modelo efectivamente detecta la forma de la dinámica inyectada 
# usando RERI del GAM (riesgo absoluto) en lugar de IOR (odds ratio)

# positivos exitosos 
pos_for_dynamics <- positives_scores[model_success == TRUE & injection_success == TRUE & !is.na(dynamic)]

# expansión de datos por etapa incluyendo RERI del GAM
pos_dynamics_expanded_reri <- pos_for_dynamics[, {
  stages <- unlist(stage)
  reri_vals <- unlist(reri_values)  # RERI del GAM
  
  # valido longitudes
  n <- min(length(stages), length(reri_vals))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      reri = reri_vals[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

# stage_name después de expandir
pos_dynamics_expanded_reri[, stage_name := niveles_nichd[stage]]

# calculo el promedio de RERI, clasificado por dinámica y etapa
dynamics_summary_reri <- pos_dynamics_expanded_reri[, .(
  mean_reri = mean(reri, na.rm = TRUE),
  sd_reri = sd(reri, na.rm = TRUE),
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

# stage_name en summary
dynamics_summary_reri[, stage_name := niveles_nichd[stage]]

# uniform como baseline para calcular la diferencia
uniform_baseline_reri <- dynamics_summary_reri[dynamic == "uniform", .(
  stage, 
  baseline_reri = mean_reri
)]

# diferencias vs uniform (delta RERI)
dynamics_diff_reri <- merge(
  dynamics_summary_reri[dynamic != "uniform"],
  uniform_baseline_reri,
  by = "stage",
  all.x = TRUE
)

dynamics_diff_reri[, reri_diff := mean_reri - baseline_reri]

print(dynamics_diff_reri[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_reri, 
  baseline = baseline_reri,
  difference = reri_diff
)])

################################################################################
# Bootstrap para IC (RERI)
################################################################################

# aplico bootstrap a todas las combinaciones
dynamics_nonuniform_reri <- unique(pos_dynamics_expanded_reri[dynamic != "uniform", dynamic])
stages_reri <- 1:7

bootstrap_results_reri <- rbindlist(pblapply(dynamics_nonuniform_reri, function(dyn) {
  rbindlist(lapply(stages_reri, function(s) {
    boot_res <- bootstrap_dynamic_diff_reri(pos_dynamics_expanded_reri, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

# merge con datos principales
dynamics_with_ci_reri <- merge(
  dynamics_diff_reri,
  bootstrap_results_reri,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci_reri[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci_reri, paste0(output_dir, "dynamics_recovery_analysis_reri.csv"))

recovery_stats_reri <- dynamics_with_ci_reri[, .(
  mean_difference = mean(reri_diff, na.rm = TRUE),
  max_difference = max(abs(reri_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats_reri)

################################################################################
# Gráfico de Delta RERI (vs uniform)
################################################################################F

# Diferencias vs uniform con intervalo de confianza para RERI

p_dynamics_diff_reri <- ggplot(
  dynamics_with_ci_reri[!is.na(mean_diff)],
  aes(x = stage_name, y = reri_diff, color = dynamic, group = dynamic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  labs(
    x = "Etapa",
    y = "Δ RERI (vs uniforme)",
    color = "Dinámica inyectada",
    fill = "Dinámica inyectada"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
  ggsave(
    paste0(output_dir, "fig_dynamics_recovery_reri.png"),
    p_dynamics_diff_reri,
    width = 12,
    height = 8,
    dpi = 300
)

print(p_dynamics_diff_reri)
