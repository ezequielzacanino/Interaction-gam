################################################################################
# Script de generación de datos semisintéticos
# Script 01_augmentation
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

# Parámetros para inyección
min_reports_triplet <- 10
n_pos <- 500
n_neg <- 2500
lambda_fc <- 0.75
dinamicas <- c("uniform","increase","decrease","plateau","inverse_plateau")
n_cores <- max(1, floor(detectCores() * 0.75)) 
z90 <- qnorm(0.95)

# Parámetros de fórmula para GAM
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                 
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML" 

output_dir <- "./augmentation_results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
source("giangreco_theme.R")
theme_set(theme_giangreco())



################################################################################
# Carga de datos
################################################################################


ade_raw_dt <- fread(ruta_ade_raw)

# valido columnas requeridas
cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")

# Si se va a usar sex, agrego a validación
if (include_sex) {
  cols_req <- c(cols_req, "sex")
}

# Proceso nichd
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

# Proceso sex si está presente
if (include_sex) {
  # Estandarizo valores
  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  
  # Convierto a factor con niveles estándar
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  
  # Estadística descriptiva
  sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
  message("\nDistribución de sexo:")
  print(sex_summary)
}

message(sprintf("Dataset %s filas", format(nrow(ade_raw_dt), big.mark = ",")))


################################################################################
# Construcción de tripletes candidatos
################################################################################

drugs_by_report <- unique(ade_raw_dt[!is.na(atc_concept_id), 
                                     .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[!is.na(meddra_concept_id), 
                                      .(safetyreportid, meddra_concept_id)])

# Genero listas por reporte
reports <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), 
                               by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), 
                                 by = safetyreportid]

# Crear reports_meta  
reports_meta <- merge(reports, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
reports_meta[is.na(drugs), drugs := list(integer(0))]
reports_meta[is.na(events), events := list(integer(0))]

# report_combo para construcción de tripletes 
report_combo <- copy(reports_meta)  

# Genero tripletes
message("Enumeración de tripletes")

triplets_list <- pblapply(seq_len(nrow(report_combo)), function(i) {
  rowi <- report_combo[i]
  make_triplets_per_report(
    rowi$drugs[[1]], 
    rowi$events[[1]], 
    rowi$safetyreportid, 
    rowi$nichd
  )
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

# Conteos por triplete único
trip_counts <- unique(triplets_dt[, .(drugA, drugB, meddra, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra)
]

# Filtrar candidatos
candidatos_pos <- trip_counts[N >= min_reports_triplet]
message("Tripletes candidatos: ", nrow(candidatos_pos))

positivos_sel <- candidatos_pos[sample(.N, n_pos)]


################################################################################
# Selección de positivos
################################################################################

# Asignación de dinámicas y fold-changes
n_per_dynamic <- ceiling(nrow(positivos_sel) / length(dinamicas))
din_for_triplets <- rep(dinamicas, length.out = nrow(positivos_sel))

fold_changes <- sample_fold_change(nrow(positivos_sel), lambda = lambda_fc)

pos_meta <- positivos_sel[, .(drugA, drugB, meddra, N)]
pos_meta[, `:=`(
  dynamic = din_for_triplets,
  fold_change = fold_changes,
  triplet_id = 1:.N
)]

# Flageo los 30 tripletes con mayor cantidad de reportes para utilizarlos después en análisis de escasez
top30_ids <- pos_meta[order(-N)][1:min(.N, 30), triplet_id]
pos_meta[, is_top30 := triplet_id %in% top30_ids]

fwrite(pos_meta, paste0(output_dir, "positive_triplets_metadata.csv"))


################################################################################
# Procesado en paralelo de tripletes positivos
################################################################################

message("Procesado de positivos ", nrow(pos_meta), " tripletes")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# objetos necesarios
clusterExport(cl, c("inject_signal", "fit_differential_gam",
                    "ade_raw_dt", "niveles_nichd", "dynamic_fun", "z90",
                    "pos_meta", "normalize_triplet_result", "spline_individuales", "include_sex",
                    "include_stage_sex", "k_spline", "nichd_spline", "bs_type", "select",
                    "output_dir"), 
              envir = environment())

clusterEvalQ(cl, {
  library(data.table)
  library(mgcv)
})

# foreach para ir viendo barra de progreso
pb <- txtProgressBar(max = nrow(pos_meta), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

positives_results <- foreach(
  i = 1:nrow(pos_meta),
  .packages = c("data.table", "mgcv"),
  .errorhandling = "pass",
  .options.snow = opts,
  .verbose = FALSE
) %dopar% {
  
  rowt <- pos_meta[i]
  
  ###########
  # Creación de dataset aumentado INDEPENDIENTE (cada triplete inyectado se introduce en una copia del dataset ORIGINAL)
  # Si pongo todos los tripletes positivos inyectados JUNTOS en el dataset original antes de ajustar contamino mucho los datos
  ###########

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
    list(success = FALSE, n_injected = 0, n_coadmin = 0, 
         message = paste("Error en inyección:", e$message))   # Chequeo por si falla
  })
  
  if (inj_result$success && rowt$is_top30) {        # Guardado de datasets aumentados para los 30 tripletes más reportados
    # nombre: dataset_idA_idB.csv 
    # Agrego triplet_id para evitar sobrescritura si el par drogas se repite con otro evento
    file_name <- paste0("dataset_sensitivity_", rowt$drugA, "_", rowt$drugB, ".csv")
    full_path <- file.path(output_dir, file_name)
    
    # guardo el dataset aumentado que está en inj_result$ade_aug
    fwrite(inj_result$ade_aug, full_path)
  }
  inj_success <- inj_result$success
  n_injected_val <- inj_result$n_injected
  n_coadmin_val <- inj_result$n_coadmin
  diag_data <- list(inj_result$diagnostics)
  inj_message <- if(!is.null(inj_result$message)) inj_result$message else NA_character_

  if (!inj_result$success) {
    result <- data.table(
      triplet_id = i,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "positive",
      N = 0,
      dynamic = rowt$dynamic,
      fold_change = rowt$fold_change,
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
      diagnostics = diag_data,
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = NA_character_,
      message = inj_message
    )
    rm(inj_result); gc(verbose = FALSE)
    return(result)
  }
  
  ###########
  # Ajuste de modelo en dataset aumentado
  ###########
  
  model_res <- tryCatch({
    fit_differential_gam(
    drugA_id = rowt$drugA,
    drugB_id = rowt$drugB,
    event_id = rowt$meddra,
    ade_data = inj_result$ade_aug,
    spline_individuales = spline_individuales,
    include_sex = include_sex,
    include_stage_sex = include_stage_sex,
    k_spline = k_spline,
    bs_type = bs_type,
    select = select,
    nichd_spline = nichd_spline
    )
  }, error = function(e) {
    list(success = FALSE, n_events = 0, n_coadmin = inj_result$n_coadmin,
         error_msg = paste("Error en modelo:", e$message))        # Chequeo por si falla
  })
  
  ###########
  # Guardo diagnósticos y resultados antes de borrar inj_result
  ###########

  model_success <- model_res$success
  n_events_val <- if(!is.null(model_res$n_events)) model_res$n_events else NA_integer_
  n_coadmin_model <- if(!is.null(model_res$n_coadmin)) model_res$n_coadmin else n_coadmin_val
  error_msg_val <- if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_
  formula_val <- if(model_success && !is.null(model_res$formula_used)) model_res$formula_used else NA_character_
  
  rm(inj_result); gc(verbose = FALSE)
  
  if (!model_res$success) {
    result <- data.table(
      triplet_id = i,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "positive",
      N = 0,
      dynamic = rowt$dynamic,
      fold_change = rowt$fold_change,
      model_success = FALSE,
      injection_success = TRUE,
      n_injected = n_injected_val,
      n_coadmin = n_coadmin_model,
      n_events = n_events_val,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      diagnostics = diag_data,
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = formula_val,
      error_msg = error_msg_val
    )
    rm(model_res); gc(verbose = FALSE)
    return(result)
  }
  
  ###########
  # resultados completos
  ###########

  result <- data.table(
    triplet_id = i,
    drugA = rowt$drugA,
    drugB = rowt$drugB,
    meddra = rowt$meddra,
    type = "positive",
    N = model_res$n_events,
    dynamic = rowt$dynamic,
    fold_change = rowt$fold_change,
    model_success = TRUE,
    injection_success = TRUE,
    n_injected = n_injected_val,
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
    diagnostics = diag_data,
    spline_individuales = spline_individuales,
    nichd_spline = nichd_spline,
    include_sex = include_sex,
    include_stage_sex = include_stage_sex,
    k_spline = k_spline,
    bs_type = bs_type,
    select = select,
    formula_used = formula_val
  )
  
  rm(model_res); gc(verbose = FALSE)
  return(result)
}

close(pb)
stopCluster(cl)

# normalizo y combino resultados

positives_results_clean <- Filter(function(x) !inherits(x, "error"), positives_results)
positives_results_normalized <- lapply(positives_results_clean, normalize_triplet_result)
positives_scores <- rbindlist(positives_results_normalized, fill = TRUE)


pos_success <- positives_scores[injection_success == TRUE]



################################################################################
# Guardado de resultados
################################################################################

# Guardo RDS 
saveRDS(positives_scores, paste0(output_dir, "positive_triplets_results.rds"))

# Guardo en versión csv
positives_scores_csv <- copy(positives_scores)

# Elimino columna diagnostics
if ("diagnostics" %in% names(positives_scores_csv)) {
  positives_scores_csv[, diagnostics := NULL]
}

# aplanado de listas
# Detección de columnas tipo lista (stage, log_ior)
list_cols <- names(positives_scores_csv)[sapply(positives_scores_csv, is.list)]


# Convierto listas a texto
for (col in list_cols) {
  positives_scores_csv[, (col) := sapply(get(col), function(x) {
    # Si es NULL o vacío, devolver NA (string) para que fwrite no falle
    if (is.null(x) || length(x) == 0) return(NA_character_)
    # Si tiene datos, unirlos con comas
    paste(x, collapse = ",")
  })]
}

# Guardo csv
fwrite(positives_scores_csv, paste0(output_dir, "positive_triplets_results.csv"))

# limpio memoria
rm(positives_results, positives_scores_csv)
gc()

message("\nPositivos exitosos: ", 
        sum(positives_scores$model_success & positives_scores$injection_success, na.rm = TRUE))



################################################################################
# Selección de datos negativos
################################################################################

# Son excluyentes de positivos
# Si coincidiese un triplete elegido para negativo con uno positivo, que justo se inyectó con un fold-change muy pequeño
# Habría demasiado solapamiento

# No hay garantía de que los negativos sean realmente negativos
# Son simplemente elegidos al azar
# La "garantía" es que la incidencia de interacciones farmacológicas sinérgicas es "baja"
# Por lo que, probabilisticamente, agarrar tripletes al azar tiene muchisimas más chances de agarrar negativos

# identificador único para tripletes positivos
pos_triplet_ids <- paste(
  pmin(pos_meta$drugA, pos_meta$drugB),
  pmax(pos_meta$drugA, pos_meta$drugB),
  pos_meta$meddra,
  sep = "_"
)

# Identifico tripletes candidatos que NO estén en positivos
candidatos_neg <- candidatos_pos[
  !paste(pmin(drugA, drugB), pmax(drugA, drugB), meddra, sep = "_") %in% pos_triplet_ids
]

message("Tripletes negativos disponibles: ", nrow(candidatos_neg))

# Muestrear negativos
n_neg_sample <- min(n_neg, nrow(candidatos_neg))
set.seed(456)
selected_negatives <- candidatos_neg[sample(.N, n_neg_sample)]
selected_negatives[, `:=`(
  type = "negativo",
  triplet_id = (nrow(pos_meta) + 1):(nrow(pos_meta) + n_neg_sample)
)]

fwrite(selected_negatives, paste0(output_dir, "negative_triplets_metadata.csv"))

message("Seleccionados ", nrow(selected_negatives), " negativos")



################################################################################
# Procesado de negativos
################################################################################



message("Procesado de negativos: ", nrow(selected_negatives), " tripletes")


# Proceso negativos en lotes para evitar que crashee 
batch_size <- 100
n_batches <- ceiling(nrow(selected_negatives) / batch_size)

negatives_scores_list <- list()

for (batch in 1:n_batches) {
  
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, nrow(selected_negatives))
  batch_indices <- start_idx:end_idx
  
  message("\nlote ", batch, " / ", n_batches, 
          " (tripletes ", start_idx, "-", end_idx, ")")
  
  # Configuración de cluster para lote
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  clusterExport(cl, c("fit_differential_gam", "ade_raw_dt", "z90",
                      "selected_negatives", "normalize_triplet_result", "spline_individuales", "include_sex",
                      "include_stage_sex", "k_spline", "nichd_spline", "bs_type", "select", "method"), 
                envir = environment())
  
  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
  })
  
  # Procesado del lote
  batch_results <- foreach(
    idx = batch_indices,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE
  ) %dopar% {
    
    rowt <- selected_negatives[idx]
    
    ###########
    # Ajuste del modelo en dataset original
    ###########
    model_res <- tryCatch({
      fit_differential_gam(
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
    }, error = function(e) {
      list(success = FALSE, n_events = 0, n_coadmin = 0,
           error_msg = paste("Error:", e$message))            # Chequeo por si falla
    })
    
    if (!model_res$success) {
      result <- data.table(
        triplet_id = rowt$triplet_id,
        drugA = rowt$drugA,
        drugB = rowt$drugB,
        meddra = rowt$meddra,
        type = "negative",
        N = rowt$N,
        model_success = FALSE,
        n_coadmin = model_res$n_coadmin,
        n_events = model_res$n_events,
        n_stages_significant = NA_integer_,
        max_ior = NA_real_,
        mean_ior = NA_real_,
        model_aic = NA_real_,
        stage = list(1:7),
        log_ior = list(rep(NA_real_, 7)),
        log_ior_lower90 = list(rep(NA_real_, 7)),
        ior_values = list(rep(NA_real_, 7)),
        formula_used = if(model_res$success) model_res$formula_used else NA_character_ ,
        message = inj_result$message
      )
      
      rm(model_res)
      gc(verbose = FALSE)
      
      return(result)
    }
    
    ###########
    # resultados completos
    ###########
    result <- data.table(
      triplet_id = rowt$triplet_id,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "negative",
      N = rowt$N,
      model_success = TRUE,
      n_coadmin = model_res$n_coadmin,
      n_events = model_res$n_events,
      n_stages_significant = model_res$n_stages_significant,
      max_ior = model_res$max_ior,
      mean_ior = model_res$mean_ior,
      model_aic = model_res$model_aic,
      stage = list(1:7),
      log_ior = list(model_res$log_ior),
      log_ior_lower90 = list(model_res$log_ior_lower90),
      ior_values = list(model_res$ior_values)
    )
    
    rm(model_res)
    gc(verbose = FALSE)
    
    return(result)
  }
  
  stopCluster(cl)
  
  # Limpieza y guardado de lote
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)
  batch_results_normalized <- lapply(batch_results_clean, normalize_triplet_result)
  batch_dt <- rbindlist(batch_results_normalized, fill = TRUE)
  
  negatives_scores_list[[batch]] <- batch_dt
  
  # Libero memoria del lote
  rm(batch_results, batch_results_clean, batch_results_normalized, batch_dt)
  gc()
  
  message("Lote completado: ", sum(negatives_scores_list[[batch]]$model_success, na.rm = TRUE), 
          " exitosos")
}

# Combinación de todos los lotes
negatives_scores <- rbindlist(negatives_scores_list, fill = TRUE)
rm(negatives_scores_list)
gc()

# Guardar negativos
negatives_scores_csv <- copy(negatives_scores)
negatives_scores_csv[, stage := sapply(stage, function(x) paste(x, collapse = ","))]
negatives_scores_csv[, log_ior := sapply(log_ior, function(x) paste(x, collapse = ","))]
negatives_scores_csv[, log_ior_lower90 := sapply(log_ior_lower90, function(x) paste(x, collapse = ","))]
negatives_scores_csv[, ior_values := sapply(ior_values, function(x) paste(x, collapse = ","))]

fwrite(negatives_scores_csv, paste0(output_dir, "negative_triplets_results.csv"))
saveRDS(negatives_scores, paste0(output_dir, "negative_triplets_results.rds"))

rm(negatives_scores_csv)
gc()

message("\nNegativos exitosos: ", sum(negatives_scores$model_success, na.rm = TRUE))



################################################################################
# Creación de pool nulo 
################################################################################

# Es por muestreo aleatorio simple, total después se permutan 


# Parámetros de muestreo
n_null_reports <- 100000  
set.seed(789)           

# Obtengo todos los reportes únicos disponibles
all_reports <- unique(ade_raw_dt$safetyreportid)


# No se excluyen reportes que puedan estar en positivos porque:
# La permutación de drogas/eventos destruye las asociaciones reales
# Maximiza el tamaño de muestra disponible para la distribución nula
selected_null_reports <- sample(all_reports, n_null_reports)

# metadata del pool nulo (safetyreportid, etapa NICHD)
null_pool_meta <- unique(ade_raw_dt[
  safetyreportid %in% selected_null_reports,
  .(safetyreportid, nichd, nichd_num)
])

# Verificación de distribución de reportes por etapa
stage_distribution <- null_pool_meta[, .N, by = nichd][order(nichd)]
message("\nDistribución de reportes por etapa NICHD:")
print(stage_distribution)

# Habría que agregar más diagnósticos para ver que los reportes seleccionados sean representativos (en etapa sobretodo)
# Al ser un dataset "grande", probabilisticamente seguro es representativo pero bueno


# Guardo metadata
fwrite(null_pool_meta, paste0(output_dir, "null_pool_reports_metadata.csv"))

message("Total de reportes: ", nrow(null_pool_meta))



################################################################################
# Comparación
################################################################################


# Filtro solo tripletes exitosos
pos_success <- positives_scores[model_success == TRUE & injection_success == TRUE]
neg_success <- negatives_scores[model_success == TRUE]

message("\nTriplets exitosos:")
message(" Positivos: ", nrow(pos_success))
message(" Negativos: ", nrow(neg_success))

if (nrow(pos_success) > 0 && nrow(neg_success) > 0) {
  
  # Expando listas para análisis
  pos_expanded <- pos_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    # unas lineas para manejar casos donde las listas puedan tener longitudes diferentes
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
  
  # estadísticas descriptivas
  message("\nEstadísticas de log-IOR:")
  message("  Positivos. Media: ", round(mean(pos_expanded$log_ior, na.rm = TRUE), 4))
  message("  Negativos.  Media: ", round(mean(neg_expanded$log_ior, na.rm = TRUE), 4))
  
  # Test estadístico para ver que negativos y positivos no son iguales
  if (nrow(pos_expanded) > 30 && nrow(neg_expanded) > 30) {
    test_result <- wilcox.test(
      pos_expanded$log_ior,
      neg_expanded$log_ior,
      alternative = "greater"
    )
    message("\nWilcoxon Test (Positivos > Negativos):")
    message("  p-value: ", format.pval(test_result$p.value, digits = 3))
  }
  
  ###########
  # Visualización
  ###########
  all_expanded <- rbind(
    pos_expanded[, .(triplet_id, type, log_ior, stage)],
    neg_expanded[, .(triplet_id, type, log_ior, stage)]
  )
  
  p1 <- ggplot(all_expanded, aes(x = type, y = log_ior, fill = type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Distribución de Log-IOR por tipo de Triplete",
      x = "Tipo",
      y = "Log-IOR"
    )
  
  ggsave(paste0(output_dir, "comparison_log_ior_distribution.png"),
         p1, width = 8, height = 6, dpi = 300)
  
}

################################################################################
# 8) RESUMEN FINAL
################################################################################

message("Resumen final")

summary_stats <- data.table(
  metric = c(
    "n_positive_total",
    "n_positive_injected",
    "n_positive_modeled",
    "n_negative_total",
    "n_negative_modeled",
    "mean_ior_positive",
    "mean_ior_negative",
    "mean_stages_sig_positive",
    "mean_stages_sig_negative"
  ),
  value = c(
    nrow(pos_meta),
    sum(positives_scores$injection_success, na.rm = TRUE),
    sum(positives_scores$model_success, na.rm = TRUE),
    nrow(selected_negatives),
    sum(negatives_scores$model_success, na.rm = TRUE),
    mean(pos_success$mean_ior, na.rm = TRUE),
    mean(neg_success$mean_ior, na.rm = TRUE),
    mean(pos_success$n_stages_significant, na.rm = TRUE),
    mean(neg_success$n_stages_significant, na.rm = TRUE)
  )
)

print(summary_stats)
fwrite(summary_stats, paste0(output_dir, "summary_statistics.csv"))




