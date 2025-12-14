################################################################################
# Script de análisis de poder estadístico
# Script 04_power_analysis.R
################################################################################

library(data.table)
library(tidyverse)
library(mgcv)
library(pROC)
library(parallel)
library(pbapply)
library(ggplot2)
library(RColorBrewer)
set.seed(9427)

################################################################################
# Configuración
################################################################################

setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)
source("giangreco_theme.R")
theme_set(theme_giangreco())

# Parámetros de análisis de poder
fold_change_levels <- seq(1.5, 5.0, by = 0.5)  # Rango de fold-changes
min_reports_levels <- seq(5, 50, by = 5)       # Rango de número de reportes
target_power <- 0.80                           # Poder objetivo (80%)
n_bootstrap <- 1000                            # Iteraciones bootstrap para CI
n_cores <- max(1, floor(detectCores() * 0.75))

# Configuración de parámetros del modelo (igual que scripts anteriores)
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                 
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML" 

# Cargar configuraciones existentes
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

# Directorios
output_dir <- paste0("./results/", suffix, "/power_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Archivos de entrada
ruta_ade_raw <- "./ade_raw.csv"
ruta_pos_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
ruta_neg_results <- paste0("./results/", suffix, "/augmentation_results/negative_triplets_results.rds")
ruta_null_thresh <- paste0("./results/", suffix, "/null_distribution_results/null_thresholds.csv")

################################################################################
# Carga de datos base
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

positives_scores <- readRDS(ruta_pos_results)
negatives_scores <- readRDS(ruta_neg_results)
null_thresholds <- fread(ruta_null_thresh)

message("Dataset cargado:")
message(sprintf("  Positivos: %d", nrow(positives_scores)))
message(sprintf("  Negativos: %d", nrow(negatives_scores)))

################################################################################
# Funciones auxiliares para análisis de poder
################################################################################

# Función para determinar si una señal es detectada considerando criterios de infinito/NA
is_signal_detected <- function(log_ior_lower90, threshold = 0, method = "GAM") {
  # Si GAM o IOR clásico dan infinito o NA, se considera NO detectado
  detected <- !is.infinite(log_ior_lower90) & !is.na(log_ior_lower90) & (log_ior_lower90 > threshold)
  return(detected)
}

# Función para calcular poder estadístico
calculate_power <- function(positive_data, negative_data, fold_change_threshold, 
                           min_reports_threshold, percentile_level = "p95") {
  
  # Filtrar datos por características
  filtered_pos <- positive_data[
    model_success == TRUE & 
    injection_success == TRUE & 
    fold_change >= fold_change_threshold & 
    n_coadmin >= min_reports_threshold &
    dynamic != "uniform"
  ]
  
  filtered_neg <- negative_data[model_success == TRUE]
  
  if (nrow(filtered_pos) == 0 || nrow(filtered_neg) == 0) {
    return(data.table(
      fold_change = fold_change_threshold,
      min_reports = min_reports_threshold,
      power_gam = NA_real_,
      power_classic = NA_real_,
      n_pos = 0,
      n_neg = 0,
      detection_rate_gam = NA_real_,
      detection_rate_classic = NA_real_
    ))
  }
  
  # Obtener umbral nulo correspondiente
  threshold_col <- paste0("threshold_", percentile_level)
  avg_threshold <- mean(null_thresholds[[threshold_col]], na.rm = TRUE)
  
  # Expandir datos positivos
  pos_expanded <- filtered_pos[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]
  
  # Expandir datos negativos  
  neg_expanded <- filtered_neg[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]
  
  # Calcular detección de señales por triplete
  pos_detected_gam <- pos_expanded[, .(
    detected_gam = any(is_signal_detected(log_ior_lower90, avg_threshold, "GAM")),
    detected_classic = any(is_signal_detected(log_ior_classic_lower90, 0, "Clásico"))
  ), by = triplet_id]
  
  neg_detected_gam <- neg_expanded[, .(
    detected_gam = any(is_signal_detected(log_ior_lower90, avg_threshold, "GAM")),
    detected_classic = any(is_signal_detected(log_ior_classic_lower90, 0, "Clásico"))
  ), by = triplet_id]
  
  # Calcular poder (tasa de verdaderos positivos)
  power_gam <- mean(pos_detected_gam$detected_gam, na.rm = TRUE)
  power_classic <- mean(pos_detected_gam$detected_classic, na.rm = TRUE)
  
  # Calcular tasas de detección
  detection_rate_gam <- mean(c(pos_detected_gam$detected_gam, !neg_detected_gam$detected_gam), na.rm = TRUE)
  detection_rate_classic <- mean(c(pos_detected_gam$detected_classic, !neg_detected_gam$detected_classic), na.rm = TRUE)
  
  return(data.table(
    fold_change = fold_change_threshold,
    min_reports = min_reports_threshold,
    power_gam = power_gam,
    power_classic = power_classic,
    n_pos = nrow(filtered_pos),
    n_neg = nrow(filtered_neg),
    detection_rate_gam = detection_rate_gam,
    detection_rate_classic = detection_rate_classic
  ))
}

################################################################################
# Cálculo de matriz de poder
################################################################################

message("\nCalculando matriz de poder...")

# Crear grid de combinaciones
power_grid <- expand.grid(
  fold_change = fold_change_levels,
  min_reports = min_reports_levels
)

# Calcular poder para cada combinación
power_results <- pbapply(power_grid, 1, function(params) {
  calculate_power(
    positive_data = positives_scores,
    negative_data = negatives_scores,
    fold_change_threshold = params[["fold_change"]],
    min_reports_threshold = params[["min_reports"]],
    percentile_level = "p95"
  )
})

power_matrix <- rbindlist(power_results)

# Guardar resultados
fwrite(power_matrix, paste0(output_dir, "power_matrix.csv"))

message(sprintf("Matriz de poder calculada para %d combinaciones", nrow(power_matrix)))

################################################################################
# Identificación de umbrales para 80% de poder
################################################################################

message("\nIdentificando umbrales para 80% de poder...")

# Encontrar combinaciones que alcanzan 80% de poder
power_80_gam <- power_matrix[power_gam >= target_power]
power_80_classic <- power_matrix[power_classic >= target_power]

if (nrow(power_80_gam) > 0) {
  min_fc_gam <- min(power_80_gam$fold_change)
  min_reports_gam <- min(power_80_gam$min_reports[power_80_gam$fold_change == min_fc_gam])
} else {
  min_fc_gam <- NA_real_
  min_reports_gam <- NA_real_
}

if (nrow(power_80_classic) > 0) {
  min_fc_classic <- min(power_80_classic$fold_change)
  min_reports_classic <- min(power_80_classic$min_reports[power_80_classic$fold_change == min_fc_classic])
} else {
  min_fc_classic <- NA_real_
  min_reports_classic <- NA_real_
}

# Guardar umbrales
thresholds_80 <- data.table(
  method = c("GAM", "IOR_Clásico"),
  min_fold_change = c(min_fc_gam, min_fc_classic),
  min_reports = c(min_reports_gam, min_reports_classic),
  target_power = target_power
)

fwrite(thresholds_80, paste0(output_dir, "power_thresholds_80.csv"))

message("Umbrales para 80% de poder:")
print(thresholds_80)

################################################################################
# Generación de heatmaps de poder
################################################################################

message("\nGenerando heatmaps de poder...")

# Heatmap para GAM
p_gam <- ggplot(power_matrix, aes(x = min_reports, y = fold_change, fill = power_gam)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("red", "yellow", "green"),
    limits = c(0, 1),
    na.value = "grey50",
    name = "Poder\nEstadístico"
  ) +
  geom_contour(aes(z = power_gam), color = "black", alpha = 0.3, bins = 5) +
  labs(
    title = "Análisis de Poder Estadístico - Modelo GAM",
    subtitle = sprintf("Línea negra: contorno de %d%% poder", target_power * 100),
    x = "Número mínimo de reportes",
    y = "Fold-change mínimo",
    fill = "Poder"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

# Agregar línea de 80% poder si existe
if (nrow(power_80_gam) > 0) {
  p_gam <- p_gam + geom_point(
    data = power_80_gam[fold_change == min_fc_gam & min_reports == min_reports_gam],
    size = 8, shape = 21, fill = "black", color = "white", stroke = 2
  )
}

ggsave(paste0(output_dir, "heatmap_power_gam.png"), p_gam, 
       width = 12, height = 8, dpi = 300)

# Heatmap para IOR Clásico
p_classic <- ggplot(power_matrix, aes(x = min_reports, y = fold_change, fill = power_classic)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("red", "yellow", "green"),
    limits = c(0, 1),
    na.value = "grey50",
    name = "Poder\nEstadístico"
  ) +
  geom_contour(aes(z = power_classic), color = "black", alpha = 0.3, bins = 5) +
  labs(
    title = "Análisis de Poder Estadístico - IOR Clásico",
    subtitle = sprintf("Línea negra: contorno de %d%% poder", target_power * 100),
    x = "Número mínimo de reportes",
    y = "Fold-change mínimo",
    fill = "Poder"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

# Agregar línea de 80% poder si existe
if (nrow(power_80_classic) > 0) {
  p_classic <- p_classic + geom_point(
    data = power_80_classic[fold_change == min_fc_classic & min_reports == min_reports_classic],
    size = 8, shape = 21, fill = "black", color = "white", stroke = 2
  )
}

ggsave(paste0(output_dir, "heatmap_power_classic.png"), p_classic, 
       width = 12, height = 8, dpi = 300)

################################################################################
# Filtrado del pool positivo según características de poder
################################################################################

message("\nFiltrando pool positivo según características de poder...")

# Aplicar filtros basados en los umbrales encontrados
if (!is.na(min_fc_gam) && !is.na(min_reports_gam)) {
  filtered_positives <- positives_scores[
    model_success == TRUE & 
    injection_success == TRUE &
    fold_change >= min_fc_gam & 
    n_coadmin >= min_reports_gam &
    dynamic != "uniform"
  ]
  
  message(sprintf("Pool positivo filtrado: %d tripletes (de %d originales)", 
                  nrow(filtered_positives), nrow(positives_scores)))
  
  # Guardar pool filtrado
  saveRDS(filtered_positives, paste0(output_dir, "filtered_positive_pool.rds"))
  
  # Metadatos del pool filtrado
  filtered_meta <- filtered_positives[, .(
    triplet_id, drugA, drugB, meddra, dynamic, fold_change, n_coadmin, n_injected
  )]
  fwrite(filtered_meta, paste0(output_dir, "filtered_positive_metadata.csv"))
  
} else {
  message("ADVERTENCIA: No se pudieron establecer umbrales para 80% de poder")
  filtered_positives <- positives_scores[
    model_success == TRUE & 
    injection_success == TRUE &
    dynamic != "uniform"
  ]
}

################################################################################
# Análisis de validación comparativa en pool filtrado
################################################################################

if (nrow(filtered_positives) > 0) {
  
  message("\nEjecutando análisis de validación comparativa en pool filtrado...")
  
  # Combinar con negativos válidos
  valid_negatives <- negatives_scores[model_success == TRUE]
  
  # Expandir datos filtrados
  pos_expanded_filtered <- filtered_positives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]
  
  neg_expanded_filtered <- valid_negatives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]
  
  # Aplicar umbrales de detección
  avg_threshold <- mean(null_thresholds$threshold_p95, na.rm = TRUE)
  
  # Clasificar señales para GAM
  pos_expanded_filtered[, signal_gam := is_signal_detected(log_ior_lower90, avg_threshold, "GAM")]
  neg_expanded_filtered[, signal_gam := is_signal_detected(log_ior_lower90, avg_threshold, "GAM")]
  
  # Clasificar señales para IOR clásico
  pos_expanded_filtered[, signal_classic := is_signal_detected(log_ior_classic_lower90, 0, "Clásico")]
  neg_expanded_filtered[, signal_classic := is_signal_detected(log_ior_classic_lower90, 0, "Clásico")]
  
  # Resumen a nivel triplete
  pos_summary <- pos_expanded_filtered[, .(
    signal_gam = any(signal_gam),
    signal_classic = any(signal_classic)
  ), by = triplet_id]
  pos_summary[, true_label := 1]
  
  neg_summary <- neg_expanded_filtered[, .(
    signal_gam = any(signal_gam),
    signal_classic = any(signal_classic)
  ), by = triplet_id]
  neg_summary[, true_label := 0]
  
  # Combinar datos para análisis
  validation_data <- rbind(pos_summary, neg_summary)
  
  ################################################################################
  # Cálculo de métricas de performance
  ################################################################################
  
  message("\nCalculando métricas de performance en pool filtrado...")
  
  # Función para calcular métricas
  calculate_validation_metrics <- function(data, method_col, method_name) {
    
    tp <- sum(data[[method_col]] == TRUE & data$true_label == 1, na.rm = TRUE)
    fn <- sum(data[[method_col]] == FALSE & data$true_label == 1, na.rm = TRUE)
    fp <- sum(data[[method_col]] == TRUE & data$true_label == 0, na.rm = TRUE)
    tn <- sum(data[[method_col]] == FALSE & data$true_label == 0, na.rm = TRUE)
    
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    ppv <- tp / (tp + fp)
    npv <- tn / (tn + fn)
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)
    
    return(data.table(
      method = method_name,
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv,
      accuracy = accuracy,
      f1_score = f1_score,
      tp = tp, fn = fn, fp = fp, tn = tn
    ))
  }
  
  # Calcular métricas para ambos métodos
  metrics_gam <- calculate_validation_metrics(validation_data, "signal_gam", "GAM")
  metrics_classic <- calculate_validation_metrics(validation_data, "signal_classic", "IOR_Clásico")
  
  metrics_combined <- rbind(metrics_gam, metrics_classic)
  
  # Guardar métricas
  fwrite(metrics_combined, paste0(output_dir, "validation_metrics_filtered.csv"))
  
  message("Métricas de validación en pool filtrado:")
  print(metrics_combined)
  
  ################################################################################
  # Métricas por etapa
  ################################################################################
  
  message("\nCalculando métricas por etapa...")
  
  # Expandir con información de etapa
  pos_expanded_with_stage <- filtered_positives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = .(triplet_id, fold_change, n_coadmin)]
  
  pos_expanded_with_stage[, true_label := 1]
  pos_expanded_with_stage[, stage_name := niveles_nichd[stage]]
  
  neg_expanded_with_stage <- valid_negatives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]
  
  neg_expanded_with_stage[, true_label := 0]
  neg_expanded_with_stage[, stage_name := niveles_nichd[stage]]
  
  # Combinar y clasificar por etapa
  all_expanded_filtered <- rbind(
    pos_expanded_with_stage,
    neg_expanded_with_stage
  )
  
  all_expanded_filtered[, signal_gam_stage := is_signal_detected(log_ior_lower90, avg_threshold, "GAM")]
  all_expanded_filtered[, signal_classic_stage := is_signal_detected(log_ior_classic_lower90, 0, "Clásico")]
  
  # Calcular métricas por etapa
  metrics_by_stage <- all_expanded_filtered[, {
    
    # GAM
    tp_gam <- sum(signal_gam_stage == TRUE & true_label == 1, na.rm = TRUE)
    fn_gam <- sum(signal_gam_stage == FALSE & true_label == 1, na.rm = TRUE)
    fp_gam <- sum(signal_gam_stage == TRUE & true_label == 0, na.rm = TRUE)
    tn_gam <- sum(signal_gam_stage == FALSE & true_label == 0, na.rm = TRUE)
    
    sens_gam <- ifelse((tp_gam + fn_gam) > 0, tp_gam / (tp_gam + fn_gam), NA_real_)
    spec_gam <- ifelse((tn_gam + fp_gam) > 0, tn_gam / (tn_gam + fp_gam), NA_real_)
    ppv_gam <- ifelse((tp_gam + fp_gam) > 0, tp_gam / (tp_gam + fp_gam), NA_real_)
    
    # IOR Clásico
    tp_classic <- sum(signal_classic_stage == TRUE & true_label == 1, na.rm = TRUE)
    fn_classic <- sum(signal_classic_stage == FALSE & true_label == 1, na.rm = TRUE)
    fp_classic <- sum(signal_classic_stage == TRUE & true_label == 0, na.rm = TRUE)
    tn_classic <- sum(signal_classic_stage == FALSE & true_label == 0, na.rm = TRUE)
    
    sens_classic <- ifelse((tp_classic + fn_classic) > 0, tp_classic / (tp_classic + fn_classic), NA_real_)
    spec_classic <- ifelse((tn_classic + fp_classic) > 0, tn_classic / (tn_classic + fp_classic), NA_real_)
    ppv_classic <- ifelse((tp_classic + fp_classic) > 0, tp_classic / (tp_classic + fp_classic), NA_real_)
    
    data.table(
      stage_name = stage_name[1],
      n_total = .N,
      n_pos = sum(true_label == 1),
      n_neg = sum(true_label == 0),
      sensitivity_gam = sens_gam,
      specificity_gam = spec_gam,
      ppv_gam = ppv_gam,
      sensitivity_classic = sens_classic,
      specificity_classic = spec_classic,
      ppv_classic = ppv_classic
    )
  }, by = stage]
  
  # Guardar métricas por etapa
  fwrite(metrics_by_stage, paste0(output_dir, "validation_metrics_by_stage_filtered.csv"))
  
  message("Métricas por etapa en pool filtrado:")
  print(metrics_by_stage)
  
  ################################################################################
  # Métricas según número de reportes
  ################################################################################
  
  message("\nCalculando métricas según número de reportes...")
  
  # Categorizar por número de reportes
  pos_expanded_with_stage[, reports_category := cut(
    n_coadmin, 
    breaks = c(0, 10, 20, 30, Inf), 
    labels = c("5-10", "11-20", "21-30", "30+"),
    include.lowest = TRUE
  )]
  
  # Calcular métricas por categoría de reportes
  metrics_by_reports <- pos_expanded_with_stage[, {
    
    # GAM
    tp_gam <- sum(signal_gam_stage == TRUE, na.rm = TRUE)
    fn_gam <- sum(signal_gam_stage == FALSE, na.rm = TRUE)
    sens_gam <- ifelse((tp_gam + fn_gam) > 0, tp_gam / (tp_gam + fn_gam), NA_real_)
    
    # IOR Clásico
    tp_classic <- sum(signal_classic_stage == TRUE, na.rm = TRUE)
    fn_classic <- sum(signal_classic_stage == FALSE, na.rm = TRUE)
    sens_classic <- ifelse((tp_classic + fn_classic) > 0, tp_classic / (tp_classic + fn_classic), NA_real_)
    
    data.table(
      reports_category = reports_category[1],
      n_triplets = .N,
      mean_fold_change = mean(fold_change, na.rm = TRUE),
      sensitivity_gam = sens_gam,
      sensitivity_classic = sens_classic
    )
  }, by = reports_category]
  
  # Guardar métricas por reportes
  fwrite(metrics_by_reports, paste0(output_dir, "validation_metrics_by_reports_filtered.csv"))
  
  message("Métricas por número de reportes en pool filtrado:")
  print(metrics_by_reports)
  
}

################################################################################
# Resumen final y reporte
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("RESUMEN DEL ANÁLISIS DE PODER")
message(paste(rep("=", 80), collapse = ""))

cat(sprintf("
ARCHIVOS GENERADOS:
- Matriz de poder: %spower_matrix.csv
- Umbrales 80%% poder: %spower_thresholds_80.csv
- Heatmaps: %sheatmap_power_*.png
- Pool filtrado: %sfiltered_positive_pool.rds
- Métricas validación: %svalidation_metrics_*.csv

UMBRALES PARA 80%% DE PODER:
- GAM: FC >= %.1f, Reportes >= %d
- IOR Clásico: FC >= %.1f, Reportes >= %d

POOL FILTRADO:
- Tripletes positivos: %d (de %d originales)
- Criterios: FC >= %.1f, Reportes >= %d, dinámico != uniform
",
output_dir, output_dir, output_dir, output_dir, output_dir,
min_fc_gam, min_reports_gam, min_fc_classic, min_reports_classic,
nrow(filtered_positives), nrow(positives_scores),
min_fc_gam, min_reports_gam))

# Crear reporte consolidado
report <- data.table(
  parameter = c(
    "Fold-change mínimo GAM", "Reportes mínimo GAM",
    "Fold-change mínimo IOR", "Reportes mínimo IOR",
    "Tripletes originales", "Tripletes filtrados",
    "Poder objetivo"
  ),
  value = c(
    min_fc_gam, min_reports_gam,
    min_fc_classic, min_reports_classic,
    nrow(positives_scores), nrow(filtered_positives),
    target_power
  )
)

fwrite(report, paste0(output_dir, "power_analysis_summary.csv"))

message("\nAnálisis de poder completado exitosamente!")
message(sprintf("Resultados guardados en: %s", output_dir))