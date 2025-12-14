################################################################################
# Script de análisis de poder estadístico
# Script 04_power_analysis
################################################################################

library(data.table)
library(tidyverse)
library(mgcv)
library(pROC)
library(boot)
library(pbapply)
library(parallel)
library(RColorBrewer)
library(gridExtra)

setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)
source("giangreco_theme.R")
theme_set(theme_giangreco())

################################################################################
# Configuración
################################################################################

# Parámetros de percentil
PERCENTILE_LEVEL <- "p95"  

# Configuración de GAM (debe coincidir con scripts anteriores)
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                 
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML" 

# Configuración del análisis de poder
POWER_TARGET <- 0.80  # Objetivo de poder estadístico

# Rangos para el análisis de poder
fold_change_range <- seq(1.1, 5, by = 0.2)  # Rangos de fold-change
n_reports_range <- seq(10, 200, by = 10)    # Rangos de número de reportes
n_simulations <- 100  # Simulaciones para estimar poder

# Directorios
ruta_ade_raw <- "./ade_raw.csv"
ruta_pos_results <- "./augmentation_results/positive_triplets_results.rds"
ruta_neg_results <- "./augmentation_results/negative_triplets_results.rds"
ruta_null_thresh <- "./null_distribution_results/null_thresholds.csv"

# Directorio de salida
output_dir <- "./power_analysis_results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Carga de datos base
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)
positives_scores <- readRDS(ruta_pos_results)
negatives_scores <- readRDS(ruta_neg_results)
null_thresholds <- fread(ruta_null_thresh)

# Configuración de percentil
percentile_config <- list(
  p90 = list(threshold_col = "threshold_p90", alpha = 0.10),
  p95 = list(threshold_col = "threshold_p95", alpha = 0.05),
  p99 = list(threshold_col = "threshold_p99", alpha = 0.01)
)

config <- percentile_config[[PERCENTILE_LEVEL]]
null_thresholds[, threshold := get(config$threshold_col)]

# Preparar datos
pos_valid <- positives_scores[model_success == TRUE & injection_success == TRUE]
neg_valid <- negatives_scores[model_success == TRUE]

# Filtrar datos uniformes para análisis de poder
pos_valid <- pos_valid[dynamic != "uniform"]

message(sprintf("Datos cargados - Positivos válidos: %d, Negativos válidos: %d", 
                nrow(pos_valid), nrow(neg_valid)))

################################################################################
# Función de cálculo de poder estadístico
################################################################################

calculate_power_for_combination <- function(fold_change, n_reports, null_thresholds, 
                                          ade_data, sim_n = 100) {
  """
  Calcula el poder estadístico para una combinación específica de fold-change y número de reportes
  """
  
  power_results <- list()
  
  # Parámetros para simulación
  dynamic_type <- "increase"  # Usamos una dinámica representativa
  
  # Tomar una muestra representativa de tripletes
  sample_triplets <- pos_valid[sample(nrow(pos_valid), min(sim_n, nrow(pos_valid)))]
  
  # Resultados para GAM y IOR clásico
  gam_signals <- c()
  classic_signals <- c()
  
  for (i in seq_len(nrow(sample_triplets))) {
    triplet <- sample_triplets[i]
    
    tryCatch({
      # Simular inyección de señal
      injection_result <- inject_signal(
        drugA_id = triplet$drugA,
        drugB_id = triplet$drugB,
        event_id = triplet$meddra,
        dynamic_type = dynamic_type,
        fold_change = fold_change,
        ade_raw_dt = ade_data
      )
      
      if (injection_result$success && injection_result$n_injected > 0) {
        
        # GAM
        gam_result <- fit_differential_gam(
          drugA_id = triplet$drugA,
          drugB_id = triplet$drugB,
          event_id = triplet$meddra,
          ade_data = injection_result$ade_aug,
          spline_individuales = spline_individuales,
          include_sex = include_sex,
          include_stage_sex = include_stage_sex,
          k_spline = k_spline,
          bs_type = bs_type,
          select = select,
          nichd_spline = nichd_spline
        )
        
        if (gam_result$success) {
          # Verificar si detecta señal según criterio doble
          signal_detection <- calculate_interaction_signal(gam_result, null_thresholds$threshold)
          gam_signals <- c(gam_signals, as.numeric(signal_detection$signal_detected))
        }
        
        # IOR clásico
        classic_result <- calculate_classic_ior(
          drugA_id = triplet$drugA,
          drugB_id = triplet$drugB,
          event_id = triplet$meddra,
          ade_data = injection_result$ade_aug
        )
        
        if (classic_result$success) {
          # Verificar detección con IOR clásico
          stages_with_signal <- sum(classic_result$results_by_stage$log_ior_classic_lower90 > 0)
          classic_signals <- c(classic_signals, as.numeric(stages_with_signal > 0))
        }
      }
      
    }, error = function(e) {
      # Continuar con siguiente simulación
      return(NULL)
    })
  }
  
  # Calcular poder
  gam_power <- if(length(gam_signals) > 0) mean(gam_signals) else 0
  classic_power <- if(length(classic_signals) > 0) mean(classic_signals) else 0
  
  return(list(
    gam_power = gam_power,
    classic_power = classic_power,
    n_simulations = length(gam_signals),
    fold_change = fold_change,
    n_reports = n_reports
  ))
}

################################################################################
# Análisis de poder por combinaciones
################################################################################

message("\nIniciando análisis de poder estadístico...")
message(sprintf("Probando %d x %d combinaciones = %d puntos", 
                length(fold_change_range), length(n_reports_range),
                length(fold_change_range) * length(n_reports_range)))

# Crear grid de combinaciones
power_grid <- CJ(fold_change = fold_change_range, n_reports = n_reports_range)

# Calcular poder para cada combinación (esto puede tomar tiempo)
pb <- txtProgressBar(max = nrow(power_grid), style = 3)
power_results <- list()

for (i in seq_len(nrow(power_grid))) {
  setTxtProgressBar(pb, i)
  
  row <- power_grid[i]
  
  power_result <- calculate_power_for_combination(
    fold_change = row$fold_change,
    n_reports = row$n_reports,
    null_thresholds = null_thresholds,
    ade_data = ade_raw_dt,
    sim_n = 20  # Reducido para tiempo de ejecución
  )
  
  power_results[[i]] <- power_result
}

close(pb)

# Convertir a data.table
power_dt <- rbindlist(power_results)

# Guardar resultados
fwrite(power_dt, paste0(output_dir, "power_analysis_results.csv"))
saveRDS(power_dt, paste0(output_dir, "power_analysis_results.rds"))

################################################################################
# Identificar umbrales para 80% de poder
################################################################################

message("\nIdentificando umbrales para 80% de poder...")

# GAM thresholds
gam_80_threshold <- power_dt[gam_power >= POWER_TARGET, .(
  min_fold_change = min(fold_change),
  min_n_reports = min(n_reports)
)]

# IOR clásico thresholds
classic_80_threshold <- power_dt[classic_power >= POWER_TARGET, .(
  min_fold_change = min(fold_change),
  min_n_reports = min(n_reports)
)]

message("Umbrales para 80% de poder:")
message("GAM:")
print(gam_80_threshold)
message("IOR Clásico:")
print(classic_80_threshold)

# Guardar umbrales
thresholds_summary <- data.table(
  method = c("GAM", "IOR_Classic"),
  power_target = POWER_TARGET,
  min_fold_change = c(gam_80_threshold$min_fold_change, classic_80_threshold$min_fold_change),
  min_n_reports = c(gam_80_threshold$min_n_reports, classic_80_threshold$min_n_reports)
)

fwrite(thresholds_summary, paste0(output_dir, "power_thresholds_80.csv"))

################################################################################
# Generar heatmaps de poder
################################################################################

message("\nGenerando heatmaps de poder...")

# Preparar datos para heatmap
power_matrix_gam <- dcast(power_dt, fold_change ~ n_reports, value.var = "gam_power")
power_matrix_classic <- dcast(power_dt, fold_change ~ n_reports, value.var = "classic_power")

# Convertir a matrices
gam_matrix <- as.matrix(power_matrix_gam[, -1])
rownames(gam_matrix) <- power_matrix_gam$fold_change
colnames(gam_matrix) <- names(power_matrix_gam)[-1]

classic_matrix <- as.matrix(power_matrix_classic[, -1])
rownames(classic_matrix) <- power_matrix_classic$fold_change
colnames(classic_matrix) <- names(power_matrix_classic)[-1]

# Heatmap para GAM
p_gam_heatmap <- ggplot(melt(gam_matrix), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1),
                     name = "Poder") +
  labs(
    title = "Análisis de Poder - Modelo GAM",
    subtitle = sprintf("Poder estadístico por fold-change y número de reportes (Target: %d%%)", 
                      POWER_TARGET * 100),
    x = "Número de reportes",
    y = "Fold-change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_contour(aes(z = value), breaks = c(0.8), color = "yellow", size = 2)

ggsave(paste0(output_dir, "heatmap_power_gam.png"), p_gam_heatmap, 
       width = 12, height = 8, dpi = 300)

# Heatmap para IOR Clásico
p_classic_heatmap <- ggplot(melt(classic_matrix), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1),
                     name = "Poder") +
  labs(
    title = "Análisis de Poder - IOR Clásico",
    subtitle = sprintf("Poder estadístico por fold-change y número de reportes (Target: %d%%)", 
                      POWER_TARGET * 100),
    x = "Número de reportes",
    y = "Fold-change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_contour(aes(z = value), breaks = c(0.8), color = "yellow", size = 2)

ggsave(paste0(output_dir, "heatmap_power_classic.png"), p_classic_heatmap, 
       width = 12, height = 8, dpi = 300)

################################################################################
# Filtrar pool positivo según características de poder
################################################################################

message("\nFiltrando pool positivo según umbrales de poder...")

# Aplicar filtros basados en umbrales de 80% poder
gam_filter <- pos_valid[
  fold_change >= gam_80_threshold$min_fold_change &
  N >= gam_80_threshold$min_n_reports
]

classic_filter <- pos_valid[
  fold_change >= classic_80_threshold$min_fold_change &
  N >= classic_80_threshold$min_n_reports
]

# Pool combinado (cumple ambos criterios)
combined_filter <- pos_valid[
  fold_change >= max(gam_80_threshold$min_fold_change, classic_80_threshold$min_fold_change) &
  N >= max(gam_80_threshold$min_n_reports, classic_80_threshold$min_n_reports)
]

message(sprintf("Pool filtrado - GAM: %d tripletes", nrow(gam_filter)))
message(sprintf("Pool filtrado - IOR Clásico: %d tripletes", nrow(classic_filter)))
message(sprintf("Pool filtrado - Combinado: %d tripletes", nrow(combined_filter)))

# Guardar pools filtrados
fwrite(gam_filter, paste0(output_dir, "positive_pool_filtered_gam.csv"))
fwrite(classic_filter, paste0(output_dir, "positive_pool_filtered_classic.csv"))
fwrite(combined_filter, paste0(output_dir, "positive_pool_filtered_combined.csv"))

################################################################################
# Análisis de validación comparativa en pool filtrado
################################################################################

message("\nEjecutando análisis de validación en pool filtrado...")

# Usar pool combinado para validación
validation_pool <- combined_filter

if (nrow(validation_pool) > 0) {
  
  # Preparar datos para validación
  pos_validation <- validation_pool
  neg_validation <- negatives_scores[model_success == TRUE]
  
  # Análisis similar al script 032_validation_param.R
  pos_expanded <- pos_validation[, {
    data.table(
      stage = unlist(stage),
      log_ior = unlist(log_ior),
      log_ior_lower90 = unlist(log_ior_lower90),
      ior_value = unlist(ior_values)
    )
  }, by = .(triplet_id, drugA, drugB, meddra, type, dynamic, fold_change)]

  neg_expanded <- neg_validation[, {
    data.table(
      stage = unlist(stage),
      log_ior = unlist(log_ior),
      log_ior_lower90 = unlist(log_ior_lower90),
      ior_value = unlist(ior_values)
    )
  }, by = .(triplet_id, drugA, drugB, meddra, type)]

  pos_expanded[, stage_name := niveles_nichd[stage]]
  neg_expanded[, stage_name := niveles_nichd[stage]]

  # Merge con umbrales
  pos_expanded <- merge(pos_expanded, null_thresholds[, .(stage, threshold)], by = "stage")
  neg_expanded <- merge(neg_expanded, null_thresholds[, .(stage, threshold)], by = "stage")

  # Aplicar criterios de señal
  pos_expanded[, `:=`(
    nominal_sig = log_ior_lower90 > 0,
    nullmodel_sig = log_ior_lower90 > threshold,
    signal_giangreco = (log_ior_lower90 > 0) & (log_ior_lower90 > threshold)
  )]

  neg_expanded[, `:=`(
    nominal_sig = log_ior_lower90 > 0,
    nullmodel_sig = log_ior_lower90 > threshold,
    signal_giangreco = (log_ior_lower90 > 0) & (log_ior_lower90 > threshold)
  )]

  ################################################################################
  # Métricas de performance
  ################################################################################

  message("\nCalculando métricas de performance...")

  # Clasificación a nivel triplete
  signal_by_triplet <- rbind(
    pos_expanded[, .(
      signal_detected = any(signal_giangreco, na.rm = TRUE),
      n_stages_sig = sum(signal_giangreco, na.rm = TRUE)
    ), by = .(triplet_id, type, dynamic, fold_change)],
    neg_expanded[, .(
      signal_detected = any(signal_giangreco, na.rm = TRUE),
      n_stages_sig = sum(signal_giangreco, na.rm = TRUE)
    ), by = .(triplet_id, type)]
  )

  # Dataset combinado
  all_data_validation <- rbind(
    pos_validation[, .(triplet_id, type, max_ior, mean_ior, dynamic, fold_change, true_label = 1)],
    neg_validation[, .(triplet_id, type, max_ior, mean_ior, dynamic = NA_character_, fold_change = NA_real_, true_label = 0)]
  )

  all_data_validation <- merge(all_data_validation, signal_by_triplet[, .(triplet_id, signal_detected)], by = "triplet_id")

  # Calcular métricas
  tp <- sum(all_data_validation$signal_detected == TRUE & all_data_validation$true_label == 1, na.rm = TRUE)
  fn <- sum(all_data_validation$signal_detected == FALSE & all_data_validation$true_label == 1, na.rm = TRUE)
  fp <- sum(all_data_validation$signal_detected == TRUE & all_data_validation$true_label == 0, na.rm = TRUE)
  tn <- sum(all_data_validation$signal_detected == FALSE & all_data_validation$true_label == 0, na.rm = TRUE)

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)

  # AUROC
  roc_obj <- roc(response = all_data_validation$true_label, predictor = all_data_validation$max_ior, quiet = TRUE)
  auroc <- auc(roc_obj)

  # Métricas globales
  metrics_global_validation <- data.table(
    metric = c("AUROC", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1_Score"),
    value = c(auroc, sensitivity, specificity, ppv, npv, accuracy, f1_score)
  )

  fwrite(metrics_global_validation, paste0(output_dir, "metrics_validation_pool.csv"))

  message("\nMétricas de Performance (Pool Filtrado):")
  print(metrics_global_validation)

  ################################################################################
  # Métricas por etapa
  ################################################################################

  # Métricas por etapa
  metrics_by_stage_validation <- pos_expanded[, {
    
    # Calcular para positivos
    pos_signals <- sum(signal_giangreco, na.rm = TRUE)
    pos_total <- .N
    
    # Calcular para negativos (usar datos neg_expanded)
    neg_data <- neg_expanded[stage == .BY$stage]
    neg_signals <- sum(neg_data$signal_giangreco, na.rm = TRUE)
    neg_total <- nrow(neg_data)
    
    # AUROC con IOR
    roc_stage <- tryCatch({
      roc(response = c(rep(1, pos_total), rep(0, neg_total)), 
          predictor = c(ior_value, neg_data$ior_value),
          direction = "<", quiet = TRUE)
    }, error = function(e) NULL)
    
    auroc_val <- if(!is.null(roc_stage)) auc(roc_stage) else NA
    
    data.table(
      stage = .BY$stage,
      stage_name = niveles_nichd[.BY$stage],
      n_pos = pos_total,
      n_neg = neg_total,
      n_signals_pos = pos_signals,
      n_signals_neg = neg_signals,
      auroc = auroc_val,
      detection_rate = pos_signals / pos_total
    )
  }, by = stage]

  fwrite(metrics_by_stage_validation, paste0(output_dir, "metrics_by_stage_validation.csv"))

  ################################################################################
  # Métricas por número de reportes
  ################################################################################

  # Análisis por número de reportes
  pos_validation[, reports_category := cut(N, breaks = c(0, 50, 100, 150, Inf), 
                                          labels = c("10-50", "51-100", "101-150", "150+"))]

  metrics_by_reports <- pos_validation[, {
    
    pos_expanded_cat <- pos_expanded[triplet_id %in% triplet_id]
    
    tp_cat <- sum(pos_expanded_cat$signal_giangreco, na.rm = TRUE)
    total_pos <- nrow(pos_expanded_cat)
    
    # Tomar muestra de negativos para comparación
    neg_sample <- negatives_scores[sample(nrow(negatives_scores), min(total_pos, nrow(negatives_scores)))]
    neg_expanded_cat <- neg_sample[, {
      data.table(
        stage = unlist(stage),
        log_ior_lower90 = unlist(log_ior_lower90),
        signal_giangreco = unlist(log_ior_lower90) > 0  # Criterio simple para negativos
      )
    }, by = triplet_id]
    
    fp_cat <- sum(neg_expanded_cat$signal_giangreco, na.rm = TRUE)
    total_neg <- nrow(neg_expanded_cat)
    
    sensitivity_cat <- tp_cat / total_pos
    specificity_cat <- 1 - (fp_cat / total_neg)
    
    data.table(
      reports_category = .BY$reports_category,
      n_triplets = .N,
      sensitivity = sensitivity_cat,
      specificity = specificity_cat,
      detection_rate = tp_cat / total_pos
    )
  }, by = reports_category]

  fwrite(metrics_by_reports, paste0(output_dir, "metrics_by_reports_validation.csv"))

  ################################################################################
  # Gráficos comparativos
  ################################################################################

  # Curva ROC para pool filtrado
  p_roc_filtered <- ggroc(roc_obj, legacy.axes = TRUE, colour = "#E31A1C", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    annotate("text", x = 0.75, y = 0.25,
             label = sprintf("AUC = %.3f\nPool Filtrado", auroc),
             size = 5, fontface = "bold", color = "#E31A1C") +
    labs(
      title = "Curva ROC - Pool Filtrado según Análisis de Poder",
      subtitle = sprintf("Umbral 80%% poder: FC≥%.2f, Reportes≥%d", 
                       max(gam_80_threshold$min_fold_change, classic_80_threshold$min_fold_change),
                       max(gam_80_threshold$min_n_reports, classic_80_threshold$min_n_reports)),
      x = "1 - Especificidad (FPR)", 
      y = "Sensibilidad (TPR)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30")
    )

  ggsave(paste0(output_dir, "fig_roc_filtered_pool.png"), p_roc_filtered, 
         width = 8, height = 7, dpi = 300)

}

################################################################################
# Resumen final
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("RESUMEN DEL ANÁLISIS DE PODER")
message(paste(rep("=", 80), collapse = ""))

cat(sprintf("
ANÁLISIS DE PODER COMPLETADO
================================

Configuración:
- Poder objetivo: %.0f%%
- Rango fold-change: %.1f - %.1f
- Rango reportes: %d - %d
- Simulaciones por punto: %d

Umbrales identificados para %d%% poder:
- GAM: FC ≥ %.2f, Reportes ≥ %d
- IOR Clásico: FC ≥ %.2f, Reportes ≥ %d

Pools filtrados:
- GAM: %d tripletes
- IOR Clásico: %d tripletes  
- Combinado: %d tripletes

Archivos generados:
- power_analysis_results.csv: Resultados completos de poder
- power_thresholds_80.csv: Umbrales para 80%% poder
- heatmap_power_gam.png: Heatmap poder GAM
- heatmap_power_classic.png: Heatmap poder IOR clásico
- metrics_validation_pool.csv: Métricas pool filtrado
- metrics_by_stage_validation.csv: Métricas por etapa
- metrics_by_reports_validation.csv: Métricas por número de reportes
- fig_roc_filtered_pool.png: ROC pool filtrado

Directorio de resultados: %s
", 
POWER_TARGET * 100, min(fold_change_range), max(fold_change_range),
min(n_reports_range), max(n_reports_range), 20,
POWER_TARGET * 100,
gam_80_threshold$min_fold_change, gam_80_threshold$min_n_reports,
classic_80_threshold$min_fold_change, classic_80_threshold$min_n_reports,
nrow(gam_filter), nrow(classic_filter), nrow(combined_filter),
output_dir))

message("\nAnálisis de poder completado exitosamente!")