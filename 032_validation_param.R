################################################################################
# Script de validación
# Script 03_validation
################################################################################

library(data.table)
library(tidyverse)
library(mgcv)
library(pROC)
library(boot)
library(pbapply)
library(parallel)


setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)
source("giangreco_theme.R")
theme_set(theme_giangreco())

################################################################################
# Configuración de parámetros
################################################################################

# umbrales de significancia
# Opciones: "p90", "p95", "p99", "p999"
PERCENTILE_LEVEL <- "p95"  

# Parámetros de fórmula para GAM
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                 
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML" 


# Mapeo de configuraciones por percentil
percentile_config <- list(
  p90 = list(
    threshold_col = "threshold_p90",
    alpha = 0.10,
    label = "Percentilo 90 (α=0.10)",
    z_score = qnorm(0.95)
  ),
  p95 = list(
    threshold_col = "threshold_p95",
    alpha = 0.05,
    label = "Percentilo 95 (α=0.05)",
    z_score = qnorm(0.975)
  ),
  p99 = list(
    threshold_col = "threshold_p99",
    alpha = 0.01,
    label = "Percentilo 99 (α=0.01)",
    z_score = qnorm(0.995)
  )
)


config <- percentile_config[[PERCENTILE_LEVEL]]

message("\n", paste(rep("=", 80), collapse = ""))
message("Percentil elegido para validación: ", config$label)
message("  Umbral: ", config$threshold_col)
message("  Alfa: ", config$alpha)
message(paste(rep("=", 80), collapse = ""))

# Archivos utilizados
ruta_ade_raw <- "./ade_raw.csv"
ruta_pos_results <- "./augmentation_results/positive_triplets_results.rds"
ruta_neg_results <- "./augmentation_results/negative_triplets_results.rds"
ruta_pos_meta <- "./augmentation_results/positive_triplets_metadata.csv"
ruta_neg_meta <- "./augmentation_results/negative_triplets_metadata.csv"

# Archivos de distribución nula
ruta_null_dist <- "./null_distribution_results/null_distribution.csv"
ruta_null_thresh <- "./null_distribution_results/null_thresholds.csv"

# Parámetros
n_bootstrap <- 1000
n_cores <- max(1, detectCores() - 1)
seed_base <- 2025

# se agrega sufijo a la carpeta de guardado para tener clasificado según los valores de los parámetros
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

# Directorio de salida con identificador de percentil
output_dir <- paste0("./results/validation_results_", PERCENTILE_LEVEL, "_", suffix, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Carga y preparado de datos
################################################################################

positives_scores <- readRDS(ruta_pos_results)
negatives_scores <- readRDS(ruta_neg_results)
pos_meta <- fread(ruta_pos_meta)
neg_meta <- fread(ruta_neg_meta)

null_distribution <- fread(ruta_null_dist)
null_thresholds <- fread(ruta_null_thresh)

# Crea columna genérica 'threshold' para código reutilizable
null_thresholds[, threshold := get(config$threshold_col)]

message("Distribución nula")
cat(sprintf("\nUmbrales %s por etapa:\n", config$label))
print(null_thresholds[, .(stage, stage_name, threshold, mean_null, sd_null, n_samples)])

# filtrado de exitosos
pos_valid <- positives_scores[model_success == TRUE & injection_success == TRUE]
neg_valid <- negatives_scores[model_success == TRUE]

# saco a uniform del análisis, porque solo sirve para comparar scores promedios de las otras dinámicas
# porque uniform por definición no inyecta reportes -->  return(rep(0, N)
pos_valid <- pos_valid[dynamic != "uniform"]

message(sprintf("\n Datos utilizados:"))
message(sprintf("  Positivos válidos: %d", nrow(pos_valid)))
message(sprintf("  Negativos válidos: %d", nrow(neg_valid)))


################################################################################
# Expansión y aplicación de criterio de señal positiva
################################################################################

# Expando por etapa
pos_expanded <- pos_valid[, {
  data.table(
    stage = unlist(stage),
    log_ior = unlist(log_ior),
    log_ior_lower90 = unlist(log_ior_lower90),
    ior_value = unlist(ior_values)
  )
}, by = .(triplet_id, drugA, drugB, meddra, type, dynamic, fold_change,
          max_ior, mean_ior, n_stages_significant)]

neg_expanded <- neg_valid[, {
  data.table(
    stage = unlist(stage),
    log_ior = unlist(log_ior),
    log_ior_lower90 = unlist(log_ior_lower90),
    ior_value = unlist(ior_values)
  )
}, by = .(triplet_id, drugA, drugB, meddra, type,
          max_ior, mean_ior, n_stages_significant)]

pos_expanded[, stage_name := niveles_nichd[stage]]
neg_expanded[, stage_name := niveles_nichd[stage]]

# merge con umbrales del percentil seleccionado
pos_expanded <- merge(
  pos_expanded,
  null_thresholds[, .(stage, threshold)],
  by = "stage", all.x = TRUE
)

neg_expanded <- merge(
  neg_expanded,
  null_thresholds[, .(stage, threshold)],
  by = "stage", all.x = TRUE
)

# aplicación de criterio de "señal positiva" (doble umbral)
# 1. Nominal: IC90% > 0
# 2. Null Model: IC90% > threshold (percentil seleccionado)
pos_expanded[, `:=`(
  nominal_sig = log_ior_lower90 > 0,
  nullmodel_sig = log_ior_lower90 > threshold,
  signal_giangreco = (log_ior_lower90 > 0) & (log_ior_lower90 > threshold)
)]

neg_expanded[, `:=`(
  nominal_sig = log_ior_lower90 > 0,
  nullmodel_sig = log_ior_lower90 > threshold,
  signal_giangreco = (log_ior_lower90 > 0) & (log_ior_lower90 > threshold),
  dynamic = NA_character_,
  fold_change = NA_real_
)]

message(sprintf("  Positivos con señal detectada (≥1 etapa): %d (%.1f%%)",
                sum(pos_expanded[, any(signal_giangreco), by = triplet_id]$V1),
                100 * mean(pos_expanded[, any(signal_giangreco), by = triplet_id]$V1)))
message(sprintf("  Negativos con señal detectada (≥1 etapa): %d (%.1f%%)",
                sum(neg_expanded[, any(signal_giangreco), by = triplet_id]$V1),
                100 * mean(neg_expanded[, any(signal_giangreco), by = triplet_id]$V1)))

# clasificación a nivel triplete
signal_by_triplet <- rbind(
  pos_expanded[, .(
    signal_detected = any(signal_giangreco, na.rm = TRUE),
    n_stages_sig = sum(signal_giangreco, na.rm = TRUE)
  ), by = .(triplet_id, type, dynamic, fold_change)],
  neg_expanded[, .(
    signal_detected = any(signal_giangreco, na.rm = TRUE),
    n_stages_sig = sum(signal_giangreco, na.rm = TRUE)
  ), by = .(triplet_id, type, dynamic, fold_change)]
)

# dataset combinado
all_data <- rbind(
  pos_valid[, .(triplet_id, type, max_ior, mean_ior, dynamic, fold_change,
                n_stages_significant, true_label = 1)],
  neg_valid[, .(triplet_id, type, max_ior, mean_ior, dynamic = NA_character_,
                fold_change = NA_real_, n_stages_significant, true_label = 0)]
)

all_data <- merge(
  all_data,
  signal_by_triplet[, .(triplet_id, signal_detected, n_stages_sig)],
  by = "triplet_id", all.x = TRUE
)

# remover infinitos y NAs (para cuando alguna formula no converge)
all_data_clean <- all_data[is.finite(max_ior) & !is.na(max_ior)]

if (nrow(all_data_clean) < nrow(all_data)) {
  message(sprintf(" Removidos %d tripletes con valores infinitos/NA",
                  nrow(all_data) - nrow(all_data_clean)))
  
  # guardo casos problemáticos para diagnóstico posterior
  problematic_cases <- all_data[!is.finite(max_ior) | is.na(max_ior)]
  fwrite(problematic_cases, paste0(output_dir, "problematic_cases.csv"))
}

all_data <- all_data_clean

all_expanded <- rbind(
  pos_expanded[, .(triplet_id, type, stage, stage_name, log_ior,
                   log_ior_lower90, ior_value, dynamic, true_label = 1,
                   threshold, nominal_sig, nullmodel_sig, signal_giangreco)],
  neg_expanded[, .(triplet_id, type, stage, stage_name, log_ior,
                   log_ior_lower90, ior_value, dynamic = NA_character_,
                   true_label = 0, threshold, nominal_sig, nullmodel_sig,
                   signal_giangreco)]
)
all_expanded[, stage_name := factor(stage_name, levels = niveles_nichd)]


################################################################################
# Métricas de rendimiento
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("MÉTRICAS DE RENDIMIENTO - GLOBAL")
message(paste(rep("=", 80), collapse = ""))

# Matriz de confusión
tp <- sum(all_data$signal_detected == TRUE & all_data$true_label == 1, na.rm = TRUE)
fn <- sum(all_data$signal_detected == FALSE & all_data$true_label == 1, na.rm = TRUE)
fp <- sum(all_data$signal_detected == TRUE & all_data$true_label == 0, na.rm = TRUE)
tn <- sum(all_data$signal_detected == FALSE & all_data$true_label == 0, na.rm = TRUE)

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)
accuracy <- (tp + tn) / (tp + tn + fp + fn)
f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)

# AUROC usando max_ior como score continuo
roc_obj <- roc(
  response = all_data$true_label,
  predictor = all_data$max_ior,
  direction = "<",
  levels = c(0, 1),
  quiet = TRUE
)

auroc_global <- auc(roc_obj)
auroc_ci <- ci.auc(roc_obj, conf.level = 1 - config$alpha)

# resultados
cat(sprintf("
================================================================================
MATRIZ DE CONFUSIÓN (%s)
================================================================================
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d

MÉTRICAS GLOBALES:
  AUROC:                    %.3f  [IC: %.3f - %.3f]
  Sensitivity (Recall):     %.3f  [%d/%d positivos detectados]
  Specificity:              %.3f  [%d/%d negativos correctos]
  PPV (Precision):          %.3f  [%d/%d predicciones+ correctas]
  NPV:                      %.3f
  Accuracy:                 %.3f
  F1-Score:                 %.3f
================================================================================
",
config$label,
tp, fn, fp, tn,
auroc_global, auroc_ci[1], auroc_ci[3],
sensitivity, tp, tp+fn,
specificity, tn, tn+fp,
ppv, tp, tp+fp,
npv, accuracy, f1_score
))

# Guardar métricas
metrics_global <- data.table(
  metric = c("AUROC", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1_Score"),
  value = c(auroc_global, sensitivity, specificity, ppv, npv, accuracy, f1_score),
  ci_lower = c(auroc_ci[1], NA, NA, NA, NA, NA, NA),
  ci_upper = c(auroc_ci[3], NA, NA, NA, NA, NA, NA)
)

fwrite(metrics_global, paste0(output_dir, "metrics_global.csv"))

################################################################################
# Métricas por etapa
################################################################################

message("MÉTRICAS POR ETAPA")

metrics_by_stage <- all_expanded[, {
  
  if (.N < 10 || length(unique(true_label)) < 2) {
    data.table(
      auroc = NA_real_, sensitivity = NA_real_,
      ppv = NA_real_, npv = NA_real_,
      auroc_ci_lower = NA_real_, auroc_ci_upper = NA_real_,
      n_signals = 0, n_total = .N
    )
  } else {
    
    # AUROC con IOR continuo
    roc_stage <- tryCatch({
      roc(response = true_label, predictor = ior_value, 
          direction = "<", levels = c(0, 1), quiet = TRUE)
    }, error = function(e) NULL)
    
    if (is.null(roc_stage)) {
      data.table(
        auroc = NA_real_, sensitivity = NA_real_,
        ppv = NA_real_, npv = NA_real_,
        auroc_ci_lower = NA_real_, auroc_ci_upper = NA_real_,
        n_signals = sum(signal_giangreco, na.rm = TRUE), n_total = .N
      )
    } else {
      
      auroc_val <- auc(roc_stage)
      auroc_ci_val <- ci.auc(roc_stage, conf.level = 1 - config$alpha)
      
      # Métricas usando doble criterio
      tp_s <- sum(signal_giangreco == TRUE & true_label == 1, na.rm = TRUE)
      fn_s <- sum(signal_giangreco == FALSE & true_label == 1, na.rm = TRUE)
      fp_s <- sum(signal_giangreco == TRUE & true_label == 0, na.rm = TRUE)
      tn_s <- sum(signal_giangreco == FALSE & true_label == 0, na.rm = TRUE)
      
      sens_s <- ifelse((tp_s + fn_s) > 0, tp_s / (tp_s + fn_s), NA_real_)
      ppv_s <- ifelse((tp_s + fp_s) > 0, tp_s / (tp_s + fp_s), NA_real_)
      npv_s <- ifelse((tn_s + fn_s) > 0, tn_s / (tn_s + fn_s), NA_real_)
      
      data.table(
        auroc = auroc_val,
        sensitivity = sens_s,
        ppv = ppv_s,
        npv = npv_s,
        auroc_ci_lower = auroc_ci_val[1],
        auroc_ci_upper = auroc_ci_val[3],
        n_signals = sum(signal_giangreco, na.rm = TRUE),
        n_total = .N
      )
    }
  }
}, by = .(stage, stage_name)]

cat("\nMÉTRICAS POR ETAPA:\n")
print(metrics_by_stage[!is.na(auroc), 
      .(stage_name, auroc, sensitivity, ppv, npv, n_signals)])

fwrite(metrics_by_stage, paste0(output_dir, "metrics_by_stage.csv"))

################################################################################
# Gráficos
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("GENERANDO FIGURAS PRINCIPALES")
message(paste(rep("=", 80), collapse = ""))

# Gráfico de distribuciones (observada vs null distribution vs negativo)
null_sample <- null_distribution[sample(.N, min(.N, 10000))]

comparison_data <- rbind(
  null_sample[, .(stage, log_value = log_lower90, source = "Null Distribution")],
  pos_expanded[, .(stage, log_value = log_ior_lower90, source = "Positives")],
  neg_expanded[, .(stage, log_value = log_ior_lower90, source = "Negatives")]
)

comparison_data[, stage_name := niveles_nichd[stage]]

p_null_vs_obs <- ggplot(comparison_data, aes(x = log_value, fill = source)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = null_thresholds,
             aes(xintercept = threshold),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ stage_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(
    "Distribución nula" = "gray60",
    "Positivos" = "#4DAF4A",
    "Negativos" = "#E41A1C"
  )) +
  labs(
    title = sprintf("Distribución nula vs señales detectadas (%s)", config$label),
    subtitle = sprintf("Linea roja = %s umbral de detección", config$label),
    x = "Log(IOR) - Rango inferior IC 90%",
    y = "Densidad",
    fill = "Source"
  ) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-10, 15))

ggsave(paste0(output_dir, "fig_CRITICAL_null_vs_observed.png"), p_null_vs_obs,
       width = 16, height = 12, dpi = 300)

# Curva ROC
message("Generando curva ROC...")

p_roc <- ggroc(roc_obj, legacy.axes = TRUE, colour = "#2C3E50", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  annotate("text", x = 0.75, y = 0.25,
           label = sprintf("AUC = %.3f\n%s", auroc_global, config$label),
           size = 5, fontface = "bold", color = "#2C3E50") +
  labs(
    title = sprintf("Curva ROC - Detección de señal (%s)", config$label),
    subtitle = sprintf("Positivos: %d  |  Negativos: %d", 
                       nrow(pos_valid), nrow(neg_valid)),
    x = "1 - Especificidad (FPR)", 
    y = "Sensibilidad (TPR)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30")
  )

ggsave(paste0(output_dir, "fig_roc_curve.png"), p_roc, 
       width = 8, height = 7, dpi = 300)

# Distribución de log-IOR por etapa

p_log_ior <- ggplot(all_expanded[type == "positive"],
                    aes(x = stage_name, y = log_ior)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_violin(fill = "#4DAF4A", alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
  labs(
    title = "Distribución de Log-IOR en etapas de desarrollo (Controles positivos)",
    subtitle = sprintf("Usando %s umbral", config$label),
    x = "Etapa del desarrollo", y = "Log(IOR)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir, "fig_log_ior_distribution.png"), p_log_ior,
       width = 12, height = 7, dpi = 300)

# Performance por etapa

metrics_by_stage_long <- melt(metrics_by_stage[!is.na(auroc)],
                              id.vars = c("stage", "stage_name"),
                              measure.vars = c("auroc", "sensitivity", "ppv", "npv"),
                              variable.name = "metric", value.name = "value")

p_perf_stage <- ggplot(metrics_by_stage_long,
                       aes(x = stage_name, y = value, color = metric, group = metric)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1", labels = c("AUROC", "Sensitivity", "PPV", "NPV")) +
  labs(
    title = sprintf("Performance de detección por etapa (%s)", config$label),
    x = "Etapa del desarrollo", y = "Valor de performance",
    color = "Métrica"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir, "fig_performance_by_stage.png"), p_perf_stage,
       width = 12, height = 7, dpi = 300)


################################################################################
# Resumen
################################################################################

# comparo criterio simple vs doble
n_nominal_only <- sum(all_expanded[, any(nominal_sig & !nullmodel_sig), by = triplet_id]$V1)
n_giangreco <- sum(all_data$signal_detected, na.rm = TRUE)

# Extraer configuración del modelo
model_config_summary <- sprintf(
  "  Spline type:              %s
  Select (penalty to zero): %s
  Stage effect:             %s
  Individual effects:       %s
  Include sex:              %s
  Stage-sex interaction:    %s
  Number of knots (k):      %d",
  ifelse(select, "TRUE", "FALSE"),
  ifelse(nichd_spline, "spline", "linear"),
  ifelse(spline_individuales, "splines", "linear"),
  ifelse(include_sex, "YES", "NO"),
  ifelse(include_stage_sex, "YES", "NO"),
  k_spline
)

executive_summary <- paste0("
================================================================================
RESUMEN DE VALIDACION (", config$label, ")
================================================================================

CRITERIO DE DETECCIÓN DE SEÑAL:
  Percentil seleccionado:   ", PERCENTILE_LEVEL, "
  1. Nominal:               IC90% log(IOR) > 0
  2. Null Model:            IC90% log(IOR) > ", config$threshold_col, "

CONFIGURACION DEL MODELO GAM:
", model_config_summary, "

DISTRIBUCION NULA:
  Total observaciones:      ", nrow(null_distribution), "
  Umbrales ", PERCENTILE_LEVEL, ":              [", 
  sprintf("%.3f", min(null_thresholds$threshold, na.rm = TRUE)), " - ",
  sprintf("%.3f", max(null_thresholds$threshold, na.rm = TRUE)), "] (rango entre etapas)

METRICAS DE RENDIMIENTO GLOBAL:
  AUROC:                    ", sprintf("%.3f", auroc_global), 
  "  [IC: ", sprintf("%.3f", auroc_ci[1]), " - ", sprintf("%.3f", auroc_ci[3]), "]
  Sensitivity (Power):      ", sprintf("%.3f", sensitivity),
  "  [", tp, "/", tp+fn, " positivos detectados]
  Specificity:              ", sprintf("%.3f", specificity),
  "  [", tn, "/", tn+fp, " negativos correctos]
  PPV (Precision):          ", sprintf("%.3f", ppv),
  "  [", tp, "/", tp+fp, " predicciones correctas]
  NPV:                      ", sprintf("%.3f", npv), "
  Accuracy:                 ", sprintf("%.3f", accuracy), "
  F1-Score:                 ", sprintf("%.3f", f1_score), "

COMPARACION CRITERIOS:
  Solo Nominal (IC90% > 0):    ", n_nominal_only + n_giangreco, " senales
  Con Criterio doble:  ", n_giangreco, " senales
  Reduccion:                   ", 
  sprintf("%.1f", 100 * (1 - n_giangreco / (n_nominal_only + n_giangreco))), "%

METRICAS POR ETAPA (Promedio de etapas validas):
  AUROC:        ", sprintf("%.3f", mean(metrics_by_stage$auroc, na.rm = TRUE)),
  "  (rango: ", sprintf("%.3f", min(metrics_by_stage$auroc, na.rm = TRUE)),
  " - ", sprintf("%.3f", max(metrics_by_stage$auroc, na.rm = TRUE)), ")
  Sensitivity:  ", sprintf("%.3f", mean(metrics_by_stage$sensitivity, na.rm = TRUE)),
  "  (rango: ", sprintf("%.3f", min(metrics_by_stage$sensitivity, na.rm = TRUE)),
  " - ", sprintf("%.3f", max(metrics_by_stage$sensitivity, na.rm = TRUE)), ")
  PPV:          ", sprintf("%.3f", mean(metrics_by_stage$ppv, na.rm = TRUE)),
  "  (rango: ", sprintf("%.3f", min(metrics_by_stage$ppv, na.rm = TRUE)),
  " - ", sprintf("%.3f", max(metrics_by_stage$ppv, na.rm = TRUE)), ")

Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
================================================================================
")

cat(executive_summary)
writeLines(executive_summary, paste0(output_dir, "EXECUTIVE_SUMMARY.txt"))

# Guardado
fwrite(signal_by_triplet, paste0(output_dir, "signal_classification.csv"))
fwrite(all_expanded, paste0(output_dir, "all_expanded_with_criteria.csv"))





