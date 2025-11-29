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

# Directorio de salida con identificador de percentil
output_dir <- paste0("./validation_results_", PERCENTILE_LEVEL, "/")
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

# Mostrar resultados
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

# ============================================================================
# SECCIÓN 4: MÉTRICAS POR ETAPA
# ============================================================================

message("\n", paste(rep("=", 80), collapse = ""))
message("MÉTRICAS POR ETAPA")
message(paste(rep("=", 80), collapse = ""))

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
      
      # Métricas usando criterio Giangreco
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

# ============================================================================
# SECCIÓN 5: FIGURAS PRINCIPALES
# ============================================================================

message("\n", paste(rep("=", 80), collapse = ""))
message("GENERANDO FIGURAS PRINCIPALES")
message(paste(rep("=", 80), collapse = ""))

# --- FIGURA CRÍTICA: NULL vs OBSERVADO ---
message("Generando figura crítica: Null vs Observado")

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

# --- FIGURA: CURVA ROC ---
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

# --- FIGURA: DISTRIBUCIÓN LOG-IOR POR ETAPA ---
message("Generando distribución log-IOR")

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

# --- FIGURA: PERFORMANCE POR ETAPA ---
message("Generando performance por etapa")

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


# --- FIGURA: PERFORMANCE POR DINÁMICA (CORREGIDO 2) ---
if ("dynamic" %in% colnames(all_data) && !all(is.na(all_data$dynamic))) {
  
  message("Generando performance por dinámica...")
  
  # 1. Aislar los controles negativos (Background)
  neg_controls <- all_data[true_label == 0]
  
  # 2. Iterar sobre cada dinámica positiva
  metrics_by_dynamic <- all_data[type == "positive" & !is.na(dynamic), {
    
    # 'curr_pos' es .SD (los datos de este grupo, SIN la columna dynamic)
    curr_pos <- .SD
    
    # CORRECCIÓN: Usar fill = TRUE porque neg_controls tiene la columna 'dynamic' 
    # y curr_pos no la tiene (por estar agrupado por ella).
    combined_dt <- rbind(curr_pos, neg_controls, fill = TRUE)
    
    if (nrow(curr_pos) < 5 || nrow(neg_controls) < 5) {
      list(auroc = NA_real_, sensitivity = NA_real_, ppv = NA_real_, npv = NA_real_)
    } else {
      # Calcular ROC con el dataset combinado
      roc_dyn <- tryCatch({
        roc(response = combined_dt$true_label, predictor = combined_dt$max_ior,
            direction = "<", levels = c(0, 1), quiet = TRUE)
      }, error = function(e) NULL)
      
      if (is.null(roc_dyn)) {
        list(auroc = NA_real_, sensitivity = NA_real_, ppv = NA_real_, npv = NA_real_)
      } else {
        auroc_val <- as.numeric(auc(roc_dyn))
        
        # 'Best' threshold para métricas
        coords_dyn <- coords(roc_dyn, "best", 
                             ret = c("sensitivity", "specificity", "ppv", "npv"),
                             transpose = FALSE)
        
        list(
          auroc = auroc_val,
          sensitivity = coords_dyn$sensitivity[1],
          ppv = coords_dyn$ppv[1],
          npv = coords_dyn$npv[1]
        )
      }
    }
  }, by = dynamic] 
  
  # Limpiar NAs antes de graficar
  metrics_by_dynamic <- metrics_by_dynamic[!is.na(auroc)]
  
  if (nrow(metrics_by_dynamic) > 0) {
    metrics_by_dynamic_long <- melt(metrics_by_dynamic, id.vars = "dynamic",
                                    variable.name = "metric", value.name = "value")
    
    p_perf_dynamic <- ggplot(metrics_by_dynamic_long,
                             aes(x = dynamic, y = value, fill = metric)) +
      geom_col(position = "dodge", alpha = 0.8, color = "black", width = 0.7) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
      scale_fill_brewer(palette = "Set2", 
                        labels = c("AUROC", "Sensibilidad", "VPP", "VPN")) +
      labs(
        title = sprintf("Performance de detección según dinámica(%s)", config$label),
        x = "Dinámica", y = "Valor de Performance",
        fill = "Métrica"
      ) +
      theme_giangreco() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p_perf_dynamic)
    
    ggsave(paste0(output_dir, "fig_performance_by_dynamic.png"), p_perf_dynamic,
           width = 10, height = 6, dpi = 300)
    
    fwrite(metrics_by_dynamic, paste0(output_dir, "metrics_by_dynamic.csv"))
  } else {
    message("⚠ No se pudieron calcular métricas por dinámica (datos insuficientes).")
  }
}

# --- FIGURA: BOOTSTRAP GLOBAL (OPCIONAL) ---
message("Generando bootstrap global")

set.seed(seed_base)

bootstrap_global_metrics <- function(data, indices) {
  d <- data[indices, ]
  
  roc_boot <- roc(response = d$true_label, predictor = d$max_ior,
                  direction = "<", levels = c(0, 1), quiet = TRUE)
  
  auroc_val <- auc(roc_boot)
  coords_boot <- coords(roc_boot, "best",
                        ret = c("sensitivity", "ppv", "npv"))
  
  c(auroc = auroc_val, sensitivity = coords_boot$sensitivity,
    ppv = coords_boot$ppv, npv = coords_boot$npv)
}

boot_global <- boot(data = all_data, statistic = bootstrap_global_metrics,
                    R = min(500, n_bootstrap), parallel = "multicore", ncpus = n_cores)

boot_ci_global <- data.table(
  metric = c("AUROC", "Sensitivity", "PPV", "NPV"),
  mean = colMeans(boot_global$t),
  ci_lower = apply(boot_global$t, 2, quantile, probs = config$alpha/2),
  ci_upper = apply(boot_global$t, 2, quantile, probs = 1 - config$alpha/2)
)

p_bootstrap <- ggplot(boot_ci_global, aes(x = metric, y = mean)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  labs(
    title = sprintf("Bootstrap Performance (Global, IC %.0f%%)", 100*(1-config$alpha)),
    subtitle = sprintf("Based on %d resamples", nrow(boot_global$t)),
    x = "Metric", y = "Value"
  )

ggsave(paste0(output_dir, "fig_bootstrap_global.png"), p_bootstrap,
       width = 8, height = 6, dpi = 300)

fwrite(boot_ci_global, paste0(output_dir, "bootstrap_ci_global.csv"))

# --- FIGURA: DISTRIBUCIÓN DE FOLD-CHANGES ---
if ("fold_change" %in% colnames(pos_meta) && !all(is.na(pos_meta$fold_change))) {
  
  message("Generando distribución de fold-changes")
  
  p_fc <- ggplot(pos_meta, aes(x = fold_change)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "steelblue", alpha = 0.7) +
    geom_density(color = "darkblue", linewidth = 1.2) +
    labs(
      title = "Distribution of Fold-Changes (Positive Controls)",
      x = "Fold Change", y = "Density"
    )
  
  ggsave(paste0(output_dir, "fig_fold_change_distribution.png"), p_fc,
         width = 8, height = 6, dpi = 300)
}

# ============================================================================
# SECCIÓN 6: RESUMEN EJECUTIVO
# ============================================================================

message("\n", paste(rep("=", 80), collapse = ""))
message("GENERANDO RESUMEN")
message(paste(rep("=", 80), collapse = ""))

# ComparaciÃ³n de criterios
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
  bs_type,
  ifelse(select, "TRUE", "FALSE"),
  ifelse(nichd_spline, "spline", "linear"),
  ifelse(spline_individuales, "splines", "linear"),
  ifelse(include_sex, "YES", "NO"),
  ifelse(include_stage_sex, "YES", "NO"),
  k_spline
)

executive_summary <- paste0("
================================================================================
RESUMEN EJECUTIVO - VALIDACION dGAM (", config$label, ")
================================================================================

CRITERIO DE SIGNIFICANCIA APLICADO (Giangreco et al. 2022):
  Percentil seleccionado:   ", PERCENTILE_LEVEL, "
  1. Nominal:               IC90% log(IOR) > 0
  2. Null Model:            IC90% log(IOR) > ", config$threshold_col, "
  -> Senal detectada si AMBOS criterios se cumplen en >=1 etapa

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
  Con Null Model (Giangreco):  ", n_giangreco, " senales
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

ARCHIVOS GENERADOS:
  - metrics_global.csv
  - metrics_by_stage.csv
  - signal_classification.csv
  - fig_CRITICAL_null_vs_observed.png
  - fig_roc_curve.png
  - fig_log_ior_distribution.png
  - fig_performance_by_stage.png
  ", ifelse(exists("p_perf_dynamic"), "- fig_performance_by_dynamic.csv\n", ""), "

Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
================================================================================
")

cat(executive_summary)
writeLines(executive_summary, paste0(output_dir, "EXECUTIVE_SUMMARY.txt"))

# Guardar clasificaciones finales
fwrite(signal_by_triplet, paste0(output_dir, "signal_classification.csv"))
fwrite(all_expanded, paste0(output_dir, "all_expanded_with_criteria.csv"))

message("\n", paste(rep("=", 80), collapse = ""))
message(" VALIDACIÓN COMPLETADA")
message(paste(rep("=", 80), collapse = ""))
message(sprintf("  Percentil usado:          %s", PERCENTILE_LEVEL))
message(sprintf("  Criterio Giangreco:       COMPLETO (doble umbral)"))
message(sprintf("  Distribución nula:        INTEGRADA"))
message(sprintf("  Resultados guardados en:  %s", output_dir))
message(paste(rep("=", 80), collapse = ""))



# ============================================================================
# SECCIÓN 4: MÉTRICAS POR ETAPA (AJUSTADO POR SEÑAL ACTIVA)
# ============================================================================

message("\n", paste(rep("=", 80), collapse = ""))
message("CALCULANDO MÉTRICAS AJUSTADAS (Solo Señales Activas)")
message(paste(rep("=", 80), collapse = ""))

# 1. Definir dónde existe realmente la señal según la dinámica
#    (Basado en las funciones tanh de 00_functions.R)
#    Increase: Sube desde etapa 4.
#    Decrease: Baja hasta etapa 4.
#    Plateau: Centro (3, 4, 5).
#    Inverse: Bordes (1, 2, 6, 7).

all_expanded[, is_active_signal := FALSE] 

# Negativos: Nunca tienen señal activa (se quedan en FALSE)
# Positivos: Depende de la dinámica
all_expanded[type == "positive" & dynamic == "uniform", is_active_signal := TRUE]
all_expanded[type == "positive" & dynamic == "increase" & stage %in% 4:7, is_active_signal := TRUE]
all_expanded[type == "positive" & dynamic == "decrease" & stage %in% 1:4, is_active_signal := TRUE]
all_expanded[type == "positive" & dynamic == "plateau" & stage %in% 3:5, is_active_signal := TRUE]
all_expanded[type == "positive" & dynamic == "inverse_plateau" & stage %in% c(1,2,6,7), is_active_signal := TRUE]

# 2. Filtrar datos para la evaluación de "Power Real"
#    Mantenemos TODOS los negativos (para especificidad/FP)
#    Mantenemos SOLO los positivos donde hay señal activa (para sensibilidad/TP)
data_for_plots <- all_expanded[true_label == 0 | (true_label == 1 & is_active_signal == TRUE)]

message(sprintf("Positivos originales: %d | Positivos activos usados: %d", 
                nrow(all_expanded[true_label==1]), nrow(data_for_plots[true_label==1])))

# 3. Función de métricas para Bootstrap
calc_perf_metrics <- function(d, indices) {
  d_boot <- d[indices, ]
  
  # Conteos básicos
  tp <- sum(d_boot$signal_giangreco == TRUE & d_boot$true_label == 1, na.rm = TRUE)
  fn <- sum(d_boot$signal_giangreco == FALSE & d_boot$true_label == 1, na.rm = TRUE)
  fp <- sum(d_boot$signal_giangreco == TRUE & d_boot$true_label == 0, na.rm = TRUE)
  tn <- sum(d_boot$signal_giangreco == FALSE & d_boot$true_label == 0, na.rm = TRUE)
  
  # Métricas
  sens <- if((tp + fn) > 0) tp / (tp + fn) else NA  # Power
  ppv  <- if((tp + fp) > 0) tp / (tp + fp) else NA
  npv  <- if((tn + fn) > 0) tn / (tn + fn) else NA
  
  # AUROC (Probabilístico)
  roc_val <- tryCatch({
    r <- roc(d_boot$true_label, d_boot$ior_value, direction="<", levels=c(0,1), quiet=TRUE)
    as.numeric(auc(r))
  }, error = function(e) NA)
  
  return(c(AUROC = roc_val, Power = sens, PPV = ppv, NPV = npv))
}

# 4. Ejecutar Bootstrap por Etapa
metric_names <- c("AUROC", "Power", "PPV", "NPV")
stage_boot_list <- list()
set.seed(seed_base)

for (s in unique(data_for_plots$stage_name)) {
  stage_subset <- data_for_plots[stage_name == s]
  
  # Validar que existan ambas clases en esta etapa tras el filtro
  if (length(unique(stage_subset$true_label)) == 2) {
    
    boot_res <- boot(data = stage_subset, statistic = calc_perf_metrics, R = 500)
    
    # Calcular medias e IC
    means <- colMeans(boot_res$t, na.rm = TRUE)
    ci_vals <- apply(boot_res$t, 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
    
    # Crear DT explícito
    dt_res <- data.table(
      stage_name = s,
      metric = metric_names,
      mean = as.numeric(means),
      lower = as.numeric(ci_vals[1, ]),
      upper = as.numeric(ci_vals[2, ])
    )
    stage_boot_list[[s]] <- dt_res
  }
}

stage_metrics_boot <- rbindlist(stage_boot_list)

# ============================================================================
# GENERAR GRÁFICOS SOLICITADOS
# ============================================================================

if (nrow(stage_metrics_boot) > 0) {
  
  # Ordenar factores para visualización correcta
  stage_metrics_boot[, stage_name := factor(stage_name, levels = niveles_nichd)]
  stage_metrics_boot[, metric := factor(metric, levels = metric_names)]
  
  # --- GRÁFICO 1: 4 PANELES POR ETAPA ---
  p_panels <- ggplot(stage_metrics_boot, aes(x = stage_name, y = mean, group = 1)) +
    # Puntos centrales
    geom_point(color = "#2c3e50", size = 2.5) +
    # Barras de error
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, color = "#2c3e50") +
    # Facetas
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    # Estética
    theme_giangreco() +
    labs(
      title = "Métricas de rendimiento según etapa del desarrollo (en señales activas)",
      x = "", 
      y = "Performance"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 12, color = "white"),
      strip.background = element_rect(fill = "#2c3e50")
    )
  
  print(p_panels)
  ggsave(paste0(output_dir, "fig_panel_stage_performance_active.png"), 
         p_panels, width = 12, height = 10, dpi = 300)
}







# Asegúrate de estar en el directorio correcto donde se guardó este archivo en el paso 03
file_path <- "./validation_results_p95/all_expanded_with_criteria.csv" 

if (!file.exists(file_path)) {
  stop("No se encuentra el archivo all_expanded_with_criteria.csv. Verifica la ruta.")
}

dt <- fread(file_path)

# 2. Calcular el Ancho del Intervalo de Confianza (Hacia el límite inferior)
# Width = Estimación Central - Límite Inferior
# Nota: Como log_ior suele ser mayor que lower90, esto da un valor positivo que representa la incertidumbre.
dt[, ci_width := log_ior - log_ior_lower90]

# 3. Definir Categorías para el Color
# Umbral de "Locura numérica": Log-IOR > 10 (aprox IOR > 22,000)
umbral_extremo <- 10

dt[, category := fcase(
  # Caso 1: Valores dentro de rango normal
  abs(log_ior) <= umbral_extremo, "Rango habitual",
  
  # Caso 2: Extremo pero NO detectado (El sistema funcionó: IC muy ancho o cruzó 0)
  abs(log_ior) > umbral_extremo & signal_giangreco == FALSE, "NO-detectado",
  
  # Caso 3: Extremo y DETECTADO (Ojo: Puede ser FP técnico o señal masiva)
  abs(log_ior) > umbral_extremo & signal_giangreco == TRUE, "Detectado"
)]

# Ordenar niveles para que la leyenda salga bonita
dt[, category := factor(category, levels = c("Rango habitual", 
                                             "NO-detectado", 
                                             "Detectado"))]

# 4. Generar el Gráfico
p_diagnostic <- ggplot(dt, aes(x = log_ior, y = ci_width, color = category)) +
  # Puntos
  geom_point(aes(shape = type), alpha = 0.6, size = 4) +
  
  # Colores manuales: 
  # Gris para normal, 
  # Azul/Verde para los rechazados (sistema funcionando), 
  # Rojo para los detectados (alerta)
  scale_color_manual(values = c(
    "Rango habitual" = "gray70",
    "NO-detectado" = "#2ecc71", # Verde: El filtro funcionó
    "Detectado" = "#e74c3c" # Rojo: Señal extrema aceptada
  )) +
  
  # Líneas de referencia para el umbral de extremo
  geom_vline(xintercept = c(-umbral_extremo, umbral_extremo), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  
  # Títulos y etiquetas
  labs(
    title = "Comportamiento de señales. Magnitud vs Incertidumbre",
    x = "Log IOR (Estimación Puntual)",
    y = "Incertidumbre (Log IOR - Límite Inferior 90%)",
    color = "Estado de Detección"
  ) +
  
  # Tema limpio
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  ) +
  
  # Zoom para no perder de vista los datos si hay un outlier infinito
  coord_cartesian(xlim = c(-30, 30), ylim = c(0, 100)) 

# 5. Mostrar y Guardar
print(p_diagnostic)

ggsave("./validation_results_p95/fig_diagnostic_extremes_classification.png", 
       p_diagnostic, width = 10, height = 7, dpi = 300)

# 6. Resumen Numérico en consola
message("\n--- Resumen de Casos Extremos ---")
print(dt[, .N, by = category])
