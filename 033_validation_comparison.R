################################################################################
# Script de validación que compara GAM vs IOR clásico
# Script 033_validation_comparison
################################################################################

library(data.table)
library(tidyverse)
library(pROC)
library(parallel)

setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)
source("giangreco_theme.R")
theme_set(theme_giangreco())


niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

################################################################################
# Configuración
################################################################################

PERCENTILE_LEVEL <- "p90"

spline_individuales <- TRUE
include_sex <- TRUE 
include_stage_sex <- FALSE 
k_spline <- 7
nichd_spline <- TRUE
bs_type <- "cs"
select <- FALSE

suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

percentile_config <- list(
  p90 = list(threshold_col = "threshold_p90", alpha = 0.10, z_score = qnorm(0.95)),
  p95 = list(threshold_col = "threshold_p95", alpha = 0.05, z_score = qnorm(0.975)),
  p99 = list(threshold_col = "threshold_p99", alpha = 0.01, z_score = qnorm(0.995))
)

config <- percentile_config[[PERCENTILE_LEVEL]]

ruta_pos_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
ruta_neg_results <- paste0("./results/", suffix, "/augmentation_results/negative_triplets_results.rds")
ruta_null_thresh <- paste0("./results/", suffix, "/null_distribution_results/null_thresholds.csv")

output_dir <- paste0("./results/", suffix, "/comparison_", PERCENTILE_LEVEL, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Carga de datos
################################################################################

positives_scores <- readRDS(ruta_pos_results)
negatives_scores <- readRDS(ruta_neg_results)
null_thresholds <- fread(ruta_null_thresh)
null_thresholds[, threshold := get(config$threshold_col)]

pos_valid <- positives_scores
neg_valid <- negatives_scores

message(sprintf("\nDatos válidos:"))
message(sprintf("  Positivos: %d", nrow(pos_valid)))
message(sprintf("  Negativos: %d", nrow(neg_valid)))

################################################################################
# Expansión y clasificación
################################################################################

pos_gam <- expand_and_classify(pos_valid, "log_ior_lower90", "GAM")
neg_gam <- expand_and_classify(neg_valid, "log_ior_lower90", "GAM")

pos_classic <- expand_and_classify(pos_valid, "log_ior_classic_lower90", "Clásico")
neg_classic <- expand_and_classify(neg_valid, "log_ior_classic_lower90", "Clásico")

pos_gam[, true_label := 1]
neg_gam[, true_label := 0]
pos_classic[, true_label := 1]
neg_classic[, true_label := 0]

all_data <- rbind(
  rbind(pos_gam, neg_gam),
  rbind(pos_classic, neg_classic),
  fill = TRUE
)

################################################################################
# Métricas por método
################################################################################

metrics_gam <- calculate_metrics(rbind(pos_gam, neg_gam), "GAM")
metrics_classic <- calculate_metrics(rbind(pos_classic, neg_classic), "Clásico")

comparison_metrics <- rbind(metrics_gam, metrics_classic)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("COMPARACIÓN DE MÉTODOS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
print(comparison_metrics[, .(method, sensitivity, specificity, ppv, accuracy, f1_score)])

fwrite(comparison_metrics, paste0(output_dir, "comparison_metrics.csv"))

################################################################################
# Curvas ROC
################################################################################

# preparo datos para GAM
scores_gam <- rbind(
  pos_valid[, .(triplet_id, max_ior, type = "positive")],
  neg_valid[, .(triplet_id, max_ior, type = "negative")]
)

scores_gam[, true_label := ifelse(type == "positive", 1, 0)]
scores_gam <- scores_gam[is.finite(max_ior) & !is.na(max_ior)]

# calculo de max_ior_classic para curva ROC
calculate_max_ior_classic <- function(dt) {
  dt[, {
    ior_vals <- unlist(ior_classic)
    max_val <- if(length(ior_vals) > 0 && !all(is.na(ior_vals))) {
      max(ior_vals, na.rm = TRUE)
    } else {
      NA_real_
    }
    list(max_ior_classic = max_val)
  }, by = triplet_id]
}

pos_classic_max <- calculate_max_ior_classic(pos_valid)
neg_classic_max <- calculate_max_ior_classic(neg_valid)

scores_classic <- rbind(
  merge(pos_valid[, .(triplet_id)], pos_classic_max, by = "triplet_id")[, type := "positive"],
  merge(neg_valid[, .(triplet_id)], neg_classic_max, by = "triplet_id")[, type := "negative"]
)

scores_classic[, true_label := ifelse(type == "positive", 1, 0)]
scores_classic <- scores_classic[is.finite(max_ior_classic) & !is.na(max_ior_classic)]

# ROC
roc_gam <- roc(
  scores_gam$true_label, 
  scores_gam$max_ior,  # valor puntual del IOR
  direction = "<", 
  levels = c(0, 1), 
  quiet = TRUE
)

roc_classic <- roc(
  scores_classic$true_label,
  scores_classic$max_ior_classic,
  direction = "<",
  levels = c(0, 1),
  quiet = TRUE
)

auroc_gam <- auc(roc_gam)
auroc_classic <- auc(roc_classic)

auroc_comparison <- data.table(
  method = c("GAM", "Clásico"),
  auroc = c(auroc_gam, auroc_classic)
)

cat("\nAUROC:\n")
print(auroc_comparison)

fwrite(auroc_comparison, paste0(output_dir, "auroc_comparison.csv"))

################################################################################
# Gráficos
################################################################################

roc_data <- rbind(
  data.table(
    method = "GAM",
    sensitivity = roc_gam$sensitivities,
    specificity = roc_gam$specificities
  ),
  data.table(
    method = "Clásico",
    sensitivity = roc_classic$sensitivities,
    specificity = roc_classic$specificities
  )
)

p_roc <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("GAM" = "#2C3E50", "Clásico" = "#E74C3C")) +
  annotate("text", x = 0.7, y = 0.3, 
           label = sprintf("AUC GAM = %.3f\nAUC Clásico = %.3f", 
                           auroc_gam, auroc_classic),
           size = 4, fontface = "bold") +
  labs(
    title = "Comparación de curvas ROC",
    subtitle = sprintf("GAM vs IOR Clásico (%s)", config$threshold_col),
    x = "1 - Especificidad (FPR)",
    y = "Sensibilidad (TPR)",
    color = "Método"
  ) +
  theme(legend.position = c(0.8, 0.2))

ggsave(paste0(output_dir, "fig_roc_comparison.png"), 
       p_roc, width = 8, height = 7, dpi = 300)

print(p_roc)

metrics_long <- melt(
  comparison_metrics[, .(method, sensitivity, specificity, ppv, accuracy, f1_score)],
  id.vars = "method",
  variable.name = "metric",
  value.name = "value"
)

p_metrics <- ggplot(metrics_long, aes(x = metric, y = value, fill = method)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("GAM" = "#2C3E50", "Clásico" = "#E74C3C")) +
  labs(
    title = "Comparación de métricas",
    subtitle = sprintf("Umbral: %s", config$threshold_col),
    x = "Métrica",
    y = "Valor",
    fill = "Método"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir, "fig_metrics_comparison.png"),
       p_metrics, width = 10, height = 6, dpi = 300)

print(p_metrics)

################################################################################
# Análisis por etapa
################################################################################

# Configuración Bootstrap
n_boot <- 2000
set.seed(9427) 

metrics_by_stage <- all_data[, {
  # preparación de vectores para ensibilidad 
  # true_label = 1 son los positivos inyectados
  hits_vec <- signal_detected[true_label == 1]
  n_pos <- length(hits_vec)
  
  # bootstrap para sensibilidad
  if (n_pos > 0) {
    # estimación puntual
    sens_point <- mean(hits_vec, na.rm = TRUE)
    
    # replicación bootstrap
    # muestreo con reemplazo del vector acierto/fallo
    boot_dist <- replicate(n_boot, {
      mean(sample(hits_vec, n_pos, replace = TRUE), na.rm = TRUE)
    })
    
    # intervalos percentiles (2.5% y 97.5%)
    sens_low <- quantile(boot_dist, 0.025, na.rm = TRUE)
    sens_high <- quantile(boot_dist, 0.975, na.rm = TRUE)
  } else {
    sens_point <- NA_real_
    sens_low <- NA_real_
    sens_high <- NA_real_
  }
  
  # resto de métricas
  # conteos
  tp <- sum(hits_vec, na.rm = TRUE) # ya filtrado true_label==1 arriba
  # especificidad necesita los negativos (true_label == 0)
  neg_vec <- signal_detected[true_label == 0]
  tn <- sum(neg_vec == FALSE, na.rm = TRUE)
  fp <- sum(neg_vec == TRUE, na.rm = TRUE)
  
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  ppv_val <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  
  # resultados
  list(
    sensitivity = as.double(sens_point),
    sens_lower = as.double(sens_low),
    sens_upper = as.double(sens_high),
    specificity = as.double(spec),
    ppv = as.double(ppv_val),
    n_signals = sum(signal_detected, na.rm = TRUE),
    n_positives_ground_truth = n_pos
  )
}, by = .(stage, method)]

metrics_by_stage[, stage_name := niveles_nichd[stage]]

metrics_by_stage[, stage_name := factor(stage_name, levels = niveles_nichd)]


fwrite(metrics_by_stage, paste0(output_dir, "metrics_by_stage.csv"))


pd <- position_dodge(width = 0.5)

p_stage <- ggplot(metrics_by_stage[!is.na(sensitivity)],
                  aes(x = stage_name, y = sensitivity, color = method, group = method)) +
  
  # IC95
  geom_errorbar(aes(ymin = sens_lower, ymax = sens_upper), 
                width = 0.4, position = pd, alpha = 0.7) +
  
  # puntos centrales
  geom_point(size = 5, position = pd) +
  
  scale_color_manual(values = c("GAM" = "#2C3E50", "Clásico" = "#E74C3C")) +
  
  labs(
    title = "Sensibilidad por etapa (IC 95% por bootstrapping)",
    subtitle = sprintf("Comparación GAM vs IOR clásico (%s)", config$threshold_col),
    x = "Etapa",
    y = "Sensibilidad",
    color = "Método"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir, "fig_sensitivity_by_stage.png"),
       p_stage, width = 12, height = 7, dpi = 300)

print(p_stage)

################################################################################
# Resumen ejecutivo
################################################################################

mejora_sens <- (metrics_gam$sensitivity - metrics_classic$sensitivity) / 
               metrics_classic$sensitivity * 100
mejora_f1 <- (metrics_gam$f1_score - metrics_classic$f1_score) / 
             metrics_classic$f1_score * 100

executive_summary <- sprintf("
================================================================================
COMPARACIÓN GAM VS IOR CLÁSICO (%s)
================================================================================

Métricas Globales:
                    GAM         Clásico     Mejora
  AUROC:            %.3f       %.3f       %+.1f%%
  Sensibilidad:     %.3f       %.3f       %+.1f%%
  Especificidad:    %.3f       %.3f       %+.1f%%
  PPV:              %.3f       %.3f       %+.1f%%
  F1-Score:         %.3f       %.3f       %+.1f%%

Matriz de Confusión - GAM:
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d

Matriz de Confusión - Clásico:
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d

Etapas con mejora (Sensibilidad): %d de %d
Fecha: %s
================================================================================
",
config$threshold_col,
auroc_gam, auroc_classic, 
(auroc_gam - auroc_classic) / auroc_classic * 100,
metrics_gam$sensitivity, metrics_classic$sensitivity, mejora_sens,
metrics_gam$specificity, metrics_classic$specificity,
(metrics_gam$specificity - metrics_classic$specificity) / metrics_classic$specificity * 100,
metrics_gam$ppv, metrics_classic$ppv,
(metrics_gam$ppv - metrics_classic$ppv) / metrics_classic$ppv * 100,
metrics_gam$f1_score, metrics_classic$f1_score, mejora_f1,
metrics_gam$tp, metrics_gam$fn,
metrics_gam$fp, metrics_gam$tn,
metrics_classic$tp, metrics_classic$fn,
metrics_classic$fp, metrics_classic$tn,
sum(metrics_by_stage[method == "GAM", sensitivity] > 
    metrics_by_stage[method == "Clásico", sensitivity], na.rm = TRUE),
7,
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

cat(executive_summary)
writeLines(executive_summary, paste0(output_dir, "resumen_comparacion.txt"))

