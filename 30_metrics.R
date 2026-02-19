################################################################################
# Script de análisis de métricas y sensibilidad 
# Script: 30_metrics.R
################################################################################

library(data.table)
library(tidyverse)
library(pROC)
library(svglite)

setwd("D:/Bioestadística/gam-farmacovigilancia")
set.seed(9427)

source("00_functions.R", local = TRUE)
source("01_theme.R", local = TRUE)

niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

################################################################################
# Configuración
################################################################################

# percentilo elegido de distribución nula
percentil <- "p95"
n_boot <- 2000

# parámetros para calculo de subset de poder
target_power <- 0.80
grid_res <- 30
tij_max <- 0.15  
n_max <- 250

# parámetros para uso de distribución nula
use_threshold_ior <- TRUE
use_threshold_reri <- TRUE

# parámetros para correr según formula guardada
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
k_spline <- 7
include_nichd <- FALSE
nichd_spline <- FALSE
bs_type <- "cs"
select <- FALSE

# parámetros para guardado según formula usada
suffix_sensitivity <- paste0(
  "sens_",
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

suffix_base <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

# clasificación de etapas con alto reporte según dinámica
stage_class <- rbind(
  data.table(nichd = niveles_nichd, dynamic = "uniform", class = 1),
  data.table(nichd = niveles_nichd, dynamic = "increase", class = c(0, 0, NA, NA, NA, 1, 1)),
  data.table(nichd = niveles_nichd, dynamic = "decrease", class = c(1, 1, NA, NA, NA, 0, 0)),
  data.table(nichd = niveles_nichd, dynamic = "plateau", class = c(0, NA, 1, 1, 1, NA, 0)),
  data.table(nichd = niveles_nichd, dynamic = "inverse_plateau", class = c(1, NA, 0, 0, 0, NA, 1))
)

# rutas
ruta_base_sensitivity <- paste0("./results/", suffix_sensitivity, "/augmentation_results/")
output_dir <- paste0("./results/", suffix_sensitivity, "/metrics_results_test/")

ruta_coadmin_pos <- paste0("./results/", suffix_sensitivity, "/augmentation_results/positive_coadmin_by_stage.csv")
ruta_coadmin_neg <- paste0("./results/", suffix_sensitivity, "/augmentation_results/negative_coadmin_by_stage.csv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# niveles de reducción
reduction_levels <- c(0, seq(10, 90, by = 10))

################################################################################
# Etiquetas para gráficos
################################################################################

# etiquetas para traducción 
dynamic_labels <- c(
  "uniform" = "Uniforme",
  "increase" = "Incremento",
  "decrease" = "Disminución",
  "plateau" = "Meseta",
  "inverse_plateau" = "Valle"
)

nichd_labels <- c(
  "term_neonatal" = "Neonato a término",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temprana",
  "late_adolescence" = "Adolescencia tardía"
)

metrics_to_plot <- c("sensitivity", "specificity", "PPV", "NPV", "F1")
metric_labels <- c(
  "sensitivity" = "Sensibilidad",
  "specificity" = "Especificidad",
  "PPV" = "Valor Predictivo Positivo",
  "NPV" = "Valor Predictivo Negativo",
  "F1" = "F1-Score"
)

method_pairs <- list(
  list(gam = "GAM-logIOR", classic = "Estratificado-IOR", label = "IOR"),
  list(gam = "GAM-RERI", classic = "Estratificado-RERI", label = "RERI"),
  list(gam = "GAM-Doble", classic = "Estratificado-Doble", label = "Doble")
)

################################################################################
# Carga de umbrales nulos
################################################################################
ruta_null_dist <- paste0("./results/", suffix_base, "/null_distribution_results/null_distribution.csv")

null_thresholds_ior <- fread(paste0("./results/", suffix_base, "/null_distribution_results/null_thresholds.csv"))
thresh_col <- paste0("threshold_", percentil)
null_thresholds_ior <- null_thresholds_ior[, .(stage, threshold_ior = get(thresh_col))]

null_thresholds_reri <- fread(paste0("./results/", suffix_base, "/null_distribution_results/null_thresholds_reri.csv"))
null_thresholds_reri <- null_thresholds_reri[, .(stage, threshold_reri = get(thresh_col))]

null_thresholds <- merge(null_thresholds_ior, null_thresholds_reri, by = "stage")

# datos de coadministración
coadmin_stage_pos <- fread(ruta_coadmin_pos)
coadmin_stage_neg <- fread(ruta_coadmin_neg)
setnames(coadmin_stage_pos, "nichd_num", "stage_num")
setnames(coadmin_stage_neg, "nichd_num", "stage_num")

################################################################################
# Función de cálculo de métricas por bootstrap
################################################################################

# Calcula métricas de clasificación con IC95% por bootstrap
# 
# Parámetros:
# dt data.table con columnas: triplet_id, detected, label
# n_boot número de replicaciones bootstrap
# agregar_por_triplete si TRUE, agrega a nivel triplete con any(). si FALSE, mantiene granularidad actual (etapa o dinámica)
# 
# Return: 
# data.table con métricas e IC95%

calcular_metricas_simple <- function(dt, n_boot = 2000, agregar_por_triplete = TRUE, score_type, score_type_auc) {
  
  # Agregación a nivel triplete si se solicita
  if (agregar_por_triplete) {
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
    # Filtrar valores NA e infinitos antes de calcular ROC  
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
      
      # Verificar que hay ambas clases en la muestra bootstrap
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
  if (agregar_por_triplete) {
    # Identificar IDs únicos para bootstrap por triplete
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
    # Bootstrap por FILA para análisis por etapa
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
  
  # Cálculo de métricas puntuales como MEDIA del bootstrap (alineado con IC)
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
  
  # Cálculo de IC percentílicos
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
  
  # Conteos del dataset original (no bootstrap) para referencia
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
# Función de expansión de datos con reducción
################################################################################

# Función que expande datos a formato largo
#
# Implementación:
# Carga datos y crea objetos según nivel de reducción
# Filtra inyecciones fallidas
# Expande a formato largo
# Filtra por etapas con alto reporte según dinamica

expandir_datos <- function(red_pct) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)  # objeto con nombre según nivel de reducción
  
  ruta_pos <- paste0(ruta_base_sensitivity, "positive_triplets_results", suffix_file, ".rds")
  ruta_neg <- paste0(ruta_base_sensitivity, "negative_triplets_results", suffix_file, ".rds")
  
  pos_raw <- readRDS(ruta_pos)
  neg_raw <- readRDS(ruta_neg)
  
  # Filtrar solo inyecciones exitosas
  pos_valid <- pos_raw[injection_success == TRUE]
  
  # Expandir a formato largo
  pos_exp <- expand_clean_all_metrics(pos_valid, 1, null_thresholds, 
                                      use_threshold_ior, use_threshold_reri)
  neg_exp <- expand_clean_all_metrics(neg_raw, 0, null_thresholds,
                                      use_threshold_ior, use_threshold_reri)
  
  # Merge con clasificación de etapas
  pos_exp <- merge(pos_exp, stage_class, by = c("nichd", "dynamic"), all.x = TRUE)
  neg_exp[, class := 0]
  
  # Merge con datos de coadministración
  pos_exp <- merge(pos_exp, 
                   coadmin_stage_pos[, .(triplet_id, stage_num, n_coadmin_stage)], 
                   by = c("triplet_id", "stage_num"), all.x = TRUE)
  neg_exp <- merge(neg_exp, 
                   coadmin_stage_neg[, .(triplet_id, stage_num, n_coadmin_stage)], 
                   by = c("triplet_id", "stage_num"), all.x = TRUE)
  
  # Dataset de alto reporte (excluir uniform)
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
# Función de detección de señal
################################################################################

aplicar_deteccion <- function(dt, method_name, detection_type, use_null) {
  
  is_gam <- grepl("GAM", method_name) # busca GAM en el string del método
  
  # Determinar columnas según método
  if (is_gam) {
    ior_col <- "gam_log_ior_lower90"
    reri_col <- "gam_reri_lower90"
    thresh_ior_col <- "threshold_ior"   # umbrales de distribución nula
    thresh_reri_col <- "threshold_reri"
  } else {
    ior_col <- "classic_log_ior_lower90"
    reri_col <- "classic_reri_lower90"
    thresh_ior_col <- NULL
    thresh_reri_col <- NULL
  }
  
  # Calcular detección según tipo
  if (detection_type == "IOR") {
    if (is_gam && use_null) {   # detección para gam con umbral de distribución nula
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)]
    } else {
      dt[, detected := !is.na(get(ior_col)) & get(ior_col) > 0]
    }
  } else if (detection_type == "RERI") {
    if (is_gam && use_null) {
      dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0 & get(reri_col) > get(thresh_reri_col)]
    } else {
      dt[, detected := !is.na(get(reri_col)) & get(reri_col) > 0]
    }
  } else { # Doble
    if (is_gam && use_null) {
      dt[, detected := (!is.na(get(ior_col)) & get(ior_col) > 0 & get(ior_col) > get(thresh_ior_col)) | # ambos criterios
                       (!is.na(get(reri_col)) & get(reri_col) > 0 & get(reri_col) > get(thresh_reri_col))]
    } else {
      dt[, detected := (!is.na(get(ior_col)) & get(ior_col) > 0) |
                       (!is.na(get(reri_col)) & get(reri_col) > 0)]
    }
  }
  
  # Reemplazar NA por FALSE
  dt[is.na(detected), detected := FALSE]
  
  return(dt)
}

################################################################################
# Carga de datos sin reducción 
################################################################################

data_baseline <- expandir_datos(0)  # carga datos sin reducción para cálculos de umbrales de poder
pos_high_base <- data_baseline$pos_high
neg_high_base <- data_baseline$neg_high

################################################################################
# Gráficos de comparación de distribuciones
################################################################################

null_distribution <- fread(ruta_null_dist)

# sampleo distribución nula para graficar
null_sample <- null_distribution[sample(.N, min(.N, 50000))]

comparison_data <- rbind(
    null_sample[, .(stage_num = stage, log_ior_lower90 = log_lower90, reri_lower90 = reri_lower90, source = "Distribución nula")],
    pos_high_base[, .(stage_num, log_ior_lower90 = gam_log_ior_lower90, reri_lower90 = gam_reri_lower90, source = "Positivos")],
    neg_high_base[, .(stage_num, log_ior_lower90 = gam_log_ior_lower90, reri_lower90 = gam_reri_lower90, source = "Negativos")]
)
comparison_data[, stage_name := factor(stage_name,
                                       levels = niveles_nichd,
                                       labels = nichd_labels)] # ordeno para las facetas

# corto data para mostrar en gráfico (esto lo debería hacer limitando ejes en el gráfico igual)
comparison_data <- comparison_data[ is.finite(log_ior_lower90)]
comparison_data <- comparison_data[ is.finite(reri_lower90)]

# Paleta de colores consistente
color_palette <- c(
  "Distribución nula" = "gray60",
  "Positivos" = "#4DAF4A",
  "Negativos" = "#E41A1C"
)

# Gráfico combinado con facets
p_null_vs_obs_ior <- ggplot(
  comparison_data, 
  aes(x = log_ior_lower90, fill = source)
) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~ stage_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(limits = c(-5, 5)) + # limito para visualizar mejor las distribuciones
  labs(
    title = sprintf("Distribución nula vs señales detectadas (%s)", percentil),
    x = "Log(IOR) - Límite inferior IC 90%",
    y = "Densidad",
    fill = "Fuente"
  )

# Guardar gráfico combinado
ggsave(
  paste0(output_dir, "fig_null_vs_observed_ior.png"), 
  p_null_vs_obs_ior,
  width = 16, 
  height = 12, 
  dpi = 300
)

print(p_null_vs_obs_ior)

# Gráfico combinado con facets
p_null_vs_obs_reri <- ggplot(
  comparison_data, 
  aes(x = reri_lower90, fill = source)
) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~ stage_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(limits = c(-0.25, 0.25)) + # limito para visualizar mejor las distribuciones
  labs(
    title = sprintf("Distribución nula vs señales detectadas (%s)", percentil),
    x = "RERI - Límite inferior IC 90%",
    y = "Densidad",
    fill = "Fuente"
  ) 
# Guardar gráfico combinado
ggsave(
  paste0(output_dir, "fig_null_vs_observed_reri.png"), 
  p_null_vs_obs_reri,
  width = 16, 
  height = 12, 
  dpi = 300
)

print(p_null_vs_obs_reri)

message(sprintf("Dataset base: %d positivos, %d negativos", nrow(pos_high_base), nrow(neg_high_base)))

################################################################################
# Cálculo de subsets de poder para cada método
################################################################################

# GAM-logIOR
message("GAM-logIOR")
power_gam_ior <- calculate_power_gam(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,
  use_threshold_ior = TRUE,
  use_threshold_reri = FALSE,
  detection = "ior"
)

# Heatmap de superficie de poder GAM
p_surface_gam_ior <- plot_power_surface(power_gam_ior, 
  target_power, detection = "Log-IOR", grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)
ggsave(
  paste0(output_dir, "fig_power_surface_gam_ior.png"), 
  p_surface_gam_ior,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_gam_ior)

# GAM-RERI
message("GAM-RERI")
power_gam_reri <- calculate_power_gam(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,
  use_threshold_ior = FALSE,
  use_threshold_reri = TRUE,
  detection = "reri"
)

# Heatmap de superficie de poder GAM-RERI
p_surface_gam_reri <- plot_power_surface(power_gam_reri, 
  target_power, detection = "RERI", grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)

ggsave(
  paste0(output_dir, "fig_power_surface_gam_reri.png"),
  p_surface_gam_reri,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_gam_reri)

# GAM-Doble
message("GAM-Doble")
power_gam_doble <- calculate_power_gam(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,
  use_threshold_ior = TRUE,
  use_threshold_reri = TRUE,
  detection = "double"
)

# Heatmap de superficie de poder GAM Doble
p_surface_gam <- plot_power_surface(power_gam_doble, 
  target_power, grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)

ggsave(
  paste0(output_dir, "fig_power_surface_gam_double.png"),
  p_surface_gam,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_gam)

# Estratificado-IOR
message("Estratificado-IOR")
power_cls_ior <- calculate_power_classic(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "ior",
  na_remove = TRUE
)

# Heatmap de superficie de poder para método Estratificado-IOR
p_surface_classic_ior <- plot_power_surface(power_cls_ior, 
  target_power, detection = "IOR", grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)

ggsave( 
  paste0(output_dir, "fig_power_surface_classic_ior.png"),
  p_surface_classic_ior,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_classic_ior)

# Estratificado-RERI
message("Estratificado-RERI")
power_cls_reri <- calculate_power_classic(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "reri",
  na_remove = TRUE
)

# Heatmap de superficie de poder para método Estratificado-RERI
p_surface_classic_reri <- plot_power_surface(power_cls_reri, 
  target_power, detection = "RERI", grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)

ggsave( 
  paste0(output_dir, "fig_power_surface_classic_reri.png"),
  p_surface_classic_reri,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_classic_reri)

# Estratificado-Doble
message("Estratificado-Doble")
power_cls_doble <- calculate_power_classic(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "double",
  na_remove = TRUE
)

# Heatmap de superficie de poder para método Estratificado-Doble
p_surface_classic <- plot_power_surface(
  power_cls_doble, 
  target_power, grid_size = grid_res,
  t_range = c(0, tij_max), n_range = c(0, n_max)
)

ggsave(
  paste0(output_dir, "fig_power_surface_classic_double.png"),
  p_surface_classic,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_surface_classic)

# Guardo IDs de tripletes para cada subset de poder
power_ids <- list(
  "GAM-logIOR" = unique(power_gam_ior$superset_pos[class == 1]$triplet_id),
  "GAM-RERI" = unique(power_gam_reri$superset_pos[class == 1]$triplet_id),
  "GAM-Doble" = unique(power_gam_doble$superset_pos[class == 1]$triplet_id),
  "Estratificado-IOR" = unique(power_cls_ior$superset_pos[class == 1]$triplet_id),
  "Estratificado-RERI" = unique(power_cls_reri$superset_pos[class == 1]$triplet_id),
  "Estratificado-Doble" = unique(power_cls_doble$superset_pos[class == 1]$triplet_id)
)

# Calcular intersección
triplets_intersection <- intersect(power_ids[["GAM-Doble"]], power_ids[["Estratificado-Doble"]])
message(sprintf("\nIntersección GAM-Doble & Estratificado-Doble: %d tripletes", length(triplets_intersection)))

# Resumen de subsets
summary_power <- data.table(
  method = names(power_ids),
  n_triplets = sapply(power_ids, length),
  t_star = c(power_gam_ior$t_star, power_gam_reri$t_star, power_gam_doble$t_star,
             power_cls_ior$t_star, power_cls_reri$t_star, power_cls_doble$t_star),
  n_star = c(power_gam_ior$n_star, power_gam_reri$n_star, power_gam_doble$n_star,
             power_cls_ior$n_star, power_cls_reri$n_star, power_cls_doble$n_star)
)
print(summary_power)

################################################################################
# Carga de datos con reducción 
################################################################################

# expando datos según niveles de reducción
datos_por_reduccion <- lapply(reduction_levels, expandir_datos)
names(datos_por_reduccion) <- as.character(reduction_levels)

################################################################################
# Definición de métodos 
################################################################################

# métodos con tipo de score mapeado 
metodos <- list(
  list(nombre = "GAM-logIOR", tipo = "IOR", es_gam = TRUE, 
       score_type = "gam_log_ior_lower90", score_type_auc = "gam_log_ior"),
  list(nombre = "GAM-RERI", tipo = "RERI", es_gam = TRUE, 
       score_type = "gam_reri_lower90", score_type_auc = "gam_reri"),
  list(nombre = "GAM-Doble", tipo = "Doble", es_gam = TRUE, 
       score_type = "gam_reri_lower90", score_type_auc = "gam_reri"),  # o podrías usar gam_reri_lower90
  list(nombre = "Estratificado-IOR", tipo = "IOR", es_gam = FALSE, 
       score_type = "classic_log_ior_lower90", score_type_auc = "classic_log_ior"),
  list(nombre = "Estratificado-RERI", tipo = "RERI", es_gam = FALSE, 
       score_type = "classic_reri_lower90", score_type_auc = "classic_reri"),
  list(nombre = "Estratificado-Doble", tipo = "Doble", es_gam = FALSE, 
       score_type = "classic_reri_lower90",  score_type_auc = "classic_reri")  # o classic_reri_lower90
)

################################################################################
# Cálculo de métricas
################################################################################

# Métricas originales sin reducción
res_global_original <- list()
res_dinamica_original <- list()
res_etapa_original <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo de niveles de reducción
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in metodos) {
    
    message(sprintf("  %s", met$nombre))
    
    # Métricas globales
    # Filtrar negativos solo a etapas de alto reporte
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    dt_global <- rbind(datos$pos_high, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    # Agregar a nivel triplete para métricas globales
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE,  # acá TRUE para evitar doble conteo
    score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`(
      method = met$nombre,
      reduction_pct = red_pct,
      dataset = "original"
    )]
    res_global_original[[length(res_global_original) + 1]] <- metrics_global
    
    # Métricas por dinámica
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # remuevo uniform (creo que redundante)
     
    for (dyn in dinamicas) {
      # Etapas altas para esta dinámica
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd] 
      pos_dyn <- datos$pos_high[nichd %in% etapas_altas & dynamic == dyn]
  
      # Negativos en las MISMAS etapas altas
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]
  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE)
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc) 
      metrics_dyn[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        dynamic = dyn,
        dataset = "original"
      )]
      res_dinamica_original[[length(res_dinamica_original) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa (NO AGREGAR POR TRIPLETE)
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      
      # Para positivos: solo etapas clasificadas como altas (redundante)
      pos_stage <- datos$pos_high[stage_num == s & class == 1]
      
      # Negativos de la misma etapa
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      # NO agregar - cada fila es una etapa específica
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        stage_num = s,
        nichd = nichd_label,
        dataset = "original"
      )]
      res_etapa_original[[length(res_etapa_original) + 1]] <- metrics_stage
    }
  }
}

metrics_global_original <- rbindlist(res_global_original, fill = TRUE)
metrics_dynamic_original <- rbindlist(res_dinamica_original, fill = TRUE)
metrics_stage_original <- rbindlist(res_etapa_original, fill = TRUE)

fwrite(metrics_global_original, paste0(output_dir, "metrics_global_original.csv"))
fwrite(metrics_dynamic_original, paste0(output_dir, "metrics_dynamic_original.csv"))
fwrite(metrics_stage_original, paste0(output_dir, "metrics_stage_original.csv"))

################################################################################
# Cálculo de métricas en subset de poder estadístico
################################################################################

res_global_filtered <- list()
res_dinamica_filtered <- list()
res_etapa_filtered <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in metodos) {
    
    message(sprintf("  %s", met$nombre))
    
    ids_filtrar <- power_ids[[met$nombre]]
    
    # Métricas globales filtradas
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    pos_fil <- datos$pos_high[triplet_id %in% ids_filtrar]
    dt_global <- rbind(pos_fil, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE, 
      score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`(
      method = met$nombre,
      reduction_pct = red_pct,
      dataset = "filtered"
    )]
    res_global_filtered[[length(res_global_filtered) + 1]] <- metrics_global
    
    # Métricas por dinámica filtradas
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # (redundante)
    
    for (dyn in dinamicas) {
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd]
      pos_dyn <- datos$pos_high[triplet_id %in% ids_filtrar & nichd %in% etapas_altas & dynamic == dyn]
      
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE) 
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_dyn[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        dynamic = dyn,
        dataset = "filtered"
      )]
      res_dinamica_filtered[[length(res_dinamica_filtered) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa filtradas
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      pos_stage <- datos$pos_high[triplet_id %in% ids_filtrar & stage_num == s & class == 1]
      
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        stage_num = s,
        nichd = nichd_label,
        dataset = "filtered"
      )]
      res_etapa_filtered[[length(res_etapa_filtered) + 1]] <- metrics_stage
    }
  }
}

metrics_global_filtered <- rbindlist(res_global_filtered, fill = TRUE)
metrics_dynamic_filtered <- rbindlist(res_dinamica_filtered, fill = TRUE)
metrics_stage_filtered <- rbindlist(res_etapa_filtered, fill = TRUE)

fwrite(metrics_global_filtered, paste0(output_dir, "metrics_global_filtered.csv"))
fwrite(metrics_dynamic_filtered, paste0(output_dir, "metrics_dynamic_filtered.csv"))
fwrite(metrics_stage_filtered, paste0(output_dir, "metrics_stage_filtered.csv"))

################################################################################
# Cálculo de métricas en subset de intersección
################################################################################

metodos_interseccion <- metodos[sapply(metodos, function(m) m$nombre %in% c("GAM-Doble", "Estratificado-Doble"))]

res_global_inter <- list()
res_dinamica_inter <- list()
res_etapa_inter <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in metodos_interseccion) {
    
    message(sprintf("  %s", met$nombre))
    
    pos_int <- datos$pos_high[triplet_id %in% triplets_intersection]
    
    # Métricas a nivel global intersección
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    pos_int <- datos$pos_high[triplet_id %in% triplets_intersection]
    dt_global <- rbind(pos_int, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE, 
      score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`(
      method = met$nombre,
      reduction_pct = red_pct,
      dataset = "intersection"
    )]
    res_global_inter[[length(res_global_inter) + 1]] <- metrics_global
    
    # Métricas por dinámica en intersección 
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # redundante
    
    for (dyn in dinamicas) {
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd]
      pos_dyn <- pos_int[nichd %in% etapas_altas & dynamic == dyn]
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]  
  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE)
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_dyn[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        dynamic = dyn,
        dataset = "intersection"
      )]
      res_dinamica_inter[[length(res_dinamica_inter) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa en intersección 
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      pos_stage <- pos_int[stage_num == s & class == 1] # redundante class == 1
      
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        stage_num = s,
        nichd = nichd_label,
        dataset = "intersection"
      )]
      res_etapa_inter[[length(res_etapa_inter) + 1]] <- metrics_stage
    }
  }
}

metrics_global_intersection <- rbindlist(res_global_inter, fill = TRUE)
metrics_dynamic_intersection <- rbindlist(res_dinamica_inter, fill = TRUE)
metrics_stage_intersection <- rbindlist(res_etapa_inter, fill = TRUE)

fwrite(metrics_global_intersection, paste0(output_dir, "metrics_global_intersection.csv"))
fwrite(metrics_dynamic_intersection, paste0(output_dir, "metrics_dynamic_intersection.csv"))
fwrite(metrics_stage_intersection, paste0(output_dir, "metrics_stage_intersection.csv"))
