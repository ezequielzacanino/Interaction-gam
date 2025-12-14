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

# Permite ejecutar el script tanto con el setwd del entorno original (Windows)
# como desde el root del repo.
project_dir <- "D:/Bioestadística/gam-farmacovigilancia"
if (dir.exists(project_dir)) setwd(project_dir)

source("00_functions.R", local = TRUE)
if (file.exists("giangreco_theme.R")) {
  source("giangreco_theme.R")
  if (exists("theme_giangreco")) theme_set(theme_giangreco())
}

# Parámetros de análisis de poder
fold_change_levels <- seq(1.5, 5.0, by = 0.5)
min_reports_levels <- seq(5, 50, by = 5)
target_power <- 0.80
n_cores <- max(1, floor(detectCores() * 0.75))

# Configuración de parámetros del modelo 
spline_individuales <- TRUE
include_sex <- TRUE
include_stage_sex <- FALSE
k_spline <- 7
nichd_spline <- TRUE
bs_type <- "cs"
select <- FALSE
method <- "fREML"

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

stopifnot(file.exists(ruta_ade_raw))
stopifnot(file.exists(ruta_pos_results))
stopifnot(file.exists(ruta_neg_results))
stopifnot(file.exists(ruta_null_thresh))

ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

positives_scores <- readRDS(ruta_pos_results)
negatives_scores <- readRDS(ruta_neg_results)
null_thresholds <- fread(ruta_null_thresh)

setDT(positives_scores)
setDT(negatives_scores)

message("Dataset cargado:")
message(sprintf("  Positivos: %d", nrow(positives_scores)))
message(sprintf("  Negativos: %d", nrow(negatives_scores)))

################################################################################
# Funciones auxiliares para análisis de poder
################################################################################

is_signal_detected <- function(log_ior_lower90, threshold = 0) {
  !is.infinite(log_ior_lower90) & !is.na(log_ior_lower90) & (log_ior_lower90 > threshold)
}

calculate_power <- function(positive_data, negative_data, fold_change_threshold,
                            min_reports_threshold, percentile_level = "p95") {
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
      n_pos = 0L,
      n_neg = 0L
    ))
  }

  threshold_col <- paste0("threshold_", percentile_level)
  avg_threshold <- mean(null_thresholds[[threshold_col]], na.rm = TRUE)

  pos_expanded <- filtered_pos[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]

  neg_expanded <- filtered_neg[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90)
    )
  }, by = triplet_id]

  pos_detected <- pos_expanded[, .(
    detected_gam = any(is_signal_detected(log_ior_lower90, avg_threshold)),
    detected_classic = any(is_signal_detected(log_ior_classic_lower90, 0))
  ), by = triplet_id]

  power_gam <- mean(pos_detected$detected_gam, na.rm = TRUE)
  power_classic <- mean(pos_detected$detected_classic, na.rm = TRUE)

  data.table(
    fold_change = fold_change_threshold,
    min_reports = min_reports_threshold,
    power_gam = power_gam,
    power_classic = power_classic,
    n_pos = nrow(filtered_pos),
    n_neg = nrow(filtered_neg)
  )
}

pick_thresholds <- function(power_dt, power_col, target_power) {
  power_dt <- copy(power_dt)
  power_dt <- power_dt[!is.na(get(power_col))]

  if (nrow(power_dt) == 0) {
    return(list(
      min_fc = NA_real_,
      min_reports = NA_real_,
      achieved_power = NA_real_,
      reached = FALSE,
      rows_target = power_dt[0]
    ))
  }

  rows_target <- power_dt[get(power_col) >= target_power]

  if (nrow(rows_target) > 0) {
    min_fc <- min(rows_target$fold_change)
    min_reports_val <- min(rows_target$min_reports[rows_target$fold_change == min_fc])
    achieved_power <- rows_target[
      fold_change == min_fc & min_reports == min_reports_val,
      get(power_col)
    ][1]

    return(list(
      min_fc = min_fc,
      min_reports = min_reports_val,
      achieved_power = achieved_power,
      reached = TRUE,
      rows_target = rows_target
    ))
  }

  best_power <- max(power_dt[[power_col]], na.rm = TRUE)
  best_rows <- power_dt[get(power_col) == best_power]

  min_fc <- min(best_rows$fold_change)
  min_reports <- min(best_rows$min_reports[best_rows$fold_change == min_fc])

  list(
    min_fc = min_fc,
    min_reports = min_reports,
    achieved_power = best_power,
    reached = FALSE,
    rows_target = best_rows
  )
}

################################################################################
# Cálculo de matriz de poder
################################################################################

power_grid <- CJ(fold_change = fold_change_levels, min_reports = min_reports_levels)

power_results <- pblapply(seq_len(nrow(power_grid)), function(i) {
  calculate_power(
    positive_data = positives_scores,
    negative_data = negatives_scores,
    fold_change_threshold = power_grid$fold_change[i],
    min_reports_threshold = power_grid$min_reports[i],
    percentile_level = "p95"
  )
})

power_matrix <- rbindlist(power_results, use.names = TRUE, fill = TRUE)

fwrite(power_matrix, paste0(output_dir, "power_matrix.csv"))
message(sprintf("Matriz de poder calculada para %d combinaciones", nrow(power_matrix)))

################################################################################
# Identificación de umbrales
################################################################################

message("\nIdentificando umbrales para 80% de poder...")

thr_gam <- pick_thresholds(power_matrix, "power_gam", target_power)
thr_classic <- pick_thresholds(power_matrix, "power_classic", target_power)

thresholds_80 <- data.table(
  method = c("GAM", "IOR_Clásico"),
  min_fold_change = c(thr_gam$min_fc, thr_classic$min_fc),
  min_reports = c(thr_gam$min_reports, thr_classic$min_reports),
  target_power = target_power,
  achieved_power = c(thr_gam$achieved_power, thr_classic$achieved_power),
  reached_target = c(thr_gam$reached, thr_classic$reached)
)

fwrite(thresholds_80, paste0(output_dir, "power_thresholds_80.csv"))

message("Umbrales (o mejor aproximación si no se alcanza el objetivo):")
print(thresholds_80)

################################################################################
# Generación de heatmaps de poder
################################################################################

message("\nGenerando heatmaps de poder...")

p_gam <- ggplot(power_matrix, aes(x = min_reports, y = fold_change, fill = power_gam)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("red", "yellow", "green"),
    limits = c(0, 1),
    na.value = "grey50",
    name = "Poder\nEstadístico"
  ) +
  geom_contour(aes(z = power_gam), breaks = target_power, color = "black", linewidth = 0.6) +
  labs(
    title = "Análisis de Poder Estadístico - Modelo GAM",
    subtitle = sprintf("Contorno negro: %.0f%% poder", target_power * 100),
    x = "Número mínimo de reportes",
    y = "Fold-change mínimo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

if (!is.na(thr_gam$min_fc) && !is.na(thr_gam$min_reports)) {
  p_gam <- p_gam +
    geom_point(
      data = power_matrix[fold_change == thr_gam$min_fc & min_reports == thr_gam$min_reports],
      size = 8, shape = 21, fill = "black", color = "white", stroke = 2
    )
}

ggsave(paste0(output_dir, "heatmap_power_gam.png"), p_gam, width = 12, height = 8, dpi = 300)

p_classic <- ggplot(power_matrix, aes(x = min_reports, y = fold_change, fill = power_classic)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("red", "yellow", "green"),
    limits = c(0, 1),
    na.value = "grey50",
    name = "Poder\nEstadístico"
  ) +
  geom_contour(aes(z = power_classic), breaks = target_power, color = "black", linewidth = 0.6) +
  labs(
    title = "Análisis de Poder Estadístico - IOR Clásico",
    subtitle = sprintf("Contorno negro: %.0f%% poder", target_power * 100),
    x = "Número mínimo de reportes",
    y = "Fold-change mínimo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

if (!is.na(thr_classic$min_fc) && !is.na(thr_classic$min_reports)) {
  p_classic <- p_classic +
    geom_point(
      data = power_matrix[fold_change == thr_classic$min_fc & min_reports == thr_classic$min_reports],
      size = 8, shape = 21, fill = "black", color = "white", stroke = 2
    )
}

ggsave(
  paste0(output_dir, "heatmap_power_classic.png"),
  p_classic,
  width = 12, height = 8, dpi = 300
)

################################################################################
# Filtrado del pool positivo según características de poder (GAM)
################################################################################

message("\nFiltrando pool positivo según características de poder...")

if (!is.na(thr_gam$min_fc) && !is.na(thr_gam$min_reports)) {
  filtered_positives <- positives_scores[
    model_success == TRUE &
      injection_success == TRUE &
      fold_change >= thr_gam$min_fc &
      n_coadmin >= thr_gam$min_reports &
      dynamic != "uniform"
  ]
} else {
  message("ADVERTENCIA: No se pudieron establecer umbrales para GAM")
  filtered_positives <- positives_scores[
    model_success == TRUE &
      injection_success == TRUE &
      dynamic != "uniform"
  ]
}

message(sprintf(
  "Pool positivo filtrado: %d tripletes (de %d originales)",
  nrow(filtered_positives), nrow(positives_scores)
))

saveRDS(filtered_positives, paste0(output_dir, "filtered_positive_pool.rds"))

filtered_meta <- filtered_positives[, .(
  triplet_id, drugA, drugB, meddra, dynamic, fold_change, n_coadmin, n_injected
)]

fwrite(filtered_meta, paste0(output_dir, "filtered_positive_metadata.csv"))

################################################################################
# Análisis de validación comparativa en pool filtrado
################################################################################

if (nrow(filtered_positives) > 0) {
  message("\nEjecutando análisis de validación comparativa en pool filtrado...")

  valid_negatives <- negatives_scores[model_success == TRUE]

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

  avg_threshold <- mean(null_thresholds$threshold_p95, na.rm = TRUE)

  pos_expanded_filtered[, signal_gam := is_signal_detected(log_ior_lower90, avg_threshold)]
  neg_expanded_filtered[, signal_gam := is_signal_detected(log_ior_lower90, avg_threshold)]

  pos_expanded_filtered[, signal_classic := is_signal_detected(log_ior_classic_lower90, 0)]
  neg_expanded_filtered[, signal_classic := is_signal_detected(log_ior_classic_lower90, 0)]

  pos_summary <- pos_expanded_filtered[, .(
    signal_gam = any(signal_gam),
    signal_classic = any(signal_classic)
  ), by = triplet_id]
  pos_summary[, true_label := 1L]

  neg_summary <- neg_expanded_filtered[, .(
    signal_gam = any(signal_gam),
    signal_classic = any(signal_classic)
  ), by = triplet_id]
  neg_summary[, true_label := 0L]

  validation_data <- rbindlist(list(pos_summary, neg_summary), use.names = TRUE, fill = TRUE)

  message("\nCalculando métricas de performance en pool filtrado...")

  calculate_validation_metrics <- function(data, method_col, method_name) {
    tp <- sum(data[[method_col]] == TRUE & data$true_label == 1, na.rm = TRUE)
    fn <- sum(data[[method_col]] == FALSE & data$true_label == 1, na.rm = TRUE)
    fp <- sum(data[[method_col]] == TRUE & data$true_label == 0, na.rm = TRUE)
    tn <- sum(data[[method_col]] == FALSE & data$true_label == 0, na.rm = TRUE)

    sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
    specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
    ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
    npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    f1_score <- ifelse(
      !is.na(ppv) && !is.na(sensitivity) && (ppv + sensitivity) > 0,
      2 * (ppv * sensitivity) / (ppv + sensitivity),
      NA_real_
    )

    data.table(
      method = method_name,
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv,
      accuracy = accuracy,
      f1_score = f1_score,
      tp = tp, fn = fn, fp = fp, tn = tn
    )
  }

  metrics_gam <- calculate_validation_metrics(validation_data, "signal_gam", "GAM")
  metrics_classic <- calculate_validation_metrics(validation_data, "signal_classic", "IOR_Clásico")

  metrics_combined <- rbindlist(list(metrics_gam, metrics_classic), use.names = TRUE, fill = TRUE)
  fwrite(metrics_combined, paste0(output_dir, "validation_metrics_filtered.csv"))

  message("Métricas de validación en pool filtrado:")
  print(metrics_combined)

  ################################################################################
  # Métricas por etapa
  ################################################################################

  message("\nCalculando métricas por etapa...")

  pos_expanded_with_stage <- filtered_positives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90),
      fold_change = fold_change[1],
      n_coadmin = n_coadmin[1]
    )
  }, by = triplet_id]

  pos_expanded_with_stage[, true_label := 1L]

  neg_expanded_with_stage <- valid_negatives[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(log_ior_lower90),
      log_ior_classic_lower90 = unlist(log_ior_classic_lower90),
      fold_change = NA_real_,
      n_coadmin = n_coadmin[1]
    )
  }, by = triplet_id]

  neg_expanded_with_stage[, true_label := 0L]

  all_expanded_filtered <- rbindlist(
    list(pos_expanded_with_stage, neg_expanded_with_stage),
    use.names = TRUE,
    fill = TRUE
  )

  all_expanded_filtered[, stage_name := niveles_nichd[stage]]

  all_expanded_filtered[, signal_gam_stage := is_signal_detected(log_ior_lower90, avg_threshold)]
  all_expanded_filtered[, signal_classic_stage := is_signal_detected(log_ior_classic_lower90, 0)]

  metrics_by_stage <- all_expanded_filtered[, {
    tp_gam <- sum(signal_gam_stage == TRUE & true_label == 1, na.rm = TRUE)
    fn_gam <- sum(signal_gam_stage == FALSE & true_label == 1, na.rm = TRUE)
    fp_gam <- sum(signal_gam_stage == TRUE & true_label == 0, na.rm = TRUE)
    tn_gam <- sum(signal_gam_stage == FALSE & true_label == 0, na.rm = TRUE)

    sens_gam <- ifelse((tp_gam + fn_gam) > 0, tp_gam / (tp_gam + fn_gam), NA_real_)
    spec_gam <- ifelse((tn_gam + fp_gam) > 0, tn_gam / (tn_gam + fp_gam), NA_real_)
    ppv_gam <- ifelse((tp_gam + fp_gam) > 0, tp_gam / (tp_gam + fp_gam), NA_real_)

    tp_classic <- sum(signal_classic_stage == TRUE & true_label == 1, na.rm = TRUE)
    fn_classic <- sum(signal_classic_stage == FALSE & true_label == 1, na.rm = TRUE)
    fp_classic <- sum(signal_classic_stage == TRUE & true_label == 0, na.rm = TRUE)
    tn_classic <- sum(signal_classic_stage == FALSE & true_label == 0, na.rm = TRUE)

    sens_classic <- ifelse(
      (tp_classic + fn_classic) > 0,
      tp_classic / (tp_classic + fn_classic),
      NA_real_
    )
    spec_classic <- ifelse(
      (tn_classic + fp_classic) > 0,
      tn_classic / (tn_classic + fp_classic),
      NA_real_
    )
    ppv_classic <- ifelse(
      (tp_classic + fp_classic) > 0,
      tp_classic / (tp_classic + fp_classic),
      NA_real_
    )

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

  fwrite(metrics_by_stage, paste0(output_dir, "validation_metrics_by_stage_filtered.csv"))
  message("Métricas por etapa en pool filtrado:")
  print(metrics_by_stage)

  ################################################################################
  # Métricas según número de reportes (solo positivos)
  ################################################################################

  message("\nCalculando métricas según número de reportes...")

  pos_expanded_with_stage[, signal_gam_stage := is_signal_detected(log_ior_lower90, avg_threshold)]
  pos_expanded_with_stage[, signal_classic_stage := is_signal_detected(log_ior_classic_lower90, 0)]

  pos_expanded_with_stage[, reports_category := cut(
    n_coadmin,
    breaks = c(0, 10, 20, 30, Inf),
    labels = c("5-10", "11-20", "21-30", "30+"),
    include.lowest = TRUE
  )]

  pos_triplet_reports <- pos_expanded_with_stage[, .(
    fold_change = fold_change[1],
    n_coadmin = n_coadmin[1],
    reports_category = reports_category[1],
    detected_gam = any(signal_gam_stage),
    detected_classic = any(signal_classic_stage)
  ), by = triplet_id]

  metrics_by_reports <- pos_triplet_reports[, .(
    n_triplets = .N,
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    sensitivity_gam = mean(detected_gam, na.rm = TRUE),
    sensitivity_classic = mean(detected_classic, na.rm = TRUE)
  ), by = reports_category]

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

cat(sprintf(
  "\nARCHIVOS GENERADOS:\n- Matriz de poder: %spower_matrix.csv\n- Umbrales: %spower_thresholds_80.csv\n- Heatmaps: %sheatmap_power_*.png\n- Pool filtrado: %sfiltered_positive_pool.rds\n- Métricas validación: %svalidation_metrics_*.csv\n\nUMBRALES (o mejor aproximación):\n- GAM: FC >= %.1f, Reportes >= %s (poder=%.3f, alcanzado=%s)\n- IOR Clásico: FC >= %.1f, Reportes >= %s (poder=%.3f, alcanzado=%s)\n\nPOOL FILTRADO (según umbrales GAM):\n- Tripletes positivos: %d (de %d originales)\n- Criterios: FC >= %.1f, Reportes >= %s, dinámico != uniform\n",
  output_dir,
  output_dir,
  output_dir,
  output_dir,
  output_dir,
  thr_gam$min_fc,
  as.character(thr_gam$min_reports),
  thr_gam$achieved_power,
  as.character(thr_gam$reached),
  thr_classic$min_fc,
  as.character(thr_classic$min_reports),
  thr_classic$achieved_power,
  as.character(thr_classic$reached),
  nrow(filtered_positives),
  nrow(positives_scores),
  thr_gam$min_fc,
  as.character(thr_gam$min_reports)
))

report <- data.table(
  parameter = c(
    "Fold-change mínimo GAM",
    "Reportes mínimo GAM",
    "Poder alcanzado GAM",
    "Alcanzó objetivo GAM",
    "Fold-change mínimo IOR",
    "Reportes mínimo IOR",
    "Poder alcanzado IOR",
    "Alcanzó objetivo IOR",
    "Tripletes originales",
    "Tripletes filtrados",
    "Poder objetivo"
  ),
  value = c(
    thr_gam$min_fc,
    thr_gam$min_reports,
    thr_gam$achieved_power,
    thr_gam$reached,
    thr_classic$min_fc,
    thr_classic$min_reports,
    thr_classic$achieved_power,
    thr_classic$reached,
    nrow(positives_scores),
    nrow(filtered_positives),
    target_power
  )
)

fwrite(report, paste0(output_dir, "power_analysis_summary.csv"))

message("\nAnálisis de poder completado exitosamente!")
message(sprintf("Resultados guardados en: %s", output_dir))
