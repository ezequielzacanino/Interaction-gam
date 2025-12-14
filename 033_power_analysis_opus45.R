################################################################################
# Script: 033_power_analysis_opus45.R
# Análisis de poder y validación comparativa: GAM vs IOR clásico
#
# Requisitos:
# - Ejecutar luego de 01_augmentation.R (y opcionalmente 02_nulldistribution_final.R)
# - Usa resultados semisintéticos (positivos inyectados + negativos) ya calculados.
#
# Objetivos:
# - Calcular potencia (sensibilidad en controles positivos) para combinaciones
#   de características mínimas: fold-change y número de reportes del triplete.
# - Encontrar umbrales (min fold-change + min reportes) con potencia >= 80%
#   para GAM y para IOR clásico.
# - Graficar potencia con heatmap.
# - Filtrar pool positivo usando umbrales encontrados y correr validación
#   comparativa (GAM vs IOR clásico).
# - Reportar métricas globales, por etapa NICHD y por número de reportes.
#
# Nota importante:
# - Si GAM o IOR clásico entregan NA o Inf, se considera NO detectado.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
})

if (file.exists("00_functions.R")) {
  source("00_functions.R", local = TRUE)
}

if (file.exists("giangreco_theme.R")) {
  source("giangreco_theme.R")
  if (exists("theme_giangreco")) theme_set(theme_giangreco())
} else {
  theme_set(theme_minimal(base_size = 12))
}

################################################################################
# Configuración
################################################################################

TARGET_POWER <- 0.80
MIN_TRIPLETS_FOR_POWER <- 25L

ruta_pos_results <- "./augmentation_results/positive_triplets_results.rds"
ruta_neg_results <- "./augmentation_results/negative_triplets_results.rds"
ruta_pos_meta <- "./augmentation_results/positive_triplets_metadata.csv"
ruta_neg_meta <- "./augmentation_results/negative_triplets_metadata.csv"

output_dir <- "./results/power_analysis_opus45/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Helpers
################################################################################

as_numeric_vec <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

signal_detected_from_lower90 <- function(lower90_vec) {
  v <- as_numeric_vec(lower90_vec)
  if (length(v) == 0) return(FALSE)
  any(is.finite(v) & v > 0)
}

safe_max_finite <- function(x) {
  v <- as_numeric_vec(x)
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  max(v)
}

calc_binary_metrics <- function(dt, true_col, pred_col) {
  true <- dt[[true_col]]
  pred <- dt[[pred_col]]

  tp <- sum(true == 1 & pred == TRUE, na.rm = TRUE)
  fn <- sum(true == 1 & pred == FALSE, na.rm = TRUE)
  fp <- sum(true == 0 & pred == TRUE, na.rm = TRUE)
  tn <- sum(true == 0 & pred == FALSE, na.rm = TRUE)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  f1 <- if (is.finite(ppv) && is.finite(sensitivity) && (ppv + sensitivity) > 0) {
    2 * (ppv * sensitivity) / (ppv + sensitivity)
  } else {
    NA_real_
  }

  data.table(
    tp = tp, fn = fn, fp = fp, tn = tn,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy,
    f1_score = f1
  )
}

safe_auroc <- function(true, score) {
  ok <- is.finite(score) & !is.na(true)
  if (sum(ok) < 10 || length(unique(true[ok])) < 2) return(list(auc = NA_real_))

  roc_obj <- tryCatch({
    roc(response = true[ok], predictor = score[ok],
        direction = "<", levels = c(0, 1), quiet = TRUE)
  }, error = function(e) NULL)

  if (is.null(roc_obj)) return(list(auc = NA_real_))
  list(auc = as.numeric(auc(roc_obj)))
}

pick_threshold_80 <- function(grid_dt, power_col, target_power = 0.8, min_n = 25L) {
  candidates <- grid_dt[
    get(power_col) >= target_power & n_eligible >= min_n & !is.na(get(power_col))
  ]

  if (nrow(candidates) == 0) {
    return(data.table(
      method = power_col,
      min_reports = NA_integer_,
      min_fold_change = NA_real_,
      n_eligible = 0L,
      power = NA_real_
    ))
  }

  # criterio: maximizar n_eligible (menos restrictivo), luego minimizar reportes y fold-change
  setorder(candidates, -n_eligible, min_reports, min_fold_change)
  best <- candidates[1]

  data.table(
    method = power_col,
    min_reports = best$min_reports,
    min_fold_change = best$min_fold_change,
    n_eligible = best$n_eligible,
    power = best[[power_col]]
  )
}

################################################################################
# Carga y preparación de datos
################################################################################

if (!file.exists(ruta_pos_results)) stop("No existe: ", ruta_pos_results)
if (!file.exists(ruta_neg_results)) stop("No existe: ", ruta_neg_results)
if (!file.exists(ruta_pos_meta)) stop("No existe: ", ruta_pos_meta)
if (!file.exists(ruta_neg_meta)) stop("No existe: ", ruta_neg_meta)

pos_scores <- as.data.table(readRDS(ruta_pos_results))
neg_scores <- as.data.table(readRDS(ruta_neg_results))

pos_meta <- fread(ruta_pos_meta)
neg_meta <- fread(ruta_neg_meta)

# Aseguro nombres consistentes (compatibilidad con versiones viejas de data.table)
if ("N" %in% names(pos_meta) && !("n_reports" %in% names(pos_meta))) {
  setnames(pos_meta, "N", "n_reports")
}
if ("N" %in% names(neg_meta) && !("n_reports" %in% names(neg_meta))) {
  setnames(neg_meta, "N", "n_reports")
}

# Merge para recuperar número de reportes del triplete (no usar N de pos_scores, porque allí es n_events)
pos_dt <- merge(
  pos_scores,
  pos_meta[, .(triplet_id, n_reports, fold_change, dynamic)],
  by = "triplet_id",
  all.x = TRUE
)

neg_dt <- merge(
  neg_scores,
  neg_meta[, .(triplet_id, n_reports)],
  by = "triplet_id",
  all.x = TRUE
)

# Filtrado de positivos válidos (se excluye uniform)
pos_dt <- pos_dt[
  injection_success == TRUE &
    !is.na(dynamic) &
    dynamic != "uniform" &
    !is.na(n_reports) &
    !is.na(fold_change)
]

# Filtrado de negativos válidos
neg_dt <- neg_dt[model_success == TRUE & !is.na(n_reports)]

message(sprintf("Positivos (inyectados y no-uniform): %d", nrow(pos_dt)))
message(sprintf("Negativos (model_success): %d", nrow(neg_dt)))

# Verificación de columnas para IOR clásico
has_classic <- all(c("classic_success", "log_ior_classic_lower90", "ior_classic") %in% names(pos_scores))
if (!has_classic) {
  warning("No se encontraron columnas de IOR clásico en los RDS. ",
          "Se calcula potencia y validación solo para GAM.")
}

################################################################################
# Detección a nivel triplete (GAM y clásico)
################################################################################

pos_dt[, gam_detected := fifelse(
  model_success == TRUE,
  sapply(log_ior_lower90, signal_detected_from_lower90),
  FALSE
)]

pos_dt[, gam_score := fifelse(
  model_success == TRUE,
  sapply(ior_values, safe_max_finite),
  NA_real_
)]

neg_dt[, gam_detected := fifelse(
  model_success == TRUE,
  sapply(log_ior_lower90, signal_detected_from_lower90),
  FALSE
)]

neg_dt[, gam_score := fifelse(
  model_success == TRUE,
  sapply(ior_values, safe_max_finite),
  NA_real_
)]

if (has_classic) {
  pos_dt[, classic_detected := fifelse(
    classic_success == TRUE,
    sapply(log_ior_classic_lower90, signal_detected_from_lower90),
    FALSE
  )]

  pos_dt[, classic_score := fifelse(
    classic_success == TRUE,
    sapply(ior_classic, safe_max_finite),
    NA_real_
  )]

  neg_dt[, classic_detected := fifelse(
    classic_success == TRUE,
    sapply(log_ior_classic_lower90, signal_detected_from_lower90),
    FALSE
  )]

  neg_dt[, classic_score := fifelse(
    classic_success == TRUE,
    sapply(ior_classic, safe_max_finite),
    NA_real_
  )]
} else {
  pos_dt[, `:=`(classic_detected = FALSE, classic_score = NA_real_)]
  neg_dt[, `:=`(classic_detected = FALSE, classic_score = NA_real_)]
}

################################################################################
# 1) Potencia por combinación mínima (fold-change, reportes)
################################################################################

# Construyo grillas compactas (basadas en cuantiles) para evitar explosionar combinaciones
report_grid <- unique(as.integer(round(
  quantile(pos_dt$n_reports, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
)))
report_grid <- sort(unique(pmax(report_grid, min(report_grid, na.rm = TRUE))))

fc_grid <- unique(round(
  quantile(pos_dt$fold_change, probs = seq(0, 1, by = 0.05), na.rm = TRUE),
  digits = 2
))
fc_grid <- sort(unique(pmax(fc_grid, 1)))

# Si quedaron demasiados puntos, reduzco
if (length(report_grid) > 40) {
  report_grid <- sort(unique(as.integer(round(quantile(pos_dt$n_reports, probs = seq(0, 1, length.out = 30), na.rm = TRUE)))))
}
if (length(fc_grid) > 40) {
  fc_grid <- sort(unique(round(quantile(pos_dt$fold_change, probs = seq(0, 1, length.out = 30), na.rm = TRUE), 2)))
}

power_grid <- CJ(
  min_reports = report_grid,
  min_fold_change = fc_grid
)

power_grid[, `:=`(n_eligible = 0L, power_gam = NA_real_, power_classic = NA_real_)]

for (i in seq_len(nrow(power_grid))) {
  mr <- power_grid$min_reports[i]
  mfc <- power_grid$min_fold_change[i]

  eligible <- pos_dt[n_reports >= mr & fold_change >= mfc]
  power_grid$n_eligible[i] <- nrow(eligible)

  if (nrow(eligible) > 0) {
    power_grid$power_gam[i] <- mean(eligible$gam_detected, na.rm = TRUE)
    power_grid$power_classic[i] <- mean(eligible$classic_detected, na.rm = TRUE)
  }
}

fwrite(power_grid, file.path(output_dir, "power_grid.csv"))

################################################################################
# 2) Umbrales (80%) por método
################################################################################

th_gam <- pick_threshold_80(power_grid, "power_gam", TARGET_POWER, MIN_TRIPLETS_FOR_POWER)

if (has_classic) {
  th_classic <- pick_threshold_80(power_grid, "power_classic", TARGET_POWER, MIN_TRIPLETS_FOR_POWER)
} else {
  th_classic <- data.table(
    method = "power_classic",
    min_reports = NA_integer_,
    min_fold_change = NA_real_,
    n_eligible = 0L,
    power = NA_real_
  )
}

thresholds <- rbind(
  th_gam[, .(method = "GAM", min_reports, min_fold_change, n_eligible, power)],
  if (has_classic) th_classic[, .(method = "IOR_classic", min_reports, min_fold_change, n_eligible, power)]
)

# Umbral conjunto para filtrar un pool único comparable (intersección de requisitos)
# Si no hay clásico, BOTH = GAM.
both_min_reports <- th_gam$min_reports
both_min_fc <- th_gam$min_fold_change

if (has_classic && !is.na(th_classic$min_reports) && !is.na(th_classic$min_fold_change)) {
  both_min_reports <- max(th_gam$min_reports, th_classic$min_reports)
  both_min_fc <- max(th_gam$min_fold_change, th_classic$min_fold_change)
}

th_both <- data.table(
  method = "BOTH",
  min_reports = both_min_reports,
  min_fold_change = both_min_fc
)

thresholds <- rbind(
  thresholds,
  th_both[, .(method, min_reports, min_fold_change, n_eligible = NA_integer_, power = NA_real_)]
)

fwrite(thresholds, file.path(output_dir, "power_thresholds_80pct.csv"))

################################################################################
# 3) Heatmap de potencia
################################################################################

power_measure_vars <- if (has_classic) c("power_gam", "power_classic") else "power_gam"

power_long <- melt(
  power_grid,
  id.vars = c("min_reports", "min_fold_change", "n_eligible"),
  measure.vars = power_measure_vars,
  variable.name = "method",
  value.name = "power"
)

power_long[, method := fifelse(method == "power_gam", "GAM", "IOR clásico")]

p_heatmap <- ggplot(power_long[n_eligible >= MIN_TRIPLETS_FOR_POWER],
                    aes(x = min_fold_change, y = min_reports, fill = power)) +
  geom_tile() +
  geom_contour(aes(z = power), breaks = TARGET_POWER, color = "white", linewidth = 0.6) +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
  facet_wrap(~method) +
  labs(
    title = "Análisis de poder (sensibilidad) según características mínimas",
    subtitle = sprintf("Señal detectada si existe al menos una etapa con IC90 log(IOR) > 0 y finito | Target = %.0f%%", 100 * TARGET_POWER),
    x = "Fold-change mínimo",
    y = "Número mínimo de reportes del triplete",
    fill = "Potencia"
  )

ggsave(
  filename = file.path(output_dir, "fig_power_heatmap.png"),
  plot = p_heatmap,
  width = 12,
  height = 6,
  dpi = 300
)

################################################################################
# 4) Filtrado de pool positivo según umbral conjunto
################################################################################

if (is.na(th_both$min_reports) || is.na(th_both$min_fold_change)) {
  stop("No se encontró un umbral conjunto BOTH. Revisar power_thresholds_80pct.csv")
}

pos_filt <- pos_dt[n_reports >= th_both$min_reports & fold_change >= th_both$min_fold_change]
neg_filt <- neg_dt[n_reports >= th_both$min_reports]

message(sprintf("\nPool filtrado (BOTH): Positivos=%d | Negativos=%d\n", nrow(pos_filt), nrow(neg_filt)))

fwrite(pos_filt, file.path(output_dir, "positives_filtered.csv"))
fwrite(neg_filt, file.path(output_dir, "negatives_filtered.csv"))

################################################################################
# 5) Validación comparativa en pool filtrado (métricas globales)
################################################################################

triplets_eval <- rbind(
  pos_filt[, .(
    triplet_id,
    true_label = 1L,
    n_reports,
    fold_change,
    gam_detected,
    classic_detected,
    gam_score,
    classic_score
  )],
  neg_filt[, .(
    triplet_id,
    true_label = 0L,
    n_reports,
    fold_change = NA_real_,
    gam_detected,
    classic_detected,
    gam_score,
    classic_score
  )],
  fill = TRUE
)

metrics_global <- rbind(
  cbind(method = "GAM", calc_binary_metrics(triplets_eval, "true_label", "gam_detected")),
  if (has_classic) cbind(method = "IOR_classic", calc_binary_metrics(triplets_eval, "true_label", "classic_detected"))
)

# AUROC (score continuo, max IOR por triplete)
auc_gam <- safe_auroc(triplets_eval$true_label, triplets_eval$gam_score)$auc
metrics_global[method == "GAM", auroc := auc_gam]

if (has_classic) {
  auc_classic <- safe_auroc(triplets_eval$true_label, triplets_eval$classic_score)$auc
  metrics_global[method == "IOR_classic", auroc := auc_classic]
}

fwrite(metrics_global, file.path(output_dir, "metrics_global_triplet.csv"))

################################################################################
# 6) Métricas por etapa (stage-level)
################################################################################

# Expansión a nivel etapa
expand_stage_dt <- function(dt, label_value, lower90_col, score_col) {
  dt[, {
    lower90 <- as_numeric_vec(get(lower90_col)[[1]])
    score <- as_numeric_vec(get(score_col)[[1]])

    # Relleno/recorte a 7 etapas
    lower90 <- c(lower90, rep(NA_real_, 7 - length(lower90)))[1:7]
    score <- c(score, rep(NA_real_, 7 - length(score)))[1:7]

    data.table(
      stage = 1:7,
      true_label = label_value,
      n_reports = n_reports,
      lower90 = lower90,
      score = score,
      detected = is.finite(lower90) & lower90 > 0
    )
  }, by = .(triplet_id)]
}

stage_gam <- rbind(
  expand_stage_dt(pos_filt, 1L, "log_ior_lower90", "ior_values"),
  expand_stage_dt(neg_filt, 0L, "log_ior_lower90", "ior_values"),
  fill = TRUE
)

stage_gam[, stage_name := if (exists("niveles_nichd")) niveles_nichd[stage] else as.character(stage)]

metrics_by_stage_gam <- stage_gam[, {
  m <- calc_binary_metrics(.SD, "true_label", "detected")
  auc <- safe_auroc(true_label, score)$auc
  cbind(m, auroc = auc, n_total = .N)
}, by = .(stage, stage_name)]
metrics_by_stage_gam[, method := "GAM"]

if (has_classic) {
  stage_classic <- rbind(
    expand_stage_dt(pos_filt, 1L, "log_ior_classic_lower90", "ior_classic"),
    expand_stage_dt(neg_filt, 0L, "log_ior_classic_lower90", "ior_classic"),
    fill = TRUE
  )

  stage_classic[, stage_name := if (exists("niveles_nichd")) niveles_nichd[stage] else as.character(stage)]

  metrics_by_stage_classic <- stage_classic[, {
    m <- calc_binary_metrics(.SD, "true_label", "detected")
    auc <- safe_auroc(true_label, score)$auc
    cbind(m, auroc = auc, n_total = .N)
  }, by = .(stage, stage_name)]
  metrics_by_stage_classic[, method := "IOR_classic"]
} else {
  metrics_by_stage_classic <- data.table()
}

metrics_by_stage <- rbind(metrics_by_stage_gam, metrics_by_stage_classic, fill = TRUE)
setcolorder(metrics_by_stage, c("method", "stage", "stage_name", "n_total", "tp", "fn", "fp", "tn",
                               "sensitivity", "specificity", "ppv", "npv", "accuracy", "f1_score", "auroc"))

fwrite(metrics_by_stage, file.path(output_dir, "metrics_by_stage.csv"))

################################################################################
# 7) Métricas por número de reportes (triplet-level)
################################################################################

# Bins por cantidad de reportes (ajustable)
report_breaks <- c(th_both$min_reports, 15, 20, 30, 50, 75, 100, 200, Inf)
report_breaks <- sort(unique(report_breaks))

triplets_eval[, report_bin := cut(
  n_reports,
  breaks = report_breaks,
  right = FALSE,
  include.lowest = TRUE
)]

metrics_by_reports <- rbind(
  triplets_eval[, {
    m <- calc_binary_metrics(.SD, "true_label", "gam_detected")
    auc <- safe_auroc(true_label, gam_score)$auc
    cbind(method = "GAM", m, auroc = auc, n_total = .N)
  }, by = .(report_bin)],
  if (has_classic) triplets_eval[, {
    m <- calc_binary_metrics(.SD, "true_label", "classic_detected")
    auc <- safe_auroc(true_label, classic_score)$auc
    cbind(method = "IOR_classic", m, auroc = auc, n_total = .N)
  }, by = .(report_bin)],
  fill = TRUE
)

fwrite(metrics_by_reports, file.path(output_dir, "metrics_by_report_bin.csv"))

################################################################################
# Resumen textual
################################################################################

summary_txt <- paste0(
  "===============================================================================\n",
  "ANÁLISIS DE PODER + VALIDACIÓN COMPARATIVA (Opus 4.5)\n",
  "===============================================================================\n\n",
  "Target de potencia: ", sprintf("%.0f%%", 100 * TARGET_POWER), "\n",
  "Mínimo tripletes para estimar potencia: ", MIN_TRIPLETS_FOR_POWER, "\n\n",
  "Umbrales (80%)\n",
  capture.output(print(thresholds)),
  "\n\n",
  "Pool filtrado (BOTH):\n",
  "  min_reports = ", th_both$min_reports, "\n",
  "  min_fold_change = ", th_both$min_fold_change, "\n",
  "  positivos = ", nrow(pos_filt), "\n",
  "  negativos = ", nrow(neg_filt), "\n\n",
  "Métricas globales (triplet-level):\n",
  capture.output(print(metrics_global)),
  "\n\n",
  "Outputs guardados en: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n",
  "===============================================================================\n"
)

writeLines(summary_txt, file.path(output_dir, "EXECUTIVE_SUMMARY.txt"))
cat(summary_txt)
