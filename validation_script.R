library(data.table)
library(tidyverse)
library(mgcv)
library(parallel)
library(pROC)
library(pbapply)


# ============================================================================
# SCRIPT DE VALIDACIÓN
# ============================================================================


# ---------- PARÁMETROS DE CONFIGURACIÓN ----------
ruta_ade_augmented <- "./ade_augmented.csv"
ruta_ground_truth <- "./ground_truth_positive.csv"
ruta_positive_meta <- "./positive_triplets_metadata.csv" 
ruta_negative_meta <- "./negative_triplets_metadata.csv" #adaptar para cuando le saque la palabra "triplets"

ruta_ground_truth_negative <- "./ground_truth_negative.csv"
n_replicas <- 500
k_triplets_sample <- 300
prioritize_by_count <- TRUE
seed_base <- 2025
min_reports_triplet <- 10
save_null_distribution_file <- "null_distribution_permuted.csv"
save_null_thresholds_file <- "null_thresholds_pc99.csv"
save_empirical_pvalues_file <- "empirical_pvalues_groundtruth.csv"


alpha_nominal <- 0.10
z90 <- qnorm(0.95)
n_cores <- max(1, detectCores() - 1)

# ---------- CARGAR POOL DE REPORTES NEGATIVOS ----------
ground_truth_negative <- fread(ruta_ground_truth_negative)
neg_report_ids <- unique(ground_truth_negative$selected_neg_reports)

message("Pool de reportes negativos: ", length(neg_report_ids))

# ---------- FUNCIONES AUXILIARES ----------

fit_differential_gam <- function(drugA_id, drugB_id, event_id, ade_data) {
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  datos_modelo <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  if (sum(datos_modelo$ea_ocurrio) < 5 || sum(datos_modelo$droga_ab) < 3) {
    return(list(success = FALSE, n_events = sum(datos_modelo$ea_ocurrio), 
                n_coadmin = sum(datos_modelo$droga_ab)))
  }
  
  tryCatch({
    modelo <- bam(
      ea_ocurrio ~ droga_a + droga_b +
        s(nichd_num, k = 7, bs = "cs") +
        s(nichd_num, bs = "cs", by = droga_ab, k = 7),
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = "fREML",
      discrete = TRUE
    )
    
    return(list(success = TRUE, model = modelo, data = datos_modelo,
                n_events = sum(datos_modelo$ea_ocurrio),
                n_coadmin = sum(datos_modelo$droga_ab)))
    
  }, error = function(e) {
    return(list(success = FALSE, error_msg = e$message,
                n_events = sum(datos_modelo$ea_ocurrio),
                n_coadmin = sum(datos_modelo$droga_ab)))
  })
}


extract_interaction_coefs <- function(model_result) {
  if (!model_result$success) return(NULL)
  
  modelo <- model_result$model
  
  grid_ab <- data.table(nichd_num = 1:7, droga_a = 1, droga_b = 1, droga_ab = 1)
  grid_00 <- data.table(nichd_num = 1:7, droga_a = 0, droga_b = 0, droga_ab = 0)
  
  pred_ab <- predict(modelo, newdata = grid_ab, type = "link", se.fit = TRUE)
  pred_00 <- predict(modelo, newdata = grid_00, type = "link", se.fit = TRUE)
  
  coef_by_stage <- pred_ab$fit - pred_00$fit
  se_by_stage <- sqrt(pred_ab$se.fit^2 + pred_00$se.fit^2)
  
  return(list(
    coef = coef_by_stage,
    se = se_by_stage,
    lower90 = coef_by_stage - z90 * se_by_stage,
    upper90 = coef_by_stage + z90 * se_by_stage
  ))
}


calculate_interaction_signal_giangreco <- function(model_result, null_thresholds) {
  if (!model_result$success) {
    return(list(
      success = FALSE, signal_detected = FALSE, n_stages_significant = 0,
      max_ior = NA, ior_values = rep(NA, 7), ior_li = rep(NA, 7), 
      ior_ls = rep(NA, 7), log_ior_lower90 = rep(NA, 7)
    ))
  }
  
  modelo <- model_result$model
  
  grid_dif <- CJ(nichd_num = 1:7, droga_a = c(0, 1), droga_b = c(0, 1))
  grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
  grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
  grid_dif[, `:=`(li = lp - z90 * se, ls = lp + z90 * se)]
  
  w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                value.var = c("lp", "se", "li", "ls"))
  
  log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
  
  Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
  Vb <- vcov(modelo, unconditional = TRUE)
  cov_link <- Xp %*% Vb %*% t(Xp)
  
  log_ior_se <- numeric(7)
  for (stage in 1:7) {
    idx_00 <- which(grid_dif$nichd_num == stage & grid_dif$droga_a == 0 & grid_dif$droga_b == 0)
    idx_01 <- which(grid_dif$nichd_num == stage & grid_dif$droga_a == 0 & grid_dif$droga_b == 1)
    idx_10 <- which(grid_dif$nichd_num == stage & grid_dif$droga_a == 1 & grid_dif$droga_b == 0)
    idx_11 <- which(grid_dif$nichd_num == stage & grid_dif$droga_a == 1 & grid_dif$droga_b == 1)
    
    cvec <- rep(0, nrow(grid_dif))
    cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)
    
    log_ior_se[stage] <- sqrt(max(as.numeric(t(cvec) %*% cov_link %*% cvec), 0))
  }
  
  log_ior_lower90 <- log_ior - z90 * log_ior_se
  log_ior_upper90 <- log_ior + z90 * log_ior_se
  
  ior_values <- exp(log_ior)
  ior_li <- exp(log_ior_lower90)
  ior_ls <- exp(log_ior_upper90)
  
  # Criterio Giangreco: nominal + null model
  nominal_significant <- log_ior_lower90 > 0
  nullmodel_significant <- sapply(1:7, function(i) {
    log_ior_lower90[i] > null_thresholds[i]
  })
  
  signal_detected <- any(nullmodel_significant)
  
  return(list(
    success = TRUE,
    signal_detected = signal_detected,
    n_stages_significant = sum(nullmodel_significant),
    max_ior = max(ior_values, na.rm = TRUE),
    mean_ior = mean(ior_values, na.rm = TRUE),
    ior_values = ior_values,
    ior_li = ior_li,
    ior_ls = ior_ls,
    log_ior = log_ior,
    log_ior_lower90 = log_ior_lower90,
    log_ior_se = log_ior_se,
    nullmodel_by_stage = nullmodel_significant,
    model_aic = AIC(modelo),
    model_deviance = deviance(modelo)
  ))
}


analyze_triplet <- function(i, triplet_row, ade_data, null_thresh) {
  drugA <- triplet_row$drugA
  drugB <- triplet_row$drugB
  event <- triplet_row$meddra
  
  model_result <- fit_differential_gam(drugA, drugB, event, ade_data)
  signal_result <- calculate_interaction_signal_giangreco(model_result, null_thresh)
  
  result <- list(
    triplet_id = i, drugA = drugA, drugB = drugB, meddra = event,
    type = triplet_row$type, model_success = model_result$success,
    n_events = model_result$n_events, n_coadmin = model_result$n_coadmin,
    signal_detected = signal_result$signal_detected,
    n_stages_significant = signal_result$n_stages_significant,
    max_ior = signal_result$max_ior,
    mean_ior = if (signal_result$success) signal_result$mean_ior else NA,
    model_aic = if (signal_result$success) signal_result$model_aic else NA,
    ior_values = if (signal_result$success) list(signal_result$ior_values) else list(rep(NA, 7)),
    ior_li = if (signal_result$success) list(signal_result$ior_li) else list(rep(NA, 7)),
    ior_ls = if (signal_result$success) list(signal_result$ior_ls) else list(rep(NA, 7)),
    log_ior_lower90 = if (signal_result$success) list(signal_result$log_ior_lower90) else list(rep(NA, 7))
  )
  
  return(result)
}

permute_events_within_strata <- function(pool_reports_meta, niveles_nichd, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pool <- copy(pool_reports_meta)
  pool[, events_perm := vector("list", .N)]
  for (stage in niveles_nichd) {
    idx <- which(pool$nichd == stage)
    if (length(idx) <= 1) {
      pool$events_perm[idx] <- pool$events[idx]
      next
    }
    perm_idx <- sample(seq_along(idx), length(idx), replace = FALSE)
    pool$events_perm[idx] <- pool$events[idx[perm_idx]]
  }
  pool[, .(safetyreportid, nichd, events_perm)]
}

make_triplets_per_report <- function(dr, ev, rid, nichd_stage) {
  if (length(dr) < 2 || length(ev) < 1) return(NULL)
  dr <- unique(dr); ev <- unique(ev)
  if (length(dr) == 2) combs <- matrix(dr, nrow = 1) else combs <- t(combn(dr, 2))
  n_combs <- nrow(combs); n_events <- length(ev)
  data.table(
    safetyreportid = rid,
    drugA = rep(combs[,1], times = n_events),
    drugB = rep(combs[,2], times = n_events),
    meddra = rep(ev, each = n_combs),
    nichd = nichd_stage
  )
}

build_triplets_from_permuted <- function(permuted_pool_dt, reports_meta_full) {
  tmp <- merge(permuted_pool_dt, reports_meta_full[, .(safetyreportid, drugs)], by = "safetyreportid", all.x = TRUE)
  trip_list <- vector("list", nrow(tmp))
  for (i in seq_len(nrow(tmp))) {
    rowi <- tmp[i]
    trip_list[[i]] <- make_triplets_per_report(rowi$drugs[[1]], rowi$events_perm[[1]], rowi$safetyreportid, rowi$nichd)
  }
  rbindlist(Filter(Negate(is.null), trip_list), use.names = TRUE)
}

# ---------- CARGAR DATOS ----------

ade_aug <- fread(ruta_ade_augmented)
ground_truth <- fread(ruta_ground_truth)
pos_meta <- fread(ruta_positive_meta)
neg_meta <- fread(ruta_negative_meta)


niveles_nichd <- c("term_neonatal","infancy","toddler","early_childhood",
                   "middle_childhood","early_adolescence","late_adolescence")
ade_aug[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_aug[, nichd_num := as.integer(nichd)]


message(sprintf("Dataset: %d filas | Ground truth: %d tripletes (%d pos, %d neg)",
                nrow(ade_aug), nrow(ground_truth),
                sum(ground_truth$type == "positive"),
                sum(ground_truth_negative$type == "negative")))

# ---------- CONSTRUIR LISTAS DE DROGAS Y EVENTOS POR REPORTE ----------

drugs_by_report <- unique(ade_aug[!is.na(atc_concept_id), .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_aug[!is.na(meddra_concept_id), .(safetyreportid, meddra_concept_id)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), by = safetyreportid]

reports_meta <- unique(ade_aug[, .(safetyreportid, nichd)])
reports_meta <- merge(reports_meta, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
reports_meta[is.na(drugs), drugs := list(integer(0))]
reports_meta[is.na(events), events := list(integer(0))]

pool_reports_meta <- reports_meta[safetyreportid %in% neg_report_ids]

# ---------- PARALELIZACIÓN POR RÉPLICA ----------
cl <- makeCluster(n_cores)
clusterExport(cl, c("pool_reports_meta", "reports_meta", "niveles_nichd",
                    "min_reports_triplet", "fit_differential_gam",
                    "calculate_interaction_signal_giangreco",
                    "permute_events_within_strata", "make_triplets_per_report",
                    "build_triplets_from_permuted", "prioritize_by_count",
                    "k_triplets_sample", "ade_aug", "z90"),
              envir = environment())
clusterEvalQ(cl, { library(data.table); library(mgcv) })

message(sprintf("Ejecutando %d réplicas permutadas en %d núcleos...", n_replicas, n_cores))

replica_results <- pblapply(seq_len(n_replicas), function(rep) {
  set.seed(2025 + rep)
  
  # 1. Permutar eventos dentro de cada nichd
  permuted_pool <- permute_events_within_strata(pool_reports_meta, niveles_nichd)
  
  # 2. Construir tripletes
  triplets_perm <- build_triplets_from_permuted(permuted_pool, reports_meta)
  if (nrow(triplets_perm) == 0) return(NULL)
  
  # 3. Contar y filtrar
  trip_counts_rep <- unique(triplets_perm[, .(drugA, drugB, meddra, safetyreportid)])[ , .N, by = .(drugA, drugB, meddra)]
  trip_counts_rep <- trip_counts_rep[N >= min_reports_triplet]
  if (nrow(trip_counts_rep) == 0) return(NULL)
  
  # 4. Submuestreo
  if (prioritize_by_count) {
    sel_k <- head(trip_counts_rep[order(-N)], k_triplets_sample)
  } else {
    sel_k <- trip_counts_rep[sample(.N, min(.N, k_triplets_sample))]
  }
  
  # 5. Ajuste GAM por triplete (paralelo dentro de réplica)
  trip_results <- lapply(seq_len(nrow(sel_k)), function(ti) {
    rowt <- sel_k[ti]
    model_res <- fit_differential_gam(rowt$drugA, rowt$drugB, rowt$meddra, ade_aug)
    if (!model_res$success) return(NULL)
    sig_res <- calculate_interaction_signal_giangreco(model_res, null_thresholds = rep(-Inf,7))
    if (!sig_res$success) return(NULL)
    data.table(
      drugA = rowt$drugA, drugB = rowt$drugB, meddra = rowt$meddra,
      stage = 1:7, log_lower90 = sig_res$log_ior_lower90
    )
  })  
  
  valid <- rbindlist(Filter(Negate(is.null), trip_results), fill = TRUE)
  if (nrow(valid) == 0) return(NULL)
  valid[, replica := rep]
  valid
}, cl = cl)

stopCluster(cl)

# ---------- CONSOLIDAR DISTRIBUCIÓN NULA ----------
null_all <- rbindlist(Filter(Negate(is.null), replica_results), fill = TRUE)
fwrite(null_all, save_null_distribution_file)
message("✓ Distribución nula guardada en ", save_null_distribution_file)

null_thresholds <- null_all[, .(threshold_p99 = quantile(log_lower90, 0.99, na.rm = TRUE)), by = stage]
fwrite(null_thresholds, save_null_thresholds_file)
message("✓ Umbrales Pc99 guardados en ", save_null_thresholds_file)

print(null_thresholds)

# ---------- CARGAR PARA NO VOLVER A CORRER ----------

# cargar null_all para no volver a ajustar toda la distribución nula
null_all <- fread("null_distribution_permuted.csv")
null_thresholds <- fread("null_thresholds_pc99.csv")

# ---------- CALCULAR P-VALORES EMPÍRICOS PARA GROUND TRUTH OBSERVADO ----------
#
# ESTO HACE LO MISMO QUE LA SECCIÓN ANÁLISIS DE TRIPLETES
#

message("\nCalculando p-valores empíricos para ground_truth...")

# Ajustar modelo real si no lo tenés ya ajustado
observed_results <- pblapply(seq_len(nrow(ground_truth)), function(i) {
  rowt <- ground_truth[i]
  model_res <- fit_differential_gam(rowt$drugA, rowt$drugB, rowt$meddra, ade_aug)
  if (!model_res$success) return(NULL)
  sig_res <- calculate_interaction_signal_giangreco(model_res, null_thresholds = rep(-Inf,7))
  if (!sig_res$success) return(NULL)
  data.table(
    triplet_id = i,
    drugA = rowt$drugA, drugB = rowt$drugB, meddra = rowt$meddra, type = rowt$type,
    stage = 1:7, log_lower90_obs = sig_res$log_ior_lower90
  )
})

observed_dt <- rbindlist(Filter(Negate(is.null), observed_results), fill = TRUE)

# Empirical p = mean(null >= observed)
pval_dt <- merge(observed_dt, null_all[, .(log_lower90_null = log_lower90), by = stage], by = "stage", allow.cartesian = TRUE)
pval_dt[, p_emp := mean(log_lower90_null >= log_lower90_obs, na.rm = TRUE), by = .(triplet_id, stage)]
pval_final <- unique(pval_dt[, .(triplet_id, drugA, drugB, meddra, type, stage, log_lower90_obs, p_emp)])

# Guardar p-values y marcar significancia (Pc99)
pval_final <- merge(pval_final, null_thresholds, by = "stage", all.x = TRUE)
pval_final[, significant_pc99 := log_lower90_obs > threshold_p99]
fwrite(pval_final, save_empirical_pvalues_file)
message("✓ P-valores empíricos guardados en ", save_empirical_pvalues_file)

# Resumen global
message("\nResumen de etapas significativas (Pc99):")
print(pval_final[, .(n_sig = sum(significant_pc99)), by = type])

# ---------- ANÁLISIS DE TRIPLETES ----------
cl <- makeCluster(n_cores)
clusterExport(cl, c("fit_differential_gam", "calculate_interaction_signal_giangreco",
                    "analyze_triplet", "z90", "niveles_nichd", "ade_aug", 
                    "ground_truth", "null_thresholds_by_stage"))
clusterEvalQ(cl, { library(data.table); library(mgcv) })


results_list <- pblapply(seq_len(nrow(ground_truth)), function(i) {
  analyze_triplet(i, ground_truth[i], ade_aug, null_thresholds_by_stage)
}, cl = cl)


stopCluster(cl)


results_dt <- rbindlist(lapply(results_list, function(x) {
  x$ior_values <- NULL; x$ior_li <- NULL; x$ior_ls <- NULL; x$log_ior_lower90 <- NULL
  as.data.table(x)
}))


ior_matrix <- do.call(rbind, lapply(results_list, function(x) x$ior_values[[1]]))
colnames(ior_matrix) <- paste0("ior_stage_", 1:7)
results_dt <- cbind(results_dt, ior_matrix)


# ---------- MÉTRICAS DE DESEMPEÑO (ESTILO GIANGRECO) ----------
positivos_res <- results_dt[type == "positive"]
negativos_res <- results_dt[type == "negative"]


tp <- sum(positivos_res$signal_detected, na.rm = TRUE)
fn <- sum(!positivos_res$signal_detected, na.rm = TRUE)
fp <- sum(negativos_res$signal_detected, na.rm = TRUE)
tn <- sum(!negativos_res$signal_detected, na.rm = TRUE)


sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
ppv <- tp / (tp + fp)
npv <- tn / (tn + fn)
accuracy <- (tp + tn) / (tp + tn + fp + fn)
f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)


cat(sprintf("
MATRIZ DE CONFUSIÓN (NULL MODEL):
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d


MÉTRICAS:
  Sensitivity (TPR):  %.3f  [%d/%d positivos detectados]
  Specificity (TNR):  %.3f  [%d/%d negativos correctos]
  PPV (Precision):    %.3f  [%d/%d predicciones+ correctas]
  NPV:                %.3f
  Accuracy:           %.3f
  F1-Score:           %.3f
",
tp, fn, fp, tn,
sensitivity, tp, tp+fn,
specificity, tn, tn+fp,
ppv, tp, tp+fp,
npv, accuracy, f1_score
))


# Curva ROC
results_for_roc <- results_dt[!is.na(max_ior)]
results_for_roc[, true_label := as.integer(type == "positive")]


if (nrow(results_for_roc) > 10 && length(unique(results_for_roc$true_label)) == 2) {
  roc_obj <- roc(response = results_for_roc$true_label,
                 predictor = results_for_roc$max_ior, direction = ">")
  auc_value <- auc(roc_obj)
  
  cat(sprintf("CURVA ROC:\n  AUC: %.3f  [IC 95%%: %.3f - %.3f]\n\n",
              auc_value, ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]))
  
  p_roc <- ggroc(roc_obj, legacy.axes = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.75, y = 0.25, label = sprintf("AUC = %.3f", auc_value), 
             size = 5, fontface = "bold") +
    labs(title = "Curva ROC - Detección de Señales (Giangreco)",
         subtitle = "Score: IOR máximo por triplete",
         x = "1 - Especificidad", y = "Sensibilidad") +
    theme_minimal()
  
  ggsave("validation_roc_curve.png", p_roc, width = 8, height = 7, dpi = 300)
} else {
  auc_value <- NA
}


# Merge con metadata
results_dt <- merge(results_dt, 
                    pos_meta[, .(drugA, drugB, meddra, N, n_injected = injected)],
                    by = c("drugA", "drugB", "meddra"), all.x = TRUE)


results_dt[, sample_size_cat := cut(N, breaks = c(0, 100, 250, 500, Inf),
                                     labels = c("50-100", "100-250", "250-500", ">500"),
                                     include.lowest = TRUE)]


perf_by_size <- results_dt[type == "positive", .(
  n = .N, n_detected = sum(signal_detected, na.rm = TRUE),
  detection_rate = mean(signal_detected, na.rm = TRUE),
  mean_max_ior = mean(max_ior, na.rm = TRUE)
), by = sample_size_cat]


cat("DESEMPEÑO POR TAMAÑO MUESTRAL:\n")
print(perf_by_size)


# ---------- VISUALIZACIONES ----------
# Distribución IOR
p_ior <- results_dt[!is.na(max_ior)] %>%
  ggplot(aes(x = type, y = max_ior, fill = type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black") +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "Distribución IOR Máximo", x = "Tipo", y = "IOR Máximo (log)") +
  scale_fill_manual(values = c("positive" = "#4DAF4A", "negative" = "#E41A1C")) +
  theme_minimal() + theme(legend.position = "none")
p_ior
ggsave("validation_ior_distribution.png", p_ior, width = 8, height = 6, dpi = 300)


# Tasa por tamaño
p_size <- ggplot(perf_by_size, aes(x = sample_size_cat, y = detection_rate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", detection_rate * 100, n_detected, n)),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
  labs(title = "Poder de Detección por Tamaño Muestral",
       x = "Reportes", y = "Tasa de Detección") +
  theme_minimal()
p_size
ggsave("validation_detection_by_sample_size.png", p_size, width = 9, height = 6, dpi = 300)


# Etapas significativas
p_stages <- results_dt[model_success == TRUE] %>%
  ggplot(aes(x = n_stages_significant, fill = type)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  labs(title = "Etapas con Señal Significativa", 
       x = "Nº Etapas NICHD", y = "Frecuencia") +
  scale_fill_manual(values = c("positive" = "#4DAF4A", "negative" = "#E41A1C")) +
  scale_x_continuous(breaks = 0:7) + theme_minimal()
p_stages
ggsave("validation_stages_distribution.png", p_stages, width = 10, height = 6, dpi = 300)


# ---------- TESTS ESTADÍSTICOS ----------
if (sum(!is.na(positivos_res$max_ior)) > 5 && sum(!is.na(negativos_res$max_ior)) > 5) {
  test_ior <- wilcox.test(positivos_res$max_ior, negativos_res$max_ior, alternative = "greater")
  
  cat(sprintf("\nMANN-WHITNEY TEST (IOR máximo):\n  W = %.2f, p = %.2e\n",
              test_ior$statistic, test_ior$p.value))
  cat(sprintf("  Medianas: Pos=%.3f, Neg=%.3f\n",
              median(positivos_res$max_ior, na.rm = TRUE),
              median(negativos_res$max_ior, na.rm = TRUE)))
}


# ---------- GUARDAR RESULTADOS ----------
fwrite(results_dt, "validation_detailed_results.csv")
fwrite(null_thresholds_dt, "validation_null_thresholds.csv")


summary_metrics <- data.table(
  metric = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1", "AUC"),
  value = c(sensitivity, specificity, ppv, npv, accuracy, f1_score, auc_value)
)
fwrite(summary_metrics, "validation_summary_metrics.csv")
fwrite(perf_by_size, "validation_performance_by_sample_size.csv")


ior_export <- results_dt[, c("triplet_id", "drugA", "drugB", "meddra", 
                              "type", "signal_detected", "n_stages_significant",
                              paste0("ior_stage_", 1:7)), with = FALSE]
fwrite(ior_export, "validation_ior_values_by_stage.csv")


save.image("validation_workspace.RData")


message("\n✓ VALIDACIÓN COMPLETADA")
message(sprintf("  Sensitivity: %.3f | PPV: %.3f | AUC: %.3f", 
                sensitivity, ppv, auc_value))


# ---------- ANÁLISIS ADICIONAL: COMPARACIÓN DE DINÁMICAS ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("ANÁLISIS EXPLORATORIO: SEÑAL POR TIPO DE DINÁMICA")
message(paste(rep("=", 70), collapse = ""))


# Merge con metadata de dinámicas
results_with_dynamic <- merge(
  results_dt,
  ground_truth[type == "positive", .(drugA, drugB, meddra, dynamic)],
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)


# Performance por dinámica verdadera (solo positivos)
perf_by_dynamic <- results_with_dynamic[type == "positive" & !is.na(dynamic), .(
  n_total = .N,
  n_detected = sum(signal_detected, na.rm = TRUE),
  detection_rate = mean(signal_detected, na.rm = TRUE),
  mean_stages_sig = mean(n_stages_significant, na.rm = TRUE),
  mean_max_ior = mean(max_ior, na.rm = TRUE),
  median_max_ior = median(max_ior, na.rm = TRUE)
), by = dynamic]


cat("\nDesempeño por Tipo de Dinámica Inyectada:\n")
cat("(Solo para interpretación exploratoria - no es objetivo primario)\n\n")
print(perf_by_dynamic[order(-detection_rate)])


fwrite(perf_by_dynamic, "validation_performance_by_dynamic.csv")
message("\n✓ Performance por dinámica: validation_performance_by_dynamic.csv")


# Gráfico exploratorio
p6 <- ggplot(perf_by_dynamic, 
             aes(x = reorder(dynamic, -detection_rate), y = detection_rate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", 
                                detection_rate * 100, n_detected, n_total)),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.15), labels = scales::percent) +
  labs(
    title = "Tasa de Detección por Tipo de Dinámica (Exploratorio)",
    subtitle = "Nota: el objetivo es detectar señal, no clasificar tipos",
    x = "Tipo de Dinámica Inyectada", y = "Tasa de Detección"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(p6)
ggsave("validation_detection_by_dynamic_exploratory.png", p6, 
       width = 10, height = 6, dpi = 300)