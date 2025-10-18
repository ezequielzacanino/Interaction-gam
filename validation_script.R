library(data.table)
library(tidyverse)
library(mgcv)
library(parallel)
library(pROC)
library(zoo)

# ============================================================================
# SCRIPT DE VALIDACIÓN: EVALUACIÓN DEL MODELO GAM DIFERENCIAL
# ============================================================================

# ---------- PARÁMETROS DE CONFIGURACIÓN ----------
ruta_ade_augmented <- "./ade_augmented.csv"
ruta_ground_truth <- "./ground_truth.csv"
ruta_positive_meta <- "./positive_triplets_metadata.csv"
ruta_negative_meta <- "./negative_triplets_metadata.csv"

# Parámetros de detección de señal
alpha_nominal <- 0.10        # nivel para IC 90%
alpha_significance <- 0.05   # p-valor para significancia
min_or_threshold <- 1.5      # IOR mínimo para considerar interacción relevante

# Parámetros de clasificación de dinámicas
n_stages <- 7
niveles_nichd <- c("term_neonatal","infancy","toddler","early_childhood",
                   "middle_childhood","early_adolescence","late_adolescence")

# Número de cores para paralelización
n_cores <- max(1, detectCores() - 1)

# ---------- FUNCIONES AUXILIARES ----------

# Clasificador de dinámicas basado en patrones de IOR
classify_dynamic_pattern <- function(ior_values, nichd_stages) {
  # ior_values: vector de IORs a lo largo de las etapas
  # Retorna: clasificación de dinámica detectada
  
  if (all(is.na(ior_values)) || length(ior_values) < 3) {
    return("unclassifiable")
  }
  
  # Normalizar valores para análisis
  ior_norm <- (ior_values - min(ior_values, na.rm = TRUE)) / 
              (max(ior_values, na.rm = TRUE) - min(ior_values, na.rm = TRUE) + 1e-10)
  
  # Calcular tendencias
  x <- seq_along(ior_norm)
  fit_linear <- lm(ior_norm ~ x)
  slope <- coef(fit_linear)[2]
  r_squared <- summary(fit_linear)$r.squared
  
  # Detectar pico/valle
  peak_idx <- which.max(ior_norm)
  valley_idx <- which.min(ior_norm)
  
  # reglas de clasificación
  if (r_squared > 0.6) {  # tendencia lineal 
    if (slope > 0.1) {
      return("increase")
    } else if (slope < -0.1) {
      return("decrease")
    } else {
      return("uniform")
    }
  } else if (peak_idx > 2 && peak_idx < (length(ior_norm) - 1)) {
    # pico en medio sería plateau
    return("plateau")
  } else if (valley_idx > 2 && valley_idx < (length(ior_norm) - 1)) {
    # valle en medio sería inverse plateau
    return("inverse_plateau")
  } else {
    # variación sin patrón claro
    sd_norm <- sd(ior_norm, na.rm = TRUE)
    if (sd_norm < 0.2) {
      return("uniform")
    } else {
      return("complex")
    }
  }
}

# Calcular métricas de similitud entre dinámicas
calculate_dynamic_similarity <- function(ior_observed, ior_expected) {
  # Correlación de Pearson entre patrones
  if (length(ior_observed) != length(ior_expected)) return(NA)
  if (all(is.na(ior_observed)) || all(is.na(ior_expected))) return(NA)
  
  cor_val <- cor(ior_observed, ior_expected, use = "pairwise.complete.obs")
  return(cor_val)
}

# Generar patrón esperado de IOR para una dinámica dada
generate_expected_ior_pattern <- function(dynamic_type, n_stages = 7, base_ior = 2) {
  x <- seq(-2, 2, length.out = n_stages)
  
  pattern <- switch(dynamic_type,
    "uniform" = rep(base_ior, n_stages),
    "increase" = base_ior * (1 + tanh(x)),
    "decrease" = base_ior * (2 - tanh(x)),
    "plateau" = {
      if (n_stages <= 4) {
        base_ior * (1 + tanh(x))
      } else {
        c(base_ior * (1 + tanh(x[1:4])), rep(base_ior * (1 + tanh(2)), n_stages - 4))
      }
    },
    "inverse_plateau" = {
      if (n_stages <= 4) {
        base_ior * (1 + tanh(x))
      } else {
        c(rep(base_ior * 0.5, n_stages - 4), base_ior * (1 + tanh(x[(n_stages-3):n_stages])))
      }
    },
    rep(base_ior, n_stages)
  )
  
  return(pattern)
}

# ---------- FUNCIÓN PRINCIPAL: AJUSTAR MODELO GAM DIFERENCIAL ----------
fit_differential_gam <- function(drugA_id, drugB_id, event_id, ade_data) {
  # Preparar datos para el triplete específico
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Dataset para modelo
  datos_modelo <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  
  # Variables
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # Verificar variabilidad mínima
  if (sum(datos_modelo$ea_ocurrio) < 5 || 
      sum(datos_modelo$droga_ab) < 3) {
    return(list(
      success = FALSE,
      reason = "insufficient_data",
      n_events = sum(datos_modelo$ea_ocurrio),
      n_coadmin = sum(datos_modelo$droga_ab)
    ))
  }
  
  # Ajustar modelo
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
    
    return(list(
      success = TRUE,
      model = modelo,
      data = datos_modelo,
      n_events = sum(datos_modelo$ea_ocurrio),
      n_coadmin = sum(datos_modelo$droga_ab)
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      reason = "model_error",
      error_msg = e$message,
      n_events = sum(datos_modelo$ea_ocurrio),
      n_coadmin = sum(datos_modelo$droga_ab)
    ))
  })
}

# ---------- CALCULAR IOR Y MÉTRICAS DEL MODELO ----------
calculate_ior_from_model <- function(model_result) {
  if (!model_result$success) {
    return(list(
      success = FALSE,
      ior_values = rep(NA, 7),
      ior_li = rep(NA, 7),
      ior_ls = rep(NA, 7),
      significant = rep(FALSE, 7)
    ))
  }
  
  modelo <- model_result$model
  datos_modelo <- model_result$data
  
  # Grid de predicción
  grid_dif <- CJ(
    nichd_num = 1:7,
    droga_a = c(0, 1),
    droga_b = c(0, 1)
  )
  grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # Predicciones
  pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
  
  grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
  
  z90 <- qnorm(0.95)
  grid_dif[, `:=`(li = lp - z90 * se, ls = lp + z90 * se)]
  
  # Calcular IOR por etapa: log(IOR) = LP_AB - LP_A - LP_B + LP_0
  w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                value.var = c("lp", "se", "li", "ls"))
  
  # Índices: 0_0, 0_1, 1_0, 1_1
  log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
  
  # Para ICs correctos, usar covarianza del lpmatrix
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
  
  log_ior_li <- log_ior - z90 * log_ior_se
  log_ior_ls <- log_ior + z90 * log_ior_se
  
  ior_values <- exp(log_ior)
  ior_li <- exp(log_ior_li)
  ior_ls <- exp(log_ior_ls)
  
  # Significancia: IC no incluye 1
  significant <- !(ior_li <= 1 & ior_ls >= 1)
  
  return(list(
    success = TRUE,
    ior_values = ior_values,
    ior_li = ior_li,
    ior_ls = ior_ls,
    log_ior = log_ior,
    log_ior_se = log_ior_se,
    significant = significant,
    n_significant = sum(significant),
    max_ior = max(ior_values),
    model_aic = AIC(modelo),
    model_deviance = deviance(modelo)
  ))
}

# ---------- CARGAR DATOS ----------
message("Cargando datos aumentados y ground truth...")
ade_aug <- fread(ruta_ade_augmented)
ground_truth <- fread(ruta_ground_truth)
pos_meta <- fread(ruta_positive_meta)
neg_meta <- fread(ruta_negative_meta)

# Normalizar nichd
ade_aug[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_aug[, nichd_num := as.integer(nichd)]

message("Dataset: ", nrow(ade_aug), " filas")
message("Ground truth: ", nrow(ground_truth), " tripletes")
message("  - Positivos: ", sum(ground_truth$type == "positive"))
message("  - Negativos: ", sum(ground_truth$type == "negative"))

# ---------- ANÁLISIS DE TRIPLETES ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("INICIANDO ANÁLISIS DE TRIPLETES")
message(paste(rep("=", 70), collapse = ""))

# Función wrapper para análisis paralelo
analyze_triplet <- function(i, triplet_row, ade_data) {
  drugA <- triplet_row$drugA
  drugB <- triplet_row$drugB
  event <- triplet_row$meddra
  
  # Ajustar modelo
  model_result <- fit_differential_gam(drugA, drugB, event, ade_data)
  
  # Calcular IOR
  ior_result <- calculate_ior_from_model(model_result)
  
  # Clasificar dinámica detectada
  if (ior_result$success) {
    detected_dynamic <- classify_dynamic_pattern(
      ior_result$ior_values, 
      niveles_nichd
    )
    
    # Si es positivo, calcular similitud con dinámica esperada
    if (triplet_row$type == "positive") {
      expected_pattern <- generate_expected_ior_pattern(
        triplet_row$dynamic, 
        n_stages = 7
      )
      similarity <- calculate_dynamic_similarity(
        ior_result$ior_values, 
        expected_pattern
      )
    } else {
      similarity <- NA
    }
  } else {
    detected_dynamic <- "failed"
    similarity <- NA
  }
  
  # Compilar resultados
  result <- list(
    triplet_id = i,
    drugA = drugA,
    drugB = drugB,
    meddra = event,
    type = triplet_row$type,
    true_dynamic = if (!is.na(triplet_row$dynamic)) triplet_row$dynamic else "none",
    detected_dynamic = detected_dynamic,
    model_success = model_result$success,
    n_events = model_result$n_events,
    n_coadmin = model_result$n_coadmin,
    signal_detected = ior_result$success && ior_result$n_significant > 0,
    n_stages_significant = if (ior_result$success) ior_result$n_significant else 0,
    max_ior = if (ior_result$success) ior_result$max_ior else NA,
    ior_values = if (ior_result$success) list(ior_result$ior_values) else list(rep(NA, 7)),
    ior_li = if (ior_result$success) list(ior_result$ior_li) else list(rep(NA, 7)),
    ior_ls = if (ior_result$success) list(ior_result$ior_ls) else list(rep(NA, 7)),
    pattern_similarity = similarity,
    model_aic = if (ior_result$success) ior_result$model_aic else NA
  )
  
  return(result)
}

# Ejecutar análisis en paralelo
message("\nAnalizando ", nrow(ground_truth), " tripletes usando ", n_cores, " cores...")
message("Esto puede tardar varios minutos...")

cl <- makeCluster(n_cores)
clusterExport(cl, c("fit_differential_gam", "calculate_ior_from_model",
                    "classify_dynamic_pattern", "calculate_dynamic_similarity",
                    "generate_expected_ior_pattern", "niveles_nichd"))
clusterEvalQ(cl, {
  library(data.table)
  library(mgcv)
})

pb <- txtProgressBar(max = nrow(ground_truth), style = 3)
results_list <- vector("list", nrow(ground_truth))

for (i in seq_len(nrow(ground_truth))) {
  results_list[[i]] <- analyze_triplet(i, ground_truth[i], ade_aug)
  setTxtProgressBar(pb, i)
}

stopCluster(cl)
close(pb)

# Convertir resultados a data.table
results_dt <- rbindlist(lapply(results_list, function(x) {
  x$ior_values <- NULL  # remover listas anidadas temporalmente
  x$ior_li <- NULL
  x$ior_ls <- NULL
  as.data.table(x)
}))

# Añadir vectores IOR por separado
ior_matrix <- do.call(rbind, lapply(results_list, function(x) x$ior_values[[1]]))
colnames(ior_matrix) <- paste0("ior_stage_", 1:7)
results_dt <- cbind(results_dt, ior_matrix)

# ---------- MÉTRICAS DE DESEMPEÑO ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("MÉTRICAS DE DESEMPEÑO DEL MODELO")
message(paste(rep("=", 70), collapse = ""))

# 1. Detección de señal (positivos vs negativos)
message("\n1. DETECCIÓN DE SEÑAL (Signal Detection)")
message(paste(rep("-", 70), collapse = ""))

positivos_res <- results_dt[type == "positive"]
negativos_res <- results_dt[type == "negative"]

# Matriz de confusión para detección
tp <- sum(positivos_res$signal_detected, na.rm = TRUE)
fn <- sum(!positivos_res$signal_detected, na.rm = TRUE)
fp <- sum(negativos_res$signal_detected, na.rm = TRUE)
tn <- sum(!negativos_res$signal_detected, na.rm = TRUE)

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
precision <- tp / (tp + fp)
f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)

cat(sprintf("
Matriz de Confusión:
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d

Métricas:
  - Sensibilidad (Recall):    %.3f  [%d/%d positivos detectados]
  - Especificidad:            %.3f  [%d/%d negativos sin señal]
  - Precisión (Precision):    %.3f
  - F1-Score:                 %.3f
  - Accuracy:                 %.3f
",
tp, fn, fp, tn,
sensitivity, tp, tp+fn,
specificity, tn, tn+fp,
precision, f1_score,
(tp + tn) / (tp + tn + fp + fn)
))

# 2. Clasificación de dinámicas (solo positivos exitosos)
message("\n2. CLASIFICACIÓN DE DINÁMICAS")
message(paste(rep("-", 70), collapse = ""))

positivos_clasificables <- positivos_res[
  model_success == TRUE & signal_detected == TRUE & detected_dynamic != "failed"
]

if (nrow(positivos_clasificables) > 0) {
  # Matriz de confusión de dinámicas
  confusion_dynamic <- table(
    Verdadera = positivos_clasificables$true_dynamic,
    Detectada = positivos_clasificables$detected_dynamic
  )
  
  cat("Matriz de Confusión de Dinámicas:\n")
  print(confusion_dynamic)
  
  # Accuracy de clasificación
  correct_classifications <- sum(
    positivos_clasificables$true_dynamic == positivos_clasificables$detected_dynamic
  )
  accuracy_dynamic <- correct_classifications / nrow(positivos_clasificables)
  
  cat(sprintf("\nAccuracy de Clasificación: %.3f [%d/%d]\n",
              accuracy_dynamic, correct_classifications, 
              nrow(positivos_clasificables)))
  
  # Similitud de patrones (correlación)
  pattern_correlations <- positivos_clasificables[
    !is.na(pattern_similarity), pattern_similarity
  ]
  
  if (length(pattern_correlations) > 0) {
    cat(sprintf("\nSimilitud de Patrones (Correlación):\n"))
    cat(sprintf("  - Media:    %.3f\n", mean(pattern_correlations)))
    cat(sprintf("  - Mediana:  %.3f\n", median(pattern_correlations)))
    cat(sprintf("  - SD:       %.3f\n", sd(pattern_correlations)))
    cat(sprintf("  - Rango:    [%.3f, %.3f]\n", 
                min(pattern_correlations), max(pattern_correlations)))
  }
} else {
  message("⚠ No hay positivos clasificables para evaluar dinámicas")
}

# 3. Análisis por dinámica verdadera
message("\n3. DESEMPEÑO POR TIPO DE DINÁMICA")
message(paste(rep("-", 70), collapse = ""))

performance_by_dynamic <- positivos_res[, .(
  n_total = .N,
  n_detected = sum(signal_detected, na.rm = TRUE),
  detection_rate = mean(signal_detected, na.rm = TRUE),
  n_classified = sum(model_success & signal_detected & detected_dynamic != "failed"),
  n_correct = sum(true_dynamic == detected_dynamic & signal_detected, na.rm = TRUE),
  mean_max_ior = mean(max_ior, na.rm = TRUE),
  mean_n_stages_sig = mean(n_stages_significant, na.rm = TRUE)
), by = true_dynamic]

performance_by_dynamic[, classification_accuracy := n_correct / n_classified]

print(performance_by_dynamic)

# 4. Efecto del tamaño muestral
message("\n4. EFECTO DEL TAMAÑO MUESTRAL")
message(paste(rep("-", 70), collapse = ""))

# Añadir info de tamaño muestral desde metadata
results_dt <- merge(
  results_dt, 
  pos_meta[, .(drugA, drugB, meddra, N, n_injected = injected)],
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)

# Categorizar por tamaño
results_dt[, sample_size_cat := cut(
  N, 
  breaks = c(0, 100, 250, 500, Inf),
  labels = c("50-100", "100-250", "250-500", ">500"),
  include.lowest = TRUE
)]

perf_by_size <- results_dt[type == "positive", .(
  n = .N,
  detection_rate = mean(signal_detected, na.rm = TRUE),
  mean_stages_sig = mean(n_stages_significant, na.rm = TRUE),
  mean_max_ior = mean(max_ior, na.rm = TRUE)
), by = sample_size_cat]

cat("\nDesempeño por Rango de Reportes:\n")
print(perf_by_size)

# 5. Causas de fallo
message("\n5. ANÁLISIS DE FALLOS")
message(paste(rep("-", 70), collapse = ""))

failed_positives <- positivos_res[signal_detected == FALSE | model_success == FALSE]

if (nrow(failed_positives) > 0) {
  cat(sprintf("Positivos no detectados: %d / %d (%.1f%%)\n\n",
              nrow(failed_positives), nrow(positivos_res),
              100 * nrow(failed_positives) / nrow(positivos_res)))
  
  cat("Distribución de causas:\n")
  cat(sprintf("  - Modelo falló:           %d\n", 
              sum(!failed_positives$model_success)))
  cat(sprintf("  - Sin señal significativa: %d\n", 
              sum(failed_positives$model_success & !failed_positives$signal_detected)))
  
  # Estadísticas de los fallidos
  cat(sprintf("\nCaracterísticas de fallidos:\n"))
  cat(sprintf("  - Media eventos:      %.1f\n", 
              mean(failed_positives$n_events, na.rm = TRUE)))
  cat(sprintf("  - Media coadmin:      %.1f\n", 
              mean(failed_positives$n_coadmin, na.rm = TRUE)))
  cat(sprintf("  - Media inyectados:   %.1f\n", 
              mean(failed_positives$n_injected, na.rm = TRUE)))
}

# ---------- VISUALIZACIONES ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("GENERANDO VISUALIZACIONES")
message(paste(rep("=", 70), collapse = ""))

# Gráfico 1: ROC-like para umbrales de IOR
if (sum(!is.na(results_dt$max_ior)) > 10) {
  results_dt[, true_positive := as.integer(type == "positive")]
  
  # Calcular TPR y FPR para diferentes umbrales
  thresholds <- seq(0.5, 5, by = 0.1)
  roc_data <- lapply(thresholds, function(thr) {
    results_dt[!is.na(max_ior), .(
      threshold = thr,
      tpr = sum(true_positive == 1 & max_ior > thr) / sum(true_positive == 1),
      fpr = sum(true_positive == 0 & max_ior > thr) / sum(true_positive == 0)
    )]
  })
  roc_dt <- rbindlist(roc_data)
  
  p1 <- ggplot(roc_dt, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1.2, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Curva ROC - Detección basada en IOR máximo",
      x = "Tasa de Falsos Positivos (FPR)",
      y = "Tasa de Verdaderos Positivos (TPR)"
    ) +
    theme_minimal(base_size = 12) +
    coord_equal()
  
  print(p1)
  ggsave("validation_roc_curve.png", p1, width = 8, height = 6, dpi = 300)
}

# Gráfico 2: Distribución de IOR por tipo
p2 <- results_dt[!is.na(max_ior)] %>%
  ggplot(aes(x = type, y = max_ior, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title = "Distribución de IOR Máximo por Tipo de Triplete",
    x = "Tipo", y = "IOR Máximo (escala log)",
    fill = "Tipo"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p2)
ggsave("validation_ior_distribution.png", p2, width = 8, height = 6, dpi = 300)

# Gráfico 3: Matriz de confusión de dinámicas
if (exists("confusion_dynamic") && nrow(positivos_clasificables) > 0) {
  conf_df <- as.data.frame(confusion_dynamic)
  
  p3 <- ggplot(conf_df, aes(x = Detectada, y = Verdadera, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue", 
                        name = "Frecuencia") +
    labs(
      title = "Matriz de Confusión - Clasificación de Dinámicas",
      x = "Dinámica Detectada", y = "Dinámica Verdadera"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p3)
  ggsave("validation_confusion_matrix.png", p3, width = 8, height = 7, dpi = 300)
}

# Gráfico 4: Tasa de detección por dinámica
p4 <- ggplot(performance_by_dynamic, 
             aes(x = reorder(true_dynamic, -detection_rate), y = detection_rate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f\n(%d/%d)", detection_rate, 
                                n_detected, n_total)),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
  labs(
    title = "Tasa de Detección por Tipo de Dinámica",
    subtitle = "Proporción de tripletes positivos con señal detectada",
    x = "Tipo de Dinámica", y = "Tasa de Detección"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)
ggsave("validation_detection_by_dynamic.png", p4, width = 10, height = 6, dpi = 300)

# Gráfico 5: Ejemplos de patrones IOR recuperados
if (nrow(positivos_clasificables) > 0) {
  # Seleccionar ejemplos representativos (uno por dinámica si es posible)
  ejemplos <- positivos_clasificables[, .SD[which.max(pattern_similarity)], 
                                      by = true_dynamic]
  
  if (nrow(ejemplos) > 0) {
    # Preparar datos para plotting
    plot_data_list <- list()
    for (i in seq_len(min(nrow(ejemplos), 5))) {
      ej <- ejemplos[i]
      ior_vals <- unlist(results_list[[ej$triplet_id]]$ior_values)
      ior_li_vals <- unlist(results_list[[ej$triplet_id]]$ior_li)
      ior_ls_vals <- unlist(results_list[[ej$triplet_id]]$ior_ls)
      expected <- generate_expected_ior_pattern(ej$true_dynamic, n_stages = 7)
      
      plot_data_list[[i]] <- data.table(
        triplet_id = i,
        dynamic = ej$true_dynamic,
        stage = 1:7,
        stage_name = factor(niveles_nichd, levels = niveles_nichd),
        observed = ior_vals,
        observed_li = ior_li_vals,
        observed_ls = ior_ls_vals,
        expected = expected,
        similarity = ej$pattern_similarity
      )
    }
    
    plot_data <- rbindlist(plot_data_list)
    
    p5 <- ggplot(plot_data, aes(x = stage)) +
      geom_ribbon(aes(ymin = observed_li, ymax = observed_ls), 
                  alpha = 0.2, fill = "steelblue") +
      geom_line(aes(y = observed, color = "Observado"), linewidth = 1) +
      geom_point(aes(y = observed, color = "Observado"), size = 2) +
      geom_line(aes(y = expected, color = "Esperado"), 
                linewidth = 1, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dotted", color = "gray40") +
      scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
      scale_color_manual(values = c("Observado" = "steelblue", 
                                     "Esperado" = "coral")) +
      facet_wrap(~ paste0(dynamic, "\n(r=", round(similarity, 2), ")"), 
                 scales = "free_y", ncol = 3) +
      labs(
        title = "Ejemplos de Patrones IOR Recuperados vs Esperados",
        subtitle = "IC 90% mostrado como área sombreada",
        x = "Etapa NICHD", y = "IOR (escala log)",
        color = "Tipo"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
      )
    
    print(p5)
    ggsave("validation_ior_patterns_examples.png", p5, 
           width = 12, height = 8, dpi = 300)
  }
}

# Gráfico 6: Relación entre tamaño muestral y detección
p6 <- results_dt[type == "positive" & !is.na(N)] %>%
  ggplot(aes(x = N, y = as.integer(signal_detected))) +
  geom_jitter(aes(color = true_dynamic), 
              width = 0, height = 0.05, alpha = 0.6, size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(breaks = c(50, 100, 250, 500, 1000)) +
  scale_y_continuous(breaks = c(0, 1), labels = c("No detectado", "Detectado")) +
  labs(
    title = "Efecto del Tamaño Muestral en la Detección",
    subtitle = "Regresión logística ajustada (curva negra)",
    x = "Número de Reportes del Triplete (escala log)",
    y = "Estado de Detección",
    color = "Dinámica\nVerdadera"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

print(p6)
ggsave("validation_sample_size_effect.png", p6, width = 10, height = 6, dpi = 300)

# ---------- GUARDAR RESULTADOS ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("GUARDANDO RESULTADOS")
message(paste(rep("=", 70), collapse = ""))

# Resultados detallados
fwrite(results_dt, "validation_detailed_results.csv")
message("✓ Resultados detallados: validation_detailed_results.csv")

# Resumen de métricas
summary_metrics <- data.table(
  metric = c("Sensitivity", "Specificity", "Precision", "F1-Score", 
             "Accuracy_Detection", "Accuracy_Classification"),
  value = c(sensitivity, specificity, precision, f1_score,
            (tp + tn) / (tp + tn + fp + fn),
            if (nrow(positivos_clasificables) > 0) accuracy_dynamic else NA)
)
fwrite(summary_metrics, "validation_summary_metrics.csv")
message("✓ Resumen de métricas: validation_summary_metrics.csv")

# Performance por dinámica
fwrite(performance_by_dynamic, "validation_performance_by_dynamic.csv")
message("✓ Performance por dinámica: validation_performance_by_dynamic.csv")

# Performance por tamaño muestral
fwrite(perf_by_size, "validation_performance_by_sample_size.csv")
message("✓ Performance por tamaño: validation_performance_by_sample_size.csv")

# Exportar IORs completos para análisis adicional
ior_export <- results_dt[, c("triplet_id", "drugA", "drugB", "meddra", 
                              "type", "true_dynamic", "detected_dynamic",
                              paste0("ior_stage_", 1:7)), with = FALSE]
fwrite(ior_export, "validation_ior_values_by_stage.csv")
message("✓ Valores IOR por etapa: validation_ior_values_by_stage.csv")

# ---------- REPORTE FINAL ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("RESUMEN EJECUTIVO")
message(paste(rep("=", 70), collapse = ""))

cat(sprintf("
VALIDACIÓN DEL MODELO GAM DIFERENCIAL
======================================

DATASET:
  • Total de tripletes analizados:     %d
  • Tripletes positivos:                %d
  • Tripletes negativos:                %d
  • Modelos exitosos:                   %d (%.1f%%)

DETECCIÓN DE SEÑAL:
  • Sensibilidad (recall):              %.3f
  • Especificidad:                      %.3f
  • Precisión:                          %.3f
  • F1-Score:                           %.3f

CLASIFICACIÓN DE DINÁMICAS:
  • Positivos clasificables:            %d
  • Clasificaciones correctas:          %s
  • Accuracy de clasificación:          %s
  • Correlación media de patrones:      %s

DINÁMICAS MÁS DIFÍCILES DE DETECTAR:
%s

RECOMENDACIONES:
%s

ARCHIVOS GENERADOS:
  - validation_detailed_results.csv
  - validation_summary_metrics.csv
  - validation_performance_by_dynamic.csv
  - validation_ior_values_by_stage.csv
  - validation_roc_curve.png
  - validation_ior_distribution.png
  - validation_confusion_matrix.png
  - validation_detection_by_dynamic.png
  - validation_ior_patterns_examples.png
  - validation_sample_size_effect.png

",
nrow(ground_truth),
nrow(positivos_res),
nrow(negativos_res),
sum(results_dt$model_success),
100 * mean(results_dt$model_success),
sensitivity,
specificity,
precision,
f1_score,
nrow(positivos_clasificables),
if (nrow(positivos_clasificables) > 0) 
  sprintf("%d / %d", correct_classifications, nrow(positivos_clasificables)) 
else "N/A",
if (nrow(positivos_clasificables) > 0) 
  sprintf("%.3f", accuracy_dynamic) 
else "N/A",
if (length(pattern_correlations) > 0) 
  sprintf("%.3f (SD = %.3f)", mean(pattern_correlations), sd(pattern_correlations))
else "N/A",
paste(performance_by_dynamic[order(detection_rate)][1:min(3, .N)], 
      collapse = "\n  • "),
if (sensitivity < 0.7) {
  "  ⚠ Sensibilidad baja: considerar aumentar n_pos o ajustar effect_size
  ⚠ Revisar tripletes con bajo N de reportes
  ⚠ Considerar reducir min_reports_triplet para más candidatos"
} else if (specificity < 0.8) {
  "  ⚠ Especificidad baja: muchos falsos positivos
  ⚠ Considerar aumentar umbral de significancia
  ⚠ Revisar negativos con señales espurias"
} else {
  "  ✓ Desempeño general satisfactorio
  ✓ El modelo detecta efectivamente las interacciones inyectadas
  ✓ Considerar aumentar complejidad para dinámicas sutiles"
}
))

# ---------- ANÁLISIS ADICIONAL: POWER ANALYSIS ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("ANÁLISIS DE POTENCIA ESTADÍSTICA")
message(paste(rep("=", 70), collapse = ""))

# Estimar poder de detección por combinación de factores
if (nrow(positivos_res) > 0) {
  # Usar results_dt[type == "positive"] en lugar de positivos_res
  power_analysis <- results_dt[type == "positive", .(
    n = .N,
    power = mean(signal_detected, na.rm = TRUE),
    mean_n_events = mean(n_events, na.rm = TRUE),
    mean_n_coadmin = mean(n_coadmin, na.rm = TRUE),
    mean_n_injected = mean(n_injected, na.rm = TRUE)
  ), by = .(true_dynamic, sample_size_cat)]
  
  cat("\nPoder de Detección por Dinámica y Tamaño Muestral:\n")
  print(power_analysis[order(true_dynamic, sample_size_cat)])
  
  fwrite(power_analysis, "validation_power_analysis.csv")
  message("\n✓ Análisis de poder: validation_power_analysis.csv")
  
  # Visualización de poder
  p7 <- ggplot(power_analysis[!is.na(sample_size_cat)], 
               aes(x = sample_size_cat, y = power, 
                   group = true_dynamic, color = true_dynamic)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Poder Estadístico por Dinámica y Tamaño Muestral",
      x = "Rango de Reportes", y = "Poder (Tasa de Detección)",
      color = "Dinámica"
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
  
  print(p7)
  ggsave("validation_power_analysis.png", p7, width = 10, height = 6, dpi = 300)
}

# ---------- TEST ESTADÍSTICO: COMPARACIÓN CON NULL ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("TEST ESTADÍSTICO: DIFERENCIA POSITIVOS vs NEGATIVOS")
message(paste(rep("=", 70), collapse = ""))

# Test de Mann-Whitney para IOR máximo
if (sum(!is.na(positivos_res$max_ior)) > 5 && 
    sum(!is.na(negativos_res$max_ior)) > 5) {
  
  test_ior <- wilcox.test(
    positivos_res$max_ior,
    negativos_res$max_ior,
    alternative = "greater"
  )
  
  cat(sprintf("
Test de Mann-Whitney para IOR Máximo:
  H0: IOR_positivos <= IOR_negativos
  Ha: IOR_positivos > IOR_negativos
  
  W = %.2f
  p-valor = %.2e
  %s
  
  Medianas:
    - Positivos: %.3f
    - Negativos: %.3f
    - Diferencia: %.3f
",
  test_ior$statistic,
  test_ior$p.value,
  if (test_ior$p.value < 0.05) "✓ RECHAZAMOS H0: Positivos tienen IOR mayor" 
  else "✗ No hay evidencia suficiente",
  median(positivos_res$max_ior, na.rm = TRUE),
  median(negativos_res$max_ior, na.rm = TRUE),
  median(positivos_res$max_ior, na.rm = TRUE) - 
    median(negativos_res$max_ior, na.rm = TRUE)
  ))
}

# Test chi-cuadrado para proporción de detección
test_prop <- prop.test(
  x = c(sum(positivos_res$signal_detected, na.rm = TRUE),
        sum(negativos_res$signal_detected, na.rm = TRUE)),
  n = c(nrow(positivos_res), nrow(negativos_res)),
  alternative = "greater"
)

cat(sprintf("
Test de Proporciones para Tasa de Detección:
  H0: prop_positivos <= prop_negativos
  Ha: prop_positivos > prop_negativos
  
  X² = %.2f
  p-valor = %.2e
  %s
  
  Proporciones:
    - Positivos: %.3f (%d/%d)
    - Negativos: %.3f (%d/%d)
",
test_prop$statistic,
test_prop$p.value,
if (test_prop$p.value < 0.05) "✓ RECHAZAMOS H0: Positivos se detectan más" 
else "✗ No hay diferencia significativa",
sum(positivos_res$signal_detected, na.rm = TRUE) / nrow(positivos_res),
sum(positivos_res$signal_detected, na.rm = TRUE), nrow(positivos_res),
sum(negativos_res$signal_detected, na.rm = TRUE) / nrow(negativos_res),
sum(negativos_res$signal_detected, na.rm = TRUE), nrow(negativos_res)
))

# ---------- FINALIZACIÓN ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("VALIDACIÓN COMPLETADA EXITOSAMENTE")
message(paste(rep("=", 70), collapse = ""))
message("\nTiempo total de ejecución: ", 
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("\nRevisa los archivos generados para análisis detallado.")
message("Considera ajustar parámetros de augmentation si el desempeño es bajo.")

# Guardar workspace para análisis posterior
save.image("validation_workspace.RData")
message("\n✓ Workspace guardado: validation_workspace.RData")
message("  (Puedes cargarlo con: load('validation_workspace.RData'))")

message("\n¡Análisis completado! 🎉\n")