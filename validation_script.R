library(data.table)
library(tidyverse)
library(mgcv)
library(parallel)
library(pROC)

# ============================================================================
# SCRIPT DE VALIDACIÓN: EVALUACIÓN ESTILO GIANGRECO
# Objetivo: Detectar PRESENCIA de señal de interacción (no clasificar tipos)
# ============================================================================

# ---------- PARÁMETROS DE CONFIGURACIÓN ----------
ruta_ade_augmented <- "./ade_augmented.csv"
ruta_ground_truth <- "./ground_truth.csv"
ruta_positive_meta <- "./positive_triplets_metadata.csv"
ruta_negative_meta <- "./negative_triplets_metadata.csv"

# Parámetros de detección de señal (estilo Giangreco)
alpha_nominal <- 0.10        # IC 90% (como en Giangreco)
z90 <- qnorm(0.95)           # z-score para IC 90%

# Criterios de señal significativa (adaptar según necesidad)
# Opción 1: Al menos 1 etapa con IC90 que no incluya IOR=1
min_stages_significant <- 1

# Opción 2: Umbral de IOR mínimo (opcional, comentar si no se usa)
# min_ior_threshold <- 1.5

# Número de cores para paralelización
n_cores <- max(1, detectCores() - 1)

# ---------- FUNCIONES AUXILIARES ----------

# Calcular señal GAM diferencial para un triplete
fit_differential_gam <- function(drugA_id, drugB_id, event_id, ade_data) {
  # Preparar datos
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  datos_modelo <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  
  # Variables de exposición y evento
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # Verificar suficientes datos
  if (sum(datos_modelo$ea_ocurrio) < 5 || sum(datos_modelo$droga_ab) < 3) {
    return(list(
      success = FALSE,
      reason = "insufficient_data",
      n_events = sum(datos_modelo$ea_ocurrio),
      n_coadmin = sum(datos_modelo$droga_ab)
    ))
  }
  
  # Ajustar modelo diferencial
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

# Calcular IOR y determinar presencia de señal
calculate_interaction_signal <- function(model_result) {
  if (!model_result$success) {
    return(list(
      success = FALSE,
      signal_detected = FALSE,
      n_stages_significant = 0,
      max_ior = NA,
      ior_values = rep(NA, 7),
      ior_li = rep(NA, 7),
      ior_ls = rep(NA, 7)
    ))
  }
  
  modelo <- model_result$model
  
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
  grid_dif[, `:=`(li = lp - z90 * se, ls = lp + z90 * se)]
  
  # Calcular IOR: log(IOR) = LP_AB - LP_A - LP_B + LP_0
  w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                value.var = c("lp", "se", "li", "ls"))
  
  log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
  
  # SE correcto usando covarianza
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
  
  # Criterio de señal: IC90 no incluye 1 (como Giangreco con "nominally significant")
  significant <- !(ior_li <= 1 & ior_ls >= 1)
  n_significant <- sum(significant)
  
  # Detectar señal si hay al menos min_stages_significant etapas significativas
  signal_detected <- n_significant >= min_stages_significant
  
  # Opcional: agregar umbral de magnitud
  # if (exists("min_ior_threshold")) {
  #   signal_detected <- signal_detected && max(ior_values, na.rm = TRUE) > min_ior_threshold
  # }
  
  return(list(
    success = TRUE,
    signal_detected = signal_detected,
    n_stages_significant = n_significant,
    max_ior = max(ior_values, na.rm = TRUE),
    mean_ior = mean(ior_values, na.rm = TRUE),
    ior_values = ior_values,
    ior_li = ior_li,
    ior_ls = ior_ls,
    log_ior = log_ior,
    log_ior_se = log_ior_se,
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
niveles_nichd <- c("term_neonatal","infancy","toddler","early_childhood",
                   "middle_childhood","early_adolescence","late_adolescence")
ade_aug[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_aug[, nichd_num := as.integer(nichd)]

message("Dataset: ", nrow(ade_aug), " filas")
message("Ground truth: ", nrow(ground_truth), " tripletes")
message("  - Positivos: ", sum(ground_truth$type == "positive"))
message("  - Negativos: ", sum(ground_truth$type == "negative"))

# ---------- ANÁLISIS DE TRIPLETES ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("INICIANDO ANÁLISIS DE DETECCIÓN DE SEÑALES")
message(paste(rep("=", 70), collapse = ""))

# Función wrapper para análisis paralelo
analyze_triplet <- function(i, triplet_row, ade_data) {
  drugA <- triplet_row$drugA
  drugB <- triplet_row$drugB
  event <- triplet_row$meddra
  
  # Ajustar modelo
  model_result <- fit_differential_gam(drugA, drugB, event, ade_data)
  
  # Detectar señal
  signal_result <- calculate_interaction_signal(model_result)
  
  # Compilar resultados
  result <- list(
    triplet_id = i,
    drugA = drugA,
    drugB = drugB,
    meddra = event,
    type = triplet_row$type,
    model_success = model_result$success,
    n_events = model_result$n_events,
    n_coadmin = model_result$n_coadmin,
    signal_detected = signal_result$signal_detected,
    n_stages_significant = signal_result$n_stages_significant,
    max_ior = signal_result$max_ior,
    mean_ior = if (signal_result$success) signal_result$mean_ior else NA,
    model_aic = if (signal_result$success) signal_result$model_aic else NA,
    ior_values = if (signal_result$success) list(signal_result$ior_values) else list(rep(NA, 7)),
    ior_li = if (signal_result$success) list(signal_result$ior_li) else list(rep(NA, 7)),
    ior_ls = if (signal_result$success) list(signal_result$ior_ls) else list(rep(NA, 7))
  )
  
  return(result)
}

# Ejecutar análisis en paralelo
message("\nAnalizando ", nrow(ground_truth), " tripletes usando ", n_cores, " cores...")

cl <- makeCluster(n_cores)
clusterExport(cl, c("fit_differential_gam", "calculate_interaction_signal",
                    "z90", "min_stages_significant", "niveles_nichd"))
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

# Convertir a data.table
results_dt <- rbindlist(lapply(results_list, function(x) {
  x$ior_values <- NULL
  x$ior_li <- NULL
  x$ior_ls <- NULL
  as.data.table(x)
}))

# Añadir vectores IOR
ior_matrix <- do.call(rbind, lapply(results_list, function(x) x$ior_values[[1]]))
colnames(ior_matrix) <- paste0("ior_stage_", 1:7)
results_dt <- cbind(results_dt, ior_matrix)

# ---------- MÉTRICAS DE DESEMPEÑO (ESTILO GIANGRECO) ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("MÉTRICAS DE DESEMPEÑO DE DETECCIÓN DE SEÑAL")
message(paste(rep("=", 70), collapse = ""))

positivos_res <- results_dt[type == "positive"]
negativos_res <- results_dt[type == "negative"]

# 1. MATRIZ DE CONFUSIÓN Y MÉTRICAS BÁSICAS
message("\n1. MATRIZ DE CONFUSIÓN")
message(paste(rep("-", 70), collapse = ""))

tp <- sum(positivos_res$signal_detected, na.rm = TRUE)
fn <- sum(!positivos_res$signal_detected, na.rm = TRUE)
fp <- sum(negativos_res$signal_detected, na.rm = TRUE)
tn <- sum(!negativos_res$signal_detected, na.rm = TRUE)

# Métricas según Giangreco (Figura 3B)
sensitivity <- tp / (tp + fn)  # Recall / TPR
specificity <- tn / (tn + fp)  # TNR
ppv <- tp / (tp + fp)          # Precision / PPV
npv <- tn / (tn + fn)          # NPV
accuracy <- (tp + tn) / (tp + tn + fp + fn)
f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)

cat(sprintf("
Matriz de Confusión:
                 Predicción
                 Señal    Sin Señal
Verdad Positivo  %4d     %4d
       Negativo  %4d     %4d

Métricas de Desempeño (estilo Giangreco):
  - Sensitivity (Recall/TPR):  %.3f  [%d/%d positivos detectados]
  - Specificity (TNR):         %.3f  [%d/%d negativos correctos]
  - PPV (Precision):           %.3f  [%d/%d predicciones positivas correctas]
  - NPV:                       %.3f  [%d/%d predicciones negativas correctas]
  - Accuracy:                  %.3f
  - F1-Score:                  %.3f
  
Interpretación:
  - TPR alto: buen poder para detectar interacciones reales
  - PPV alto: pocas falsas alarmas en señales detectadas
  - NPV alto: confianza en ausencia de señal cuando no se detecta
",
tp, fn, fp, tn,
sensitivity, tp, tp+fn,
specificity, tn, tn+fp,
ppv, tp, tp+fp,
npv, tn, tn+fn,
accuracy,
f1_score
))

# 2. CURVA ROC (como Giangreco Figura 3C)
message("\n2. ANÁLISIS ROC")
message(paste(rep("-", 70), collapse = ""))

# Usar max_ior como score continuo
results_for_roc <- results_dt[!is.na(max_ior)]
results_for_roc[, true_label := as.integer(type == "positive")]

if (nrow(results_for_roc) > 10 && 
    length(unique(results_for_roc$true_label)) == 2) {
  
  roc_obj <- roc(
    response = results_for_roc$true_label,
    predictor = results_for_roc$max_ior,
    direction = ">"
  )
  
  auc_value <- auc(roc_obj)
  
  cat(sprintf("
Curva ROC (usando IOR máximo como score):
  - AUC:  %.3f
  - IC 95%%: [%.3f, %.3f]
  
", 
  auc_value,
  ci.auc(roc_obj)[1],
  ci.auc(roc_obj)[3]
  ))
  
  # Encontrar umbral óptimo (Youden's index)
  coords_roc <- coords(roc_obj, "best", best.method = "youden")
  cat(sprintf("Umbral óptimo (Youden): IOR = %.2f\n", coords_roc$threshold))
  cat(sprintf("  En este umbral: Sens=%.3f, Spec=%.3f\n", 
              coords_roc$sensitivity, coords_roc$specificity))
  
} else {
  message("⚠ Insuficientes datos para análisis ROC")
  auc_value <- NA
}

# 3. ANÁLISIS POR TAMAÑO MUESTRAL (Giangreco evalúa esto)
message("\n3. EFECTO DEL TAMAÑO MUESTRAL")
message(paste(rep("-", 70), collapse = ""))

# Merge con metadata para obtener N
results_dt <- merge(
  results_dt, 
  pos_meta[, .(drugA, drugB, meddra, N, n_injected = injected)],
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)

# Categorizar tamaño
results_dt[, sample_size_cat := cut(
  N, 
  breaks = c(0, 100, 250, 500, Inf),
  labels = c("50-100", "100-250", "250-500", ">500"),
  include.lowest = TRUE
)]

perf_by_size <- results_dt[type == "positive", .(
  n = .N,
  n_detected = sum(signal_detected, na.rm = TRUE),
  detection_rate = mean(signal_detected, na.rm = TRUE),
  mean_stages_sig = mean(n_stages_significant, na.rm = TRUE),
  mean_max_ior = mean(max_ior, na.rm = TRUE)
), by = sample_size_cat]

cat("\nDesempeño por Rango de Reportes (Positivos):\n")
print(perf_by_size)

# 4. ANÁLISIS DE FALLOS
message("\n4. ANÁLISIS DE FALLOS")
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
  
  cat(sprintf("\nCaracterísticas de fallidos:\n"))
  cat(sprintf("  - Media eventos:      %.1f\n", 
              mean(failed_positives$n_events, na.rm = TRUE)))
  cat(sprintf("  - Media coadmin:      %.1f\n", 
              mean(failed_positives$n_coadmin, na.rm = TRUE)))
  cat(sprintf("  - Media inyectados:   %.1f\n", 
              mean(failed_positives$n_injected, na.rm = TRUE)))
}

false_positives <- negativos_res[signal_detected == TRUE]
if (nrow(false_positives) > 0) {
  cat(sprintf("\nNegativos con señal espuria: %d / %d (%.1f%%)\n",
              nrow(false_positives), nrow(negativos_res),
              100 * nrow(false_positives) / nrow(negativos_res)))
  cat(sprintf("  - Media eventos:      %.1f\n", 
              mean(false_positives$n_events, na.rm = TRUE)))
  cat(sprintf("  - Media coadmin:      %.1f\n", 
              mean(false_positives$n_coadmin, na.rm = TRUE)))
}

# ---------- VISUALIZACIONES ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("GENERANDO VISUALIZACIONES")
message(paste(rep("=", 70), collapse = ""))

# Gráfico 1: Curva ROC
if (!is.na(auc_value)) {
  p1 <- ggroc(roc_obj, legacy.axes = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.75, y = 0.25, 
             label = sprintf("AUC = %.3f", auc_value), 
             size = 5, fontface = "bold") +
    labs(
      title = "Curva ROC - Detección de Señales de Interacción",
      subtitle = "Score: IOR máximo por triplete",
      x = "Tasa de Falsos Positivos (1 - Especificidad)",
      y = "Tasa de Verdaderos Positivos (Sensibilidad)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p1)
  ggsave("validation_roc_curve.png", p1, width = 8, height = 7, dpi = 300)
}

# Gráfico 2: Distribución de IOR por tipo
p2 <- results_dt[!is.na(max_ior)] %>%
  ggplot(aes(x = type, y = max_ior, fill = type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "white", color = "black") +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "Distribución de IOR Máximo por Tipo de Triplete",
    subtitle = "Rombo blanco = media; líneas = cuartiles",
    x = "Tipo", y = "IOR Máximo (escala log)",
    fill = "Tipo"
  ) +
  scale_fill_manual(values = c("positive" = "#4DAF4A", "negative" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p2)
ggsave("validation_ior_distribution.png", p2, width = 8, height = 6, dpi = 300)

# Gráfico 3: Tasa de detección por tamaño muestral
p3 <- ggplot(perf_by_size, aes(x = sample_size_cat, y = detection_rate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", 
                                detection_rate * 100, n_detected, n)),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
  labs(
    title = "Poder de Detección por Tamaño Muestral",
    subtitle = "Proporción de tripletes positivos con señal detectada",
    x = "Número de Reportes", y = "Tasa de Detección"
  ) +
  theme_minimal(base_size = 12)

print(p3)
ggsave("validation_detection_by_sample_size.png", p3, width = 9, height = 6, dpi = 300)

# Gráfico 4: Relación continua N vs detección
p4 <- results_dt[type == "positive" & !is.na(N)] %>%
  ggplot(aes(x = N, y = as.integer(signal_detected))) +
  geom_jitter(aes(color = signal_detected), 
              width = 0, height = 0.05, alpha = 0.5, size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, color = "black", linewidth = 1.2) +
  scale_x_log10(breaks = c(50, 100, 250, 500, 1000, 2000)) +
  scale_y_continuous(breaks = c(0, 1), labels = c("No detectado", "Detectado")) +
  scale_color_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C"),
                     guide = "none") +
  labs(
    title = "Efecto del Tamaño Muestral en la Detección",
    subtitle = "Regresión logística ajustada (IC 95%)",
    x = "Número de Reportes del Triplete (escala log)",
    y = "Estado de Detección"
  ) +
  theme_minimal(base_size = 12)

print(p4)
ggsave("validation_sample_size_effect.png", p4, width = 10, height = 6, dpi = 300)

# Gráfico 5: Número de etapas significativas
p5 <- results_dt[model_success == TRUE] %>%
  ggplot(aes(x = n_stages_significant, fill = type)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  labs(
    title = "Distribución de Etapas con Señal Significativa",
    subtitle = "Por tipo de triplete (solo modelos exitosos)",
    x = "Número de Etapas NICHD Significativas", 
    y = "Frecuencia",
    fill = "Tipo"
  ) +
  scale_fill_manual(values = c("positive" = "#4DAF4A", "negative" = "#E41A1C")) +
  scale_x_continuous(breaks = 0:7) +
  theme_minimal(base_size = 12)

print(p5)
ggsave("validation_stages_distribution.png", p5, width = 10, height = 6, dpi = 300)

# ---------- TESTS ESTADÍSTICOS ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("TESTS ESTADÍSTICOS")
message(paste(rep("=", 70), collapse = ""))

# Test 1: Mann-Whitney para IOR máximo
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
  
  W = %.2f, p-valor = %.2e
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

# Test 2: Chi-cuadrado para proporción de detección
test_prop <- prop.test(
  x = c(tp, fp),
  n = c(tp + fn, fp + tn),
  alternative = "greater"
)

cat(sprintf("
Test de Proporciones para Tasa de Detección:
  H0: prop_positivos <= prop_negativos
  Ha: prop_positivos > prop_negativos
  
  X² = %.2f, p-valor = %.2e
  %s
  
  Proporciones:
    - Positivos: %.3f (%d/%d)
    - Negativos: %.3f (%d/%d)
",
test_prop$statistic,
test_prop$p.value,
if (test_prop$p.value < 0.05) "✓ RECHAZAMOS H0: Positivos se detectan más" 
else "✗ No hay diferencia significativa",
tp / (tp + fn), tp, tp + fn,
fp / (fp + tn), fp, fp + tn
))

# ---------- GUARDAR RESULTADOS ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("GUARDANDO RESULTADOS")
message(paste(rep("=", 70), collapse = ""))

# Resultados detallados
fwrite(results_dt, "validation_detailed_results.csv")
message("✓ Resultados detallados: validation_detailed_results.csv")

# Resumen de métricas
summary_metrics <- data.table(
  metric = c("Sensitivity_TPR", "Specificity_TNR", "PPV_Precision", 
             "NPV", "Accuracy", "F1_Score", "AUC"),
  value = c(sensitivity, specificity, ppv, npv, accuracy, f1_score, auc_value),
  interpretation = c(
    sprintf("%d/%d positivos detectados", tp, tp+fn),
    sprintf("%d/%d negativos correctos", tn, tn+fp),
    sprintf("%d/%d predicciones+ correctas", tp, tp+fp),
    sprintf("%d/%d predicciones- correctas", tn, tn+fn),
    sprintf("%d/%d tripletes correctos", tp+tn, tp+tn+fp+fn),
    "Media armónica Prec-Recall",
    "Área bajo curva ROC"
  )
)
fwrite(summary_metrics, "validation_summary_metrics.csv")
message("✓ Resumen de métricas: validation_summary_metrics.csv")

# Performance por tamaño muestral
fwrite(perf_by_size, "validation_performance_by_sample_size.csv")
message("✓ Performance por tamaño: validation_performance_by_sample_size.csv")

# Exportar IORs completos
ior_export <- results_dt[, c("triplet_id", "drugA", "drugB", "meddra", 
                              "type", "signal_detected", "n_stages_significant",
                              paste0("ior_stage_", 1:7)), with = FALSE]
fwrite(ior_export, "validation_ior_values_by_stage.csv")
message("✓ Valores IOR por etapa: validation_ior_values_by_stage.csv")

# Matriz de confusión
confusion_matrix <- data.table(
  Truth = c("Positive", "Positive", "Negative", "Negative"),
  Prediction = c("Signal", "No Signal", "Signal", "No Signal"),
  Count = c(tp, fn, fp, tn),
  Proportion = c(tp/(tp+fn), fn/(tp+fn), fp/(fp+tn), tn/(fp+tn))
)
fwrite(confusion_matrix, "validation_confusion_matrix.csv")
message("✓ Matriz de confusión: validation_confusion_matrix.csv")

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

# ---------- ANÁLISIS DE PODER ESTADÍSTICO ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("ANÁLISIS DE PODER ESTADÍSTICO")
message(paste(rep("=", 70), collapse = ""))

# Power analysis: combinación de dinámica y tamaño
power_analysis <- results_with_dynamic[
  type == "positive" & !is.na(dynamic) & !is.na(sample_size_cat), 
  .(
    n = .N,
    power = mean(signal_detected, na.rm = TRUE),
    mean_n_events = mean(n_events, na.rm = TRUE),
    mean_n_coadmin = mean(n_coadmin, na.rm = TRUE),
    mean_n_injected = mean(n_injected, na.rm = TRUE),
    mean_max_ior = mean(max_ior, na.rm = TRUE)
  ), 
  by = .(dynamic, sample_size_cat)
]

cat("\nPoder de Detección por Dinámica y Tamaño Muestral:\n")
print(power_analysis[order(dynamic, sample_size_cat)])

fwrite(power_analysis, "validation_power_analysis.csv")
message("\n✓ Análisis de poder: validation_power_analysis.csv")

# Visualización de poder
if (nrow(power_analysis) > 0) {
  p7 <- ggplot(power_analysis[!is.na(sample_size_cat)], 
               aes(x = sample_size_cat, y = power, 
                   group = dynamic, color = dynamic)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_text(aes(label = n), vjust = -1, size = 2.5, color = "black") +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
    labs(
      title = "Poder Estadístico por Dinámica y Tamaño Muestral",
      subtitle = "Números = cantidad de tripletes en cada celda",
      x = "Rango de Reportes", y = "Poder (Tasa de Detección)",
      color = "Dinámica"
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
  
  print(p7)
  ggsave("validation_power_analysis.png", p7, width = 11, height = 6, dpi = 300)
}

# ---------- EJEMPLOS DE PATRONES IOR ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("SELECCIONANDO EJEMPLOS REPRESENTATIVOS")
message(paste(rep("=", 70), collapse = ""))

# Seleccionar ejemplos: mejores detecciones por tipo
ejemplos_positivos <- results_with_dynamic[
  type == "positive" & signal_detected == TRUE & !is.na(dynamic),
  .SD[which.max(n_stages_significant)],
  by = dynamic
]

ejemplos_falsos_pos <- results_dt[
  type == "negative" & signal_detected == TRUE,
  .SD[which.max(max_ior)]
]

ejemplos_falsos_neg <- results_with_dynamic[
  type == "positive" & signal_detected == FALSE & !is.na(dynamic),
  .SD[which.max(max_ior)]
]

if (nrow(ejemplos_positivos) > 0) {
  cat(sprintf("\nEjemplos de Verdaderos Positivos detectados: %d\n", 
              nrow(ejemplos_positivos)))
  print(ejemplos_positivos[, .(dynamic, n_stages_significant, max_ior, N)])
}

if (nrow(ejemplos_falsos_pos) > 0) {
  cat(sprintf("\nEjemplo de Falso Positivo: max_ior = %.2f\n", 
              ejemplos_falsos_pos$max_ior))
}

if (nrow(ejemplos_falsos_neg) > 0) {
  cat(sprintf("\nEjemplo de Falso Negativo: max_ior = %.2f\n", 
              ejemplos_falsos_neg$max_ior))
}

# Graficar ejemplos
if (nrow(ejemplos_positivos) > 0) {
  plot_data_list <- list()
  
  for (i in seq_len(min(nrow(ejemplos_positivos), 6))) {
    ej <- ejemplos_positivos[i]
    ior_vals <- unlist(results_list[[ej$triplet_id]]$ior_values)
    ior_li_vals <- unlist(results_list[[ej$triplet_id]]$ior_li)
    ior_ls_vals <- unlist(results_list[[ej$triplet_id]]$ior_ls)
    
    plot_data_list[[i]] <- data.table(
      example_id = i,
      dynamic = ej$dynamic,
      stage = 1:7,
      stage_name = factor(niveles_nichd, levels = niveles_nichd),
      ior = ior_vals,
      ior_li = ior_li_vals,
      ior_ls = ior_ls_vals,
      n_stages_sig = ej$n_stages_significant,
      max_ior = ej$max_ior
    )
  }
  
  plot_data <- rbindlist(plot_data_list)
  
  p8 <- ggplot(plot_data, aes(x = stage)) +
    geom_ribbon(aes(ymin = ior_li, ymax = ior_ls), 
                alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = ior), linewidth = 1.2, color = "steelblue") +
    geom_point(aes(y = ior), size = 2.5, color = "steelblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
    scale_x_continuous(breaks = 1:7, labels = niveles_nichd) +
    facet_wrap(~ paste0(dynamic, "\n(", n_stages_sig, " etapas sig, max=", 
                        round(max_ior, 1), ")"), 
               scales = "free_y", ncol = 3) +
    labs(
      title = "Ejemplos de Señales Detectadas (Verdaderos Positivos)",
      subtitle = "Área sombreada = IC 90%; línea roja = sin interacción (IOR=1)",
      x = "Etapa NICHD", y = "IOR (escala log)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 9)
    )
  
  print(p8)
  ggsave("validation_example_patterns.png", p8, 
         width = 12, height = 8, dpi = 300)
}

# ---------- REPORTE FINAL ----------
message("\n", paste(rep("=", 70), collapse = ""))
message("RESUMEN EJECUTIVO - VALIDACIÓN ESTILO GIANGRECO")
message(paste(rep("=", 70), collapse = ""))

cat(sprintf("
╔════════════════════════════════════════════════════════════════════╗
║           VALIDACIÓN DEL MODELO GAM DIFERENCIAL                    ║
║              Detección de Señales de Interacción                   ║
╚════════════════════════════════════════════════════════════════════╝

DATASET:
  • Total de tripletes analizados:     %d
  • Tripletes positivos (con señal):   %d
  • Tripletes negativos (sin señal):   %d
  • Modelos exitosos:                  %d (%.1f%%)

MÉTRICAS PRINCIPALES (estilo Giangreco et al. 2022):
  • Sensitivity (TPR):                 %.3f
  • Specificity (TNR):                 %.3f
  • PPV (Precision):                   %.3f
  • NPV:                               %.3f
  • AUC:                               %.3f
  • F1-Score:                          %.3f

INTERPRETACIÓN:
  %s
  %s
  %s

EFECTO DEL TAMAÑO MUESTRAL:
  • Rango 50-100:   %.1f%% detección (%d tripletes)
  • Rango 100-250:  %.1f%% detección (%d tripletes)
  • Rango 250-500:  %.1f%% detección (%d tripletes)
  • Rango >500:     %.1f%% detección (%d tripletes)

TIPOS DE ERROR:
  • Falsos Negativos (no detectados): %d / %d (%.1f%%)
  • Falsos Positivos (señal espuria): %d / %d (%.1f%%)

ARCHIVOS GENERADOS:
  ✓ validation_detailed_results.csv
  ✓ validation_summary_metrics.csv
  ✓ validation_confusion_matrix.csv
  ✓ validation_ior_values_by_stage.csv
  ✓ validation_performance_by_sample_size.csv
  ✓ validation_performance_by_dynamic.csv (exploratorio)
  ✓ validation_power_analysis.csv
  
  ✓ validation_roc_curve.png
  ✓ validation_ior_distribution.png
  ✓ validation_detection_by_sample_size.png
  ✓ validation_sample_size_effect.png
  ✓ validation_stages_distribution.png
  ✓ validation_detection_by_dynamic_exploratory.png
  ✓ validation_power_analysis.png
  ✓ validation_example_patterns.png

",
nrow(ground_truth),
nrow(positivos_res),
nrow(negativos_res),
sum(results_dt$model_success),
100 * mean(results_dt$model_success),
sensitivity,
specificity,
ppv,
npv,
auc_value,
f1_score,
if (sensitivity >= 0.8) "✓ EXCELENTE sensibilidad: el método detecta >80% de las interacciones reales"
else if (sensitivity >= 0.6) "○ BUENA sensibilidad: detecta mayoría de interacciones, pero hay margen de mejora"
else "✗ BAJA sensibilidad: muchas interacciones reales no son detectadas",
if (ppv >= 0.8) "✓ ALTA precisión: pocas falsas alarmas, las señales detectadas son confiables"
else if (ppv >= 0.6) "○ PRECISIÓN MODERADA: algunas falsas alarmas, validación adicional recomendada"
else "✗ BAJA precisión: muchas señales detectadas son falsas alarmas",
if (auc_value >= 0.85) "✓ EXCELENTE discriminación: AUC >0.85 indica separación clara positivos/negativos"
else if (auc_value >= 0.7) "○ BUENA discriminación: AUC >0.7 indica método útil"
else "✗ DISCRIMINACIÓN INSUFICIENTE: AUC <0.7, revisar parámetros",
if (nrow(perf_by_size) >= 1) perf_by_size[1]$detection_rate * 100 else NA,
if (nrow(perf_by_size) >= 1) perf_by_size[1]$n else 0,
if (nrow(perf_by_size) >= 2) perf_by_size[2]$detection_rate * 100 else NA,
if (nrow(perf_by_size) >= 2) perf_by_size[2]$n else 0,
if (nrow(perf_by_size) >= 3) perf_by_size[3]$detection_rate * 100 else NA,
if (nrow(perf_by_size) >= 3) perf_by_size[3]$n else 0,
if (nrow(perf_by_size) >= 4) perf_by_size[4]$detection_rate * 100 else NA,
if (nrow(perf_by_size) >= 4) perf_by_size[4]$n else 0,
fn, tp + fn, 100 * fn / (tp + fn),
fp, fp + tn, 100 * fp / (fp + tn)
))

# Recomendaciones basadas en resultados
cat("\nRECOMENDACIONES:\n")
if (sensitivity < 0.7) {
  cat("  ⚠ Sensibilidad baja:\n")
  cat("    - Considerar aumentar n_pos o effect_size en augmentation\n")
  cat("    - Revisar min_reports_triplet (quizá demasiado alto)\n")
  cat("    - Evaluar si min_stages_significant es muy estricto\n\n")
}
if (ppv < 0.7) {
  cat("  ⚠ Precisión baja (muchos falsos positivos):\n")
  cat("    - Considerar criterios más estrictos (ej: min_stages_significant > 1)\n")
  cat("    - Agregar umbral de magnitud mínima (min_ior_threshold)\n")
  cat("    - Revisar negativos: ¿son realmente no-interacciones?\n\n")
}
if (sensitivity >= 0.7 && ppv >= 0.7) {
  cat("  ✓ Desempeño equilibrado y satisfactorio\n")
  cat("  ✓ El método es adecuado para screening de interacciones\n")
  cat("  ✓ Las señales detectadas requieren validación clínica\n\n")
}

cat("COMPARACIÓN CON GIANGRECO ET AL. (2022):\n")
cat("  • Giangreco reporta Precision ~0.56, Recall ~0.31 para dGAMs\n")
cat("  • Nuestro método:\n")
cat(sprintf("    - Precision (PPV):  %.3f\n", ppv))
cat(sprintf("    - Recall (TPR):     %.3f\n", sensitivity))
cat("  • Contexto: Giangreco detecta efectos ontogénicos de drogas individuales\n")
cat("    mientras que aquí detectamos interacciones entre pares de drogas\n\n")

# ---------- GUARDAR WORKSPACE ----------
save.image("validation_workspace.RData")
message("✓ Workspace guardado: validation_workspace.RData")
message("  (Cargar con: load('validation_workspace.RData'))")

message("\n╔════════════════════════════════════════════════════════════════════╗")
message("║                   ¡VALIDACIÓN COMPLETADA! 🎉                       ║")
message("╚════════════════════════════════════════════════════════════════════╝")
message("\nTiempo de finalización: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("\nRevisa los archivos CSV y PNG generados para análisis detallado.\n")