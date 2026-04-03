################################################################################
# Script de diagnóstico de supuestos del modelo GAM
# Script: 50_supuestos.R
#
# Objetivo: evaluar si los estimados GAM-IOR y GAM-RERI (y sus IC90) son
# válidos para cada triplete, identificando las características de reporte
# que determinan la viabilidad de la estimación.
#
# Estructura:
#   PARTE 1 — Criterios de usabilidad estructural (post-ajuste)
#     1a. Convergencia PIRLS
#     1b. Coeficientes finitos
#     1c. Estimabilidad del contraste (4 celdas con observaciones)
#     1d. n_events_coadmin > 0
#
#   PARTE 2 — Informatividad del triplete (pre-ajuste, desde datos)
#     2a. n_coadmin: reportes A-B totales
#     2b. n_events_coadmin: reportes A-B-E
#     2c. n_stages_with_AB: etapas NICHD con ≥1 reporte A-B
#     2d. n_stages_AB_event: etapas con ≥1 reporte A-B-E  ← señal del spline
#     2e. n_stages_bilateral_AB: etapas donde AB tiene ambos outcomes (0 y 1)
#     2f. prop_event_in_AB: tasa bruta de evento en grupo AB
#     2g. min_coadmin_stage: mínimo de reportes AB en etapas con datos
#
#   PARTE 3 — Calidad del spline de interacción (post-ajuste)
#     3a. EDF ratio  [0 = spline colapsado, 1 = saturado]
#     3b. p_interaction del spline
#
#   PARTE 4 — Calibración local en grupo coadministración (post-ajuste)
#     Compara tasas observadas vs predichas dentro del grupo droga_ab == 1
#     por etapa NICHD. MAE y max |residuo| agrupados.
#     Este es el chequeo de ajuste relevante: IOR/RERI depende del ajuste
#     en ese subgrupo, no del dataset completo.
#
# Supuestos fijos por diseño (no se evalúan):
#   - Familia binomial(logit): outcome 0/1 por construcción
#   - k = 7 knots: saturado por las 7 etapas NICHD
#   - Base "cs": shrinkage hacia cero (regularización incorporada)
#   - Independencia entre reportes: limitación estructural de FAERS
#
# Nota sobre DHARMa:
#   No se usa en este script. El KS global sobre ~264k filas Bernoulli
#   tiene poder casi infinito: cualquier desviación mínima del modelo
#   produce p < 0.05 independientemente de la validez práctica. El test
#   de dispersión es inapropiado para observaciones Bernoulli individuales
#   (la varianza está fijada por p_i(1-p_i) por definición). La calibración
#   local en el grupo AB (Parte 4) provee información más directa y útil.
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuración
################################################################################

ruta_candidatos <- "./results/sics/augmentation_results/positive_triplets_results.rds"
output_dir      <- "./results/diagnostics/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_sample    <- 100   # tripletes a evaluar
min_coadmin <- 5     # mínimo reportes A-B para incluir un triplete

# Parámetros del modelo — deben coincidir con el ajuste principal
spline_individuales <- TRUE
include_nichd       <- FALSE
nichd_spline        <- FALSE
bs_type             <- "cs"
k_spline            <- 7
select              <- FALSE
include_sex         <- FALSE
include_stage_sex   <- FALSE
method              <- "fREML"

################################################################################
# PARTE 2 — Función de informatividad del triplete (sin ajustar modelo)
################################################################################

# Computa las características de datos del triplete que determinan si el
# contraste IOR/RERI es identificable. Se puede llamar antes de ajustar.
#
# Argumentos clave:
#   r_ab: vector de safetyreportid con co-administración A-B
#   r_e:  vector de safetyreportid con el evento E
#   dm:   data.table con columnas (safetyreportid, nichd_num, droga_a, droga_b,
#          droga_ab, ea_ocurrio)

triplet_informatividad <- function(r_ab, r_e, dm) {

  n_coadmin        <- length(r_ab)
  n_events_coadmin <- length(intersect(r_ab, r_e))

  # Subconjunto AB
  dm_ab <- dm[droga_ab == 1]

  # Por etapa en grupo AB
  stage_ab <- dm_ab[, .(
    n_ab       = .N,
    n_ab_event = sum(ea_ocurrio),
    bilateral  = as.integer(var(ea_ocurrio) > 0 & .N >= 2)
  ), by = nichd_num]

  n_stages_with_AB       <- nrow(stage_ab)
  n_stages_AB_event      <- sum(stage_ab$n_ab_event > 0)
  n_stages_bilateral_AB  <- sum(stage_ab$bilateral,  na.rm = TRUE)
  prop_event_in_AB       <- if (n_coadmin > 0) n_events_coadmin / n_coadmin else NA_real_
  min_coadmin_stage      <- if (n_stages_with_AB > 0) min(stage_ab$n_ab) else NA_integer_

  list(
    n_coadmin            = n_coadmin,
    n_events_coadmin     = n_events_coadmin,
    n_stages_with_AB     = n_stages_with_AB,
    n_stages_AB_event    = n_stages_AB_event,     # clave: señal del spline
    n_stages_bilateral_AB = n_stages_bilateral_AB, # clave: separación local
    prop_event_in_AB     = prop_event_in_AB,
    min_coadmin_stage    = min_coadmin_stage
  )
}

################################################################################
# Función principal: ajuste + diagnósticos Partes 1, 3 y 4
################################################################################

# Ajusta el GAM para un triplete y extrae todos los diagnósticos.
# Devuelve una lista con: modelo ajustado, dm, y tabla de resultados.

fit_gam_diagnostics <- function(drugA_id, drugB_id, event_id, ade_data,
                                spline_individuales = TRUE,
                                include_nichd       = FALSE,
                                nichd_spline        = FALSE,
                                bs_type             = "cs",
                                k_spline            = 7,
                                select              = FALSE,
                                include_sex         = FALSE,
                                include_stage_sex   = FALSE,
                                method              = "fREML") {

  triplet_label <- paste(drugA_id, drugB_id, event_id, sep = "_")

  # ── Preparar datos ──────────────────────────────────────────────────────────
  r_a  <- unique(ade_data[atc_concept_id == drugA_id,    safetyreportid])
  r_b  <- unique(ade_data[atc_concept_id == drugB_id,    safetyreportid])
  r_e  <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  r_ab <- intersect(r_a, r_b)

  if (length(r_ab) < 2) {
    return(list(modelo = NULL, dm = NULL, result = data.table(
      triplet          = triplet_label,
      usable           = FALSE,
      usability_reason = "insufficient_coadmin")))
  }

  dm <- unique(ade_data[, .(safetyreportid, nichd_num)])
  dm[, `:=`(
    ea_ocurrio = as.integer(safetyreportid %in% r_e),
    droga_a    = as.integer(safetyreportid %in% r_a),
    droga_b    = as.integer(safetyreportid %in% r_b)
  )]
  dm[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]

  # ── PARTE 2: Informatividad (antes de ajustar) ───────────────────────────────
  info <- triplet_informatividad(r_ab, r_e, dm)

  # ── PARTE 1a: n_events_coadmin > 0 ─────────────────────────────────────────
  if (info$n_events_coadmin == 0) {
    return(list(modelo = NULL, dm = NULL, result = data.table(
      triplet          = triplet_label,
      usable           = FALSE,
      usability_reason = "zero_events_in_coadmin",
      as.data.table(info))))
  }

  # ── PARTE 1b: n_stages_AB_event >= 1 ────────────────────────────────────────
  # Sin ninguna etapa donde AB y el evento coexistan, el spline de interacción
  # no tiene señal local: los estimados IOR/RERI serán pura extrapolación.
  if (info$n_stages_AB_event == 0) {
    return(list(modelo = NULL, dm = NULL, result = data.table(
      triplet          = triplet_label,
      usable           = FALSE,
      usability_reason = "no_stage_with_AB_and_event",
      as.data.table(info))))
  }

  # ── PARTE 1c: estimabilidad del contraste IOR (4 celdas) ────────────────────
  cell_counts   <- dm[, .(n = .N), by = .(droga_a, droga_b)]
  missing_cells <- character(0)
  for (da in c(0, 1)) for (db in c(0, 1)) {
    if (nrow(cell_counts[droga_a == da & droga_b == db]) == 0)
      missing_cells <- c(missing_cells, paste0("(", da, ",", db, ")"))
  }
  cells_estimable <- length(missing_cells) == 0

  if (!cells_estimable) {
    return(list(modelo = NULL, dm = NULL, result = data.table(
      triplet          = triplet_label,
      usable           = FALSE,
      usability_reason = paste("missing_cells:", paste(missing_cells, collapse = ";")),
      as.data.table(info))))
  }

  # ── Fórmula (idéntica a fit_gam en el pipeline principal) ───────────────────
  fp <- "ea_ocurrio ~ "
  if (!spline_individuales) {
    fp <- paste0(fp, "droga_a + droga_b + ")
  } else {
    fp <- paste0(fp,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", k_spline, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", k_spline, bs_type))
  }
  if (include_nichd) {
    fp <- if (nichd_spline) {
      paste0(fp, sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline, bs_type))
    } else paste0(fp, "nichd_num + ")
  }
  fp <- paste0(fp,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline, bs_type))

  # ── Ajuste ──────────────────────────────────────────────────────────────────
  modelo <- tryCatch(
    bam(as.formula(fp), data = dm, family = binomial(link = "logit"),
        method = method, select = select, discrete = TRUE, nthreads = 1),
    error = function(e) NULL
  )

  if (is.null(modelo)) {
    return(list(modelo = NULL, dm = NULL, result = data.table(
      triplet          = triplet_label,
      usable           = FALSE,
      usability_reason = "fit_failed",
      as.data.table(info))))
  }

  # ── PARTE 1d: convergencia ──────────────────────────────────────────────────
  converged <- isTRUE(modelo$converged)

  # ── PARTE 1e: coeficientes finitos ──────────────────────────────────────────
  coefs        <- coef(modelo)
  coefs_finite <- all(is.finite(coefs))
  n_nonfinite  <- sum(!is.finite(coefs))

  usable <- converged & cells_estimable & (info$n_events_coadmin > 0) &
            coefs_finite & (info$n_stages_AB_event >= 1)

  usability_reason <- if (usable) "OK" else paste(
    c(if (!converged)                      "no_convergence",
      if (!coefs_finite)                   paste0("nonfinite_coefs(n=", n_nonfinite, ")"),
      if (!cells_estimable)                "missing_cells",
      if (info$n_stages_AB_event < 1)      "no_stage_with_AB_event"),
    collapse = ";"
  )

  # ── PARTE 3: EDF del spline de interacción ───────────────────────────────────
  # EDF ~ 0 (shrinkage total): el GAM no detectó variación ontogénica en AB
  # EDF = k-1:  spline usando toda la flexibilidad disponible
  edf_all         <- summary(modelo)$s.table
  idx_int         <- grep("droga_ab", rownames(edf_all))
  edf_interaction <- if (length(idx_int) > 0) edf_all[idx_int, "edf"]     else NA_real_
  edf_ratio       <- edf_interaction / (k_spline - 1)
  p_interaction   <- if (length(idx_int) > 0) edf_all[idx_int, "p-value"] else NA_real_

  # ── PARTE 4: Calibración local en grupo coadministración ────────────────────
  # Compara tasas observadas vs predichas dentro de droga_ab == 1 por etapa.
  # Este es el subgrupo relevante para IOR/RERI: si el modelo no ajusta bien
  # las probabilidades en el grupo AB, los estimados de interacción serán
  # sesgados independientemente del ajuste global.
  dm[, pred_prob := predict(modelo, type = "response")]

  calib_ab <- dm[droga_ab == 1, .(
    obs_rate  = mean(ea_ocurrio),
    pred_rate = mean(pred_prob),
    n_ab      = .N
  ), by = nichd_num]
  calib_ab[, resid_ab := obs_rate - pred_rate]

  mae_ab     <- mean(abs(calib_ab$resid_ab), na.rm = TRUE)
  max_resid_ab <- max(abs(calib_ab$resid_ab), na.rm = TRUE)

  # Calibración global (referencia)
  calib_global <- dm[, .(
    obs_rate  = mean(ea_ocurrio),
    pred_rate = mean(pred_prob)
  ), by = nichd_num]
  calib_global[, resid := obs_rate - pred_rate]
  mae_global <- mean(abs(calib_global$resid), na.rm = TRUE)

  list(
    modelo = modelo,
    dm     = dm,
    result = data.table(
      triplet                = triplet_label,
      drugA                  = drugA_id,
      drugB                  = drugB_id,
      event                  = event_id,

      # Parte 1: Usabilidad estructural
      usable                 = usable,
      usability_reason       = usability_reason,
      converged              = converged,
      coefs_finite           = coefs_finite,
      n_nonfinite_coefs      = n_nonfinite,
      cells_estimable        = cells_estimable,

      # Parte 2: Informatividad del triplete
      n_coadmin              = info$n_coadmin,
      n_events_coadmin       = info$n_events_coadmin,
      n_stages_with_AB       = info$n_stages_with_AB,
      n_stages_AB_event      = info$n_stages_AB_event,    # clave para IOR/RERI
      n_stages_bilateral_AB  = info$n_stages_bilateral_AB, # clave para separación
      prop_event_in_AB       = info$prop_event_in_AB,
      min_coadmin_stage      = info$min_coadmin_stage,

      # Parte 3: Calidad del spline de interacción
      edf_interaction        = edf_interaction,
      edf_max                = k_spline - 1L,
      edf_ratio              = edf_ratio,
      p_interaction          = p_interaction,

      # Parte 4: Calibración local (grupo AB)
      mae_ab                 = mae_ab,        # calibración donde importa
      max_resid_ab           = max_resid_ab,
      mae_global             = mae_global     # referencia
    )
  )
}

################################################################################
# Carga de datos
################################################################################

message("Cargando datos...")
ade_raw_dt <- fread(ruta_ade_raw)

cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")
ade_raw_dt <- ade_raw_dt[, ..cols_req]

if (!"nichd_num" %in% names(ade_raw_dt)) {
  ade_raw_dt[, nichd_num := match(nichd, niveles_nichd)]
}
ade_raw_dt <- ade_raw_dt[!is.na(nichd_num)]

message("Cargando tripletes candidatos...")
candidatos     <- readRDS(ruta_candidatos)
candidatos_ok  <- candidatos[injection_success == TRUE & n_coadmin >= min_coadmin]
idx_sample     <- sample(nrow(candidatos_ok), min(n_sample, nrow(candidatos_ok)))
triplets_sample <- candidatos_ok[idx_sample]

message(sprintf("Evaluando %d tripletes...", nrow(triplets_sample)))

################################################################################
# Loop principal
################################################################################

all_results <- pblapply(seq_len(nrow(triplets_sample)), function(i) {

  row <- triplets_sample[i]

  tryCatch(
    fit_gam_diagnostics(
      drugA_id            = row$drugA,
      drugB_id            = row$drugB,
      event_id            = row$meddra,
      ade_data            = ade_raw_dt,
      spline_individuales = spline_individuales,
      include_nichd       = include_nichd,
      nichd_spline        = nichd_spline,
      bs_type             = bs_type,
      k_spline            = k_spline,
      select              = select,
      include_sex         = include_sex,
      include_stage_sex   = include_stage_sex,
      method              = method
    )$result,
    error = function(e) data.table(
      triplet          = paste(row$drugA, row$drugB, row$meddra, sep = "_"),
      usable           = FALSE,
      usability_reason = paste("error:", e$message)
    )
  )
})

diag_results <- rbindlist(all_results, fill = TRUE)

fwrite(diag_results, file.path(output_dir, "gam_diagnostics.csv"))
message(sprintf("Resultados guardados en %s", output_dir))

################################################################################
# Resúmenes numéricos
################################################################################

cat("\n================================================\n")
cat("DIAGNÓSTICO DE SUPUESTOS GAM (IOR/RERI)\n")
cat("================================================\n\n")

cat("── PARTE 1: Criterios de usabilidad estructural ──\n")
cat(sprintf("  Usables:                      %d / %d (%.1f%%)\n",
            sum(diag_results$usable, na.rm = TRUE), nrow(diag_results),
            100 * mean(diag_results$usable, na.rm = TRUE)))
cat(sprintf("  Convergencia:                 %d / %d\n",
            sum(diag_results$converged, na.rm = TRUE), nrow(diag_results)))
cat(sprintf("  Coeficientes finitos:         %d / %d\n",
            sum(diag_results$coefs_finite, na.rm = TRUE), nrow(diag_results)))
cat(sprintf("  4 celdas estimables:          %d / %d\n",
            sum(diag_results$cells_estimable, na.rm = TRUE), nrow(diag_results)))
cat(sprintf("  n_events_coadmin > 0:         %d / %d\n",
            sum(diag_results$n_events_coadmin > 0, na.rm = TRUE), nrow(diag_results)))
cat(sprintf("  n_stages_AB_event >= 1:       %d / %d\n",
            sum(diag_results$n_stages_AB_event >= 1, na.rm = TRUE), nrow(diag_results)))

usable_dt <- diag_results[usable == TRUE]
n_usable  <- nrow(usable_dt)

cat(sprintf("\n── PARTE 2: Informatividad del triplete  (n usables = %d) ──\n", n_usable))
cat(sprintf("  n_coadmin — mediana:          %.0f  (rango: %d–%d)\n",
            median(usable_dt$n_coadmin, na.rm = TRUE),
            min(usable_dt$n_coadmin,    na.rm = TRUE),
            max(usable_dt$n_coadmin,    na.rm = TRUE)))
cat(sprintf("  n_events_coadmin — mediana:   %.0f\n",
            median(usable_dt$n_events_coadmin, na.rm = TRUE)))
cat(sprintf("  n_stages_with_AB — mediana:   %.1f / 7\n",
            median(usable_dt$n_stages_with_AB, na.rm = TRUE)))
cat(sprintf("  n_stages_AB_event — mediana:  %.1f / 7\n",
            median(usable_dt$n_stages_AB_event, na.rm = TRUE)))
cat(sprintf("  n_stages_bilateral_AB — med.: %.1f / 7\n",
            median(usable_dt$n_stages_bilateral_AB, na.rm = TRUE)))
cat(sprintf("  Tripletes con n_stages_AB_event == 1:  %d (%.1f%%)\n",
            sum(usable_dt$n_stages_AB_event == 1, na.rm = TRUE),
            100 * mean(usable_dt$n_stages_AB_event == 1, na.rm = TRUE)))

cat("\n── PARTE 3: Calidad del spline de interacción ──\n")
cat(sprintf("  EDF mediana:                  %.3f  (máx teórico: %d)\n",
            median(usable_dt$edf_interaction, na.rm = TRUE), k_spline - 1))
cat(sprintf("  EDF colapsado (ratio < 0.05): %d (%.1f%%)\n",
            sum(usable_dt$edf_ratio < 0.05, na.rm = TRUE),
            100 * mean(usable_dt$edf_ratio < 0.05, na.rm = TRUE)))
cat(sprintf("  EDF no-lineal (ratio >= 0.2): %d (%.1f%%)\n",
            sum(usable_dt$edf_ratio >= 0.2, na.rm = TRUE),
            100 * mean(usable_dt$edf_ratio >= 0.2, na.rm = TRUE)))
cat(sprintf("  p_interaction < 0.05:         %d (%.1f%%)\n",
            sum(usable_dt$p_interaction < 0.05, na.rm = TRUE),
            100 * mean(usable_dt$p_interaction < 0.05, na.rm = TRUE)))

cat("\n── PARTE 4: Calibración local en grupo coadministración ──\n")
cat(sprintf("  MAE_AB mediana:               %.5f\n",
            median(usable_dt$mae_ab, na.rm = TRUE)))
cat(sprintf("  MAE_AB máximo:                %.5f\n",
            max(usable_dt$mae_ab, na.rm = TRUE)))
cat(sprintf("  MAE_global mediana:           %.5f\n",
            median(usable_dt$mae_global, na.rm = TRUE)))
cat(sprintf("  Tripletes MAE_AB > 0.10:      %d (%.1f%%)\n",
            sum(usable_dt$mae_ab > 0.10, na.rm = TRUE),
            100 * mean(usable_dt$mae_ab > 0.10, na.rm = TRUE)))

################################################################################
# Tabla de tripletes flaggeados
#
# Criterios:
#   - No usables: fallo en criterios estructurales
#   - Spline colapsado (edf_ratio < 0.05): estimados IOR/RERI son ~constantes
#     en todas las etapas; no hay dinámica ontogénica detectada
#   - Calibración pobre en AB (mae_ab > 0.10): el modelo mal-especifica las
#     probabilidades en el subgrupo que determina el contraste de interacción
#   - n_stages_AB_event == 1: señal basada en una sola etapa; alta dependencia
#     del suavizado hacia etapas adyacentes (mayor incertidumbre ontogénica)
################################################################################

flaggeados <- diag_results[
  usable == FALSE |
  (!is.na(edf_ratio)      & edf_ratio < 0.05)       |
  (!is.na(mae_ab)         & mae_ab    > 0.10)        |
  (!is.na(n_stages_AB_event) & n_stages_AB_event == 1 & usable == TRUE)
][order(usable, mae_ab)]

fwrite(
  flaggeados[, .(triplet, usable, usability_reason,
                 n_coadmin, n_events_coadmin,
                 n_stages_with_AB, n_stages_AB_event, n_stages_bilateral_AB,
                 prop_event_in_AB, min_coadmin_stage,
                 edf_ratio, p_interaction,
                 mae_ab, max_resid_ab, mae_global)],
  file.path(output_dir, "gam_diagnostics_flagged.csv")
)

cat(sprintf("\nTripletes flaggeados: %d\n", nrow(flaggeados)))
cat(sprintf("  No usables:                  %d\n",
            sum(!flaggeados$usable, na.rm = TRUE)))
cat(sprintf("  Spline colapsado:            %d\n",
            sum(!is.na(flaggeados$edf_ratio) & flaggeados$edf_ratio < 0.05, na.rm = TRUE)))
cat(sprintf("  Calibración pobre en AB:     %d\n",
            sum(!is.na(flaggeados$mae_ab) & flaggeados$mae_ab > 0.10, na.rm = TRUE)))
cat(sprintf("  n_stages_AB_event == 1:      %d\n",
            sum(!is.na(flaggeados$n_stages_AB_event) & flaggeados$n_stages_AB_event == 1, na.rm = TRUE)))
cat(sprintf("\nGuardado en: %s\n", file.path(output_dir, "gam_diagnostics_flagged.csv")))

################################################################################
# Visualizaciones
################################################################################

theme_set(theme_minimal(base_size = 12))
library(patchwork)

# P1: Distribución de n_stages_AB_event (el predictor más relevante)
p1 <- usable_dt |>
  mutate(stages_cat = factor(pmin(n_stages_AB_event, 5),
                             labels = c("1","2","3","4","5+"))) |>
  ggplot(aes(x = stages_cat, fill = stages_cat)) +
  geom_bar(color = "white", alpha = 0.85) +
  scale_fill_brewer(palette = "Blues") +
  labs(title   = "Etapas NICHD con reportes A-B-E",
       subtitle = "Determina la señal disponible para el spline de interacción",
       x = "n etapas con AB+evento", y = "N tripletes") +
  theme(legend.position = "none")

# P2: EDF ratio vs n_stages_AB_event
p2 <- usable_dt |>
  filter(!is.na(edf_ratio), !is.na(n_stages_AB_event)) |>
  ggplot(aes(x = n_stages_AB_event, y = edf_ratio,
             size = n_events_coadmin, color = n_events_coadmin)) +
  geom_jitter(alpha = 0.6, width = 0.15) +
  scale_color_viridis_c(trans = "log10") +
  scale_size_continuous(range = c(1, 5)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title  = "Etapas AB+evento vs EDF ratio del spline",
       x = "n etapas con AB+evento", y = "EDF ratio",
       color = "n eventos\ncoadmin", size = "n eventos\ncoadmin")

# P3: MAE_AB vs EDF ratio
p3 <- usable_dt |>
  filter(!is.na(mae_ab), !is.na(edf_ratio)) |>
  ggplot(aes(x = edf_ratio, y = mae_ab)) +
  geom_point(alpha = 0.5, color = "#2c7bb6") +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "red") +
  labs(title   = "EDF ratio vs MAE en grupo AB",
       subtitle = "Línea roja = umbral orientativo 0.10",
       x = "EDF ratio  [0 = colapsado]", y = "MAE obs vs pred (grupo AB)")

# P4: n_coadmin vs n_stages_AB_event
p4 <- usable_dt |>
  filter(!is.na(n_coadmin), !is.na(n_stages_AB_event)) |>
  ggplot(aes(x = n_coadmin, y = n_stages_AB_event)) +
  geom_jitter(alpha = 0.5, color = "#1a9641", height = 0.15) +
  scale_x_log10() +
  geom_vline(xintercept = min_coadmin, linetype = "dashed", color = "red") +
  labs(title    = "N reportes A-B vs etapas con señal",
       subtitle = "Línea roja = umbral mínimo de inclusión",
       x = "N reportes A-B (escala log)", y = "N etapas con AB+evento")

fig <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title   = "Diagnóstico de supuestos GAM (IOR / RERI)",
    caption = sprintf("n = %d tripletes  |  usables: %d (%.0f%%)",
                      nrow(diag_results),
                      n_usable,
                      100 * n_usable / nrow(diag_results))
  )

ggsave(file.path(output_dir, "gam_diagnostics.png"),
       fig, width = 14, height = 10, dpi = 150)

message("Diagnóstico completo.")