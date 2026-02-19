################################################################################
# Script de funciones
# Script: 00_functions.R
# para usar, source("00_functions.R", local = TRUE)
################################################################################

# cargado de librerias con pacman 
library(pacman)
pacman::p_load(tidyverse, mgcv, MASS, data.table, parallel, akima, doRNG)

# ordenado de niveles para el resto de los scripts
niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

Z90 <- qnorm(0.95)  # Cuantil 90% para intervalos de confianza
source("01_theme.R", local = TRUE)

################################################################################
# Función para construir tripletes
################################################################################

# Genera tripletes (drugA, drugB, event) para un reporte individual
#
# Parámetros:
# drug: ids de drogas en el reporte (atc_concept_id)
# event: ids de eventos en el reporte (meddra_concept_id)
# report_id: id del reporte (safetyreportid)
# nichd_stage: etapa NICHD numérica del reporte
# 
# Return:
# data.table con columnas: safetyreportid, drugA, drugB, meddra, nichd_num
# 
# Implementación:
# Genera todas las combinaciones de pares de drogas (si >= 2 drogas)
# Cruza cada par con cada evento del reporte
# Aplica orden: drugA <= drugB (para evitar duplicados por orden inverso)
# Devuelve NULL si el reporte NO tiene >= 2 drogas o eventos

make_triplets <- function(drug, event, report_id, nichd_stage) {
  
  if (length(drug) < 2 || length(event) < 1) return(NULL)
  
  drug <- unique(drug)
  event <- unique(event)
  
  # combinaciones de pares de drogas
  if (length(drug) == 2) {
    combination <- matrix(c(min(drug), max(drug)), nrow = 1)
  } else {
    combination <- t(combn(drug, 2))
    # Orden: min primero
    combination <- t(apply(combination, 1, function(x) c(min(x), max(x))))
  }
  
  n_combination <- nrow(combination)
  n_events <- length(event)
  
  data.table(
    safetyreportid = report_id,
    drugA = rep(combination[,1], times = n_events),
    drugB = rep(combination[,2], times = n_events),
    meddra = rep(event, each = n_combination),
    nichd_num = nichd_stage
  )
}

################################################################################
# Función de tamaño de efecto (fold-change)
################################################################################

# Muestrea fold-changes siguiendo distribución exponencial negativa
#
# parámetros:
# n: Número de fold-changes a generar
# lambda: parámetro de tasa para distribución exponencial (0.75)
# 
# return: vector numérico con fold-changes >= 1
# 
# Implementación:
# - FC ~ 1 + exp(λ = 0.75)
# - rango típico: [1, 10] con sesgo a la derecha

fold_change <- function(n, lambda = 0.75) {
  1 + rexp(n, rate = lambda)
}

################################################################################
# Función auxiliar: calcular coadministración por etapa para un triplete
################################################################################

# Calcula la cantidad de reportes A+B por etapa para cada triplete
# Permite calcular superset para método clásico
# 
# return:
# n_coadmin_stage: número de reportes A+B para evento en etapa
#
# Implementación:
# Identifica reportes A+B
# Realiza conteo de reportes por etapa
# Completa etapas sin reportes con 0

coadmin_by_stage <- function(drugA, drugB, meddra, ade_data) {
  
  # reportes con cada droga
  reports_A <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  reports_B <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  
  # Reportes con ambas drogas (coadministración)
  reports_AB <- intersect(reports_A, reports_B)
  
  if (length(reports_AB) == 0) {
    # Si no hay coadministración, retornar estructura vacía
    return(data.table(
      nichd_num = 1:7,
      nichd = niveles_nichd,
      n_coadmin_stage = 0L
    ))
  }
  
  # etapa NICHD de cada reporte con coadministración
  stage_counts <- unique(ade_data[
    safetyreportid %in% reports_AB,
    .(safetyreportid, nichd, nichd_num)
  ])[, .N, by = .(nichd, nichd_num)]
  
  # ceros para etapas sin reportes
  full_stages <- data.table(
    nichd_num = 1:7,
    nichd = niveles_nichd
  )
  
  stage_counts <- merge(
    full_stages,
    stage_counts,
    by = c("nichd_num", "nichd"),
    all.x = TRUE
  )
  
  stage_counts[is.na(N), N := 0L]
  setnames(stage_counts, "N", "n_coadmin_stage")
  
  return(stage_counts[order(nichd_num)])
}

################################################################################
# Función de generación de dinámicas
################################################################################

# Genera patrones dinámicos normalizados para inyección de señales
#
# parámetros:
# type: Tipo de dinámica: "uniform", "increase", "decrease", "plateau", "inverse_plateau"
# N: Número de etapas (7 para NICHD)
# 
# return: vector numérico normalizado en [-1, 1] con el patrón temporal
# 
# formas de dinámicas
# uniform: señal constante (0 en todas las etapas)
# increase: crecimiento monotónico de -1 a +1
# decrease: decrecimiento monotónico de +1 a -1
# plateau: pico en etapas centrales (forma de campana)
# inverse_plateau: valle en etapas centrales (forma de U invertida)

generate_dynamic <- function(type, N = 7) {
  type <- as.character(type)
  
  if (type == "uniform") {
    return(rep(0, N))
  }
  
  if (type == "increase") {
    return(tanh(seq(-pi, pi, length.out = N)))
  }
  
  if (type == "decrease") {
    return(-tanh(seq(-pi, pi, length.out = N)))
  }
  
  if (type == "plateau") {
    return(c(
      tanh(seq(-pi, pi, length.out = floor(N/2))),
      tanh(seq(pi, -pi, length.out = ceiling(N/2)))
    ))
  }
  
  if (type == "inverse_plateau") {
    return(c(
      tanh(seq(pi, -pi, length.out = floor(N/2))),
      tanh(seq(-pi, pi, length.out = ceiling(N/2)))
    ))
  }
}

################################################################################
# Función de inyección de señales
################################################################################

# Inyecta señal de interacción droga-droga
#
# parámetros:
# drugA_id: id droga A (ATC concept_id)
# drugB_id: id droga B (ATC concept_id)
# event_id: id del evento adverso (MedDRA concept_id)
# dynamic_type: tipo de dinámica (llama generate_dynamic)
# fold_change: magnitud del efecto (llama fold_change)
# ade_raw_dt: data.table de dataset original
# 
# return Lista con: success, n_injected, n_coadmin, ade_aug, diagnostics
# 
# Implementación:
# 
# 1 Tasa base(e_j):
# Valor base para aplicar fold-change
# El valor base debe ser calculado teniendo en cuenta componentes individuales 
# El valor base altera la consistencia de la inyección en distintos escenarios
# Script "simulacion_inyeccion" muestra consistencia en distintos escenarios
# Los métodos que asumen "independencia" de los eventos solo se usan como proxy
# no suponen independencia real 
#
# 2 Fold_change (FC):
# Factor multiplicativo del efecto
# 
# 3 Dinámica f(j):
# Función normalizada en [0, 1] por etapa
# f(j) = llama a generate_dynamic(type, N=7)  <--- N siempre 7 por etapas nichd 
# 
# 4 Probabilidad dinámica: (formula final de probabilidad de reporte)
# p_dynamic(j) = e_j × FC + f(j)
# 
# 5 Simulación:
# Para cada reporte con drugA + drugB
# Y_new ~ Bernoulli(p_dynamic(stage_j))
# Inyecta SOLO en eventos que no existen (si ya hay evento, se deja, trata de modelar dinámica sobre existentes)

inject_signal <- function(drugA_id, drugB_id, event_id, 
                          dynamic_type, fold_change, 
                          ade_raw_dt) {
  
  if (drugA_id > drugB_id) {   # agrego ordenamiento (ya se hace en otras funciones pero por las dudas)
    temp <- drugA_id
    drugA_id <- drugB_id
    drugB_id <- temp
  }
  
  # copia independiente (para que datos inyectados no vayan corrompiendo el dataset original)
  ade_aug <- copy(ade_raw_dt)
  ade_aug[, simulated_event := FALSE]
  
  # 1- reportes con ambas drogas (coadministración)
  reports_A <- unique(ade_aug[atc_concept_id == drugA_id, safetyreportid])
  reports_B <- unique(ade_aug[atc_concept_id == drugB_id, safetyreportid])
  reports_AB <- intersect(reports_A, reports_B)
  
  if (length(reports_AB) <= 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      ade_aug = NULL,
      message = sprintf(
        "coadministración insuficiente: %d reportes",
        length(reports_AB),
      ),
      diagnostics = list(
        reason = "insufficient_coadmin",
        n_coadmin = length(reports_AB),
        drugA = drugA_id,
        drugB = drugB_id,
        event = event_id
      )
    ))
  }
  
  # 2- reportes objetivo (solo coadministración)
  target_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_AB, 
    .(safetyreportid, nichd, nichd_num)
  ])
  
  # Identifica si el evento YA existe en cada reporte
  event_in_report <- unique(ade_raw_dt[
    meddra_concept_id == event_id, 
    safetyreportid
  ])
  
  target_reports[, e_old := as.integer(safetyreportid %in% event_in_report)]
  
  # 3- Calculo tasa base (e_j)

  # Calculo tasas base individuales
  reports_A_clean <- setdiff(reports_A, reports_AB)
  reports_B_clean <- setdiff(reports_B, reports_AB)
  p_baseA <- mean(reports_A_clean %in% event_in_report)
  p_baseB <- mean(reports_B_clean %in% event_in_report)
  p_base0 <- length(event_in_report) / length(unique(ade_raw_dt$safetyreportid))

  # tener en cuenta que no son probabilidades, son odds, pero la extrapolación es válida si p <0.1  

  # Existen diversos métodos para considerar tasa de base
  # Todos se comportan distinto según el escenario

  # - Método aditivo:
  # consistente, genera IOR altos cuando riesgos individuales bajos + riesgo basal alto
  # e_j = P(evento | A ∪ B) asumiendo independencia
  # fórmula: P(A ∪ B) = P(A) + P(B) - P(A) × P(B)
  # p_baseA + p_baseB - (p_baseA * p_baseB)

  # - Método multiplicativo
  # consistente, solo genera IOR altos cuando riesgo basal bajo + riesgos individuales altos
  # tiene en cuenta todos los componentes (riesgos individuales y global)
  # e_j = (p_A × p_B) / p_0 
  # fórmula: P(A) × P(B) / P(0)
  # (p_baseA * p_baseB) / p_base0
  e_j <- p_baseA + p_baseB - (p_baseA * p_baseB)

  # t_ij = fold_change * e_j (tamaño de efecto)
  t_ij <- fold_change * e_j

  # 4- Probabilidades por etapa 
  # bprobs = rep(tij, N)
  # dy = tanh(...) * tij  (la dinámica se escala por tij)
  # rprobs = bprobs + dy
  
  N <- 7
  bprobs <- rep(t_ij, N)
  
  if (dynamic_type == "uniform") {
    dy <- rep(0, N)
  } else if (dynamic_type == "increase") {
    dy <- tanh(seq(-pi, pi, length.out = N)) * t_ij
  } else if (dynamic_type == "decrease") {
    dy <- -tanh(seq(-pi, pi, length.out = N)) * t_ij
  } else if (dynamic_type == "plateau") {
    dy <- c(
      tanh(seq(-pi, pi, length.out = floor(N/2))) * t_ij,
      tanh(seq(pi, -pi, length.out = ceiling(N/2))) * t_ij
    )
  } else if (dynamic_type == "inverse_plateau") {
    dy <- c(
      tanh(seq(pi, -pi, length.out = floor(N/2))) * t_ij,
      tanh(seq(-pi, pi, length.out = ceiling(N/2))) * t_ij
    )
  }
  
  rprobs <- bprobs + dy

  # clippeo de probabilidades
  rprobs <- pmax(pmin(rprobs, 0.999), 0.001)

  # 5- tabla de probabilidades por etapa
  stage_probs <- data.table(
    nichd_num = 1:N,
    bprobs = bprobs,
    dy = dy,
    p_dynamic = rprobs
  )
  
  # 6- Merge y generación de eventos nuevos
  target_reports <- merge(
    target_reports, 
    stage_probs[, .(nichd_num, p_dynamic)], 
    by = "nichd_num", 
    all.x = TRUE
  )
  
  # Simulación Bernoulli por reporte
  target_reports[, e_new := rbinom(.N, 1, p_dynamic)]
  
  # Combino con eventos existentes
  target_reports[, e_final := pmax(e_old, e_new)]
  
  # 7- marco reportes a inyectar
  # Identifico reportes que DEBERÍAN tener el evento (simulado)
  reports_to_mark <- target_reports[e_old == 0 & e_final == 1, safetyreportid]
  
  # validación de al menos 1 evento inyectado
  if (length(reports_to_mark) == 0) {
    return(list(
      success = FALSE,
      injection_success = FALSE,
      n_injected = 0,
      n_coadmin = length(reports_AB),
      ade_aug = NULL,
      message = sprintf(
        "Inyección fallida: 0 eventos generados (prob. media = %.4f, max = %.4f)",
        mean(target_reports$p_dynamic),
        max(target_reports$p_dynamic)
      ),
      diagnostics = list(
        reason = "zero_events_injected",
        low_probability_injection = TRUE,
        e_j = e_j,
        t_ij = t_ij,
        fold_change = fold_change,
        dynamic_type = dynamic_type,
        mean_p_dynamic = mean(target_reports$p_dynamic),
        max_p_dynamic = max(target_reports$p_dynamic),
        min_p_dynamic = min(target_reports$p_dynamic),
        n_eligible = nrow(target_reports[e_old == 0]),
        n_already_with_event = sum(target_reports$e_old),
        stage_probs = stage_probs
      )
    ))
  }
  
  # 8- marco reportes que ahora tienen el evento simulado (inyectados)
  ade_aug[
    safetyreportid %in% reports_to_mark,
    `:=`(
      simulated_event = TRUE,
      simulated_meddra = event_id,
      simulated_drugA = drugA_id,
      simulated_drugB = drugB_id
    )
  ]

  injection_rate <- length(reports_to_mark) / nrow(target_reports[e_old == 0])
  
  # 9- Diagnósticos
  diagnostics <- list(
    e_j = e_j,
    t_ij = t_ij,
    fold_change = fold_change,
    dynamic_type = dynamic_type,
    stage_probs = stage_probs,
    mean_p_dynamic = mean(target_reports$p_dynamic),
    max_p_dynamic = max(target_reports$p_dynamic),
    min_p_dynamic = min(target_reports$p_dynamic),
    n_eligible = nrow(target_reports),
    n_already_with_event = sum(target_reports$e_old),
    n_without_event = nrow(target_reports[e_old == 0]),
    n_new_events = length(reports_to_mark),
    injection_rate = length(reports_to_mark) / nrow(target_reports[e_old == 0]),
    injection_by_stage = target_reports[
      safetyreportid %in% reports_to_mark, 
      .N, 
      by = nichd_num
    ]
  )
  
  return(list(
    success = TRUE,
    injection_success = TRUE,
    n_injected = length(reports_to_mark),
    n_coadmin = length(reports_AB),
    ade_aug = ade_aug,
    message = sprintf(
      "Inyección exitosa: %d eventos en %d reportes (tasa: %.2f%%)",
      length(reports_to_mark),
      length(reports_AB),
      injection_rate * 100
    ),
    diagnostics = diagnostics
  ))
}

################################################################################
# Función auxiliar: calcular conteos básicos
################################################################################

# Calcula conteos básicos (eventos, coadministraciones) para un par droga-evento
#
# Parámetros:
#   ade_data: data.table aumentado con reportes (columnas atc_concept_id, meddra_concept_id, safetyreportid).
#   drugA: id de droga A (atc_concept_id).
#   drugB: id de droga B (atc_concept_id).
#   meddra: id del evento (meddra_concept_id).
#
# Return:
#   lista con: n_events, n_events_coadmin y n_coadmin

calc_basic_counts <- function(ade_data, drugA, drugB, meddra) {
  r_a <- unique(ade_data[atc_concept_id == drugA, safetyreportid])
  r_b <- unique(ade_data[atc_concept_id == drugB, safetyreportid])
  r_coadmin <- intersect(r_a, r_b)
  r_ea <- unique(ade_data[meddra_concept_id == meddra, safetyreportid])
  if("simulated_event" %in% names(ade_data)) {
    r_ea_sim <- unique(ade_data[simulated_event == TRUE & simulated_meddra == meddra, safetyreportid])
    r_ea <- union(r_ea, r_ea_sim)
  }
  list(
    n_events = length(r_ea),
    n_events_coadmin = length(intersect(r_coadmin, r_ea)),
    n_coadmin = length(r_coadmin)
  )
}

################################################################################
# Función para ajuste de GAM 
################################################################################

# Ajusta modelo GAM para interacción droga-droga 
#
# parametrizado para ir corriendo pruebas 
#
# parámetros:
# drugA_id: id Droga A
# drugB_id: id Droga B
# event_id: id del evento adverso
# ade_data: data.table con dataset original
# 
# include_nichd: Si TRUE, agrega efecto base de stage como covariable
# nichd_spline: Si TRUE, usa spline para efecto base de stage.
#               Si FALSE, usa coeficiente lineal (default: TRUE)
# spline_individuales: Si TRUE, usa splines para efectos individuales (suaviza riesgos basales drogaA y B por separado)
# bs_type: elección de tipo de spline: "cs", "tp", "cr"
# select: Si TRUE, permite penalización hasta cero (Para ver si algún coeficiente no aporta)
# include_sex: Si TRUE, incluye sexo como covariable
# include_stage_sex: Si TRUE, incluye interacción stage-sex
# k_spline: número de knots para splines (debería ser siempre 7 por nichd)
# method; método de ajuste GAM (dejar "fREML")
# 
# Return: 
# Lista con: success, n_events, n_coadmin, log_ior, reri, etc.

fit_gam <- function(drugA_id, drugB_id, event_id, ade_data,
                                 nichd_spline = TRUE,
                                 include_nichd = TRUE,
                                 spline_individuales = FALSE,
                                 bs_type = "cs",
                                 select = FALSE,
                                 include_sex = FALSE,
                                 include_stage_sex = FALSE,
                                 k_spline = 7,
                                 method = "fREML") {
  ###########
  # 1- identifico reportes 
  ###########
  
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  reportes_coadmin <- intersect(reportes_droga_a, reportes_droga_b)
  
  # agrego reportes simulados "flaggeados" en función de inyección
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  # Eventos previos 
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])

  # Combinar reales + simulados
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  n_events_total <- length(reportes_ea)  # total eventos (con o sin fármaco)
  n_coadmin <- length(reportes_coadmin)  # total reportes A+B (con o sin evento)
  n_events_coadmin <- length(intersect(reportes_coadmin, reportes_ea))  # eventos A+B

  ###########
  # 2- Construcción de dataset para ajustar
  ###########
  
  # Columnas base necesarias
  cols_necesarias <- c("safetyreportid", "nichd", "nichd_num")
  
  # sex
  if (include_sex) {
    cols_necesarias <- c(cols_necesarias, "sex")
  }
  
  datos_modelo <- unique(ade_data[, ..cols_necesarias])
  
  # variables de exposición
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # si sex está presente
  if (include_sex) {
    # estandarizo valores
    datos_modelo[, sex := toupper(trimws(sex))]
    datos_modelo[sex == "M", sex := "MALE"]
    datos_modelo[sex == "F", sex := "FEMALE"]
    # convierto a factor con niveles estándar
    datos_modelo[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  }
  
  ###########
  # 4- Construcción de fórmula con parámetros
  ###########
  
  # respuesta
  formula_parts <- "ea_ocurrio ~ "
  
  # opción A: efectos individuales lineales
  if (!spline_individuales) {
    formula_parts <- paste0(formula_parts, "droga_a + droga_b + ")
  } else {
    # opción B: efectos individuales con splines
    formula_parts <- paste0(
      formula_parts,
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_a) + ", 
              k_spline, bs_type),
      sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_b) + ", 
              k_spline, bs_type)
    )
  }

  # Efecto de nichd - spline o lineal (si include_nichd = TRUE)
  if (include_nichd) {
    if (nichd_spline) {
    # spline de base
      formula_parts <- paste0(
        formula_parts,
        sprintf("s(nichd_num, k = %d, bs = '%s') + ", k_spline, bs_type)
      )
    } else {
    # efecto lineal de nichd
      formula_parts <- paste0(formula_parts, "nichd_num + ")
    }
  }

  # spline de interacción (este no modificar)
  formula_parts <- paste0(
    formula_parts,
    sprintf("s(nichd_num, k = %d, bs = '%s', by = droga_ab)", k_spline, bs_type)
  )
  
  # si sexo TRUE
  if (include_sex) {
    if (include_stage_sex) {
      # Si sexo con spline por nichd
      formula_parts <- paste0(
        formula_parts,
        sprintf(" + s(nichd_num, k = %d, bs = '%s', by = sex)", 
                k_spline, bs_type)
      )
    } else {
      # Si solo sex lineal
      formula_parts <- paste0(formula_parts, " + sex")
    }
  }
  
  # Formula final
  formula_final <- as.formula(formula_parts)
  
  ###########
  # 5- Ajuste de modelo
  ###########
  
  tryCatch({
    
    modelo <- bam(
      formula = formula_final,
      data = datos_modelo,
      family = binomial(link = "logit"),
      method = method,
      select = select,    
      discrete = TRUE,
      nthreads = 1      # si no pongo esto, me da problemas en 10_augmentation con el uso de paralelización
    )
    
    ###########
    # 6- Calculo de log-IOR por nichd
    ###########
    
    # Grid de predicción (todas las combinaciones)
    grid_dif <- CJ(
      nichd_num = 1:7, 
      droga_a = c(0, 1), 
      droga_b = c(0, 1)
    )
    grid_dif[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
    
    # si sex en la formula, elegir nivel de referencia
    if (include_sex) {
      # MALE de referencia
      grid_dif[, sex := factor("MALE", levels = c("MALE", "FEMALE"))]
    }
    
    # predicciones
    pred_dif <- predict(modelo, newdata = grid_dif, type = "link", se.fit = TRUE)
    grid_dif[, `:=`(lp = pred_dif$fit, se = pred_dif$se.fit)]
    
    # Pivoteo para contraste
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                  value.var = c("lp", "se"))
    
    # Calculo final:  log(IOR) = log(OR₁₁) - log(OR₁₀) - log(OR₀₁) + log(OR₀₀)
    # Esto es lo mismo que hacer log( OR₁₁ / OR₁₀ . OR₀₁)
    log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
    
    ###########
    # 7- Calculo de standard error de log-IOR con matriz de CoVAR 
    ###########
    
    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE)  # unconditional = TRUE aplica corrección para incluir incertidumbre del suavizado
    cov_link <- Xp %*% Vb %*% t(Xp)
    
    log_ior_se <- numeric(7)
    for (stage in 1:7) {
      idx_00 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 0)
      idx_01 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 0 & grid_dif$droga_b == 1)
      idx_10 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 0)
      idx_11 <- which(grid_dif$nichd_num == stage & 
                        grid_dif$droga_a == 1 & grid_dif$droga_b == 1)
      
      # vector de contraste
      cvec <- rep(0, nrow(grid_dif))
      cvec[c(idx_11, idx_10, idx_01, idx_00)] <- c(1, -1, -1, 1)
      
      # SE = sqrt(c' Σ c)
      log_ior_se[stage] <- sqrt(max(
        as.numeric(t(cvec) %*% cov_link %*% cvec), 
        0
      ))
    }
    
    ###########
    # 8- Calculo de IC y métricas
    ###########
    
    # z90 para IC
    if (!exists("Z90")) {
      Z90 <- qnorm(0.95)
    }
    
    log_ior_lower90 <- log_ior - Z90 * log_ior_se
    log_ior_upper90 <- log_ior + Z90 * log_ior_se
    ior_values <- exp(log_ior)
    
    n_stages_significant <- sum(log_ior_lower90 > 0)
    max_ior <- max(ior_values)
    mean_ior <- mean(ior_values)

    ###########
    # 9- Cálculo de RERI por etapa (riesgo absoluto) con IC90
    ###########
    
    # etapas
    stages <- sort(unique(datos_modelo$nichd_num))
    
    # newdata con las 4 combinaciones A/B por etapa
    nd_reri <- rbindlist(lapply(stages, function(s) {
      data.table(
        nichd_num = s,
        droga_a   = c(0, 1, 0, 1),
        droga_b   = c(0, 0, 1, 1),
        droga_ab  = c(0, 0, 0, 1)
      )
    }), use.names = TRUE)
    
    # Covariables adicionales (si existen)
    if (include_sex) {
      nd_reri[, sex := factor(levels(datos_modelo$sex)[1], 
                             levels = levels(datos_modelo$sex))]
    }
    
    if (include_nichd && !nichd_spline) {
      nd_reri[, nichd := factor(niveles_nichd[nichd_num],
                                levels = niveles_nichd,
                                ordered = TRUE)]
    }
    
    # Predicción en escala de riesgo
    pred_reri <- predict(modelo, newdata = nd_reri, type = "link", se.fit = TRUE)
    nd_reri[, `:=`(
      eta = pred_reri$fit,
      se  = pred_reri$se.fit
    )]
    
    # Bootstrap paramétrico --> a diferencia del método delta, permite capturar no linealidad del método
    # matriz de diseño y coeficientes
    X_reri <- predict(modelo, newdata = nd_reri, type = "lpmatrix")
    beta_hat <- coef(modelo)
    V_beta <- vcov(modelo, unconditional = TRUE)
    
    # n simulaciones
    B <- 2000
    
    # simulación de coeficientes desde su distribución conjunta
    # beta_sim ~ MVN(beta_hat, V_beta)
    beta_sims <- mvrnorm(n = B, mu = beta_hat, Sigma = V_beta)
    
    # calculo predicciones para cada set de los coeficientes simulados
    p_sims <- matrix(NA, nrow = nrow(nd_reri), ncol = B)
    
    for (b in 1:B) {
      eta_b <- X_reri %*% beta_sims[b, ]
      p_sims[, b] <- plogis(eta_b)
    }
    
    # helper para RERI
    calc_reri <- function(p) {
      # p es un vector de 4 probabilidades por etapa: [p00, p10, p01, p11]
      p11 <- p[4]; p10 <- p[2]; p01 <- p[3]; p00 <- p[1]
      p11 - p10 - p01 + p00
    }
    
    # Cálculo por etapa
    reri_dt <- nd_reri[, {
      idx <- .I
      p_mat <- p_sims[idx, , drop = FALSE]
      reri_sim <- apply(p_mat, 2, calc_reri)
      
      data.table(
        RERI = mean(reri_sim),
        RERI_lower90 = quantile(reri_sim, 0.05),
        RERI_upper90 = quantile(reri_sim, 0.95)
      )
    }, by = nichd_num]
    
    # vectores de retorno
    reri_values <- reri_dt$RERI
    reri_lower90 <- reri_dt$RERI_lower90
    reri_upper90 <- reri_dt$RERI_upper90
    
    ###########
    # 10- Resultados
    ###########

    return(list(
      success = TRUE,
      n_events = n_events_total,
      n_coadmin = n_coadmin,
      n_events_coadmin = n_events_coadmin,
      log_ior = log_ior,
      log_ior_lower90 = log_ior_lower90,
      log_ior_upper90 = log_ior_upper90,
      log_ior_se = log_ior_se,
      ior_values = ior_values,
      n_stages_significant = n_stages_significant,
      max_ior = max_ior,
      mean_ior = mean_ior,
      reri_values = reri_values,
      reri_lower90 = reri_lower90,
      reri_upper90 = reri_upper90,
      n_stages_reri_significant = sum(reri_lower90 > 0),
      model_aic = AIC(modelo),
      model_deviance = deviance(modelo),
      formula_used = formula_parts,  # Guardo fórmula usada
      nichd_spline = nichd_spline,
      include_nichd = include_nichd,
      spline_individuales = spline_individuales,
      bs_type = bs_type,              
      select = select,                
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE, 
      n_events = n_events_total, 
      n_coadmin = n_coadmin,
      error_msg = e$message,
      formula_attempted = formula_parts
    ))
  })
}

################################################################################
# Cálculo de IOR clásico 
################################################################################

# Calcula IOR clásico usando tablas 2x2 por etapa
#
# parámetros:
# drugA_id: id droga A (ATC concept_id)
# drugB_id: id droga B (ATC concept_id)
# event_id: id del evento adverso (MedDRA concept_id)
# ade_data: data.table con dataset
# 
# return: Lista con stage, ior_classic, ior_classic_lower90, ior_classic_upper90
# 
# Implementación:
# Para cada etapa j:
# tabla 2x2:
# a: evento + coadmin
# b: sin evento + coadmin
# c: evento + sin coadmin
# d: sin evento + sin coadmin
# OR_11 = (a/b) / (c/d) = ad/bc
# OR_10 = reportes con A solo
# OR_01 = reportes con B solo
# OR_00 = reportes sin A ni B
# IOR = (OR_11 × OR_00) / (OR_10 × OR_01)
# calcula IC90 por método de Woolf (escala log)

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Identificación de reportes
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reportes con el evento (reales)
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Reportes con evento simulado 
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  
  # reales + simulados
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # Dataset único por reporte con exposiciones
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Cálculo por etapa
  resultados_por_etapa <- datos_unicos[, {
    
    # Tabla 2x2 para cada combinación de exposición
    # Grupo 11: A + B (coadministración)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_no_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 0)
    
    # Grupo 10: solo A
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_no_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 0)
    
    # Grupo 01: solo B
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_no_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 0)
    
    # Grupo 00: ni A ni B
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_no_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 0)
    
    # calculo de OR para cada grupo
    or_11 <- (n_11_evento / n_11_no_evento) / (n_00_evento / n_00_no_evento)
    or_10 <- (n_10_evento / n_10_no_evento) / (n_00_evento / n_00_no_evento)
    or_01 <- (n_01_evento / n_01_no_evento) / (n_00_evento / n_00_no_evento)
    or_00 <- 1  # por definición
    
    # como por definición OR₀₀ = 1, es lo mismo que hacer IOR = OR₁₁ / OR₁₀ . OR₀₁
    # calculo de IOR
    ior_val <- (or_11 * or_00) / (or_10 * or_01)
    log_ior <- log(ior_val)
    
    # Varianza en escala log (método de Woolf)
    # var(log(IOR)) = 1/a + 1/b + 1/c + 1/d para cada OR involucrado
    var_log_or_11 <- (1/n_11_evento + 1/n_11_no_evento + 
                      1/n_00_evento + 1/n_00_no_evento)
    var_log_or_10 <- (1/n_10_evento + 1/n_10_no_evento + 
                      1/n_00_evento + 1/n_00_no_evento)
    var_log_or_01 <- (1/n_01_evento + 1/n_01_no_evento + 
                      1/n_00_evento + 1/n_00_no_evento)
    
    # varianza del log(IOR) = var(log(OR_11)) + var(log(OR_10)) + var(log(OR_01))
    var_log_ior <- var_log_or_11 + var_log_or_10 + var_log_or_01
    se_log_ior <- sqrt(var_log_ior)
    
    # IC90 en escala log
    z90 <- qnorm(0.95)
    log_ior_lower90 <- log_ior - z90 * se_log_ior
    log_ior_upper90 <- log_ior + z90 * se_log_ior
    
    # IOR en escala original
    ior_lower90 <- exp(log_ior_lower90)
    ior_upper90 <- exp(log_ior_upper90)
    
    data.table(
      stage = nichd_num[1],
      ior_classic = ior_val,
      log_ior_classic = log_ior,
      ior_classic_lower90 = ior_lower90,
      ior_classic_upper90 = ior_upper90,
      log_ior_classic_lower90 = log_ior_lower90,
      log_ior_classic_upper90 = log_ior_upper90,
      se_log_ior_classic = se_log_ior,
      # diagnósticos
      n_11_evento = n_11_evento,
      n_11_total = n_11_evento + n_11_no_evento
    )
    
  }, by = nichd_num]
  
  # ordeno por etapa
  setorder(resultados_por_etapa, nichd_num)
  
  return(list(
    success = TRUE,
    results_by_stage = resultados_por_etapa
  ))
}

################################################################################
# Cálculo de RERI clásico 
################################################################################

# Calcula RERI clásico usando tablas 2x2 por etapa
#
# Parámetros:
# drugA_id: id droga A (ATC concept_id)
# drugB_id: id droga B (ATC concept_id)
# event_id: id del evento adverso (MedDRA concept_id)
# ade_data: data.table con dataset (puede ser aumentado)
# 
# Return: Lista con success, results_by_stage (stage, RERI, RERI_lower90, RERI_upper90)
# 
# Implementación:
# Para cada etapa j:
# Construye tabla 2x2:
# R11: riesgo con A+B
# R10: riesgo con solo A
# R01: riesgo con solo B
# R00: riesgo sin A ni B (referencia)
# RERI = R11 - R10 - R01 + R00
# IC90% usando método delta (propagación de varianzas)
# Varianza de cada riesgo: Var(R) = p(1-p)/n
# 5Varianza de RERI: suma de varianzas (asumiendo independencia)

calculate_classic_reri <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Identificación de reportes
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  
  # Reportes con el evento (reales)
  reportes_ea_real <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
  # Reportes con evento simulado (si existen)
  reportes_ea_sim <- if("simulated_event" %in% names(ade_data)) {
    unique(ade_data[
      simulated_event == TRUE & simulated_meddra == event_id,
      safetyreportid
    ])
  } else {
    integer(0)
  }
  
  # Combinar reales + simulados
  reportes_ea <- union(reportes_ea_real, reportes_ea_sim)
  
  # Dataset único por reporte con exposiciones
  datos_unicos <- unique(ade_data[, .(safetyreportid, nichd, nichd_num)])
  datos_unicos[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_unicos[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_unicos[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  
  # Cálculo por etapa
  resultados_por_etapa <- datos_unicos[, {
    
    # Conteos para cada combinación de exposición
    # Grupo 11: A + B (coadministración)
    n_11_evento <- sum(droga_a == 1 & droga_b == 1 & ea_ocurrio == 1)
    n_11_total <- sum(droga_a == 1 & droga_b == 1)
    
    # Grupo 10: solo A
    n_10_evento <- sum(droga_a == 1 & droga_b == 0 & ea_ocurrio == 1)
    n_10_total <- sum(droga_a == 1 & droga_b == 0)
    
    # Grupo 01: solo B
    n_01_evento <- sum(droga_a == 0 & droga_b == 1 & ea_ocurrio == 1)
    n_01_total <- sum(droga_a == 0 & droga_b == 1)
    
    # Grupo 00: ni A ni B (referencia)
    n_00_evento <- sum(droga_a == 0 & droga_b == 0 & ea_ocurrio == 1)
    n_00_total <- sum(droga_a == 0 & droga_b == 0)
    
    # Validación: todos los grupos deben tener datos
    if (n_11_total == 0 || n_10_total == 0 || n_01_total == 0 || n_00_total == 0) {
      return(data.table(
        stage = nichd_num[1],
        RERI_classic = NA_real_,
        RERI_classic_lower90 = NA_real_,
        RERI_classic_upper90 = NA_real_,
        RERI_classic_se = NA_real_,
        # Datos diagnósticos
        n_11_evento = n_11_evento,
        n_11_total = n_11_total,
        n_10_evento = n_10_evento,
        n_10_total = n_10_total,
        n_01_evento = n_01_evento,
        n_01_total = n_01_total,
        n_00_evento = n_00_evento,
        n_00_total = n_00_total,
        insufficient_data = TRUE
      ))
    }
    
    # Cálculo de riesgos (proporciones)
    R11 <- n_11_evento / n_11_total
    R10 <- n_10_evento / n_10_total
    R01 <- n_01_evento / n_01_total
    R00 <- n_00_evento / n_00_total
    
    # Cálculo de RERI
    # RERI = R11 - R10 - R01 + R00
    reri_val <- R11 - R10 - R01 + R00
    
    # Varianza de cada riesgo (distribución binomial)
    # Var(R) = p(1-p)/n
    # Para evitar división por cero, usar corrección de continuidad si p=0 o p=1
    var_R11 <- ifelse(R11 > 0 & R11 < 1, 
                      R11 * (1 - R11) / n_11_total,
                      0.25 / n_11_total)  # Máxima varianza posible
    
    var_R10 <- ifelse(R10 > 0 & R10 < 1,
                      R10 * (1 - R10) / n_10_total,
                      0.25 / n_10_total)
    
    var_R01 <- ifelse(R01 > 0 & R01 < 1,
                      R01 * (1 - R01) / n_01_total,
                      0.25 / n_01_total)
    
    var_R00 <- ifelse(R00 > 0 & R00 < 1,
                      R00 * (1 - R00) / n_00_total,
                      0.25 / n_00_total)
    
    # Varianza de RERI (método delta)
    # Como RERI = R11 - R10 - R01 + R00
    # Var(RERI) = Var(R11) + Var(R10) + Var(R01) + Var(R00)
    # (asumiendo independencia entre grupos)
    var_reri <- var_R11 + var_R10 + var_R01 + var_R00
    se_reri <- sqrt(var_reri)
    
    # IC90% 
    z90 <- qnorm(0.95)
    reri_lower90 <- reri_val - z90 * se_reri
    reri_upper90 <- reri_val + z90 * se_reri
    
    data.table(
      stage = nichd_num[1],
      RERI_classic = reri_val,
      RERI_classic_lower90 = reri_lower90,
      RERI_classic_upper90 = reri_upper90,
      RERI_classic_se = se_reri,
      # Riesgos individuales para diagnóstico
      R11 = R11,
      R10 = R10,
      R01 = R01,
      R00 = R00,
      # Conteos diagnósticos
      n_11_evento = n_11_evento,
      n_11_total = n_11_total,
      n_10_evento = n_10_evento,
      n_10_total = n_10_total,
      n_01_evento = n_01_evento,
      n_01_total = n_01_total,
      n_00_evento = n_00_evento,
      n_00_total = n_00_total,
      insufficient_data = FALSE
    )
    
  }, by = nichd_num]
  
  # Ordenar por etapa
  setorder(resultados_por_etapa, nichd_num)
  
  # Verificar si hay datos suficientes en al menos una etapa
  if (all(is.na(resultados_por_etapa$RERI_classic))) {
    return(list(
      success = FALSE,
      message = "Datos insuficientes en todas las etapas",
      results_by_stage = resultados_por_etapa
    ))
  }
  
  return(list(
    success = TRUE,
    results_by_stage = resultados_por_etapa
  ))
}

################################################################################
# Funcion de normalizado
################################################################################

# Normaliza resultado de triplete con listas vacías por defecto
#
# Return:
# Lista normalizada con vectores por defecto en lugar de NULL
# 
# Implementación:
# Previene errores de rbindlist cuando algunos campos de lista son NULL
# Aparentemente necesario cuando uso paralelización por si falla la convergencia de algunos tripletes

normalize_triplet_result <- function(result) {
  
  # campo tipo lista de columnas que deben existir
  list_fields <- c("stage", "log_ior", "log_ior_lower90", "ior_values",
                  "log_ior_classic", "log_ior_classic_lower90", "ior_classic")
  
  for (field in list_fields) {
    if (is.null(result[[field]])) {
      result[[field]] <- list(numeric(0))
    }
  }
  
  # Si las listas están vacías, llenar con valores por defecto
  if (length(result$stage[[1]]) == 0) {
    result$stage <- list(1:7)
  }
  
  if (length(result$log_ior[[1]]) == 0) {
    result$log_ior <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_lower90[[1]]) == 0) {
    result$log_ior_lower90 <- list(rep(NA_real_, 7))
  }
  
  if (length(result$ior_values[[1]]) == 0) {
    result$ior_values <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_classic[[1]]) == 0) {
    result$log_ior_classic <- list(rep(NA_real_, 7))
  }
  
  if (length(result$log_ior_classic_lower90[[1]]) == 0) {
    result$log_ior_classic_lower90 <- list(rep(NA_real_, 7))
  }
  
  if (length(result$ior_classic[[1]]) == 0) {
    result$ior_classic <- list(rep(NA_real_, 7))
  }
  
  return(result)
}

################################################################################
# Función de bootstrap por dinámica y etapa
################################################################################

# Realiza bootstrap de diferencia de log-IOR entre dinámica y base (uniform)
#
# Parámetros:
# data: data.table con columnas `dynamic`, `stage` y `log_ior`
# dynamic_type: nombre de la dinámica 
# stage_num: número de etapa NICHD.
# n_boot: número de réplicas bootstrap (por defecto 100)
#
# Return:
# data.table con estadísticos de bootstrap para la diferencia (media, sd, IC)

bootstrap_dynamic_diff <- function(data, dynamic_type, stage_num, n_boot = 100) {
  
  # Datos para la dinámica objetivo
  target_data <- data[dynamic == dynamic_type & stage == stage_num, log_ior]
  
  # Datos para uniform (baseline)
  uniform_data <- data[dynamic == "uniform" & stage == stage_num, log_ior]
  
  if (length(target_data) < 3 || length(uniform_data) < 3) {
    return(data.table(
      mean_diff = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  # Bootstrap
  boot_diffs <- replicate(n_boot, {
    target_sample <- sample(target_data, replace = TRUE)
    uniform_sample <- sample(uniform_data, replace = TRUE)
    mean(target_sample) - mean(uniform_sample)
  })
  data.table(
    mean_diff = mean(boot_diffs, na.rm = TRUE),
    ci_lower = quantile(boot_diffs, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_diffs, 0.975, na.rm = TRUE)
  )
}


################################################################################
# Función para procesar un solo triplete positivo con sensibilidad
################################################################################

# Procesa un triplete individual con todos los niveles de reducción
# Diseñada para usar en paralelo por lotes (consume mucha memoria)
#
# Parámetros:
# idx: id del triplete en pos_meta
# pos_meta: data.table con metadatos de tripletes positivos 
# ade_raw_dt: data.table con dataset original 
# reduction_levels: vector numérico con porcentajes de reducción a aplicar
# parámetros de configuración para la formula GAM
# base_seed: semilla base para reproducibilidad 
#
# Return:
#  data.table combinado con resultados para todos los niveles de reducción
#
# Implementación:
# inject_signal() para crear dataset aumentado independiente
# Para cada nivel de reducción --> reduce dataset, ajusta GAM, calcula IOR/RERI clásicos
# fit_reduced_model() como wrapper de ajuste

process_single_positive <- function(idx, pos_meta, ade_raw_dt, reduction_levels, 
                                    spline_individuales, include_sex, include_stage_sex,
                                    k_spline, bs_type, select, nichd_spline, z90, base_seed = 9427) {
  # Seed único para el triplete (mismo esquema que 10_augmentation)
  set.seed(9427 + idx)
  
  rowt <- pos_meta[idx]
  rowt$type <- "positive"
  
  # Creación de dataset aumentado INDEPENDIENTE
  inj_result <- tryCatch({
    inject_signal(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      dynamic_type = rowt$dynamic,
      fold_change = rowt$fold_change,
      ade_raw_dt = ade_raw_dt
    )
  }, error = function(e) {
    list(
      success = FALSE, 
      injection_success = FALSE,
      n_injected = 0, 
      n_coadmin = 0,
      ade_aug = NULL,
      message = paste("Error en inyección:", e$message),
      diagnostics = list(reason = "exception", error = e$message)
    )
  })

  inj_success <- inj_result$success
  n_injected_val <- inj_result$n_injected
  n_coadmin_val <- inj_result$n_coadmin
  diag_data <- list(inj_result$diagnostics)
  inj_message <- if(!is.null(inj_result$message)) inj_result$message else NA_character_
  
  t_ij_val <- if(inj_success && !is.null(inj_result$diagnostics$t_ij)) {
    inj_result$diagnostics$t_ij
  } else {
    NA_real_
  }

  rowt$t_ij <- t_ij_val

  if (!inj_success) {
    # Resultado de fallo sin sensibilidad
    base_result <- data.table(
      triplet_id = idx,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = "positive",
      reduction_pct = 0,
      N = rowt$N,
      dynamic = rowt$dynamic,
      fold_change = rowt$fold_change,
      t_ij = t_ij_val,
      model_success = FALSE,
      injection_success = FALSE,
      n_injected = n_injected_val,
      n_coadmin = n_coadmin_val,
      n_events = NA_integer_,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = FALSE,
      log_ior_classic = list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = list(rep(NA_real_, 7)),
      ior_classic = list(rep(NA_real_, 7)),
      reri_classic_success = FALSE,
      reri_values = list(rep(NA_real_, 7)),
      reri_lower90 = list(rep(NA_real_, 7)),
      reri_upper90 = list(rep(NA_real_, 7)),
      n_stages_reri_significant = NA_integer_,
      RERI_classic = list(rep(NA_real_, 7)),
      RERI_classic_lower90 = list(rep(NA_real_, 7)),
      RERI_classic_upper90 = list(rep(NA_real_, 7)),
      RERI_classic_se = list(rep(NA_real_, 7)),
      diagnostics = diag_data,
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = NA_character_,
      error_msg = inj_message
    )
    rm(inj_result); gc(verbose = FALSE)
    return(base_result)
  }

  # Loop de sensibilidad
  all_results <- list()
  
  # Resultado base (0% reducción)
  base_result <- fit_reduced_model(inj_result$ade_aug, rowt, 0)
  base_result$n_injected <- n_injected_val
  base_result$injection_success <- TRUE
  base_result$diagnostics <- diag_data
  
  all_results[[1]] <- base_result
  
  # Loop sobre niveles de reducción
  for (red_pct in reduction_levels) {
    # Reducir dataset
    ade_reduced <- reduce_dataset_by_stage(
      inj_result$ade_aug, 
      red_pct,
      seed = base_seed + idx)  # Propagación de seed para mayor reproductibilidad
    
    # Ajustar modelo en dataset reducido
    reduced_result <- fit_reduced_model(ade_reduced, rowt, red_pct)
    reduced_result$n_injected <- n_injected_val
    reduced_result$injection_success <- TRUE
    reduced_result$diagnostics <- diag_data
    
    all_results[[length(all_results) + 1]] <- reduced_result
    
    rm(ade_reduced); gc(verbose = FALSE)
  }
  
  rm(inj_result); gc(verbose = FALSE)
  
  # Combinar todos los resultados de sensibilidad
  combined_results <- rbindlist(all_results, fill = TRUE)
  
  return(combined_results)
}

################################################################################
# Función auxiliar: Reducir dataset por etapas
################################################################################

# Reduce un dataset aumentado removiendo un porcentaje aleatorio de filas por etapa
# 
# Parámetros:
# ade_aug: data.table aumentado
# reduction_pct: porcentaje a remover (ej: 10 para 10%)
# nichd_col: columna con la etapa NICHD
# seed: implementado para hacer reproductible 
#
# Return:
# data.table reducido

reduce_dataset_by_stage <- function(ade_aug, reduction_pct, nichd_col = "nichd", seed = NULL) {
  
  # Para cada etapa, remover el porcentaje especificado
  stages <- unique(ade_aug[[nichd_col]])
  
  reduced_rows <- lapply(stages, function(stage) {
    stage_data <- ade_aug[get(nichd_col) == stage]
    n_rows <- nrow(stage_data)
    n_keep <- ceiling(n_rows * (1 - reduction_pct / 100))
    
    if (n_keep < 1 && n_rows > 0) {
      n_keep <- 1  # Mantener al menos una fila si hay datos
    }
    
    if (n_rows > 0) {
      # Seed determinístico basado en stage + reduction_pct + dimensiones
      if (!is.null(seed)) {
        local_seed <- seed + as.integer(stage) * 100 + as.integer(reduction_pct) * 10
        set.seed(local_seed)
      }
      stage_data[sample(.N, min(n_keep, n_rows))]
    } else {
      stage_data
    }
  })
  
  rbindlist(reduced_rows, use.names = TRUE)
}

################################################################################
# Función de ajuste del modelo en dataset reducido
################################################################################

# Ajusta modelo GAM e IOR/RERI clásicos en un dataset reducido
# Wrapper de funciones para ajuste y reducción del dataset 
#
# Parámetros:
# ade_reduced: dataset reducido
# rowt: metadata del triplete (drugA, drugB, meddra, etc.)
# reduction_pct: porcentaje de reducción aplicado
# 
# Return:
# data.table con resultados del ajuste

fit_reduced_model <- function(ade_reduced, rowt, reduction_pct) {
  
  counts_reduced <- calc_basic_counts(ade_reduced, rowt$drugA, rowt$drugB, rowt$meddra)
  
  # Ajuste GAM
  model_res <- tryCatch({
    fit_gam(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced,
      spline_individuales = spline_individuales,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      nichd_spline = nichd_spline
    )
  }, error = function(e) {
    list(
      success = FALSE, 
      n_events_total = counts_reduced$n_events_total, 
      n_coadmin = counts_reduced$n_coadmin,
      error_msg = paste("Error en modelo reducido:", e$message)
    )
  })
  
  # Cálculo IOR clásico
  classic_res <- tryCatch({
    calculate_classic_ior(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
  }, error = function(e) {
    list(success = FALSE)
  })
  
  # Cálculo RERI clásico
  classic_reri <- tryCatch({
    calculate_classic_reri(
      drugA_id = rowt$drugA,
      drugB_id = rowt$drugB,
      event_id = rowt$meddra,
      ade_data = ade_reduced
    )
  }, error = function(e) {
    list(success = FALSE)
  })
  
  # Preparar resultado
  if (!model_res$success) {
    result <- data.table(
      triplet_id = rowt$triplet_id,
      drugA = rowt$drugA,
      drugB = rowt$drugB,
      meddra = rowt$meddra,
      type = rowt$type,
      reduction_pct = reduction_pct,
      N = counts_reduced$n_events_coadmin,
      dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
      fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
      t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
      model_success = FALSE,
      injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
      n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
      n_coadmin = counts_reduced$n_coadmin,
      n_events = counts_reduced$n_events_total,
      n_stages_significant = NA_integer_,
      max_ior = NA_real_,
      mean_ior = NA_real_,
      model_aic = NA_real_,
      stage = list(1:7),
      log_ior = list(rep(NA_real_, 7)),
      log_ior_lower90 = list(rep(NA_real_, 7)),
      ior_values = list(rep(NA_real_, 7)),
      classic_success = classic_res$success,
      log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
      log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
      ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
      reri_classic_success = classic_reri$success,
      reri_values = list(rep(NA_real_, 7)),
      reri_lower90 = list(rep(NA_real_, 7)),
      reri_upper90 = list(rep(NA_real_, 7)),
      n_stages_reri_significant = NA_integer_,
      RERI_classic = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic) else list(rep(NA_real_, 7)),
      RERI_classic_lower90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_lower90) else list(rep(NA_real_, 7)),
      RERI_classic_upper90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_upper90) else list(rep(NA_real_, 7)),
      RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7)),
      diagnostics = list(list(error = model_res$error_msg)),
      spline_individuales = spline_individuales,
      nichd_spline = nichd_spline,
      include_sex = include_sex,
      include_stage_sex = include_stage_sex,
      k_spline = k_spline,
      bs_type = bs_type,
      select = select,
      formula_used = if(!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
      error_msg = if(!is.null(model_res$error_msg)) model_res$error_msg else NA_character_
    )
    return(result)
  }
  
  # Resultado exitoso
  result <- data.table(
    triplet_id = rowt$triplet_id,
    drugA = rowt$drugA,
    drugB = rowt$drugB,
    meddra = rowt$meddra,
    type = rowt$type,
    reduction_pct = reduction_pct,
    N = model_res$n_events_coadmin,
    dynamic = if(!is.null(rowt$dynamic)) rowt$dynamic else NA_character_,
    fold_change = if(!is.null(rowt$fold_change)) rowt$fold_change else NA_real_,
    t_ij = if(!is.null(rowt$t_ij)) rowt$t_ij else NA_real_,
    model_success = TRUE,
    injection_success = if(!is.null(rowt$injection_success)) rowt$injection_success else NA,
    n_injected = if(!is.null(rowt$n_injected)) rowt$n_injected else NA_integer_,
    n_coadmin = model_res$n_coadmin,
    n_events = model_res$n_events_total,
    n_stages_significant = model_res$n_stages_significant,
    max_ior = model_res$max_ior,
    mean_ior = model_res$mean_ior,
    model_aic = model_res$model_aic,
    stage = list(1:7),
    log_ior = list(model_res$log_ior),
    log_ior_lower90 = list(model_res$log_ior_lower90),
    ior_values = list(model_res$ior_values),
    classic_success = classic_res$success,
    log_ior_classic = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
    log_ior_classic_lower90 = if(classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
    ior_classic = if(classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
    reri_classic_success = classic_reri$success,
    reri_values = list(model_res$reri_values),
    reri_lower90 = list(model_res$reri_lower90),
    reri_upper90 = list(model_res$reri_upper90),
    n_stages_reri_significant = model_res$n_stages_reri_significant,
    RERI_classic = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic) else list(rep(NA_real_, 7)),
    RERI_classic_lower90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_lower90) else list(rep(NA_real_, 7)),
    RERI_classic_upper90 = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_upper90) else list(rep(NA_real_, 7)),
    RERI_classic_se = if(classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7)),
    diagnostics = list(list()),
    spline_individuales = spline_individuales,
    nichd_spline = nichd_spline,
    include_sex = include_sex,
    include_stage_sex = include_stage_sex,
    k_spline = k_spline,
    bs_type = bs_type,
    select = select,
    formula_used = model_res$formula_used
  )
  
  return(result)
}

################################################################################
# Formula de permutación para distribución nula
################################################################################

# Permuta eventos entre reportes
# 
# Parámetros:
# perm_events: si TRUE, permuta eventos (meddra_concept_id)
# perm_drugs: si TRUE, permuta drogas también entre reportes
#
# Utiliza pool de reportes seleccionados en 10_augmentation (pool_meta)
# Rompe asociación droga-evento para tener un ground truth negativo

permute_pool <- function(pool_meta, niveles_nichd, 
                         perm_events = TRUE, 
                         perm_drugs = FALSE, 
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  pool_copy <- copy(pool_meta)
  pool_copy[, drugs_perm := vector("list", .N)]
  pool_copy[, events_perm := vector("list", .N)]
  
  for (stage in niveles_nichd) {
    idx <- which(pool_copy$nichd == stage)
    
    if (length(idx) <= 1) {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
      next
    }
    
    if (perm_drugs) {
      perm_idx_drugs <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx[perm_idx_drugs]]
    } else {
      pool_copy$drugs_perm[idx] <- pool_copy$drugs[idx]
    }
    
    if (perm_events) {
      perm_idx_events <- sample(seq_along(idx), length(idx), replace = FALSE)
      pool_copy$events_perm[idx] <- pool_copy$events[idx[perm_idx_events]]
    } else {
      pool_copy$events_perm[idx] <- pool_copy$events[idx]
    }
  }
  
  pool_copy[, .(safetyreportid, nichd, nichd_num, drugs_perm, events_perm)]
}

# Función para volver a meter los permutados en el dataset original antes de ajustar

reintroduce_permuted_reports <- function(ade_original, permuted_pool) {
  
  pool_report_ids <- unique(permuted_pool$safetyreportid)
  ade_without_pool <- ade_original[!safetyreportid %in% pool_report_ids]
  
  permuted_rows <- permuted_pool[, {
    drugs_vec <- drugs_perm[[1]]
    events_vec <- events_perm[[1]]
    
    if (length(drugs_vec) > 0 && length(events_vec) > 0) {
      CJ(safetyreportid = safetyreportid,
         atc_concept_id = drugs_vec,
         meddra_concept_id = events_vec,
         nichd = nichd,
         nichd_num = nichd_num)
    } else {
      data.table()
    }
  }, by = safetyreportid]
  
  rbindlist(list(ade_without_pool, permuted_rows), use.names = TRUE, fill = TRUE)
}

################################################################################
# Función de expansión
################################################################################

# Expande resultados de tripletes a formato largo con todas las métricas
#
# Parámetros:
# dt: data.table con resultados de tripletes (con columnas tipo lista)
# label_val: etiqueta para clasificación (1 = positivo, 0 = negativo)
# null_thresholds_dt: data.table con umbrales de la distribución nula por etapa
# use_threshold_ior: si TRUE, aplica umbral de IOR de la distribución nula
# use_threshold_reri: si TRUE, aplica umbral de RERI de la distribución nula
# 
# Return:
# data.table expandido con una fila por triplete-etapa
# 
# Implementación:
# Descompone columnas tipo lista (stage, log_ior, etc.) en filas individuales
# Une con umbrales de la distribución nula
# Mantiene metadatos del triplete (dynamic, t_ij, n_coadmin)

expand_clean_all_metrics <- function(dt, label_val, null_thresholds_dt,
                                     use_threshold_ior = TRUE, use_threshold_reri = TRUE) {

  has_dynamic <- "dynamic" %in% names(dt)
  by_cols <- "triplet_id"
  if ("dynamic" %in% names(dt)) {
    by_cols <- c(by_cols, "dynamic", "t_ij")
  }
  if ("n_coadmin" %in% names(dt)) {  
    by_cols <- c(by_cols, "n_coadmin")
  }

  expanded <- dt[, {
    stages <- unlist(stage)
    
    # GAM
    gam_log_ior <- unlist(log_ior)
    gam_log_ior_lower90 <- unlist(log_ior_lower90)
    gam_reri <- unlist(reri_values)
    gam_reri_lower90 <- unlist(reri_lower90)
    
    # Estratificado
    cls_log_ior <- unlist(log_ior_classic)
    cls_log_ior_lower90 <- unlist(log_ior_classic_lower90)
    cls_reri <- unlist(RERI_classic)
    cls_reri_lower90 <- unlist(RERI_classic_lower90)
    
    n <- min(length(stages), length(gam_log_ior), length(gam_log_ior_lower90),
             length(gam_reri), length(gam_reri_lower90),
             length(cls_log_ior), length(cls_log_ior_lower90),
             length(cls_reri), length(cls_reri_lower90))
    
    if (n > 0) {
      data.table(
        stage_num = stages[1:n],
        # GAM
        gam_log_ior = gam_log_ior[1:n],
        gam_log_ior_lower90 = gam_log_ior_lower90[1:n],
        gam_reri = gam_reri[1:n],
        gam_reri_lower90 = gam_reri_lower90[1:n],
        # Estratificado
        classic_log_ior = cls_log_ior[1:n],
        classic_log_ior_lower90 = cls_log_ior_lower90[1:n],
        classic_reri = cls_reri[1:n],
        classic_reri_lower90 = cls_reri_lower90[1:n]
      )
    }
  }, by = by_cols]
  
  if (!has_dynamic) {
    expanded[, `:=`(dynamic = "control", t_ij = 0)]
  }
  
  expanded[, nichd := niveles_nichd[stage_num]]
  expanded[, label := label_val]
  
  # Merge con thresholds
  expanded <- merge(expanded, null_thresholds_dt, 
                   by.x = "stage_num", by.y = "stage", all.x = TRUE)
  
  # Guardado de parámetros de uso
  expanded[, `:=`(
    use_threshold_ior = use_threshold_ior,
    use_threshold_reri = use_threshold_reri
  )]
  
  return(expanded)
}

################################################################################
# Función de cálculo de poder estadístico GAM
################################################################################

# Calcula poder estadístico filtrando por tamaño de efecto y cantidad de reportes de coadmin
#
# parámetros:
# data_pos: data.table con tripletes positivos expandidos
# data_neg: data.table con tripletes negativos expandidos (no utilizado)
# target_power: nivel de poder objetivo (0.80)
# metric_n: nombre de la columna a usar como filtro de cantidad ("n_coadmin" o "n_events")
# grid_resolution: número de pasos para la grilla de búsqueda (30x30 por defecto) 
# use_threshold_ior: si TRUE, usa umbral de la distribución nula
# use_threshold_reri: si TRUE, usa umbral de la distribución nula
# detection: modo de detección ("log-ior", "reri", o "double")
#
# return: 
# power_surface: data.table con t_threshold, n_threshold, power, len
# t_star: umbral óptimo de t_ij (tamaño de poder en formula de inyección)
# n_star: umbral óptimo de n_coadmin
# superset_pos: tripletes positivos que cumplen los filtros óptimos
# achieved_power: poder alcanzado en el punto óptimo
#
# Implementación:
# Genera una grilla 2D (según t_ij × n_coadmin) usando cuantiles
# Para cada punto de grilla, calcula el poder estadístico (TP / (TP + FN))
# Identifica configuración óptima que alcanza target_power con máxima retención
# Agrega a nivel triplete: detección si CUALQUIER etapa detecta

calculate_power_gam <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,  
  use_threshold_ior = TRUE,    # Parámetro explícito
  use_threshold_reri = TRUE,   # Parámetro explícito
  detection = "double"  ) {    # "ior", "reri", o "double"
  
    library(data.table)

    ###########  
    # 1. Preparación de datos positivos
    ###########
  
    pos_clean <- copy(data_pos)
  
    # Verificar existencia de columnas
    if (!metric_n %in% names(pos_clean)) {
      stop(sprintf("La columna métrica '%s' no existe en los datos.", metric_n))
    }
  
    # Limpieza básica
    pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
    # Merge thresholds si se pasan y faltan las columnas específicas
    if (!is.null(null_thresholds)) {
      # Verificar si faltan las columnas de threshold que vamos a usar
      need_merge <- FALSE
      if (use_threshold_ior && !"threshold_ior" %in% names(pos_clean)) need_merge <- TRUE
      if (use_threshold_reri && !"threshold_reri" %in% names(pos_clean)) need_merge <- TRUE
      if (need_merge) {
        # Asegurar que null_thresholds tenga la columna 'stage'
        if (!"stage" %in% names(null_thresholds)) {
          stop("null_thresholds debe tener una columna llamada 'stage'")
        }
        pos_clean <- merge(
          pos_clean,
          null_thresholds,
          by.x = "stage_num",
          by.y = "stage",
          all.x = TRUE
        )
      }
    }
  
  ###########
  # 2. Grilla 2D (t_ij x n_coadmin)
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Cálculo de poder
  ###########
  
  # Pre-cálculo de detección según parámetro "detection"
  
  # Criterio IOR
  if (detection %in% c("ior", "double")) {
    if (use_threshold_ior) {
      if (!"threshold_ior" %in% names(pos_clean)) {
        stop("columna threshold_ior no existe en los datos ni en null_thresholds")
      }
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & 
          gam_log_ior_lower90 > 0 & 
          gam_log_ior_lower90 > threshold_ior
      )]
    } else {
      pos_clean[, ior_detected := (
        !is.na(gam_log_ior_lower90) & gam_log_ior_lower90 > 0
      )]
    }
  }
  
  # Criterio RERI
  if (detection %in% c("reri", "double")) {
    if (use_threshold_reri) {
      if (!"threshold_reri" %in% names(pos_clean)) {
        stop("columna threshold_reri no existe en los datos ni en null_thresholds")
      }
      pos_clean[, reri_detected := (
        !is.na(gam_reri_lower90) & 
          gam_reri_lower90 > 0 & 
          gam_reri_lower90 > threshold_reri
      )]
    } else {
      pos_clean[, reri_detected := (
        !is.na(gam_reri_lower90) & gam_reri_lower90 > 0
      )]
    }
  }
  
  # Detección final según el modo elegido
  if (detection == "ior") {
    pos_clean[, is_detected := ior_detected]
  } else if (detection == "reri") {
    pos_clean[, is_detected := reri_detected]
  } else {  # "double"
    pos_clean[, is_detected := ior_detected | reri_detected]
  }

  # Agregación a nivel triplete: SI cualquier etapa detecta, triplete detectado
  triplet_detection_gam <- pos_clean[, .(
    triplet_detected = any(is_detected, na.rm = TRUE),
    t_ij_triplet = unique(t_ij)[1],
    n_coadmin_triplet = unique(get(metric_n))[1]
  ), by = triplet_id]
  
  # Iterar sobre grilla
  power_surface <- search_grid[, {
    
    # Filtro a nivel triplete: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      triplet_detection_gam$t_ij_triplet >= t_threshold & 
      triplet_detection_gam$n_coadmin_triplet >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(triplet_detection_gam$triplet_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identificación de punto óptimo
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("ADVERTENCIA (GAM): No se alcanzó poder objetivo (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # Mensaje informativo según modo de detección
  detection_label <- switch(detection,
    "ior" = "solo IOR",
    "reri" = "solo RERI", 
    "double" = "IOR O RERI"
  )
  
  message(sprintf("\nUMBRALES ÓPTIMOS (GAM - %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Poder alcanzado: %.1f%%", achieved_power * 100))
  message(sprintf("  Tripletes retenidos: %d / %d (%.1f%%)",
                  n_retained,
                  uniqueN(triplet_detection_gam$triplet_id),
                  100 * n_retained / uniqueN(triplet_detection_gam$triplet_id)))
  
  ###########
  # 5. Construcción de supersets
  ###########
  
  # Superset positivo: TODAS las observaciones de tripletes que cumplen filtros
  triplets_passed <- triplet_detection_gam[
    t_ij_triplet >= t_star & n_coadmin_triplet >= n_star,
    triplet_id
  ]
  
  superset_pos <- pos_clean[triplet_id %in% triplets_passed]
  
  power_surface[, method := paste0("GAM-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "gam",
    detection_mode = detection,
    metric_n_used = metric_n
  ))
}

################################################################################
# Función de cálculo de poder estadístico para método Estratificado
################################################################################

# Calcula poder estadístico filtrando por tamaño de efecto y cantidad de reportes de coadmin
#
# Como los métodos clasicos tienen más a arrojar NAs o IC con Inf, parametrizo manejo de NA
#
# parámetros:
# data_pos: data.table con tripletes positivos expandidos
# data_neg: data.table con tripletes negativos expandidos
# target_power: nivel de poder objetivo (0.80)
# metric_n: nombre de la columna a usar como filtro de cantidad ("n_coadmin" o "n_events")
# grid_resolution: número de pasos para la grilla de búsqueda (30x30 por defecto)
# detection: para hacer superset de IOR, RERI o ambos 
# na_remove: elegir manejo de NA. 
#  TRUE: superset incluye TODOS los tripletes con criterios de detección mínima
#  FALSE: superset final excluye NA aunque cumplan criterios de umbral mínimo
#
# return: 
# power_surface: data.table con t_threshold, n_threshold, power, len
# t_star: umbral óptimo de t_ij
# n_star: umbral óptimo de n_coadmin
# superset_pos: tripletes positivos que cumplen los filtros óptimos
# achieved_power: poder alcanzado en el punto óptimo
#
# Implementación:
# Ver calculate_power_gam
# cambios: calcula el poder a nivel de etapa

calculate_power_classic <- function(
  data_pos,
  target_power = 0.80,
  null_thresholds = NULL,
  metric_n = "n_coadmin",
  grid_resolution = 30,
  detection = "double",  # "ior", "reri", o "double"
  na_remove = "TRUE") {
    
  ###########
  # 1. Preparación de datos positivos
  ###########
  detection <- match.arg(detection, choices = c("ior", "reri", "double"))

  pos_clean <- copy(data_pos)
  
  n_triplets_original_total <- uniqueN(pos_clean$triplet_id)
  n_obs_original_total <- nrow(pos_clean)

  # Limpieza de NA según criterio de detección
  
  if (na_remove) {
    n_obs_before_na <- nrow(pos_clean)
    n_triplets_before_na <- uniqueN(pos_clean$triplet_id)
    if (detection == "ior") {
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90)]
    } else if (detection == "reri") {
      pos_clean <- pos_clean[!is.na(classic_reri_lower90)]
    } else {  # double
      pos_clean <- pos_clean[!is.na(classic_log_ior_lower90) & !is.na(classic_reri_lower90)]
    }
  }
  # registro los que se eliminaron
  n_obs_after_na <- nrow(pos_clean)
  n_triplets_after_na <- uniqueN(pos_clean$triplet_id)
  n_triplets_lost_na <- n_triplets_before_na - n_triplets_after_na
 
  # Limpieza básica
  pos_clean <- pos_clean[is.finite(t_ij) & is.finite(get(metric_n))]
  
  if (!is.null(null_thresholds) && !"threshold" %in% names(pos_clean)) {
    pos_clean <- merge(
      pos_clean,
      null_thresholds,
      by.x = "stage_num",
      by.y = "stage",
      all.x = TRUE
    )
  }
  
  ###########
  # 2. Grilla 2D
  ###########
  
  probs_grid <- seq(0, 0.95, length.out = grid_resolution)
  
  t_vals <- unique(quantile(pos_clean$t_ij, probs = probs_grid, na.rm = TRUE))
  n_vals <- unique(quantile(pos_clean[[metric_n]], probs = probs_grid, na.rm = TRUE))
  
  t_vals <- sort(unique(c(min(pos_clean$t_ij), t_vals)))
  n_vals <- sort(unique(c(min(pos_clean[[metric_n]]), n_vals)))
  
  search_grid <- CJ(t_threshold = t_vals, n_threshold = n_vals)
  
  ###########
  # 3. Cálculo de poder 
  ###########
  
  # cada fila (triplet_id + stage) es una observación independiente
  # porque el método está estratificado por etapa
  
  # Pre-cálculo de detección según parámetro 'detection'
  
  if (detection == "ior") {
    # Solo IOR: IC90 > 0
    pos_clean[, is_detected := (
      !is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0
    )]
    
  } else if (detection == "reri") {
    # Solo RERI: IC90 > 0
    pos_clean[, is_detected := (
      !is.na(classic_reri_lower90) & classic_reri_lower90 > 0
    )]
    
  } else {  # detection == "double"
    # Criterio doble: IOR O RERI, IC90 > 0
    pos_clean[, is_detected := (
      (!is.na(classic_log_ior_lower90) & classic_log_ior_lower90 > 0) |
      (!is.na(classic_reri_lower90) & classic_reri_lower90 > 0)
    )]
  }
    
  # guardo copia con detecciones
  pos_all_with_detection <- copy(pos_clean)

  # se calcula a nivel observación para cada fila: t_ij, n_coadmin, detección
  # No se agrega a nivel triplete
  
  # Iterar sobre grilla
  power_surface <- search_grid[, {
    
    # Filtro a nivel de observación: t_ij >= t_thresh AND n >= n_thresh
    idx_subset <- which(
      pos_clean$t_ij >= t_threshold & 
        pos_clean[[metric_n]] >= n_threshold
    )
    
    n_total <- length(idx_subset)
    
    if (n_total < 5) {
      list(tp = 0L, len = 0L, power = NA_real_)
    } else {
      n_tp <- sum(pos_clean$is_detected[idx_subset], na.rm = TRUE)
      list(
        tp = n_tp,
        len = n_total,
        power = n_tp / n_total
      )
    }
  }, by = .(t_threshold, n_threshold)]
  
  power_surface <- power_surface[!is.na(power)]
  
  ###########
  # 4. Identificación de punto óptimo
  ###########
  
  valid_configs <- power_surface[power >= target_power]
  
  if (nrow(valid_configs) == 0) {
    message(sprintf("No se alcanzó poder objetivo (%.0f%%). Max: %.1f%%",
                    target_power*100, max(power_surface$power, na.rm=TRUE)*100))
    best_config <- power_surface[which.max(power)]
  } else {
    setorder(valid_configs, -len, power)
    best_config <- valid_configs[1]
  }
  
  t_star <- best_config$t_threshold
  n_star <- best_config$n_threshold
  achieved_power <- best_config$power
  n_retained <- best_config$len
  
  # modo de detección para mensaje
  detection_label <- switch(detection,
    "ior" = "solo IOR",
    "reri" = "solo RERI", 
    "double" = "IOR O RERI"
  )

  ###########
  # 5. Construcción de supersets
  ###########
  
  # Superset positivo: observaciones que cumplen filtros
  superset_pos <- pos_all_with_detection[t_ij >= t_star & get(metric_n) >= n_star]
  
  # Métricas del superset
  n_retained_total <- nrow(superset_pos)
  n_detected_total <- sum(superset_pos$is_detected, na.rm = TRUE)
  power_total <- ifelse(n_retained_total > 0, n_detected_total / n_retained_total, 0)

  # conteo de tripletes únicos que cumplen criterio de retención 
  n_triplets_retained <- uniqueN(superset_pos$triplet_id)
  n_triplets_total <- uniqueN(pos_all_with_detection$triplet_id)

  n_triplets_retained_vs_original <- n_triplets_retained
  pct_retained_vs_original <- 100 * n_triplets_retained_vs_original / n_triplets_original_total

  message(sprintf("\nUMBRALES ÓPTIMOS (Método Estratificado, %s):", detection_label))
  message(sprintf("  t_ij >= %.4f", t_star))
  message(sprintf("  %s >= %.1f", metric_n, n_star))
  message(sprintf("  Poder alcanzado: %.1f%%", achieved_power * 100))
  message(sprintf("  Tripletes retenidos: %d / %d (%.1f%%)",
                  n_triplets_retained, n_triplets_total,
                  100 * n_triplets_retained / n_triplets_total))
  message(sprintf("  Poder en muestra completa: %.1f%%", power_total * 100))

  # Lineas para ver retención verdadera, contando NA que se excluyeron al principio de la función
  message(sprintf("  Tripletes retenidos (vs total original): %d / %d (%.1f%%) [INCLUYE NA]",
                  n_triplets_retained_vs_original, n_triplets_original_total,
                  pct_retained_vs_original))
  
  power_surface[, method := paste0("Estratificado-", toupper(detection))]
  
  return(list(
    power_surface = power_surface,
    t_star = t_star,
    n_star = n_star,
    superset_pos = superset_pos,
    achieved_power = achieved_power,
    criterion_type = "classic",
    metric_n_used = metric_n,
    detection_mode = detection,
    na_remove = na_remove,
    achieved_power_total = power_total,  # poder en muestra completa
    n_retained_grid = n_retained,  # retenidos en grid
    n_retained_total = n_retained_total,  # retenidos totales
    n_detected_total = n_detected_total,  # detectados totales
    n_triplets_retained = n_triplets_retained,
    n_triplets_total = n_triplets_total
  ))
}

################################################################################
# Función para visualizar superficie de poder con grilla completa
################################################################################

# Genera heatmap de superficie de poder con interpolación a grilla regular
#
# Parámetros:
# power_result: lista retornada por calculate_power_gam() o calculate_power_classic()
# target_power: nivel de poder objetivo (para referencia visual)
# detection: tipo de detección ("IOR", "RERI", "double")
# t_range: rango del eje X (tamaño de efecto)
# n_range: rango del eje Y (reportes de coadministración)
# grid_size: tamaño de la grilla interpolada (default: 50x50 para suavidad)
# 
# Return:
# Objeto ggplot con superficie interpolada completa
#
# Implementación:
# Remueve NAs e Inf
# Interpola a grilla regular usando akima::interp (para datos faltantes)
# Genera heatmap con geom_raster
# retorna estadísticas del punto óptimo

plot_power_surface <- function(
  power_result, 
  target_power = 0.80, 
  detection = "double",
  t_range = c(0, 0.5),      
  n_range = c(0, 300),
  grid_size = 30) 
  {
  
  library(ggplot2)
  library(data.table)
  library(scales)
  library(akima)
  
  ###########
  # 1- Validación y preparación de datos
  ###########
  
  surface <- as.data.table(power_result$power_surface)
  
  # parámetros 
  t_star <- power_result$t_star
  n_star <- power_result$n_star
  achieved_power <- power_result$achieved_power
  criterion_type <- power_result$criterion_type
  metric_n <- power_result$metric_n_used
  
  ###########
  # 2- Limpieza de datos originales
  ###########
  
  # Remover NA y valores no finitos
  surface_clean <- surface[is.finite(t_threshold) & is.finite(n_threshold) & is.finite(power)]
  
  ###########
  # 3- Interpolación a grilla regular completa
  ###########
  
  # Crear secuencias regulares para la grilla
  t_seq <- seq(t_range[1], t_range[2], length.out = grid_size)
  n_seq <- seq(n_range[1], n_range[2], length.out = grid_size)
  
  # Interpolación bilineal con akima
  interp_result <- tryCatch({
    interp(
      x = surface_clean$t_threshold,
      y = surface_clean$n_threshold,
      z = surface_clean$power,
      xo = t_seq,
      yo = n_seq,
      linear = TRUE,     # Interpolación lineal 
      extrap = FALSE     # No extrapolar fuera de datos
    )
  }, error = function(e) {
    warning(sprintf("Error en interpolación: %s", e$message))
    return(NULL)
  })

  # Convertir resultado de interpolación a data.table
  surface_interp <- data.table(
    expand.grid(
      t_threshold = interp_result$x,
      n_threshold = interp_result$y
    )
  )
    
  # Aplanar matriz a vector
  surface_interp[, power := as.vector(interp_result$z)]
    
  # Remover NA generados por extrap = FALSE
  surface_plot <- surface_interp[!is.na(power)]
    
  # Clipear poder a [0, 1] por si hay overshooting de interpolación
  surface_plot[, power := pmin(pmax(power, 0), 1)]

  message(sprintf("  Grilla final: %d celdas", nrow(surface_plot)))
  
  ###########
  # 5- Etiquetas de método
  ###########

  method_label <- fifelse(
    criterion_type == "gam",
    "GAM",
    "Estratificado"
  )
  
  detection_label <- switch(
    detection,
    "IOR" = "IOR",
    "RERI" = "RERI",
    "DOUBLE" = "IOR/RERI",
    "DOBLE" = "IOR/RERI",
    detection
  )
  
  ###########
  # 6- Construcción del gráfico
  ###########
  
  p <- ggplot(surface_plot, aes(x = t_threshold, y = n_threshold, fill = power)) +
    
    # geom_raster() para grillas regulares
    geom_raster(interpolate = FALSE) +  # interpolate=TRUE suaviza visualmente
    
    # Escala de color
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1), 
      name = "Poder\n(TPR)",
      option = "plasma",
      na.value = "white",
      guide = guide_colorbar(label.theme = element_text(angle = 45, vjust = 1, hjust = 1)) # Celdas sin datos en blanco
    ) +
    
    # Escalas de ejes
    scale_x_continuous(
      limits = t_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 5)
    ) +
    
    scale_y_continuous(
      limits = n_range,
      expand = c(0, 0),
      breaks = pretty_breaks(n = 6)
    ) +
    
    # Etiquetas
    labs(
      title = sprintf(
        "Superficie de Poder - Método %s - Detección %s",
        method_label,
        detection_label
      ),
      subtitle = sprintf(
        "Óptimo: t≥%.4f, N≥%.0f | Poder alcanzado: %.1f%% (objetivo: %.0f%%)",
        t_star,
        n_star,
        achieved_power * 100,
        target_power * 100
      ),
      x = expression("Tamaño de Efecto " * (t[italic(ij)])),
      y = sprintf("número de reportes A-B"),
      caption = sprintf(
        "Grilla: %dx%d celdas interpoladas desde %d puntos originales",
        grid_size,
        grid_size,
        nrow(surface_clean)
      )
    ) 
  return(p)
}






