################################################################################
# Script de funciones
# Script: 00_functions.R
# para usar, source("00_functions.R", local = TRUE)
################################################################################

# ordenado de niveles para el resto de los scripts

niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

Z90 <- qnorm(0.95)  # Cuantil 90% para intervalos de confianza


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
# 
dynamic_fun <- function(type, N = 7) {
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
# - FC ~ 1 + Exp(λ = 0.75)
# - rango típico: [1, 10] con sesgo a la derecha

sample_fold_change <- function(n, lambda = 0.75) {
  1 + rexp(n, rate = lambda)
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
# dynamic_type: tipo de dinámica (llama dynamic_fun)
# fold_change: magnitud del efecto (llama sample_fold_change)
# ade_raw_dt: data.table de dataset original
# 
# return Lista con: success, n_injected, n_coadmin, ade_aug, diagnostics
# 
# Implementación:
# 
# 1- Tasa base(t_ij):
#  Probabilidad poblacional del evento
#  Ahora calculada con intersección de droga A y B (distinto a Giangreco)
#  Independiente de coadministración
# 
# 2- Fold_change (FC):
#  Factor multiplicativo del efecto
# 
# 3- Dinámica f(j):
#   Función normalizada en [0, 1] por etapa
#   f(j) = llama a dynamic_fun(type, N=7)  <--- N siempre 7 por etapas nichd 
# 
# 4- Probabilidad dinámica: (formula final de probabilidad de reporte)
#   p_dynamic(j) = t_ij × FC + f(j)
# 
# 5- Simulación:
#   Para cada reporte con drugA + drugB
#   Y_new ~ Bernoulli(p_dynamic(stage_j))
#   Inyecta SOLO en eventos que no existen (si ya hay evento, se deja, trata de modelar dinámica sobre existentes)

inject_signal <- function(drugA_id, drugB_id, event_id, 
                          dynamic_type, fold_change, 
                          ade_raw_dt) {
  
  # copia independiente (para que datos inyectados no vayan corrompiendo el dataset original)
  ade_aug <- copy(ade_raw_dt)
  ade_aug[, simulated_flag := FALSE]
  
  # 1- reportes con ambas drogas (coadministración)
  reports_A <- unique(ade_aug[atc_concept_id == drugA_id, safetyreportid])
  reports_B <- unique(ade_aug[atc_concept_id == drugB_id, safetyreportid])
  reports_both <- intersect(reports_A, reports_B)
  
  if (length(reports_both) == 0) {
    return(list(
      success = FALSE, n_injected = 0, n_coadmin = 0, 
      ade_aug = NULL, 
      message = "No coadministración encontrada",
      diagnostics = NULL
    ))
  }
  
  # 2- reportes objetivo (solo coadministración)
  target_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_both, 
    .(safetyreportid, nichd, nichd_num)
  ])
  
  # Identifica si el evento YA existe en cada reporte
  event_in_report <- unique(ade_raw_dt[
    meddra_concept_id == event_id, 
    safetyreportid
  ])
  
  target_reports[, e_old := as.integer(safetyreportid %in% event_in_report)]
  
  # 3- Calculo tasa base (t_ij)
  # cambio por probabilidad de la unión de eventos independientes
  # Calcular tasas base individuales
  p_baseA <- mean(reports_A %in% event_in_report)
  p_baseB <- mean(reports_B %in% event_in_report)

  # e_j = P(evento | A ∪ B) asumiendo independencia
  # fórmula: P(A ∪ B) = P(A) + P(B) - P(A) × P(B)
  e_j <- p_baseA + p_baseB - (p_baseA * p_baseB)

  # t_ij = fold_change * e_j 
  t_ij <- fold_change * e_j
  
  # CHEQUEO: Validación de tasa base mínima  (revisar si saco o ajusto)
  if (e_j < 0.001) {
    return(list(
      success = FALSE, 
      n_injected = 0, 
      n_coadmin = length(reports_both),
      ade_aug = NULL,
      message = sprintf("Tasa base muy baja: e_j = %.4f", e_j),
      diagnostics = list(e_j = e_j, t_ij = t_ij, fold_change = fold_change)
    ))
  }

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
  
  # CHEQUEO: asegurar que las probabilidades estén en [0, 1]
  rprobs <- pmin(pmax(rprobs, 0), 1)
  
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
  
  if (length(reports_to_mark) == 0) {
    return(list(
      success = TRUE,  
      n_injected = 0,
      n_coadmin = length(reports_both),
      ade_aug = ade_aug,  
      message = "Probabilidad de inyección baja. 0 eventos nuevos",
      diagnostics = list(
        low_probability_injection = TRUE,
        e_j = e_j,
        t_ij = t_ij,
        fold_change = fold_change,
        mean_p_dynamic = mean(target_reports$p_dynamic),
        max_p_dynamic = max(target_reports$p_dynamic),
        n_eligible = nrow(target_reports[e_old == 0]),
        dynamic_type = dynamic_type,
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
  
  # 9- CHEQUEO: diagnósticos
  diagnostics <- list(
    # Parámetros fundamentales
    e_j = e_j,
    t_ij = t_ij,
    fold_change = fold_change,
    dynamic_type = dynamic_type,
    
    # Perfil dinámico
    stage_probs = stage_probs,
    
    # Estadísticas de probabilidades
    mean_p_dynamic = mean(target_reports$p_dynamic),
    max_p_dynamic = max(target_reports$p_dynamic),
    min_p_dynamic = min(target_reports$p_dynamic),
    
    # Reportes
    n_eligible = nrow(target_reports),
    n_already_with_event = sum(target_reports$e_old),
    n_without_event = nrow(target_reports[e_old == 0]),
    
    # Resultados de inyección
    n_new_events = length(reports_to_mark),
    injection_rate = length(reports_to_mark) / nrow(target_reports[e_old == 0]),
    
    # Distribución por etapa
    injection_by_stage = target_reports[
      safetyreportid %in% reports_to_mark, 
      .N, 
      by = nichd_num
  )
  
  return(list(
    success = TRUE,
    n_injected = nrow(reports_to_mark),
    n_coadmin = length(reports_both),
    ade_aug = ade_aug,
    diagnostics = diagnostics
  ))
}


################################################################################
# Función para ajuste de GAM 
################################################################################

# Ajusta modelo GAM para interacción droga-droga 
#
# parametrizado para ir corriendo pruebas 

# parámetros:
# drugA_id: id Droga A
# drugB_id: id Droga B
# event_id: id del evento adverso
# ade_data: data.table con dataset original
# 
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
# Lista con: success, n_events, n_coadmin, log_ior, etc.


fit_differential_gam <- function(drugA_id, drugB_id, event_id, ade_data,
                                 nichd_spline = TRUE,
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
  
  # agrego reportes simulados "flaggeados" en función de inyección
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
  
  ###########
  # 2- Construcción de dataset para ajustar
  ###########
  
  # Columnas base necesarias
  cols_necesarias <- c("safetyreportid", "nichd", "nichd_num")
  
  # sex
  if (include_sex) {
    cols_necesarias <- c(cols_necesarias, "sex")
  }
  
  # datos únicos por reporte
  datos_modelo <- unique(ade_data[, ..cols_necesarias])
  
  # variable de exposición 
  datos_modelo[, ea_ocurrio := as.integer(safetyreportid %in% reportes_ea)]
  datos_modelo[, droga_a := as.integer(safetyreportid %in% reportes_droga_a)]
  datos_modelo[, droga_b := as.integer(safetyreportid %in% reportes_droga_b)]
  datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]
  
  # si sex está presente
  if (include_sex) {
    # Estandarizar valores 
    datos_modelo[, sex := toupper(trimws(sex))]
    datos_modelo[sex == "M", sex := "MALE"]
    datos_modelo[sex == "F", sex := "FEMALE"]
    
    # convertir a factor con niveles estándar
    datos_modelo[, sex := factor(sex, levels = c("MALE", "FEMALE"))]

    }
  
  
  n_events <- sum(datos_modelo$ea_ocurrio)
  n_coadmin <- sum(datos_modelo$droga_ab)
  
  ###########
  # 3- Requisitos mínimos (cantidad de reportes en conjunto y eventos)
  ###########
  
  if (n_events < 5 || n_coadmin < 3) {
    return(list(
      success = FALSE,
      n_events = n_events,
      n_coadmin = n_coadmin,
      message = sprintf("Datos insuficientes: %d eventos, %d coadmin", 
                        n_events, n_coadmin)
    ))
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
  
  # Efecto de nichd - spline o lineal
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
      nthreads = 1              # si no pongo esto, me da problemas en 01_augmentation con el uso de paralelización
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
    
    # Pivotear para contraste
    w_lp <- dcast(grid_dif, nichd_num ~ droga_a + droga_b, 
                  value.var = c("lp", "se"))
    
    # Calculo final:  log(IOR) = log(OR₁₁) - log(OR₁₀) - log(OR₀₁) + log(OR₀₀)
    log_ior <- w_lp$lp_1_1 - w_lp$lp_1_0 - w_lp$lp_0_1 + w_lp$lp_0_0
    
    ###########
    # 7- Calculo de standard error de log-IOR con matriz de CoVAR 
    ###########
    
    Xp <- predict(modelo, newdata = grid_dif, type = "lpmatrix")
    Vb <- vcov(modelo, unconditional = TRUE)
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
    # 9- Resultados
    ###########
    
    return(list(
      success = TRUE,
      n_events = n_events,
      n_coadmin = n_coadmin,
      log_ior = log_ior,
      log_ior_lower90 = log_ior_lower90,
      log_ior_upper90 = log_ior_upper90,
      log_ior_se = log_ior_se,
      ior_values = ior_values,
      n_stages_significant = n_stages_significant,
      max_ior = max_ior,
      mean_ior = mean_ior,
      model_aic = AIC(modelo),
      model_deviance = deviance(modelo),
      formula_used = formula_parts,  # Guardo fórmula usada
      nichd_spline = nichd_spline,
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
      n_events = n_events, 
      n_coadmin = n_coadmin,
      error_msg = e$message,
      formula_attempted = formula_parts
    ))
  })
}


################################################################################
# Función de descripción del modelo usado
################################################################################

# Genera descripción del modelo ajustado
#
# Usa model_result (lista resultado de fit_differential_gam)
# Da string con descripción del modelo

describe_model_config <- function(model_result) {
  parts <- c()
  
  # tipo de spline
  parts <- c(parts, sprintf("Spline: %s", model_result$bs_type))
  
  # select
  if (model_result$select) {
    parts <- c(parts, "con penalización hasta cero")
  }
  
  # NICHD
  if (model_result$nichd_spline) {
    parts <- c(parts, "stage: spline")
  } else {
    parts <- c(parts, "stage: lineal")
  }
  
  # Splines para drogas individuales
  if (model_result$spline_individuales) {
    parts <- c(parts, "efectos individuales: splines")
  } else {
    parts <- c(parts, "efectos individuales: lineales")
  }
  
  # sex
  if (model_result$include_sex) {
    if (model_result$include_stage_sex) {
      parts <- c(parts, "sex: con interacción stage")
    } else {
      parts <- c(parts, "sex: efecto principal")
    }
  }
  
  paste(parts, collapse = " | ")
}

################################################################################
# Función para construir tripletes
################################################################################

# Genera tripletes (drugA, drugB, event) para un reporte individual
#
# Parámetros:
# dr: vector ids de drogas en el reporte (atc_concept_id)
# ev: vector ids de eventos en el reporte (meddra_concept_id)
# rid: id del reporte (safetyreportid)
# nichd_stage: etapa NICHD numérica del reporte
# 
# Return:
# data.table con columnas: safetyreportid, drugA, drugB, meddra, nichd_num
# 
# Implementación
# Genera todas las combinaciones de pares de drogas (si >= 2 drogas)
# Cruza cada par con cada evento del reporte
# Aplica orden: drugA <= drugB (para evitar duplicados por orden inverso)
# Devuelve NULL si el reporte NO tiene >= 2 drogas o eventos

make_triplets_per_report <- function(dr, ev, rid, nichd_stage) {
  
  if (length(dr) < 2 || length(ev) < 1) return(NULL)
  
  dr <- unique(dr)
  ev <- unique(ev)
  
  # Genero combinaciones de pares de drogas
  if (length(dr) == 2) {
    combs <- matrix(dr, nrow = 1)
  } else {
    combs <- t(combn(dr, 2))
  }
  
  # Orden: min primero
  combs <- t(apply(combs, 1, function(x) c(min(x), max(x))))
  
  n_combs <- nrow(combs)
  n_events <- length(ev)
  
  data.table(
    safetyreportid = rid,
    drugA = rep(combs[,1], times = n_events),
    drugB = rep(combs[,2], times = n_events),
    meddra = rep(ev, each = n_combs),
    nichd_num = nichd_stage
  )
}


################################################################################
# Formula de permutación para distribución nula
################################################################################

# Permuta EVENTOS entre reportes
# 
# Parámetros:
# perm_events: si TRUE, permuta eventos (meddra_concept_id)
# perm_drugs: si TRUE, permuta drogas también entre reportes
#
# Utiliza pool de reportes seleccionados en 01_augmentation (pool_meta)
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
# Funciones de utilidad
################################################################################

# Normaliza resultado de triplete con listas vacías por defecto
#
# Return:
# Lista normalizada con vectores por defecto en lugar de NULL
# 
# Previene errores de rbindlist cuando algunos campos de lista son NULL
# Aparentemente necesario cuando uso paralelización por si falla la convergencia de algunos tripletes


normalize_triplet_result <- function(result) {
  
  # campo tipo lista de columnas que deben existir
  list_fields <- c("stage", "log_ior", "log_ior_lower90", "ior_values")
  
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
  
  return(result)
}


################################################################################
# Funciones de expansión y clasificación
################################################################################
# Expande datos de formato ancho (1 fila = 1 triplete) a largo (7 filas = 7 etapas)
# clasifica señales según criterios de detección por método
#
# Parámetros:
# dt: data.table con tripletes en formato ancho (columnas tipo lista)
# method_col: nombre de columna con log_ior_lower90 a usar ("log_ior_lower90" o "log_ior_classic_lower90")
# type_label: identificador del método ("GAM" o "Clásico")
# 
# Return:
# data.table expandido con clasificaciones de señal por etapa y método
# 
# Implementación:
# -Expande listas por etapa (stage, log_ior_lower90)
# -Aplica umbrales de distribución nula
# -Clasifica señales con criterio diferencial por método:
#  GAM: doble criterio (nominal + null model)
#  IOR clásico: solo criterio nominal (IC90% > 0)


expand_and_classify <- function(dt, method_col, type_label) {
  
  by_cols <- c("triplet_id", "type")
  if ("dynamic" %in% names(dt)) by_cols <- c(by_cols, "dynamic")
  if ("fold_change" %in% names(dt)) by_cols <- c(by_cols, "fold_change")
  
  expanded <- dt[, {
    data.table(
      stage = unlist(stage),
      log_ior_lower90 = unlist(get(method_col))
    )
  }, by = by_cols]
  
  if (!"dynamic" %in% names(expanded)) {
    expanded[, dynamic := NA_character_]
  }
  if (!"fold_change" %in% names(expanded)) {
    expanded[, fold_change := NA_real_]
  }
  
  expanded <- merge(
    expanded,
    null_thresholds[, .(stage, threshold)],
    by = "stage", all.x = TRUE
  )
  
  expanded[, `:=`(
    nominal_sig = log_ior_lower90 > 0,
    nullmodel_sig = log_ior_lower90 > threshold,
    
    # condicional según el método
    signal_detected = if (type_label == "Clásico") {
      log_ior_lower90 > 0  # Solo criterio nominal para el IOR clásico
    } else {
      (log_ior_lower90 > 0) & (log_ior_lower90 > threshold) # doble umbral para GAM
    },
    
    method = type_label
  )]
  
  return(expanded)
}
                   
################################################################################
# Funciones de métricas de validación
################################################################################
# Calcula matriz de confusión y métricas derivadas
#
# usa expand_and_classify (data.table "expanded" con: triplet_id, signal_detected)
# Da data.table con método y métricas (sensibilidad, especificidad, ppv, npv, f1, precisión)

calculate_metrics <- function(expanded_data, method_name) {
  
  by_triplet <- expanded_data[, .(
    signal_detected = any(signal_detected, na.rm = TRUE),
    n_stages_sig = sum(signal_detected, na.rm = TRUE)
  ), by = .(triplet_id, true_label)]
  
  tp <- sum(by_triplet$signal_detected == TRUE & by_triplet$true_label == 1)
  fn <- sum(by_triplet$signal_detected == FALSE & by_triplet$true_label == 1)
  fp <- sum(by_triplet$signal_detected == TRUE & by_triplet$true_label == 0)
  tn <- sum(by_triplet$signal_detected == FALSE & by_triplet$true_label == 0)
  
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  f1_score <- 2 * (ppv * sensitivity) / (ppv + sensitivity)
  
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


################################################################################
# Cálculo de IOR clásico (sin modelo)
################################################################################

# Calcula IOR clásico usando tablas 2x2 por etapa
#
# parámetros:
# drugA_id: id droga A (ATC concept_id)
# drugB_id: id droga B (ATC concept_id)
# event_id: id del evento adverso (MedDRA concept_id)
# ade_data: data.table con dataset (puede ser aumentado)
# 
# return: Lista con stage, ior_classic, ior_classic_lower90, ior_classic_upper90
# 
# Implementación:
# Para cada etapa j:
# 1. Construye tabla 2x2:
#    - a: evento + coadmin
#    - b: sin evento + coadmin
#    - c: evento + sin coadmin
#    - d: sin evento + sin coadmin
# 2. OR_11 = (a/b) / (c/d) = ad/bc
# 3. OR_10 = reportes con A solo
# 4. OR_01 = reportes con B solo
# 5. OR_00 = reportes sin A ni B
# 6. IOR = (OR_11 × OR_00) / (OR_10 × OR_01)
# 7. IC90% usando método de Woolf (log scale)

calculate_classic_ior <- function(drugA_id, drugB_id, event_id, ade_data) {
  
  # Identificación de reportes
  reportes_droga_a <- unique(ade_data[atc_concept_id == drugA_id, safetyreportid])
  reportes_droga_b <- unique(ade_data[atc_concept_id == drugB_id, safetyreportid])
  reportes_ea <- unique(ade_data[meddra_concept_id == event_id, safetyreportid])
  
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
      # Datos diagnósticos
      n_11_evento = n_11_evento,
      n_11_total = n_11_evento + n_11_no_evento
    )
    
  }, by = nichd_num]
  
  # Ordenar por etapa
  setorder(resultados_por_etapa, nichd_num)
  
  return(list(
    success = TRUE,
    results_by_stage = resultados_por_etapa
  ))
}







