library(data.table)
library(tidyverse)
set.seed(42)

# ---------- Parámetros (ajusta según necesites) ----------
ruta_ade_raw <- "./ade_raw.csv"
min_reports_triplet <- 10
n_pos <- 1000
n_neg <- 10000
lambda_fc <- 0.75
dinamicas <- c("uniform","increase","decrease","plateau","inverse_plateau")
efecto_base <- 0.3
# ------------------------------------------------------------

# ---------- FUNCIONES AUXILIARES  ----------

dynamic_fun <- function(type, stages) {
  type <- as.character(type)
  stages <- as.character(stages)
  n <- length(stages)
  if (n <= 0) stop("dynamic_fun: 'stages' vacío o NULL")
  
  scaled_x <- seq(-2, 2, length.out = n)
  
  switch(type,
         "uniform" = rep(0, n),
         
         "increase" = tanh(scaled_x),
         
         "decrease" = -tanh(scaled_x),
         
         "plateau" = {
           vals <- -tanh(scaled_x)^2
           vals <- (vals - min(vals)) / (max(vals) - min(vals)) * 2 - 1  # reescala a [-1, 1]
           vals
         },
         
         "inverse_plateau" = {
           # Siempre bajo al inicio, luego crece
           # Punto de crecimiento: 30% del recorrido
           growth_point <- ceiling(n * 0.3)
           vals <- abs(tanh(scaled_x))
           vals
         },
         
         stop("dinámica desconocida: ", type)
  )
}

# fold-change sampling (negexp truncated a [1,10] como Giangreco)
sample_fold_change <- function(n, lambda = lambda_fc, max_fc = 10) {
  x <- rexp(n, rate = lambda)
  fc <- 1 + x
  fc <- pmin(fc, max_fc)
  return(fc)
}

# ------------------------------------------------------------
# testeo de dinámicas
# ------------------------------------------------------------

test_dynamics <- function() {
  niveles_test <- c("term_neonatal","infancy","toddler","early_childhood",
                    "middle_childhood","early_adolescence","late_adolescence")
  
  dynamics_list <- c("uniform","increase","decrease","plateau","inverse_plateau")
  
  plot_data <- data.table()
  
  for (dyn in dynamics_list) {
    vals <- dynamic_fun(dyn, niveles_test)
    plot_data <- rbind(plot_data, data.table(
      dynamic = dyn,
      stage = factor(niveles_test, levels = niveles_test),
      stage_num = 1:7,
      value = vals
    ))
  }
  
  p <- ggplot(plot_data, aes(x = stage_num, y = value, color = dynamic, group = dynamic)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = 1:7, labels = niveles_test) +
    labs(
      title = "Verificación de Dinámicas de Interacción",
      x = "Etapa NICHD", 
      y = "Valor de f(stage)",
      color = "Dinámica"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~dynamic, ncol = 3, scales = "free_y")
  
  print(p)
  ggsave("dynamic_patterns_verification.png", p, width = 12, height = 8)
  
  # También imprimir valores
  print(dcast(plot_data, stage ~ dynamic, value.var = "value"))
}

# Ejecutar antes del augmentation
test_dynamics()



# ------------------------------------------------------------
# 1) CARGAR DATOS
# ------------------------------------------------------------
message("Cargando datos...")
ade_raw_dt <- fread(ruta_ade_raw)

# Validar columnas requeridas
cols_req <- c("safetyreportid","atc_concept_id","meddra_concept_id","nichd")
if (!all(cols_req %in% colnames(ade_raw_dt))) {
  stop("Faltan columnas requeridas: ", 
       paste(setdiff(cols_req, colnames(ade_raw_dt)), collapse = ", "))
}

# Normalizar nichd a factor ordenado
niveles_nichd <- c("term_neonatal","infancy","toddler","early_childhood",
                   "middle_childhood","early_adolescence","late_adolescence")
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

message("Dataset cargado: ", nrow(ade_raw_dt), " filas, ", 
        length(unique(ade_raw_dt$safetyreportid)), " reportes únicos")

# ------------------------------------------------------------
# 2) CONSTRUIR LISTA DE TRIPLETES CANDIDATOS
# ------------------------------------------------------------
message("Construyendo tripletes candidatos...")

# Tablas de drogas y eventos por reporte
drugs_by_report <- unique(ade_raw_dt[!is.na(atc_concept_id), 
                                     .(safetyreportid, atc_concept_id)])
setkey(drugs_by_report, safetyreportid)

events_by_report <- unique(ade_raw_dt[!is.na(meddra_concept_id), 
                                      .(safetyreportid, meddra_concept_id)])
setkey(events_by_report, safetyreportid)

# Preparar listas por reporte
reports <- unique(ade_raw_dt[, .(safetyreportid, nichd)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), 
                               by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), 
                                 by = safetyreportid]

report_combo <- merge(reports, drugs_list, by = "safetyreportid", all.x = TRUE)
report_combo <- merge(report_combo, events_list, by = "safetyreportid", all.x = TRUE)
report_combo[is.na(drugs), drugs := list(integer(0))]
report_combo[is.na(events), events := list(integer(0))]

# Función optimizada para generar tripletes
make_triplets_per_report <- function(dr, ev, rid, nichd_stage) {
  if (length(dr) < 2 || length(ev) < 1) return(NULL)
  
  dr <- unique(dr)
  ev <- unique(ev)
  
  # Todas las combinaciones de 2 drogas (sin orden)
  if (length(dr) == 2) {
    combs <- matrix(dr, nrow = 1)
  } else {
    combs <- t(combn(dr, 2))
  }
  
  # Expandir por cada evento
  n_combs <- nrow(combs)
  n_events <- length(ev)
  
  out <- data.table(
    safetyreportid = rid,
    drugA = rep(combs[,1], times = n_events),
    drugB = rep(combs[,2], times = n_events),
    meddra = rep(ev, each = n_combs),
    nichd = nichd_stage
  )
  
  return(out)
}

# Generar tripletes (con barra de progreso)
message("Enumerando tripletes observados...")
triplets_list <- vector("list", nrow(report_combo))

pb <- txtProgressBar(max = nrow(report_combo), style = 3)
for (i in seq_len(nrow(report_combo))) {
  rowi <- report_combo[i]
  triplets_list[[i]] <- make_triplets_per_report(
    rowi$drugs[[1]], 
    rowi$events[[1]], 
    rowi$safetyreportid, 
    rowi$nichd
  )
  if (i %% 10000 == 0) setTxtProgressBar(pb, i)
}
close(pb)

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

if (nrow(triplets_dt) == 0) {
  stop("No se encontraron reportes con >=2 drogas y >=1 evento")
}

message("Tripletes generados: ", nrow(triplets_dt))

# Conteos por triplete único
trip_counts <- unique(triplets_dt[, .(drugA, drugB, meddra, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra)
]

# Filtrar candidatos
candidatos_pos <- trip_counts[N >= min_reports_triplet]
message("Tripletes candidatos (>= ", min_reports_triplet, " reportes): ", 
        nrow(candidatos_pos))

# ------------------------------------------------------------
# 3) SELECCIONAR POSITIVOS
# ------------------------------------------------------------
if (nrow(candidatos_pos) < n_pos) {
  stop("Insuficientes candidatos positivos. Disponibles: ", nrow(candidatos), 
       ", requeridos: ", n_pos + n_neg)
}

set.seed(123)
candidatos_pos <- candidatos_pos[sample(.N)]

positivos_sel <- candidatos_pos[1:n_pos, .(drugA, drugB, meddra, N)]
positivos_sel[, type := "positivo"]

message("Seleccionados: ", n_pos, " positivos, ")

# ------------------------------------------------------------
# 4) SELECCIONAR NEGATIVOS
# ------------------------------------------------------------
message("\nSeleccionando ", n_neg, " reportes NEGATIVOS aleatorios...")

# Identificar drogas y eventos en positivos (para excluir)
drugs_in_pos <- unique(c(positivos_sel$drugA, positivos_sel$drugB))
events_in_pos <- unique(positivos_sel$meddra)

message("Excluyendo: ", length(drugs_in_pos), " drogas y ", 
        length(events_in_pos), " eventos de positivos")

# Identificar reportes que NO contienen ninguna droga o evento de positivos
reports_with_pos_drugs <- unique(ade_raw_dt[atc_concept_id %in% drugs_in_pos, safetyreportid])
reports_with_pos_events <- unique(ade_raw_dt[meddra_concept_id %in% events_in_pos, safetyreportid])
reports_to_exclude <- unique(c(reports_with_pos_drugs, reports_with_pos_events))

# Reportes candidatos para negativos
all_reports <- unique(ade_raw_dt$safetyreportid)
candidate_reports <- setdiff(all_reports, reports_to_exclude)

message("Reportes disponibles como negativos: ", length(candidate_reports), 
        " / ", length(all_reports))

if (length(candidate_reports) < n_neg) {
  warning("Solo hay ", length(candidate_reports), " reportes elegibles, pero se pidieron ", n_neg)
  n_neg <- length(candidate_reports)
}

# Samplear reportes negativos
set.seed(456)
selected_neg_reports <- sample(candidate_reports, n_neg, replace = FALSE)

# Crear tabla de reportes negativos (sin agrupar en tripletes)
negativos_sel <- data.table(
  safetyreportid = selected_neg_reports,
  type = "negativo"
)

message("\nSeleccionados: ", nrow(negativos_sel), " reportes negativos")
message("  - Sin overlap de drogas/eventos con positivos")

# ------------------------------------------------------------
# 4) CREAR COPIA AUMENTADA
# ------------------------------------------------------------
message("Creando copia del dataset para augmentation...")
ade_aug <- copy(ade_raw_dt)
ade_aug[, simulated_flag := FALSE]  # marcar datos originales

# ------------------------------------------------------------
# 5) FUNCIÓN DE INYECCIÓN MEJORADA
# ------------------------------------------------------------
inject_dynamic_triplet <- function(drugA_id, drugB_id, event_id, 
                                   dynamic_type, effect_size = efecto_base) {
  
  # Encontrar reportes con ambas drogas
  reports_A <- unique(ade_aug[atc_concept_id == drugA_id, safetyreportid])
  reports_B <- unique(ade_aug[atc_concept_id == drugB_id, safetyreportid])
  reports_both <- intersect(reports_A, reports_B)
  
  if (length(reports_both) == 0) {
    warning("No hay reportes con ambas drogas: ", drugA_id, " y ", drugB_id)
    return(0)
  }
  
  # Identificar reportes que ya tienen el evento
  event_in_report <- unique(
    ade_raw_dt[meddra_concept_id == event_id, safetyreportid]
  )
  
  # Crear tabla de trabajo con reportes objetivo
  dt_reports <- unique(ade_raw_dt[
    safetyreportid %in% reports_both, 
    .(safetyreportid, nichd)
  ])
  dt_reports[, event_present := as.integer(safetyreportid %in% event_in_report)]
  
  # Probabilidad base por etapa
  probs_base <- dt_reports[, .(p_base = mean(event_present, na.rm = TRUE)), 
                            by = nichd]
  probs_base[, p_base := pmax(p_base, 1e-4)]  # evitar ceros
  
  # Sample fold-change
  fc <- sample_fold_change(1)
  
  # Obtener dinámica normalizada
  f_vec <- dynamic_fun(dynamic_type, niveles_nichd)
  f_min <- min(f_vec)
  f_max <- max(f_vec)
  
  if (f_max - f_min < 1e-10) {
    f_scaled <- rep(0.5, length(f_vec))  # uniform case
  } else {
    f_scaled <- (f_vec - f_min) / (f_max - f_min)
  }
  
  # Merge y calcular nuevas probabilidades
  merged_r <- merge(dt_reports, probs_base, by = "nichd", all.x = TRUE)
  merged_r[, nichd_idx := as.integer(nichd)]
  merged_r[, f_stage := f_scaled[nichd_idx]]
  
  # Fórmula: p_new = p_base * (1 + (fc - 1) * f_stage * effect_size)
  merged_r[, p_new := p_base * (1 + (fc - 1) * f_stage * effect_size)]
  merged_r[, p_new := pmin(p_new, 0.95)]  # cap máximo
  
  # NO FILTRAR - inyectar en todos los reportes
  reports_to_inject <- merged_r
  
  # Probabilidad ajustada: mayor efecto si ya no tiene el evento
  reports_to_inject[, p_inject := p_new * (1 + 2 * (1 - event_present))]
  reports_to_inject[, add_event := rbinom(.N, 1, p_inject)]
  
  # Crear nuevas filas SOLO para quienes no lo tenían y add_event == 1
  new_rows <- reports_to_inject[event_present == 0 & add_event == 1, .(
    safetyreportid,
    atc_concept_id = NA_integer_,
    meddra_concept_id = event_id,
    nichd,
    nichd_num = as.integer(nichd),
    simulated_flag = TRUE
  )]
  
  if (nrow(new_rows) == 0) return(0)
  
  # Añadir al dataset aumentado
  ade_aug <<- rbindlist(list(ade_aug, new_rows), use.names = TRUE, fill = TRUE)
  
  return(nrow(new_rows))
}

# ------------------------------------------------------------
# 6) APLICAR INYECCIÓN
# ------------------------------------------------------------
set.seed(9427)
din_for_triplets <- sample(dinamicas, n_pos, replace = TRUE)

pos_meta <- positivos_sel[, .(drugA, drugB, meddra, N)]
pos_meta[, dynamic := din_for_triplets]

message("\nInyectando dinámicas en ", n_pos, " tripletes positivos...")
message("Distribución de dinámicas: ")
print(table(din_for_triplets))

injected_counts <- integer(n_pos)
failed_injections <- integer(0)

pb <- txtProgressBar(max = n_pos, style = 3)
for (i in seq_len(n_pos)) {
  rowi <- pos_meta[i]
  
  tryCatch({
    injected_counts[i] <- inject_dynamic_triplet(
      rowi$drugA, 
      rowi$drugB, 
      rowi$meddra, 
      rowi$dynamic
    )
  }, error = function(e) {
    message("\nError en triplete ", i, ": ", e$message)
    failed_injections <<- c(failed_injections, i)
    injected_counts[i] <<- 0
  })
  
  setTxtProgressBar(pb, i)
}
close(pb)

# ------------------------------------------------------------
# 7) RESUMEN DE RESULTADOS
# ------------------------------------------------------------
message("\n" , paste(rep("=", 60), collapse = ""))
message("RESUMEN DE AUGMENTATION")
message(paste(rep("=", 60), collapse = ""))
message("Dataset original: ", nrow(ade_raw_dt), " filas")
message("Dataset aumentado: ", nrow(ade_aug), " filas")
message("Nuevas filas simuladas: ", sum(ade_aug$simulated_flag))
message("\nInyecciones por triplete:")
message("  - Exitosas: ", sum(injected_counts > 0))
message("  - Sin cambios: ", sum(injected_counts == 0))
message("  - Fallidas: ", length(failed_injections))
message("  - Total eventos inyectados: ", sum(injected_counts))
message("  - Promedio por triplete exitoso: ", 
        round(mean(injected_counts[injected_counts > 0]), 1))

pos_metadata <- cbind(pos_meta, injected = injected_counts)
pos_metadata[, success := injected > 0]

message("\nDistribución de inyecciones por dinámica:")
print(pos_metadata[, .(
  n_tripletes = .N,
  n_exitosos = sum(success),
  total_eventos = sum(injected),
  promedio = round(mean(injected), 1)
), by = dynamic])

message("\nNegativos: ", nrow(negativos_sel), " tripletes de reportes independientes")
message("  (Reportes sin overlap de drogas/eventos con positivos)")

# ------------------------------------------------------------
# 9) GUARDAR RESULTADOS
# ------------------------------------------------------------
message("\nGuardando resultados...")

fwrite(ade_aug, "ade_augmented.csv")
message("  ✓ Dataset aumentado: ade_augmented.csv")

fwrite(pos_metadata, "positive_triplets_metadata.csv")
message("  ✓ Metadata positivos: positive_triplets_metadata.csv")

fwrite(negativos_sel, "negative_triplets_metadata.csv") # después sacar la palabra "triplets" porque ya no son triplets
message("  ✓ Metadata negativos: negative_triplets_metadata.csv")

ground_truth_positive <- rbind(
  pos_metadata[, .(drugA, drugB, meddra, type = "positive", 
                   dynamic, n_injected = injected)])

ground_truth_negative <- rbind(
    negativos_sel[, .(selected_neg_reports, type = "negative", 
                    dynamic = NA_character_, n_injected = 0L)]
)

fwrite(ground_truth_positive, "ground_truth_positive.csv")

fwrite(ground_truth_negative, "ground_truth_negative.csv")

message("\n✓ Augmentation completado exitosamente!")
message(paste(rep("=", 60), collapse = ""))
