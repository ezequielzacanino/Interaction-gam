################################################################################
# Script de visualización de gráficos facetados
# Script: 41_graphs.R
################################################################################

source("00_functions.R", local = TRUE)

output_dir <- paste0("./results/", suffix, "/metrics_results/")
fig_output_dir <- paste0(output_dir, "facet_figures/")

# Crea directorio de salida si no existe
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Preprocesamiento
################################################################################

# Traducción de etapas a español
nichd_labels <- c(
  "term_neonatal" = "Neonato a term.",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temp.",
  "late_adolescence" = "Adolescencia tardía"
)

# Etiquetas de métricas
metric_labels <- c(
  "AUC" = "AUC",
  "sensitivity" = "Sensibilidad",
  "specificity" = "Especificidad",
  "PPV" = "VPP",
  "NPV" = "VPN",
  "F1" = "F1-Score"
)

# Etiquetas de dinámicas
dynamic_labels <- c(
  "uniform" = "Uniforme",
  "increase" = "Aumento",
  "plateau" = "Meseta",
  "decrease" = "Disminución",
  "inverse_plateau" = "Valle"
)

# Pares de métodos a comparar
method_pairs <- list(
  list(
    name = "IOR",
    gam = "GAM-logIOR",
    classic = "Estratificado-IOR",
    gam_label = "GAM-IOR",
    classic_label = "IOR-Estratificado"
  ),
  list(
    name = "RERI",
    gam = "GAM-RERI",
    classic = "Estratificado-RERI",
    gam_label = "GAM-RERI",
    classic_label = "RERI-Estratificado"
  )
)

# versiones de datasets a procesar
dataset_versions <- c("original", "filtered", "intersection")

# Helper: título legible según versión
version_title_label <- function(version) {
  switch(version,
    "original" = "Dataset Original",
    "filtered" = "Dataset Filtrado por Poder",
    "intersection" = "Dataset de Intersección",
    version
  )
}

# colores y linetypes según par de métodos
make_scales <- function(pair) {
  list(
    color = setNames(c("#16A085", "#C0392B"), c(pair$gam, pair$classic)),
    linetype = setNames(c("solid", "dashed"), c(pair$gam, pair$classic)),
    labels = setNames(c(pair$gam_label, pair$classic_label), c(pair$gam, pair$classic))
  )
}

################################################################################
# Gráfico de dinámicas
################################################################################

# Gráfico de funciones tangenciales por dinámica
dt_dyn_plot <- rbindlist(lapply(
  c("increase", "decrease", "plateau", "inverse_plateau"),
  function(d) data.table(
    stage = 1:7,
    valor = generate_dynamic(d),
    dinamica = factor(
      c(increase = "Creciente", decrease = "Decreciente",
        plateau = "Meseta", inverse_plateau = "Valle")[d],
      levels = c("Creciente", "Decreciente", "Meseta", "Valle")
    )
  )
))

dynamics <- ggplot(dt_dyn_plot, aes(x = stage, y = valor, color = dinamica)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~dinamica, nrow = 2) +
  scale_color_manual(values = c(
    "Creciente" = "#16A085",
    "Decreciente" = "#C0392B",
    "Meseta" = "#2980B9",
    "Valle" = "#8E44AD"
  ), guide = "none") +
  scale_x_continuous(breaks = 1:7, labels = nichd_labels, name = NULL) +
  scale_y_continuous(name = "Peso relativo en probabilidad de reporte",
                     limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +
  labs(title = NULL) +
  coord_cartesian(clip = "off") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 25)
  )
ggsave(paste0(output_dir, "dynamics.png"), dynamics, width = 15, height = 15, dpi = 300)


################################################################################
# carga de datos
################################################################################

metrics_global <- list(
  original = fread(paste0(output_dir, "metrics_global_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_global_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_global_intersection.csv"))
)

metrics_dynamic <- list(
  original = fread(paste0(output_dir, "metrics_dynamic_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_dynamic_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_dynamic_intersection.csv"))
)

metrics_stage <- list(
  original = fread(paste0(output_dir, "metrics_stage_original.csv")),
  filtered = fread(paste0(output_dir, "metrics_stage_filtered.csv")),
  intersection = fread(paste0(output_dir, "metrics_stage_intersection.csv"))
)

# proceso columna method_type + factor nichd/dynamic
for (v in dataset_versions) {
  metrics_global[[v]][, method_type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]

  dt <- metrics_dynamic[[v]]
  if ("dinamica" %in% names(dt)) setnames(dt, "dinamica", "dynamic")
  dt[, dynamic := factor(dynamic)]
  dt[, method_type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  metrics_dynamic[[v]] <- dt

  metrics_stage[[v]][, nichd := factor(nichd, levels = niveles_nichd)]
  metrics_stage[[v]][, method_type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
}

################################################################################
# Preparación de datos para facetado
################################################################################

# Prepara datos en formato largo para facetado
#
# Convierte las métricas de formato ancho a largo
#
# data.table con métricas por etapa
# lista con información del par de métodos (gam, classic)
# data.table en formato largo listo para ggplot2

prepare_facet_data <- function(dt, pair) {
  # Filtrar solo los métodos del par actual
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  # Métricas a incluir
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")
  # Convertir a formato largo
  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    
    # veo que las columnas existan
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) {
      return(NULL)
    }
    
    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      nichd = dt_filtered$nichd,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)
  # Crear etiqueta combinada de método
  dt_long[, method_label := ifelse(
    method == pair$gam,
    pair$gam_label,
    pair$classic_label
  )]
  # Ordenar métricas
  dt_long[, metric_label := factor(
    metric_label,
    levels = metric_labels
  )]
  return(dt_long)
}

# Convierte métricas de formato ancho a largo — por dinámica
prepare_facet_data_dynamic <- function(dt, pair) {
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")

  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) return(NULL)
    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      dynamic = dt_filtered$dynamic,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)

  dt_long[, method_label := ifelse(method == pair$gam, pair$gam_label, pair$classic_label)]
  dt_long[, metric_label := factor(metric_label, levels = metric_labels)]
  return(dt_long)
}

# Convierte métricas de formato ancho a largo globales (sin estratificación)
prepare_facet_data_global <- function(dt, pair) {
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")

  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    if (!(lower_col %in% names(dt_filtered)) || !(upper_col %in% names(dt_filtered))) return(NULL)
    data.table(
      method = dt_filtered$method,
      method_type = dt_filtered$method_type,
      reduction_pct = dt_filtered$reduction_pct,
      metric = m,
      metric_label = metric_labels[m],
      value = dt_filtered[[m]],
      lower = dt_filtered[[lower_col]],
      upper = dt_filtered[[upper_col]]
    )
  }), use.names = TRUE)

  dt_long[, method_label := ifelse(method == pair$gam, pair$gam_label, pair$classic_label)]
  dt_long[, metric_label := factor(metric_label, levels = metric_labels)]
  return(dt_long)
}

################################################################################
# Función de generación de gráficos facetados
################################################################################

# Genera gráfico facetado de métricas por etapa
#
# Filas: métricas (sensibilidad, especificidad, PPV, NPV, F1)
# Columnas: etapas NICHD
# Eje X: reducción del dataset 
#
# pair: Información del par de métodos
# version: nombre de la versión del dataset

plot_facet_metrics <- function(dt_long, pair, version) {
  
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("Sin datos para %s - %s", pair$name, version))
    return(NULL)
  }
  
  # Título según versión
  version_title <- switch(version,
    "original" = "Dataset Original",
    "filtered" = "Dataset Filtrado por Poder",
    "intersection" = "Dataset de Intersección",
    version
  )
  
  color_values <- setNames(
    c("#16A085", "#C0392B"), 
    c(pair$gam, pair$classic)
  )
  linetype_values <- setNames(
    c("solid", "dashed"), 
    c(pair$gam, pair$classic)
  )

  p <- ggplot(dt_long, aes(
    x = reduction_pct,
    y = value,
    color = method,
    group = method
  )) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 4,
      alpha = 0.8,
      position = position_dodge(width = 5)
    ) +
    facet_grid(
      metric_label ~ nichd,
      labeller = labeller(nichd = nichd_labels),
      scales = "free_y"
    ) +
    scale_color_manual(
      name = "Método",
      values = color_values,
      labels = setNames(c(pair$gam_label, pair$classic_label), c(pair$gam, pair$classic))
    ) +
    scale_linetype_manual(
      name = "Método",
      values = linetype_values,
      labels = setNames(c(pair$gam_label, pair$classic_label), c(pair$gam, pair$classic))
    ) +
    scale_x_continuous(
      breaks = seq(0, 90, by = 10),
      name = "Reducción del Dataset (%)"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      name = "Valor de la Métrica"
    ) +
    # Títulos
    labs(
      title = sprintf(
        "Métricas por Etapa - %s vs %s",
        pair$gam_label, pair$classic_label
      ),
      subtitle = version_title
    ) +
    # Tema
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),
      strip.background.x = element_rect(fill = "gray90", color = "black"),
      strip.background.y = element_rect(fill = "gray95", color = "black"),
      strip.text.x = element_text(face = "bold", size = 14, hjust = 0.5),
      strip.text.y = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 9, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 14),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
  
  return(p)
}

################################################################################
# Función principal de generación
################################################################################

# Genera todos los gráficos facetados para todas las combinaciones
# Itera sobre todas las versiones de datasets y pares de métodos

plot_facet_metrics <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("Sin datos para %s - %s (etapa)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_grid(metric_label ~ nichd, labeller = labeller(nichd = nichd_labels), scales = "free_y") +
    scale_color_manual(name = "Método", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Método", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Reducción del Dataset (%)") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), name = "Valor de la Métrica") +
    labs( title = sprintf("Métricas por Etapa - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version))
}

# Gráfico facetado: métricas x dinámicas
plot_facet_metrics_dynamic <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("Sin datos para %s - %s (dinámica)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_grid(metric_label ~ dynamic, labeller = labeller(dynamic = dynamic_labels), scales = "free_y") +
    scale_color_manual(name = "Método", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Método", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Reducción del Dataset (%)") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), name = "Valor de la Métrica") +
    labs( title = sprintf("Métricas por Dinámica - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version)) 
}

# Gráfico facetado: métricas globales (facet_wrap por métrica)
plot_facet_metrics_global <- function(dt_long, pair, version) {
  if (is.null(dt_long) || nrow(dt_long) == 0) {
    message(sprintf("Sin datos para %s - %s (global)", pair$name, version))
    return(NULL)
  }

  sc <- make_scales(pair)

  ggplot(dt_long, aes(x = reduction_pct, y = value, color = method, group = method)) +
    geom_point(size = 2.5, position = position_dodge(width = 5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.8, position = position_dodge(width = 5)) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    scale_color_manual(name = "Método", values = sc$color, labels = sc$labels) +
    scale_linetype_manual(name = "Método", values = sc$linetype, labels = sc$labels) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), name = "Reducción del Dataset (%)") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), name = "Valor de la Métrica") +
    labs(
      title = sprintf("Métricas Globales - %s vs %s", pair$gam_label, pair$classic_label),
      subtitle = version_title_label(version)
    ) 
}

################################################################################
# Guardado
################################################################################

save_plot <- function(p, file_suffix, width, height) {
  png_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".png")
  svg_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".svg")
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(svg_path, p, width = width, height = height, device = svglite)
}

# Gráficos facetados por etapa
generate_all_facet_plots <- function() {
  plots_generated <- 0
  for (version in dataset_versions) {
    dt_stage <- metrics_stage[[version]]
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data(dt_stage, pair)
      p <- plot_facet_metrics(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s", tolower(pair$name), version), width = 16, height = 12)
        plots_generated <- plots_generated + 1
      }
    }
  }
}

# Gráficos facetados por dinámica
generate_all_facet_plots_dynamic <- function() {
  plots_generated <- 0
  for (version in dataset_versions) {
    dt_dynamic <- metrics_dynamic[[version]]
    if (is.null(dt_dynamic)) next
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data_dynamic(dt_dynamic, pair)
      p <- plot_facet_metrics_dynamic(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s_dynamic", tolower(pair$name), version), width = 16, height = 12)
        plots_generated <- plots_generated + 1
      }
    }
  }
}

# Gráficos facetados globales (métricas globales sin estratificación)
generate_all_facet_plots_global <- function() {
  plots_generated <- 0
  for (version in dataset_versions) {
    dt_global <- metrics_global[[version]]
    if (is.null(dt_global)) next
    for (pair in method_pairs) {
      dt_long <- prepare_facet_data_global(dt_global, pair)
      p <- plot_facet_metrics_global(dt_long, pair, version)
      if (!is.null(p)) {
        save_plot(p, sprintf("%s_%s_global", tolower(pair$name), version), width = 12, height = 10)
        plots_generated <- plots_generated + 1
      }
    }
  }
}

# generación de gráficos
generate_all_facet_plots()
generate_all_facet_plots_dynamic()
generate_all_facet_plots_global()

################################################################################
# Generación de gráficos para tripletes positivos
################################################################################

# Rutas de datos para tripletes positivos
ruta_positive_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_info <- "./drug.csv"
output_dir_positive <- paste0(fig_output_dir, "positive_triplets/")

# positivos detectados en dataset original
ruta_network_graphs <- paste0("./results/", suffix, "/network/network_triplets_for_graphs.rds")

# verifica existencia de archivos requeridos
check_positive_files <- function() {
  archivos_requeridos <- c(ruta_positive_results, ruta_ade_raw, ruta_drug_info)
  nombres <- c("positive_triplets_results.rds", "ade_raw.csv", "drug.csv")
  return(TRUE)
}

################################################################################
# Funciones de carga y preprocesamiento
################################################################################

# Carga y preprocesa datos de tripletes positivos
# 
# Parámetros:
# ruta_results Ruta al archivo RDS con resultados de augmentation
# ruta_ade Ruta al archivo CSV con datos ADE originales
#
# return Lista con data.table de resultados y data.table de ADE procesado
#
# Filtra solo casos base (reduction_pct == 0) y exitosos

carga_datos_positive <- function(ruta_results, ruta_ade) {
  # Cargar resultados
  dt <- readRDS(ruta_results)
  
  # Filtrar solo base y exitosos
  dt <- dt[reduction_pct == 0 & model_success == TRUE & injection_success == TRUE]
  
  # Cargar ADE raw
  ade_raw_dt <- fread(ruta_ade)
  
  # Procesar sexo si está habilitado
  if (include_sex) {
    cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd", "sex")
    ade_raw_dt[, sex := toupper(trimws(sex))]
    ade_raw_dt[sex == "M", sex := "MALE"]
    ade_raw_dt[sex == "F", sex := "FEMALE"]
    ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
    
    sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
    message("  Distribución de sexo:")
    print(sex_summary)
  } else {
    cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")
  }
  
  message(sprintf("  Dataset ADE: %s filas", format(nrow(ade_raw_dt), big.mark = ",")))
  
  return(list(results = dt, ade_raw = ade_raw_dt))
}

# Crea tabla de traducción de IDs de drogas canónicos
# 
# Parámetros
# ruta_drug_info Ruta al archivo CSV con información de drogas
#
# return 
# data.table con mapeo atc_concept_id -> canonical_id
#
# Agrupa drogas por nombre base (antes de ;,)

translation_table <- function(ruta_drug_info) {
  drug_info_original <- fread(ruta_drug_info)
  drug_info_original[, atc_concept_id := as.character(atc_concept_id)]
  drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]
  
  # Crear mapeo canónico
  canonical_map <- drug_info_original[, .(
    canonical_id = min(atc_concept_id)
  ), by = base_name]
  
  translation_table <- merge(
    drug_info_original[, .(atc_concept_id, base_name)],
    canonical_map,
    by = "base_name"
  )
  
  message(sprintf("IDs originales: %d", uniqueN(translation_table$atc_concept_id)))
  message(sprintf("IDs únicos: %d", uniqueN(translation_table$canonical_id)))
  
  return(translation_table[, .(atc_concept_id, canonical_id)])
}

# Preprocesa datos aplicando traducción de ids y creando lookup
# 
# Parámetros:
# ade_raw_dt data.table con datos ADE originales
# translation_table data.table con mapeo de IDs
# 
# return data.table procesado con key establecida

translate_ade <- function(ade_raw_dt, translation_table) {
  # Convertir y merge
  ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]
  
  ade_raw_dt <- merge(
    ade_raw_dt,
    translation_table[, .(atc_concept_id, canonical_id)],
    by = "atc_concept_id",
    all.x = TRUE
  )
  
  # Aplicar traducción
  ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
  ade_raw_dt[, canonical_id := NULL]
  
  # Eliminar duplicados
  nrow_antes <- nrow(ade_raw_dt)
  cols_unicos <- c("safetyreportid", "atc_concept_id", "meddra_concept_id")
  if (include_sex) cols_unicos <- c(cols_unicos, "sex")
  
  ade_raw_dt <- unique(ade_raw_dt, by = cols_unicos)
  message(sprintf("  Duplicados eliminados: %d", nrow_antes - nrow(ade_raw_dt)))
  
  # Crear factor ordenado para nichd y variable numérica
  ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
  ade_raw_dt[, nichd_num := as.integer(nichd)]
  
  # Establecer key para lookups eficientes
  setkey(ade_raw_dt, atc_concept_id, meddra_concept_id, nichd_num)
  
  return(ade_raw_dt)
}

################################################################################
# cálculo de conteos
################################################################################

# Calcula conteos de reportes por etapa para un triplete específico
# 
# Parámetros
# drug_a: id droga A
# drug_b: id droga B  
# meddra_event: id del evento
# 
# return 
# data.table con conteos por etapa (n_a, n_b, n_ab, n_evento, n_evento_ab)
# 
# Esto es increiblemente ineficiente para el pipeline, debería estar hecho en 10_augmentation
# Pero eso implicaría volver a correr 10_augmentation (2 días y medio)

calculate_triplet_counts <- function(drug_a, drug_b, meddra_event, ade_dt) {
  # ids únicos por exposición
  ids_a <- unique(ade_dt[atc_concept_id == drug_a, safetyreportid])
  ids_b <- unique(ade_dt[atc_concept_id == drug_b, safetyreportid])
  ids_event <- unique(ade_dt[meddra_concept_id == meddra_event, safetyreportid])
  
  # intersecciones
  ids_ab <- intersect(ids_a, ids_b)
  ids_a_only <- setdiff(ids_a, ids_b)
  ids_b_only <- setdiff(ids_b, ids_a)
  
  # Función para contar por etapa
  count_per_stage <- function(ids_subset) {
    if (length(ids_subset) == 0) {
      return(data.table(nichd_num = 1:7, n = 0L))
    }
    counts <- ade_dt[safetyreportid %in% ids_subset, .(n = uniqueN(safetyreportid)), by = nichd_num]
    # Completa etapas faltantes
    etapas_completas <- data.table(nichd_num = 1:7)
    counts <- merge(etapas_completas, counts, by = "nichd_num", all.x = TRUE)
    counts[is.na(n), n := 0L]
    return(counts[order(nichd_num)])
  }
  
  # Calcula conteos
  n_a_dt <- count_per_stage(ids_a_only)
  n_b_dt <- count_per_stage(ids_b_only)
  n_ab_dt <- count_per_stage(ids_ab)
  n_event_dt <- count_per_stage(ids_event)
  n_event_ab_dt <- count_per_stage(intersect(ids_event, ids_ab))
  
  # Combino resultado
  result <- data.table(
    nichd_num = 1:7,
    n_a = n_a_dt$n,
    n_b = n_b_dt$n,
    n_ab = n_ab_dt$n,
    n_evento = n_event_dt$n,
    n_evento_ab = n_event_ab_dt$n
  )
  
  return(result)
}

################################################################################
# expansión de datos
################################################################################

# Expande un triplete con sus métricas y conteos
# 
# Parámetros
# row: fila individual de resultados
# de_dt: data.table con datos preprocesados
# 
# return
# data.table expandido con métricas y conteos
# 
# Extrae métricas de listas y calcula conteos desde ade_raw
# El problema es que los resultados de inyección no están en ade_raw
# Entonces se extraen de los diagnostics del rsd 

expand_triplets_counts <- function(row, ade_dt, precomputed_counts = NULL) {
  # Extraer etapas
  stages <- unlist(row$stage)
  
  # Métricas GAM
  gam_log_ior <- unlist(row$log_ior)
  gam_log_ior_lower <- unlist(row$log_ior_lower90)
  gam_log_ior_upper <- gam_log_ior + (gam_log_ior - gam_log_ior_lower)
  
  gam_reri <- unlist(row$reri_values)
  gam_reri_lower <- unlist(row$reri_lower90)
  gam_reri_upper <- unlist(row$reri_upper90)
  
  # Métricas Clásico
  cls_log_ior <- unlist(row$log_ior_classic)
  cls_log_ior_lower <- unlist(row$log_ior_classic_lower90)
  cls_log_ior_upper <- cls_log_ior + (cls_log_ior - cls_log_ior_lower)
  
  cls_reri <- unlist(row$RERI_classic)
  cls_reri_lower <- unlist(row$RERI_classic_lower90)
  cls_reri_upper <- unlist(row$RERI_classic_upper90)
  
  # bloque para admitir también los positivos de 40_network del dataset real
  if (!is.null(precomputed_counts)) {
    counts_dt <- precomputed_counts[triplet_id == row$triplet_id]
    if (nrow(counts_dt) == 0)
      counts_dt <- calculate_triplet_counts(row$drugA, row$drugB, row$meddra, ade_dt)
  } else if (!is.null(row$counts_by_stage[[1]])) {
    # Conteos embebidos directamente en la fila del objeto de red
    counts_dt <- row$counts_by_stage[[1]]
  } else {
    counts_dt <- calculate_triplet_counts(row$drugA, row$drugB, row$meddra, ade_dt)
  }
  
  # diagnostics[[1]]$injection_by_stage tiene nichd_num y N (inyectados por etapa)
    inj_by_stage <- tryCatch({
    diag <- row$diagnostics[[1]]
    if (!is.null(diag) && !is.null(diag$injection_by_stage)) {
      ibs <- merge(data.table(nichd_num = 1:7), diag$injection_by_stage, by = "nichd_num", all.x = TRUE)
      ibs[is.na(N), N := 0L]
      ibs[order(nichd_num)]$N
    } else { rep(0L, 7) }
  }, error = function(e) rep(0L, 7))
  
  # Sumar inyectados a n_evento_ab y n_evento
  counts_dt[, n_evento_ab := n_evento_ab + inj_by_stage]
  counts_dt[, n_evento := n_evento + inj_by_stage]

  # Verificar longitudes consistentes
  n <- min(length(stages), length(gam_log_ior), nrow(counts_dt))
  if (n == 0) return(NULL)
  
  # Combinar todo
  result <- data.table(
    triplet_id = row$triplet_id,
    dynamic = row$dynamic,
    stage_num = stages[1:n],
    # Métricas GAM
    gam_log_ior = gam_log_ior[1:n],
    gam_log_ior_lower = gam_log_ior_lower[1:n],
    gam_log_ior_upper = gam_log_ior_upper[1:n],
    gam_reri = gam_reri[1:n],
    gam_reri_lower = gam_reri_lower[1:n],
    gam_reri_upper = gam_reri_upper[1:n],
    # Métricas Clásico
    cls_log_ior = cls_log_ior[1:n],
    cls_log_ior_lower = cls_log_ior_lower[1:n],
    cls_log_ior_upper = cls_log_ior_upper[1:n],
    cls_reri = cls_reri[1:n],
    cls_reri_lower = cls_reri_lower[1:n],
    cls_reri_upper = cls_reri_upper[1:n],
    # Conteos
    n_a = counts_dt$n_a[1:n],
    n_b = counts_dt$n_b[1:n],
    n_ab = counts_dt$n_ab[1:n],
    n_evento = counts_dt$n_evento[1:n],
    n_evento_ab = counts_dt$n_evento_ab[1:n]
  )
  
  return(result)
}

################################################################################
# visualización
################################################################################

# Genera gráfico de métricas con conteos para un triplete
# 
# Parámetros:
# plot_dt: datos expandidos del triplete
# metric_col: columna con valores de métrica
# lower_col: columna con límite inferior IC
# upper_col: columna con límite superior IC
# y_label: etiqueta para eje Y principal
# file_suffix: Sufijo para nombre de archivo
# max_count: máximo esperado para escala logarítmica secundaria

graph_metrics_counts <- function(plot_dt, metric_col, lower_col, upper_col,
                                  y_label, file_suffix, y_limit = 10, max_count = 5000,
                                  plot_title = NULL, plot_subtitle = NULL) {
  # datos en formato largo
  counts_long <- melt(
    plot_dt,
    id.vars = "stage_num",
    measure.vars = c("n_a", "n_b", "n_ab", "n_evento", "n_evento_ab"),
    variable.name = "metric",
    value.name = "count"
  )
  
  # Etiquetas para leyenda
  metric_labels <- c(
    n_a = "A", n_b = "B", n_ab = "A-B",
    n_evento = "Evento", n_evento_ab = "A-B-Evento"
  )
  counts_long[, metric_label := factor(metric_labels[as.character(metric)], levels = unname(metric_labels))]
  
  # calculo posiciones X para barras (con dodging manual)
  bar_w <- 0.12
  counts_long[, metric_idx := as.integer(metric)]
  counts_long[, x_center := stage_num + (metric_idx - 3) * bar_w]
  counts_long[, xmin := x_center - bar_w/2]
  counts_long[, xmax := x_center + bar_w/2]
  
  # Transformación logarítmica para eje secundario
  # Mapeo: 1 -> -y_limit, max_count -> y_limit
  transforma_a_log <- function(count) {
    count_adj <- pmax(count, 1)
    -y_limit + (log10(count_adj) / log10(max_count)) * (2 * y_limit)
  }
  
  transforma_desde_log <- function(y) {
    ratio <- (y + y_limit) / (2 * y_limit)
    10^(ratio * log10(max_count))
  }
  
  counts_long[, ymax := transforma_a_log(count)]
  counts_long[, ymin := -y_limit]
  
  # gráfico
  p <- ggplot() +
    # Barras de conteos
    geom_rect(
      data = counts_long,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = metric_label),
      alpha = 0.6, color = NA
    ) +
    # Línea de cero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # IC 
    geom_ribbon(
      data = plot_dt,
      aes(x = stage_num, ymin = .data[[lower_col]], ymax = .data[[upper_col]]),
      alpha = 0.3, fill = "#2c3e50"
    ) +
    # Línea del métrico
    geom_line(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", linewidth = 1.2
    ) +
    # puntos de métrica
    geom_point(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", size = 4, fill = "white", shape = 21, stroke = 1.5
    ) +
    # Escalas
    scale_x_continuous(
      breaks = 1:7,
      labels = nichd_labels,
      name = NULL
    ) +
    scale_y_continuous(
      name = y_label,
      limits = c(-y_limit, y_limit),
      breaks = {
        # breaks visibles independientemente del rango
        step <- 10^floor(log10(y_limit / 4))
        step <- step * if (y_limit / step > 8) 2 else 1
        seq(-y_limit, y_limit, by = step)
      },
      sec.axis = sec_axis(
        trans = transforma_desde_log,
        name = "N° Reportes (log)",
        breaks = c(1, 10, 100, 1000, 10000),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)
      )
    ) +
    scale_fill_brewer(palette = "Set2", name = "") +
    # Etiquetas
    labs(
      title = if (!is.null(plot_title)) plot_title
      else paste0(plot_dt$dynamic[1], " | Triplete: ", plot_dt$triplet_id[1]),
    subtitle = if (!is.null(plot_subtitle)) plot_subtitle else "A + B"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y.right = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8)
    )
  return(p)
}

################################################################################
# Generación de gráficos para tripletes positivos
################################################################################

# Guarda un gráfico de triplete positivo.
#
# Parámetros:
# p: objeto ggplot
# triplet_id: id numérico del triplete
# dynamic: etiqueta de dinámica (usada en el nombre de archivo)
# suffix: sufijo identificador del método/métrica (e.g. "GAM_LogIOR")
# dir: directorio de salida; por defecto output_dir_positive

save_positive_graph <- function(p, triplet_id, dynamic, suffix,
                                 dir = output_dir_positive) {
  safe_name <- sprintf("Triplete_%d_%s_%s", triplet_id, dynamic, suffix)
  safe_name <- gsub("[^a-zA-Z0-9._-]", "_", safe_name)
  ggsave(
    filename = file.path(dir, paste0(safe_name, ".png")),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# Genera gráficos para tripletes positivos del pipeline de augmentation.
#
# fuente : positive_triplets_results.rds  (reduction_pct == 0, línea de base)
# conteos: calculados on-the-fly desde ade_raw.csv con calculate_triplet_counts()
# métricas: GAM-LogIOR, GAM-RERI, Clásico-LogIOR, Clásico-RERI

generate_positive_graphs_from_results <- function() {
  dir.create(output_dir_positive, showWarnings = FALSE, recursive = TRUE)

  datos <- carga_datos_positive(ruta_positive_results, ruta_ade_raw)
  dt <- datos$results
  trans_table <- translation_table(ruta_drug_info)
  ade_processed <- translate_ade(datos$ade_raw, trans_table)

  pb <- txtProgressBar(min = 0, max = nrow(dt), style = 3)
  generated_graphs <- 0

  for (i in seq_len(nrow(dt))) {
    tryCatch({
      row <- dt[i]

      ptitle <- paste0(row$dynamic, " | Triplete: ", row$triplet_id)
      psubtitle <- "A + B"

      plot_dt <- expand_triplets_counts(row, ade_dt = ade_processed, precomputed_counts = NULL)
      if (is.null(plot_dt) || nrow(plot_dt) == 0) {
        message(sprintf("\n  Salteando triplete %d: sin datos expandibles", i))
        next
      }

      if (any(is.finite(plot_dt$gam_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "gam_log_ior", "gam_log_ior_lower", "gam_log_ior_upper",
                                   "Log-IOR (GAM, IC 90%)", "GAM_LogIOR",
                                   plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "GAM_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$gam_reri_lower))) {
        reri_vals_gam <- c(plot_dt$gam_reri, plot_dt$gam_reri_lower, plot_dt$gam_reri_upper)
        reri_max_abs_gam <- max(abs(reri_vals_gam[is.finite(reri_vals_gam)]), na.rm = TRUE)
        y_limit_gam_reri <- max(ceiling(reri_max_abs_gam * 1.2), 2)
        p <- graph_metrics_counts(plot_dt, "gam_reri", "gam_reri_lower", "gam_reri_upper",
                             "RERI (GAM, IC 90%)", "GAM_RERI",
                             y_limit = y_limit_gam_reri,
                             plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "GAM_RERI")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "cls_log_ior", "cls_log_ior_lower", "cls_log_ior_upper",
                                   "Log-IOR (Estratificado, IC 90%)", "Classic_LogIOR",
                                   plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "Classic_LogIOR")
        generated_graphs <- generated_graphs + 1
      }

      if (any(is.finite(plot_dt$cls_reri_lower))) {
        reri_vals_cls <- c(plot_dt$cls_reri, plot_dt$cls_reri_lower, plot_dt$cls_reri_upper)
        reri_max_abs_cls <- max(abs(reri_vals_cls[is.finite(reri_vals_cls)]), na.rm = TRUE)
        y_limit_cls_reri <- max(ceiling(reri_max_abs_cls * 1.2), 2)
        p <- graph_metrics_counts(plot_dt, "cls_reri", "cls_reri_lower", "cls_reri_upper",
                             "RERI (Estratificado, IC 90%)", "Classic_RERI",
                             y_limit = y_limit_cls_reri,
                             plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "Classic_RERI")
        generated_graphs <- generated_graphs + 1
      }

    }, error = function(e) {
      message(sprintf("\n  Error en triplete %d: %s", i, e$message))
    })

    if (i %% 50 == 0) gc()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message(sprintf("\n  Gráficos generados: %d", generated_graphs))
  return(invisible(generated_graphs))
}

# Genera gráficos para tripletes del dataset ORIGINAL detectados como positivos en 40_network.
#
# fuente: network_triplets_for_graphs.rds
#  ya filtrado a positivos (triplet_id %in% positives$triplet_id en 40_network.R)
#  incluye drugA_name, drugB_name, meddra_name resueltos por 40_network.R
#  incluye counts_by_stage preembebidos (no se recalculan aquí)
# título : "<drugA_name> + <drugB_name> | Triplete: <id>"   /   "<meddra_name>"
# Métrica: solo GAM-LogIOR (RERI y métricas clásicas son NULL en este pipeline; se omiten)
# Salida: subcarpeta "network/" dentro de output_dir_positive

generate_positive_graphs_from_network <- function() {
  if (!file.exists(ruta_network_graphs)) {
    message("Archivo de red no encontrado: ", ruta_network_graphs)
    return(invisible(0L))
  }

  output_dir_net <- file.path(output_dir_positive, "network")
  dir.create(output_dir_net, showWarnings = FALSE, recursive = TRUE)

  dt <- readRDS(ruta_network_graphs)
  message(sprintf("  Tripletes del dataset original a graficar: %d", nrow(dt)))

  pb <- txtProgressBar(min = 0, max = nrow(dt), style = 3)
  generated_graphs <- 0

  for (i in seq_len(nrow(dt))) {
    tryCatch({
      row <- dt[i]

      # drugA_name, drugB_name y meddra_name están garantizados por 40_network.R
      ptitle <- paste0(row$drugA_name, " + ", row$drugB_name,
                          " | Triplete: ", row$triplet_id)
      psubtitle <- row$meddra_name

      # Los conteos vienen en counts_by_stage; ade_dt no es necesario
      plot_dt <- expand_triplets_counts(row, ade_dt = NULL, precomputed_counts = NULL)
      if (is.null(plot_dt) || nrow(plot_dt) == 0) {
        message(sprintf("\n  Salteando triplete %d: sin datos expandibles", i))
        next
      }

      if (any(is.finite(plot_dt$gam_log_ior_lower))) {
        p <- graph_metrics_counts(plot_dt, "gam_log_ior", "gam_log_ior_lower", "gam_log_ior_upper",
                                   "Log-IOR (GAM, IC 90%)", "GAM_LogIOR",
                                   plot_title = ptitle, plot_subtitle = psubtitle)
        save_positive_graph(p, plot_dt$triplet_id[1], plot_dt$dynamic[1], "GAM_LogIOR",
                             dir = output_dir_net)
        generated_graphs <- generated_graphs + 1
      }

    }, error = function(e) {
      message(sprintf("\n  Error en triplete %d: %s", i, e$message))
    })

    if (i %% 50 == 0) gc()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message(sprintf("\n  Gráficos generados: %d", generated_graphs))
  return(invisible(generated_graphs))
}

# Ejecución
generate_positive_graphs_from_network()
generate_positive_graphs_from_results()
