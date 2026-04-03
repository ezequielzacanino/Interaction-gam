################################################################################
# Análisis Descriptivo del Dataset de Farmacovigilancia
# Script: descriptive_analysis.R
################################################################################

library(data.table)
library(tidyverse)
library(scales)

set.seed(9427)
setwd("D:/Bioestadística/gam-farmacovigilancia")

################################################################################
# Configuración
################################################################################

ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_gene <- "./drug_gene.csv"
ruta_drug_info <- "./drug.csv" 
min_reports_threshold <- 10  # Umbral mínimo para análisis de cobertura

output_dir <- "./results/descriptive_analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source("00_functions.R", local = TRUE)
source("01_theme.R")
theme_set(theme_base())

################################################################################
# Carga y preprocesamiento
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CARGA DE DATOS")
message(paste(rep("=", 80), collapse = ""))

ade_raw_dt <- fread(ruta_ade_raw)

################################################################################
# Preprocesamiento
################################################################################

# diccionario de drogas
drug_info_original <- fread(ruta_drug_info)

# chequeo de variable
drug_info_original[, atc_concept_id := as.character(atc_concept_id)]

# limpieza de nombre
drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]

# mapeo de nombres
# si coincide nombre de base (ej propranolol; systemic = propranolol; topic)
# se unifica atc_concept_id
canonical_map <- drug_info_original[, .(
  canonical_id = min(atc_concept_id)
), by = base_name]

# merge de id
translation_table <- merge(
  drug_info_original[, .(atc_concept_id, base_name)],
  canonical_map,
  by = "base_name"
)

# resultados de unificación
cat(sprintf("IDs originales: %d\n", uniqueN(translation_table$atc_concept_id)))
cat(sprintf("IDs únicos: %d\n", uniqueN(translation_table$canonical_id)))

###########
# Unificación de ids en dataset 
###########

ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

# merge con tabla de ids canónicos
ade_raw_dt <- merge(
  ade_raw_dt, 
  translation_table[, .(atc_concept_id, canonical_id)], 
  by = "atc_concept_id", 
  all.x = TRUE
)

# Si atc_concept_id tiene otro id canónico, reemplaza
ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL] # limpio columna auxiliar

# Unifico ids crea filas duplicadas
# antes --> 1 evento con propranolol podía tener 2 filas
# ej: id 001 atc_concept_id 123 (propranolol; systemic), meddra id 456
#     id 001 atc_concept_id 124 (propranolol; topical), meddra id 456
# con ids canónicos, ambas filas se vuelven idénticas e introducen colinealidad, borro
nrow_before <- nrow(ade_raw_dt)

# unique de las columnas clave
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

# Proceso nichd
ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]


# Distribución por etapa
stage_dist <- ade_raw_dt[, .(n_reports = uniqueN(safetyreportid)), by = nichd]
setorder(stage_dist, nichd)
message("\nDistribución de reportes por etapa:")
print(stage_dist)

message(sprintf("  Dataset: %s filas", format(nrow(ade_raw_dt), big.mark = ",")))
message(sprintf("  Reportes únicos: %s", format(uniqueN(ade_raw_dt$safetyreportid), big.mark = ",")))
message(sprintf("  Drogas únicas: %s", format(uniqueN(ade_raw_dt$atc_concept_id), big.mark = ",")))
message(sprintf("  Eventos únicos: %s", format(uniqueN(ade_raw_dt$meddra_concept_id), big.mark = ",")))

################################################################################
# Construcción de tripletes
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CONSTRUCCIÓN DE TRIPLETES")
message(paste(rep("=", 80), collapse = ""))

# Metadata por reporte
drugs_by_report <- unique(ade_raw_dt[!is.na(atc_concept_id), 
                                     .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[!is.na(meddra_concept_id), 
                                      .(safetyreportid, meddra_concept_id)])

reports_meta <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), 
                               by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), 
                                 by = safetyreportid]

reports_meta <- merge(reports_meta, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
reports_meta[is.na(drugs), drugs := list(integer(0))]
reports_meta[is.na(events), events := list(integer(0))]

# Generación de tripletes
triplets_list <- lapply(seq_len(nrow(reports_meta)), function(i) {
  rowi <- reports_meta[i]
  make_triplets(
    drug = rowi$drugs[[1]], 
    event = rowi$events[[1]], 
    report_id = rowi$safetyreportid, 
    nichd_stage = rowi$nichd_num
  )
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

# Conteos de tripletes
triplet_counts <- unique(triplets_dt[, .(drugA, drugB, meddra, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra)
]

message(sprintf("  Tripletes únicos: %s", format(uniqueN(triplet_counts), big.mark = ",")))
message(sprintf("  Media reportes por triplete: %.1f", mean(triplet_counts$N)))
message(sprintf("  Mediana reportes por triplete: %.0f", median(triplet_counts$N)))

################################################################################
# Construcción de dobletes (droga simple - evento)
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CONSTRUCCIÓN DE DOBLETES")
message(paste(rep("=", 80), collapse = ""))

# Dobletes: cada droga individual con cada evento
doublet_dt <- unique(ade_raw_dt[, .(safetyreportid, atc_concept_id, meddra_concept_id)])

doublet_counts <- doublet_dt[, .N, by = .(atc_concept_id, meddra_concept_id)]

message(sprintf("  Dobletes únicos: %s", format(nrow(doublet_counts), big.mark = ",")))
message(sprintf("  Media reportes por doblete: %.1f", mean(doublet_counts$N)))
message(sprintf("  Mediana reportes por doblete: %.0f", median(doublet_counts$N)))

################################################################################
# Análisis de coadministraciones
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("ANÁLISIS DE COADMINISTRACIONES")
message(paste(rep("=", 80), collapse = ""))

# Pares de drogas por reporte
drug_pairs_dt <- drugs_by_report[, {
  drug_list <- unique(atc_concept_id)
  if (length(drug_list) >= 2) {
    pairs <- t(combn(drug_list, 2))
    data.table(
      drugA = pmin(pairs[,1], pairs[,2]),
      drugB = pmax(pairs[,1], pairs[,2])
    )
  } else {
    data.table()
  }
}, by = safetyreportid]

# Conteos de coadministraciones
coadmin_counts <- drug_pairs_dt[, .N, by = .(drugA, drugB)]

message(sprintf("  Coadministraciones únicas: %s", format(nrow(coadmin_counts), big.mark = ",")))
message(sprintf("  Media reportes por coadmin: %.1f", mean(coadmin_counts$N)))
message(sprintf("  Mediana reportes por coadmin: %.0f", median(coadmin_counts$N)))

# Coadministraciones por etapa
coadmin_by_stage <- merge(
  drug_pairs_dt,
  reports_meta[, .(safetyreportid, nichd)],
  by = "safetyreportid"
)

coadmin_stage_summary <- coadmin_by_stage[, 
  .(n_unique_coadmin = uniqueN(paste(drugA, drugB, sep = "_")),
    n_reports = .N), 
  by = nichd
]
setorder(coadmin_stage_summary, nichd)

message("\nCoadministraciones por etapa:")
print(coadmin_stage_summary)

################################################################################
# Análisis de sparsity
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("ANÁLISIS DE SPARSITY")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Sparsity de tripletes
###########

n_drugs <- uniqueN(ade_raw_dt$atc_concept_id)
n_events <- uniqueN(ade_raw_dt$meddra_concept_id)

# Combinaciones posibles de tripletes (pares de drogas × eventos)
n_possible_drug_pairs <- choose(n_drugs, 2)
n_possible_triplets <- n_possible_drug_pairs * n_events

# Tripletes observados
n_observed_triplets <- nrow(triplet_counts)

sparsity_triplets <- 1 - (n_observed_triplets / n_possible_triplets)

message(sprintf("\nSparsity de tripletes:"))
message(sprintf("  Combinaciones posibles: %s", format(n_possible_triplets, big.mark = ",")))
message(sprintf("  Tripletes observados: %s", format(n_observed_triplets, big.mark = ",")))
message(sprintf("  Sparsity: %.6f (%.2f%% de celdas vacías)", 
                sparsity_triplets, 100 * sparsity_triplets))

###########
# 2- Sparsity de dobletes
###########

n_possible_doublets <- n_drugs * n_events
n_observed_doublets <- nrow(doublet_counts)
sparsity_doublets <- 1 - (n_observed_doublets / n_possible_doublets)

message(sprintf("\nSparsity de dobletes:"))
message(sprintf("  Combinaciones posibles: %s", format(n_possible_doublets, big.mark = ",")))
message(sprintf("  Dobletes observados: %s", format(n_observed_doublets, big.mark = ",")))
message(sprintf("  Sparsity: %.6f (%.2f%% de celdas vacías)", 
                sparsity_doublets, 100 * sparsity_doublets))

###########
# 3- Cobertura de coadministraciones
###########

n_observed_coadmin <- nrow(coadmin_counts)
coverage_coadmin <- n_observed_coadmin / n_possible_drug_pairs

message(sprintf("\nCobertura de coadministraciones:"))
message(sprintf("  Pares posibles: %s", format(n_possible_drug_pairs, big.mark = ",")))
message(sprintf("  Pares observados: %s", format(n_observed_coadmin, big.mark = ",")))
message(sprintf("  Cobertura: %.6f (%.2f%%)", coverage_coadmin, 100 * coverage_coadmin))

###########
# 4- Distribución de reportes por droga
###########

reports_per_drug <- drugs_by_report[, .(n_reports = .N), by = atc_concept_id]

drug_summary <- data.table(
  metric = c("Media", "Mediana", "Min", "Max", "Q25", "Q75"),
  value = c(
    mean(reports_per_drug$n_reports),
    median(reports_per_drug$n_reports),
    min(reports_per_drug$n_reports),
    max(reports_per_drug$n_reports),
    quantile(reports_per_drug$n_reports, 0.25),
    quantile(reports_per_drug$n_reports, 0.75)
  )
)

message("\nReportes por droga:")
print(drug_summary)

###########
# 5- Distribución de reportes por evento
###########

reports_per_event <- events_by_report[, .(n_reports = .N), by = meddra_concept_id]

event_summary <- data.table(
  metric = c("Media", "Mediana", "Min", "Max", "Q25", "Q75"),
  value = c(
    mean(reports_per_event$n_reports),
    median(reports_per_event$n_reports),
    min(reports_per_event$n_reports),
    max(reports_per_event$n_reports),
    quantile(reports_per_event$n_reports, 0.25),
    quantile(reports_per_event$n_reports, 0.75)
  )
)

message("\nReportes por evento:")
print(event_summary)

###########
# 6- Drogas y eventos con suficiente cobertura
###########

n_drugs_sufficient <- sum(reports_per_drug$n_reports >= min_reports_threshold)
n_events_sufficient <- sum(reports_per_event$n_reports >= min_reports_threshold)

message(sprintf("\nEntidades con >= %d reportes:", min_reports_threshold))
message(sprintf("  Drogas: %d / %d (%.1f%%)",
                n_drugs_sufficient, n_drugs,
                100 * n_drugs_sufficient / n_drugs))
message(sprintf("  Eventos: %d / %d (%.1f%%)",
                n_events_sufficient, n_events,
                100 * n_events_sufficient / n_events))

################################################################################
# Comparación: Reportes de droga única vs Coadministración
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("COMPARACIÓN: DROGA ÚNICA VS COADMINISTRACIÓN")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Estadísticas de drogas individuales
###########

# Ya calculado anteriormente: reports_per_drug
drug_stats <- data.table(
  entity = "Droga única",
  n_unique = nrow(reports_per_drug),
  mean_reports = mean(reports_per_drug$n_reports),
  median_reports = median(reports_per_drug$n_reports),
  sd_reports = sd(reports_per_drug$n_reports),
  min_reports = min(reports_per_drug$n_reports),
  max_reports = max(reports_per_drug$n_reports),
  q25_reports = quantile(reports_per_drug$n_reports, 0.25),
  q75_reports = quantile(reports_per_drug$n_reports, 0.75),
  q95_reports = quantile(reports_per_drug$n_reports, 0.95),
  q99_reports = quantile(reports_per_drug$n_reports, 0.99)
)

message("\nEstadísticas de reportes por droga individual:")
message(sprintf("  N drogas únicas: %s", format(drug_stats$n_unique, big.mark = ",")))
message(sprintf("  Media reportes/droga: %.1f", drug_stats$mean_reports))
message(sprintf("  Mediana reportes/droga: %.0f", drug_stats$median_reports))
message(sprintf("  DE reportes/droga: %.1f", drug_stats$sd_reports))
message(sprintf("  Rango: [%d, %s]", 
                drug_stats$min_reports, 
                format(drug_stats$max_reports, big.mark = ",")))
message(sprintf("  Q25-Q75: [%.0f, %.0f]", 
                drug_stats$q25_reports, 
                drug_stats$q75_reports))

###########
# 2- Estadísticas de coadministraciones
###########

# Ya calculado anteriormente: coadmin_counts
coadmin_stats <- data.table(
  entity = "Coadministración (A+B)",
  n_unique = nrow(coadmin_counts),
  mean_reports = mean(coadmin_counts$N),
  median_reports = median(coadmin_counts$N),
  sd_reports = sd(coadmin_counts$N),
  min_reports = min(coadmin_counts$N),
  max_reports = max(coadmin_counts$N),
  q25_reports = quantile(coadmin_counts$N, 0.25),
  q75_reports = quantile(coadmin_counts$N, 0.75),
  q95_reports = quantile(coadmin_counts$N, 0.95),
  q99_reports = quantile(coadmin_counts$N, 0.99)
)

message("\nEstadísticas de reportes por coadministración:")
message(sprintf("  N coadministraciones únicas: %s", 
                format(coadmin_stats$n_unique, big.mark = ",")))
message(sprintf("  Media reportes/coadmin: %.1f", coadmin_stats$mean_reports))
message(sprintf("  Mediana reportes/coadmin: %.0f", coadmin_stats$median_reports))
message(sprintf("  DE reportes/coadmin: %.1f", coadmin_stats$sd_reports))
message(sprintf("  Rango: [%d, %s]", 
                coadmin_stats$min_reports, 
                format(coadmin_stats$max_reports, big.mark = ",")))
message(sprintf("  Q25-Q75: [%.0f, %.0f]", 
                coadmin_stats$q25_reports, 
                coadmin_stats$q75_reports))

###########
# 3- Comparación estadística
###########

# Ratio de medias
ratio_means <- drug_stats$mean_reports / coadmin_stats$mean_reports
ratio_medians <- drug_stats$median_reports / coadmin_stats$median_reports

message("\nComparación droga única vs coadministración:")
message(sprintf("  Ratio medias (única/coadmin): %.2f", ratio_means))
message(sprintf("  Ratio medianas (única/coadmin): %.2f", ratio_medians))
message(sprintf("  Diferencia de medias: %.1f reportes", 
                drug_stats$mean_reports - coadmin_stats$mean_reports))

# Test de Wilcoxon (distribuciones no paramétricas)
wilcox_test <- wilcox.test(
  reports_per_drug$n_reports,
  coadmin_counts$N,
  alternative = "two.sided"
)

message(sprintf("  Test de Wilcoxon p-value: %s", 
                format.pval(wilcox_test$p.value, digits = 3)))

###########
# 4- Tabla comparativa completa
###########

comparison_table <- rbind(drug_stats, coadmin_stats)

fwrite(comparison_table, paste0(output_dir, "comparison_drug_vs_coadmin.csv"))

print(comparison_table)

###########
# 5- Distribución de frecuencias
###########

# Categorías de frecuencia para análisis
freq_breaks <- c(1, 2, 5, 10, 20, 50, 100, 500, 1000, Inf)
freq_labels <- c("1", "2-4", "5-9", "10-19", "20-49", "50-99", 
                 "100-499", "500-999", "1000+")

reports_per_drug[, freq_category := cut(
  n_reports, 
  breaks = freq_breaks, 
  labels = freq_labels,
  right = FALSE
)]

coadmin_counts[, freq_category := cut(
  N, 
  breaks = freq_breaks, 
  labels = freq_labels,
  right = FALSE
)]

drug_freq_dist <- reports_per_drug[, .(
  entity = "Droga única",
  n_entities = .N,
  pct_entities = 100 * .N / nrow(reports_per_drug)
), by = freq_category]

coadmin_freq_dist <- coadmin_counts[, .(
  entity = "Coadministración",
  n_entities = .N,
  pct_entities = 100 * .N / nrow(coadmin_counts)
), by = freq_category]

freq_comparison <- rbind(drug_freq_dist, coadmin_freq_dist)
setorder(freq_comparison, freq_category, entity)

fwrite(freq_comparison, paste0(output_dir, "frequency_distribution_comparison.csv"))

message("\nDistribución de frecuencias:")
print(freq_comparison[, .(entity, freq_category, n_entities, 
                          pct = sprintf("%.1f%%", pct_entities))])

################################################################################
# Visualizaciones comparativas
################################################################################

###########
# 1- Distribución lado a lado
###########

# Preparar datos para visualización
drug_plot <- reports_per_drug[, .(
  entity = "Droga única",
  n_reports = n_reports
)]

coadmin_plot <- coadmin_counts[, .(
  entity = "Coadministración",
  n_reports = N
)]

combined_data <- rbind(drug_plot, coadmin_plot)


# 1. Preparar datos de Droga Simple por Etapa
single_counts_by_stage <- ade_raw_dt[!is.na(atc_concept_id), 
  .(n_reports = uniqueN(safetyreportid)), 
  by = .(nichd, atc_concept_id)
][, .(nichd, n_reports, entity = "Droga única")]

# 2. Preparar datos de Coadministración por Etapa
# Usamos el objeto coadmin_by_stage creado anteriormente
coadmin_counts_by_stage <- coadmin_by_stage[, 
  .(n_reports = .N), 
  by = .(nichd, drugA, drugB)
][, .(nichd, n_reports, entity = "Coadministración")]

# 3. Combinar datasets
combined_stage_data <- rbind(single_counts_by_stage, coadmin_counts_by_stage)
combined_stage_data[, entity := factor(entity, levels = c("Droga única", "Coadministración"))]

# 4. Muestreo estratificado para visualización (si hay demasiados datos)
# Tomamos hasta 5000 puntos por cada combinación de etapa y entidad para no saturar el renderizado
set.seed(9427)
plot_data_stage <- combined_stage_data[, {
  if (.N > 50000) .SD[sample(.N, 50000)] else .SD
}, by = .(nichd, entity)]

# 5. Generar Plot modificado
p_comparison <- ggplot(plot_data_stage, aes(x = nichd, y = n_reports, fill = entity)) +
  geom_violin(
    position = position_dodge(width = 0.8), 
    alpha = 0.6, 
    trim = FALSE, 
    scale = "width"
  ) +
  scale_y_continuous(limits = c(0,1000)) +
  scale_fill_manual(
    values = c("Droga única" = "#C0392B", "Coadministración" = "#16A085")
  ) +
  annotation_logticks(sides = "l") +
  labs(
    title = "Distribución de reportes por Etapa: Droga única vs Coadministración (A + B)",
    x = "Etapa del desarrollo",
    y = "Número de reportes por entidad",
    fill = "Tipo"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(paste0(output_dir, "fig_drug_vs_coadmin_distribution_by_stage.png"),
       p_comparison, width = 12, height = 8, dpi = 300)

print(p_comparison)

###########
# 2- Densidad superpuesta
###########

p_density <- ggplot(plot_data_stage, aes(x = n_reports, fill = entity, color = entity)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_x_log10(
    labels = comma,
    breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
  ) +
  scale_fill_manual(
    values = c("Droga única" = "#C0392B", "Coadministración" = "#16A085")
  ) +
  scale_color_manual(
    values = c("Droga única" = "#C0392B", "Coadministración" = "#C0392B")
  ) +
  annotation_logticks(sides = "b") +
  labs(
    title = "Densidad de distribución de reportes",
    subtitle = sprintf(
      "Mediana droga única: %.0f | Mediana coadmin: %.0f ",
      drug_stats$median_reports,
      coadmin_stats$median_reports
    ),
    x = "Número de reportes (log)",
    y = "Densidad",
    fill = "Entidad",
    color = "Entidad"
  ) +
  theme(legend.position = "bottom")

ggsave(paste0(output_dir, "fig_drug_vs_coadmin_density.png"),
       p_density, width = 10, height = 7, dpi = 300)

print(p_density)

###########
# 3- Distribución de frecuencias (barras)
###########

p_frequency <- ggplot(freq_comparison, 
                      aes(x = freq_category, y = pct_entities, fill = entity)) +
  geom_col(position = "dodge", alpha = 0.8, color = "black") +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_entities)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(
    values = c("Droga única" = "#C0392B", "Coadministración" = "#16A085")
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Distribución de frecuencias: Droga única vs Coadministración",
    subtitle = "Porcentaje de entidades por categoría de número de reportes",
    x = "Número de reportes",
    y = "Porcentaje de entidades",
    fill = "Tipo de entidad"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(paste0(output_dir, "fig_frequency_distribution_comparison.png"),
       p_frequency, width = 12, height = 8, dpi = 300)

print(p_frequency)

###########
# 4- Curvas acumuladas (ECDF)
###########

p_ecdf <- ggplot(plot_data_stage, aes(x = n_reports, color = entity)) +
  stat_ecdf(linewidth = 1.2) +
  scale_x_log10(
    labels = comma,
    breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
  ) +
  scale_color_manual(
    values = c("Droga única" = "#C0392B", "Coadministración" = "#16A085")
  ) +
  annotation_logticks(sides = "b") +
  labs(
    title = "Función de distribución acumulada empírica (ECDF)",
    subtitle = "Proporción de entidades con ≤ N reportes",
    x = "Número de reportes (escala log)",
    y = "Proporción acumulada",
    color = "Entidad"
  ) +
  theme(legend.position = "bottom")

ggsave(paste0(output_dir, "fig_drug_vs_coadmin_ecdf.png"),
       p_ecdf, width = 10, height = 7, dpi = 300)

print(p_ecdf)

################################################################################
# Guardado de tablas
################################################################################

fwrite(triplet_counts, paste0(output_dir, "triplet_counts.csv"))
fwrite(doublet_counts, paste0(output_dir, "doublet_counts.csv"))
fwrite(coadmin_counts, paste0(output_dir, "coadmin_counts.csv"))
fwrite(coadmin_stage_summary, paste0(output_dir, "coadmin_by_stage.csv"))

sparsity_summary <- data.table(
  entity = c("Tripletes", "Dobletes", "Coadmin"),
  n_possible = c(n_possible_triplets, n_possible_doublets, n_possible_drug_pairs),
  n_observed = c(n_observed_triplets, n_observed_doublets, n_observed_coadmin),
  sparsity = c(sparsity_triplets, sparsity_doublets, 1 - coverage_coadmin),
  coverage_pct = c(
    100 * (1 - sparsity_triplets),
    100 * (1 - sparsity_doublets),
    100 * coverage_coadmin
  )
)
fwrite(sparsity_summary, paste0(output_dir, "sparsity_summary.csv"))

################################################################################
# Visualizaciones
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("GENERACIÓN DE GRÁFICOS")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Tripletes vs Dobletes
###########

# Preparación de datos
triplet_plot_data <- triplet_counts[, .(
  type = "Triplete",
  id = paste(drugA, drugB, meddra, sep = "_"),
  n_reports = N
)]

doublet_plot_data <- doublet_counts[, .(
  type = "Doblete",
  id = paste(atc_concept_id, meddra_concept_id, sep = "_"),
  n_reports = N
)]

combined_plot_data <- rbind(triplet_plot_data, doublet_plot_data)

max_points <- 1000000      # cantidad máxima de puntos a graficar por grupo
v_repel    <- 0.01      # jitter vertical (repel)

# Muestreo para visualización (máximo 5000 puntos por tipo)
plot_sample <- combined_plot_data[, {
  if (.N > max_points) .SD[sample(.N, max_points)] else .SD
}, by = type]

p1 <- ggplot(plot_sample, aes(x = type, y = n_reports, color = type)) +
  geom_jitter(
    alpha = 0.3,
    width = 0.5,
    height = v_repel,
    size = 0.5
  ) +
  geom_violin(alpha = 0.2, color = NA, aes(fill = type)) +
  scale_y_continuous(
  limits = c(1, 500),
  breaks = seq(0, 500, by = 100),
  labels = comma) +
  scale_color_manual(values = c("Triplete" = "#2C3E50", "Doblete" = "#E74C3C")) +
  scale_fill_manual(values = c("Triplete" = "#2C3E50", "Doblete" = "#E74C3C")) +
  labs(
    title = "Distribución de reportes: Dobletes (A + Evento) vs Tripletes (A + B + Evento) ",
    x = NULL,
    y = "Número de reportes"
  ) +
  theme(legend.position = "none")

ggsave(
  paste0(output_dir, "fig_triplets_vs_doublets.png"),
  p1, width = 10, height = 7, dpi = 300
)
print(p1)

###########
# 2- Reportes por etapa: Droga Simple vs Coadministración
###########

# 1. Preparar datos de Droga Simple (Ocurrencias únicas de Droga-Reporte por etapa)
single_drug_by_stage <- unique(ade_raw_dt[!is.na(atc_concept_id), 
                                          .(safetyreportid, atc_concept_id, nichd)])

single_stage_summary <- single_drug_by_stage[, .(n_reports = .N), by = nichd]
single_stage_summary[, type := "Droga Simple"]

# 2. Preparar datos de Coadministración (Usando el resumen calculado anteriormente)
# Nota: n_reports en coadmin_stage_summary representa ocurrencias de Pares-Reporte
coadmin_plot_data <- coadmin_stage_summary[, .(nichd, n_reports, type = "Coadministración")]

# 3. Combinar datasets
stage_comparison_data <- rbind(single_stage_summary, coadmin_plot_data)
stage_comparison_data[, type := factor(type, levels = c("Droga Simple", "Coadministración"))]

# 4. Generar Plot P2 modificado
p2 <- ggplot(stage_comparison_data, aes(x = nichd, y = n_reports, fill = type)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9, color = "black") +
  geom_text(
    aes(label = comma(n_reports)), 
    position = position_dodge(width = 0.9), 
    vjust = -0.5, 
    size = 3,
    check_overlap = TRUE
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(
    values = c("Droga Simple" = "#3498DB", "Coadministración" = "#16A085")
  ) +
  labs(
    title = "Número de reportes por etapa: Droga Simple vs Coadministración (A + B)",
    x = "Etapa del desarrollo",
    y = "Número de ocurrencias",
    fill = "Tipo de evidencia"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(paste0(output_dir, "fig_coadmin_vs_single_by_stage.png"),
       p2, width = 12, height = 7, dpi = 300)

print(p2)

###########
# 3- Distribución de tripletes por número de reportes
###########

triplet_dist <- triplet_counts[, .(n_triplets = .N), by = N]
setorder(triplet_dist, N)

p3 <- ggplot(triplet_dist, aes(x = N, y = n_triplets)) +
  geom_point(color = "#2C3E50", size = 2) +
  scale_x_log10(labels = comma, breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000)) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Distribución de tripletes por número de reportes",
    subtitle = sprintf("%s tripletes únicos (A + B + Evento )", 
                      format(nrow(triplet_counts), big.mark = ",")),
    x = "Número de reportes por triplete",
    y = "Número de tripletes"
  ) +
  annotation_logticks()

ggsave(paste0(output_dir, "fig_triplet_distribution.png"),
       p3, width = 10, height = 7, dpi = 300)

print(p3)

###########
# 4- Distribución de coadministraciones por número de reportes
###########

coadmin_dist <- coadmin_counts[, .(n_coadmin = .N), by = N]
setorder(coadmin_dist, N)

p4 <- ggplot(coadmin_dist, aes(x = N, y = n_coadmin)) +
  geom_point(color = "#E74C3C", size = 2) +
  scale_x_log10(labels = comma, breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000)) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Distribución de co-administraciones (A + B) por número de reportes",
    subtitle = sprintf("%s coadministraciones únicas", 
                      format(nrow(coadmin_counts), big.mark = ",")),
    x = "Número de reportes por coadministración",
    y = "Número de coadministraciones"
  ) +
  annotation_logticks()

ggsave(paste0(output_dir, "fig_coadmin_distribution.png"),
       p4, width = 10, height = 7, dpi = 300)

print(p4)


###########
# 5- Top drogas y eventos
###########

top_n <- 20

top_drugs <- reports_per_drug[order(-n_reports)][1:top_n]
top_drugs[, rank := .I]

p5 <- ggplot(top_drugs, aes(x = reorder(atc_concept_id, n_reports), y = n_reports)) +
  geom_col(fill = "#2C3E50", alpha = 0.8, color = "black") +
  geom_text(aes(label = comma(n_reports)), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = sprintf("Top %d drogas por número de reportes", top_n),
    x = "ATC Concept ID",
    y = "Número de reportes"
  )

ggsave(paste0(output_dir, "fig_top_drugs.png"),
       p5, width = 10, height = 10, dpi = 300)

print(p5)

top_events <- reports_per_event[order(-n_reports)][1:top_n]
top_events[, rank := .I]

p6 <- ggplot(top_events, aes(x = reorder(meddra_concept_id, n_reports), y = n_reports)) +
  geom_col(fill = "#E74C3C", alpha = 0.8, color = "black") +
  geom_text(aes(label = comma(n_reports)), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = sprintf("Top %d eventos por número de reportes", top_n),
    x = "MedDRA Concept ID",
    y = "Número de reportes"
  )

ggsave(paste0(output_dir, "fig_top_events.png"),
       p6, width = 10, height = 10, dpi = 300)

print(p6)
message("  ✓ Gráficos de top entidades guardados")

################################################################################
# Resumen ejecutivo
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("RESUMEN EJECUTIVO")
message(paste(rep("=", 80), collapse = ""))

executive_summary <- sprintf("
================================================================================
ANÁLISIS DESCRIPTIVO DEL DATASET DE FARMACOVIGILANCIA
================================================================================

DIMENSIONES DEL DATASET:
  Filas totales: %s
  Reportes únicos: %s
  Drogas únicas: %s
  Eventos únicos: %s

TRIPLETES (DrugA - DrugB - Evento):
  Tripletes únicos: %s
  Combinaciones posibles: %s
  Sparsity: %.6f (%.2f%%%% vacío)
  Media reportes/triplete: %.1f
  Mediana reportes/triplete: %.0f

DOBLETES (Drug - Evento):
  Dobletes únicos: %s
  Combinaciones posibles: %s
  Sparsity: %.6f (%.2f%%%% vacío)
  Media reportes/doblete: %.1f
  Mediana reportes/doblete: %.0f

COADMINISTRACIONES (DrugA - DrugB):
  Pares únicos: %s
  Pares posibles: %s
  Cobertura: %.2f%%%%
  Media reportes/par: %.1f
  Mediana reportes/par: %.0f

COBERTURA POR ETAPA NICHD:
%s

ENTIDADES CON >= %d REPORTES:
  Drogas: %d / %d (%.1f%%%%)
  Eventos: %d / %d (%.1f%%%%)

CONCLUSIONES:
  • El dataset presenta alta sparsity (%.2f%%%% para tripletes)
  • Solo %.2f%%%% de coadministraciones posibles están observadas
  • Distribuciones fuertemente sesgadas (mediana << media)
  • Necesidad de filtros mínimos para análisis estadísticos robustos

Fecha: %s
================================================================================
",
format(nrow(ade_raw_dt), big.mark = ","),
format(uniqueN(ade_raw_dt$safetyreportid), big.mark = ","),
format(n_drugs, big.mark = ","),
format(n_events, big.mark = ","),
format(n_observed_triplets, big.mark = ","),
format(n_possible_triplets, big.mark = ","),
sparsity_triplets,
100 * sparsity_triplets,
mean(triplet_counts$N),
median(triplet_counts$N),
format(n_observed_doublets, big.mark = ","),
format(n_possible_doublets, big.mark = ","),
sparsity_doublets,
100 * sparsity_doublets,
mean(doublet_counts$N),
median(doublet_counts$N),
format(n_observed_coadmin, big.mark = ","),
format(n_possible_drug_pairs, big.mark = ","),
100 * coverage_coadmin,
mean(coadmin_counts$N),
median(coadmin_counts$N),
paste(capture.output(print(coadmin_stage_summary[, .(nichd, n_reports, n_unique_coadmin)])), collapse = "\n  "),
min_reports_threshold,
n_drugs_sufficient, n_drugs,
100 * n_drugs_sufficient / n_drugs,
n_events_sufficient, n_events,
100 * n_events_sufficient / n_events,
100 * sparsity_triplets,
100 * coverage_coadmin,
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

cat(executive_summary)
writeLines(executive_summary, paste0(output_dir, "resumen_ejecutivo.txt"))

message("\n", paste(rep("=", 80), collapse = ""))
message("ANÁLISIS COMPLETADO")
message(paste(rep("=", 80), collapse = ""))

