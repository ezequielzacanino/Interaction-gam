################################################################################
# Script de análisis descriptivo del dataset
# Script: 03_descriptive.R
################################################################################

################################################################################
# 1. CONFIGURACIÓN Y LIBRERÍAS
################################################################################

library(data.table)
library(tidyverse)
library(ggplot2)
library(svglite)
library(scales)
library(patchwork)

setwd("D:/Bioestadística/gam-farmacovigilancia")

source("00_functions.R", local = TRUE)
source("01_theme.R", local = TRUE)
theme_set(theme_base())
set.seed(9427)

################################################################################
# 2. PARÁMETROS CONFIGURABLES
################################################################################

min_reports_triplet  <- 2
min_nichd_with_rep   <- 2
top_n_drugs          <- 20   # Top N drogas más frecuentes en gráficos
top_n_events         <- 20   # Top N eventos más frecuentes en gráficos

dpi_figuras   <- 300
ancho_figura  <- 12
alto_figura   <- 8

suffix <- paste0("desc_", min_reports_triplet, "rep_", min_nichd_with_rep, "stage")

ruta_ade_raw   <- "./ade_raw.csv"
ruta_drug_info <- "./drug.csv"
output_dir     <- paste0("./results/", suffix, "/descriptive_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# HELPERS
################################################################################

# Etiquetas reutilizables para las 7 etapas NICHD
nichd_labels <- c(
  "term_neonatal"    = "Neonato",
  "infancy"          = "Lactante",
  "toddler"          = "Deambulador",
  "early_childhood"  = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence"= "Adol. Temprana",
  "late_adolescence" = "Adol. Tardía"
)

# Guarda figura en PNG + SVG con un solo llamado
save_fig <- function(p, name, w = ancho_figura, h = alto_figura) {
  base <- paste0(output_dir, name)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = dpi_figuras)
  ggsave(paste0(base, ".svg"), p, width = w, height = h, device = svglite)
  invisible(p)
}

################################################################################
# 3. CARGA Y PREPARACIÓN DE DATOS
################################################################################

drug_info_original <- fread(ruta_drug_info)
drug_info_original[, atc_concept_id := as.character(atc_concept_id)]
drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]

canonical_map <- drug_info_original[, .(canonical_id = min(atc_concept_id)), by = base_name]

translation_table <- merge(
  drug_info_original[, .(atc_concept_id, base_name)],
  canonical_map,
  by = "base_name"
)

cat(sprintf("IDs originales: %d\n", uniqueN(translation_table$atc_concept_id)))
cat(sprintf("IDs canónicos:  %d\n", uniqueN(translation_table$canonical_id)))

ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

ade_raw_dt <- merge(
  ade_raw_dt,
  translation_table[, .(atc_concept_id, canonical_id)],
  by = "atc_concept_id", all.x = TRUE
)
ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))
ade_raw_dt[, nichd     := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

################################################################################
# 4. MÉTRICAS BÁSICAS DEL DATASET
################################################################################

message("\n--- MÉTRICAS BÁSICAS ---")

n_filas_totales    <- nrow(ade_raw_dt)
n_reportes_totales <- uniqueN(ade_raw_dt$safetyreportid)
n_drogas_unicas    <- uniqueN(ade_raw_dt$atc_concept_id, na.rm = TRUE)
n_eventos_unicos   <- uniqueN(ade_raw_dt$meddra_concept_id, na.rm = TRUE)

message(sprintf("Filas totales:      %s", format(n_filas_totales,    big.mark = ",")))
message(sprintf("Reportes únicos:    %s", format(n_reportes_totales, big.mark = ",")))
message(sprintf("Drogas únicas:      %s", format(n_drogas_unicas,    big.mark = ",")))
message(sprintf("Eventos únicos:     %s", format(n_eventos_unicos,   big.mark = ",")))

# 4.1 Reportes por etapa NICHD
reportes_por_etapa <- ade_raw_dt[, .(
  n_reportes = uniqueN(safetyreportid),
  n_filas    = .N
), by = .(nichd_num, nichd)]
setorder(reportes_por_etapa, nichd_num)

message("\nReportes por etapa:")
print(reportes_por_etapa)

# 4.2 Complejidad por reporte (drogas y eventos)
complejidad_report <- merge(
  ade_raw_dt[!is.na(atc_concept_id),  .(n_drogas  = .N), by = safetyreportid],
  ade_raw_dt[!is.na(meddra_concept_id), .(n_eventos = .N), by = safetyreportid],
  by = "safetyreportid", all = TRUE
)
complejidad_report[is.na(n_drogas),  n_drogas  := 0L]
complejidad_report[is.na(n_eventos), n_eventos := 0L]

stats_drogas  <- complejidad_report[, .(media = mean(n_drogas),  mediana = median(n_drogas),
                                        q95 = quantile(n_drogas, 0.95),  max = max(n_drogas))]
stats_eventos <- complejidad_report[, .(media = mean(n_eventos), mediana = median(n_eventos),
                                        q95 = quantile(n_eventos, 0.95), max = max(n_eventos))]

message(sprintf("Drogas/reporte  – Media: %.2f | Mediana: %.1f | p95: %.1f | Max: %d",
                stats_drogas$media,  stats_drogas$mediana,  stats_drogas$q95,  stats_drogas$max))
message(sprintf("Eventos/reporte – Media: %.2f | Mediana: %.1f | p95: %.1f | Max: %d",
                stats_eventos$media, stats_eventos$mediana, stats_eventos$q95, stats_eventos$max))

################################################################################
# 5. CONSTRUCCIÓN Y ANÁLISIS DE TRIPLETES
################################################################################

message("\n--- ANÁLISIS DE TRIPLETES ---")

drugs_by_report  <- unique(ade_raw_dt[!is.na(atc_concept_id),   .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[!is.na(meddra_concept_id), .(safetyreportid, meddra_concept_id)])
reports          <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])

drugs_list  <- drugs_by_report[,  .(drugs  = list(unique(atc_concept_id))),  by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), by = safetyreportid]

reports_meta <- merge(reports, drugs_list,  by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)
reports_meta[is.na(drugs),  drugs  := list(integer(0))]
reports_meta[is.na(events), events := list(integer(0))]

message("Construyendo tripletes...")
triplets_list <- lapply(seq_len(nrow(reports_meta)), function(i) {
  rowi <- reports_meta[i]
  make_triplets(rowi$drugs[[1]], rowi$events[[1]], rowi$safetyreportid, rowi$nichd_num)
})
triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc(verbose = FALSE)

# 5.1 Totales
n_triplets_total <- nrow(triplets_dt)
message(sprintf("Tripletes totales (sin filtrar): %s", format(n_triplets_total, big.mark = ",")))

# 5.2 Reportes únicos por triplete
triplet_counts <- unique(triplets_dt[, .(drugA, drugB, meddra, safetyreportid)])[,
  .N, by = .(drugA, drugB, meddra)]

media_reportes_por_triplete   <- mean(triplet_counts$N)
mediana_reportes_por_triplete <- median(triplet_counts$N)

message(sprintf("Reportes/triplete – Media: %.2f | Mediana: %.1f",
                media_reportes_por_triplete, mediana_reportes_por_triplete))

# 5.3 Filtrado de candidatos
trip_counts_by_stage <- unique(triplets_dt[, .(drugA, drugB, meddra, nichd_num, safetyreportid)])[,
  .N, by = .(drugA, drugB, meddra, nichd_num)]

trip_summary <- trip_counts_by_stage[, .(
  N        = sum(N),
  n_stages = uniqueN(nichd_num)
), by = .(drugA, drugB, meddra)]

candidatos_filtrados <- trip_summary[N >= min_reports_triplet & n_stages >= min_nichd_with_rep]
n_candidatos_filtrados <- nrow(candidatos_filtrados)

message(sprintf("Candidatos (≥%d rep, ≥%d etapas): %s (%.2f%%)",
                min_reports_triplet, min_nichd_with_rep,
                format(n_candidatos_filtrados, big.mark = ","),
                100 * n_candidatos_filtrados / nrow(triplet_counts)))

# 5.4 Estadísticas de candidatos (sin duplicación)
candidatos_filtrados[, media_reportes_por_etapa := N / n_stages]
media_etapas_candidatos          <- mean(candidatos_filtrados$n_stages)
media_reportes_por_etapa_total   <- mean(candidatos_filtrados$media_reportes_por_etapa)

message(sprintf("Etapas/triplete (media):    %.2f / 7", media_etapas_candidatos))
message(sprintf("Reportes/etapa (media):     %.2f",    media_reportes_por_etapa_total))

################################################################################
# 6. DISTRIBUCIONES DETALLADAS
################################################################################

message("\n--- DISTRIBUCIONES ---")

# Distribución de tripletes por número de etapas
dist_etapas <- trip_summary[, .(
  n_triplets  = .N,
  pct_triplets = 100 * .N / nrow(trip_summary)
), by = n_stages][order(n_stages)]
print(dist_etapas)

# Distribución de reportes totales por candidato (rangos)
dist_reportes <- candidatos_filtrados[, .(
  n_triplets  = .N,
  pct_triplets = 100 * .N / n_candidatos_filtrados
), by = .(rango = cut(N,
                      breaks = c(0, 2, 5, 10, 20, 50, 100, Inf),
                      labels = c("2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100")))]
print(dist_reportes)

# Cobertura por etapa en candidatos
cand_key <- paste(candidatos_filtrados$drugA, candidatos_filtrados$drugB, candidatos_filtrados$meddra)
dist_por_etapa_candidatos <- trip_counts_by_stage[
  paste(drugA, drugB, meddra) %in% cand_key,
  .(n_triplets_con_etapa = uniqueN(paste(drugA, drugB, meddra)),
    total_reportes = sum(N),
    media_reportes = mean(N)),
  by = nichd_num
]
dist_por_etapa_candidatos[, nichd := niveles_nichd[nichd_num]]
setorder(dist_por_etapa_candidatos, nichd_num)
print(dist_por_etapa_candidatos[, .(nichd, n_triplets_con_etapa, total_reportes, media_reportes)])

# Top drogas y eventos
top_drugs  <- drugs_by_report[,  .N, by = atc_concept_id][order(-N)][1:top_n_drugs]
top_events <- events_by_report[, .N, by = meddra_concept_id][order(-N)][1:top_n_events]

# Coadministración: número de pares de drogas únicos
drug_pairs <- unique(triplets_dt[, .(drugA, drugB)])[, .N, by = .(drugA, drugB)]
n_drug_pairs <- nrow(drug_pairs)
message(sprintf("Pares de drogas únicos con coadministración: %s", format(n_drug_pairs, big.mark = ",")))

################################################################################
# 7. GRÁFICOS
################################################################################

message("\n--- GENERANDO GRÁFICOS ---")

# ── 7.1 Reportes por etapa NICHD ─────────────────────────────────────────────
p_reportes_etapa <- ggplot(reportes_por_etapa, aes(x = nichd, y = n_reportes)) +
  geom_col(fill = "#16A085", alpha = 0.85) +
  geom_text(aes(label = comma(n_reportes)), vjust = -0.5, size = 3.5) +
  scale_x_discrete(labels = nichd_labels) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Reportes por Etapa NICHD",
    subtitle = sprintf("Total: %s reportes únicos", comma(n_reportes_totales)),
    x = "Etapa NICHD", y = "Reportes únicos"
  )

save_fig(p_reportes_etapa, "fig_reportes_por_etapa")

# ── 7.2 Distribución de tripletes por número de etapas ───────────────────────
p_dist_etapas <- ggplot(dist_etapas, aes(x = factor(n_stages), y = n_triplets)) +
  geom_col(fill = "#E67E22", alpha = 0.85) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", comma(n_triplets), pct_triplets)),
            vjust = 1.4, size = 3, color = "white", fontface = "bold") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Tripletes por Número de Etapas con Datos",
    subtitle = sprintf("Total: %s | Media: %.2f etapas",
                       comma(n_triplets_total), mean(trip_summary$n_stages)),
    x = "Etapas NICHD con reportes", y = "Cantidad de tripletes"
  )

save_fig(p_dist_etapas, "fig_dist_etapas_tripletes")

# ── 7.3 Total vs candidatos filtrados ────────────────────────────────────────
datos_comparacion <- data.table(
  categoria = factor(c("Tripletes totales", "Candidatos filtrados"),
                     levels = c("Tripletes totales", "Candidatos filtrados")),
  cantidad  = c(nrow(triplet_counts), n_candidatos_filtrados),
  fill_col  = c("#C0392B", "#27AE60")
)

p_comparacion <- ggplot(datos_comparacion, aes(x = categoria, y = cantidad, fill = fill_col)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(cantidad)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Filtrado de Tripletes Candidatos",
    subtitle = sprintf("Criterio: ≥%d reportes AND ≥%d etapas | Retención: %.2f%%",
                       min_reports_triplet, min_nichd_with_rep,
                       100 * n_candidatos_filtrados / nrow(triplet_counts)),
    x = NULL, y = "Cantidad de tripletes"
  )

save_fig(p_comparacion, "fig_comparacion_filtro_tripletes", w = 10)

# ── 7.4 Distribución de reportes por triplete (candidatos, escala log) ───────
p_dist_reportes <- ggplot(candidatos_filtrados, aes(x = N)) +
  geom_histogram(bins = 40, fill = "#8E44AD", alpha = 0.75, color = "white") +
  geom_vline(xintercept = mediana_reportes_por_triplete,
             linetype = "dashed", color = "#C0392B", linewidth = 0.9) +
  annotate("text",
           x     = mediana_reportes_por_triplete * 1.8,
           y     = Inf, vjust = 2,
           label = sprintf("Mediana: %.1f", mediana_reportes_por_triplete),
           color = "#C0392B", fontface = "bold", size = 4) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Distribución de Reportes por Triplete (candidatos)",
    subtitle = sprintf("n = %s tripletes | Media: %.1f | Mediana: %.1f | Max: %d",
                       comma(n_candidatos_filtrados),
                       media_reportes_por_triplete,
                       mediana_reportes_por_triplete,
                       max(candidatos_filtrados$N)),
    x = "Reportes por triplete (escala log₁₀)", y = "Frecuencia"
  )

save_fig(p_dist_reportes, "fig_dist_reportes_por_triplete")

# ── 7.5 Cobertura de etapas en candidatos ────────────────────────────────────
p_cobertura_etapas <- ggplot(dist_por_etapa_candidatos,
                             aes(x = reorder(nichd, nichd_num), y = n_triplets_con_etapa)) +
  geom_col(fill = "#16A085", alpha = 0.85) +
  geom_text(aes(label = comma(n_triplets_con_etapa)), vjust = -0.5, size = 3.5) +
  scale_x_discrete(labels = nichd_labels) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Cobertura por Etapa en Candidatos Filtrados",
    subtitle = sprintf("n = %s candidatos | Media: %.2f etapas/triplete",
                       comma(n_candidatos_filtrados), media_etapas_candidatos),
    x = "Etapa", y = "n tripletes"
  )

save_fig(p_cobertura_etapas, "fig_cobertura_por_etapa")

# ── 7.6 Heatmap: reportes/etapa por triplete (muestra 500) ───────────────────
set.seed(9427)
muestra_idx <- candidatos_filtrados[sample(.N, min(500, .N)), .(drugA, drugB, meddra)]
muestra_idx[, triplet_id := .I]  # id ordinal para el eje Y

heatmap_long <- merge(
  muestra_idx,
  trip_counts_by_stage,
  by = c("drugA", "drugB", "meddra")
)[, .(triplet_id, nichd_num, N)]

# Completar con 0 en etapas sin reporte
full_grid <- CJ(triplet_id = muestra_idx$triplet_id, nichd_num = 1:7)
heatmap_long <- merge(full_grid, heatmap_long, by = c("triplet_id", "nichd_num"), all.x = TRUE)
heatmap_long[is.na(N), N := 0L]
heatmap_long[, nichd := factor(niveles_nichd[nichd_num], levels = niveles_nichd)]

p_heatmap <- ggplot(heatmap_long, aes(x = nichd, y = factor(triplet_id), fill = log1p(N))) +
  geom_tile(color = NA) +
  scale_x_discrete(labels = nichd_labels) +
  scale_fill_gradient(low = "#EBF5FB", high = "#154360",
                      name = "log(N+1)", labels = function(x) round(expm1(x), 1)) +
  labs(
    title    = "Distribución de Reportes por Etapa (muestra de tripletes)",
    subtitle = sprintf("Muestra aleatoria de %d candidatos filtrados", nrow(muestra_idx)),
    x = "Etapa NICHD", y = "Triplete (índice)"
  ) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

save_fig(p_heatmap, "fig_heatmap_reportes_triplete", h = 10)

# ── 7.7 Top drogas más frecuentes ────────────────────────────────────────────
# Intenta adjuntar nombre si está disponible en drug_info_original
drug_names <- unique(drug_info_original[, .(atc_concept_id, atc_concept_name)])
top_drugs_plot <- merge(top_drugs, drug_names, by = "atc_concept_id", all.x = TRUE)
top_drugs_plot[, label := ifelse(!is.na(atc_concept_name),
                                 substr(atc_concept_name, 1, 30),
                                 atc_concept_id)]

p_top_drugs <- ggplot(top_drugs_plot, aes(x = reorder(label, N), y = N)) +
  geom_col(fill = "#1A5276", alpha = 0.85) +
  geom_text(aes(label = comma(N)), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = sprintf("Top %d Drogas más Frecuentes", top_n_drugs),
    subtitle = sprintf("Por número de reportes únicos | Total drogas: %s", comma(n_drogas_unicas)),
    x = NULL, y = "Número de reportes"
  )

save_fig(p_top_drugs, "fig_top_drogas", h = 7)

# ── 7.8 Top eventos adversos más frecuentes ──────────────────────────────────
p_top_events <- ggplot(top_events, aes(x = reorder(as.character(meddra_concept_id), N), y = N)) +
  geom_col(fill = "#78281F", alpha = 0.85) +
  geom_text(aes(label = comma(N)), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = sprintf("Top %d Eventos Adversos más Frecuentes (MedDRA)", top_n_events),
    subtitle = sprintf("Por número de reportes únicos | Total eventos: %s", comma(n_eventos_unicos)),
    x = "MedDRA concept_id", y = "Número de reportes"
  )

save_fig(p_top_events, "fig_top_eventos", h = 7)

# ── 7.9 Distribución de drogas y eventos por reporte (panel) ─────────────────
p_drogas_hist <- ggplot(complejidad_report, aes(x = n_drogas)) +
  geom_histogram(bins = 30, fill = "#1A5276", alpha = 0.8, color = "white") +
  geom_vline(xintercept = stats_drogas$mediana,
             linetype = "dashed", color = "#C0392B", linewidth = 0.9) +
  scale_x_continuous(limits = c(0, quantile(complejidad_report$n_drogas, 0.99))) +
  scale_y_continuous(labels = comma) +
  labs(title = "Drogas por reporte",
       subtitle = sprintf("Mediana: %.0f | Media: %.1f", stats_drogas$mediana, stats_drogas$media),
       x = "Número de drogas", y = "Frecuencia")

p_eventos_hist <- ggplot(complejidad_report, aes(x = n_eventos)) +
  geom_histogram(bins = 30, fill = "#78281F", alpha = 0.8, color = "white") +
  geom_vline(xintercept = stats_eventos$mediana,
             linetype = "dashed", color = "#C0392B", linewidth = 0.9) +
  scale_x_continuous(limits = c(0, quantile(complejidad_report$n_eventos, 0.99))) +
  scale_y_continuous(labels = comma) +
  labs(title = "Eventos por reporte",
       subtitle = sprintf("Mediana: %.0f | Media: %.1f", stats_eventos$mediana, stats_eventos$media),
       x = "Número de eventos adversos", y = "Frecuencia")

p_complejidad <- p_drogas_hist + p_eventos_hist +
  plot_annotation(
    title    = "Complejidad de los Reportes",
    subtitle = sprintf("n = %s reportes únicos", comma(n_reportes_totales))
  )

save_fig(p_complejidad, "fig_complejidad_reportes", w = 14, h = 6)

################################################################################
# 8. EXPORTACIÓN DE RESULTADOS
################################################################################

message("\n--- EXPORTANDO RESULTADOS ---")

resumen_completo <- data.table(
  parametro = c(
    "n_filas_totales", "n_reportes_totales",
    "n_drogas_unicas", "n_eventos_unicos",
    "n_reportes_term_neonatal", "n_reportes_infancy",
    "n_reportes_toddler", "n_reportes_early_childhood",
    "n_reportes_middle_childhood", "n_reportes_early_adolescence",
    "n_reportes_late_adolescence",
    "media_drogas_por_reporte", "mediana_drogas_por_reporte",
    "media_eventos_por_reporte", "mediana_eventos_por_reporte",
    "n_triplets_total",
    "media_reportes_por_triplete", "mediana_reportes_por_triplete",
    "n_candidatos_filtrados", "pct_candidatos_vs_total",
    "media_etapas_candidatos", "media_reportes_por_etapa_candidatos",
    "n_pares_drogas_coadmin",
    "min_reports_filtro", "min_etapas_filtro"
  ),
  valor = c(
    n_filas_totales, n_reportes_totales,
    n_drogas_unicas, n_eventos_unicos,
    reportes_por_etapa[nichd == "term_neonatal",     n_reportes],
    reportes_por_etapa[nichd == "infancy",           n_reportes],
    reportes_por_etapa[nichd == "toddler",           n_reportes],
    reportes_por_etapa[nichd == "early_childhood",   n_reportes],
    reportes_por_etapa[nichd == "middle_childhood",  n_reportes],
    reportes_por_etapa[nichd == "early_adolescence", n_reportes],
    reportes_por_etapa[nichd == "late_adolescence",  n_reportes],
    stats_drogas$media,  stats_drogas$mediana,
    stats_eventos$media, stats_eventos$mediana,
    n_triplets_total,
    media_reportes_por_triplete, mediana_reportes_por_triplete,
    n_candidatos_filtrados,
    100 * n_candidatos_filtrados / nrow(triplet_counts),
    media_etapas_candidatos, media_reportes_por_etapa_total,
    n_drug_pairs,
    min_reports_triplet, min_nichd_with_rep
  )
)

fwrite(resumen_completo,              paste0(output_dir, "resumen_metricas.csv"))
fwrite(dist_etapas,                   paste0(output_dir, "distribucion_etapas.csv"))
fwrite(dist_reportes,                 paste0(output_dir, "distribucion_reportes_candidatos.csv"))
fwrite(dist_por_etapa_candidatos,     paste0(output_dir, "cobertura_por_etapa.csv"))
fwrite(reportes_por_etapa,            paste0(output_dir, "reportes_por_etapa.csv"))
fwrite(head(candidatos_filtrados[order(-N)], 100), paste0(output_dir, "top100_candidatos.csv"))

message("Archivos CSV exportados.")

################################################################################
# 9. REPORTE FINAL EN CONSOLA
################################################################################

cat(sprintf("
%s
RESUMEN EJECUTIVO
%s

DATASET BASE
  Filas totales:              %s
  Reportes únicos:            %s
  Drogas únicas:              %s
  Eventos únicos:             %s

COMPLEJIDAD POR REPORTE
  Drogas/reporte  – Media: %.2f | Mediana: %.0f | Max: %d
  Eventos/reporte – Media: %.2f | Mediana: %.0f | Max: %d

DISTRIBUCIÓN POR ETAPA NICHD\n",
strrep("=", 80), strrep("=", 80),
comma(n_filas_totales), comma(n_reportes_totales),
comma(n_drogas_unicas), comma(n_eventos_unicos),
stats_drogas$media,  stats_drogas$mediana,  stats_drogas$max,
stats_eventos$media, stats_eventos$mediana, stats_eventos$max
))

for (i in seq_len(nrow(reportes_por_etapa))) {
  cat(sprintf("  %-20s: %s (%.1f%%)\n",
              reportes_por_etapa$nichd[i],
              comma(reportes_por_etapa$n_reportes[i]),
              100 * reportes_por_etapa$n_reportes[i] / n_reportes_totales))
}

cat(sprintf("
ANÁLISIS DE TRIPLETES
  Tripletes totales:          %s
  Pares drogas (coadmin):     %s
  Reportes/triplete – Media: %.2f | Mediana: %.1f

FILTRADO DE CANDIDATOS (≥%d reportes, ≥%d etapas)
  Candidatos seleccionados:   %s (%.2f%% del total)
  Etapas/triplete (media):    %.2f / 7
  Reportes/etapa (media):     %.2f

SALIDA: %s
%s\n",
comma(n_triplets_total), comma(n_drug_pairs),
media_reportes_por_triplete, mediana_reportes_por_triplete,
min_reports_triplet, min_nichd_with_rep,
comma(n_candidatos_filtrados),
100 * n_candidatos_filtrados / nrow(triplet_counts),
media_etapas_candidatos, media_reportes_por_etapa_total,
output_dir,
strrep("=", 80)
))

message("Gráficos generados:")
for (f in c("fig_reportes_por_etapa", "fig_dist_etapas_tripletes",
            "fig_comparacion_filtro_tripletes", "fig_dist_reportes_por_triplete",
            "fig_cobertura_por_etapa", "fig_heatmap_reportes_triplete",
            "fig_top_drogas", "fig_top_eventos", "fig_complejidad_reportes")) {
  message(sprintf("  - %s.png / .svg", f))
}

message(strrep("=", 80))
message("ANÁLISIS DESCRIPTIVO COMPLETADO")
message(strrep("=", 80))

gc(verbose = FALSE)