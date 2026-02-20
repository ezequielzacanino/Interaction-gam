################################################################################
# Script de análisis de métricas y sensibilidad 
# Script: 30_metrics.R
################################################################################

library(data.table)
library(tidyverse)
library(pROC)
library(svglite)

setwd("D:/Bioestadística/gam-farmacovigilancia")
set.seed(9427)

source("00_functions.R", local = TRUE)
source("01_theme.R", local = TRUE)

niveles_nichd <- c(
  "term_neonatal", "infancy", "toddler", "early_childhood",
  "middle_childhood", "early_adolescence", "late_adolescence"
)

################################################################################
# Configuración
################################################################################

# percentilo elegido de distribución nula
percentil <- "p95"
n_boot <- 2000

# parámetros para calculo de subset de poder
target_power <- 0.80
grid_res <- 30
tij_max <- 0.15  
n_max <- 250

# parámetros para uso de distribución nula
use_threshold_ior <- TRUE
use_threshold_reri <- TRUE

# parámetros para correr según formula guardada
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
k_spline <- 7
include_nichd <- FALSE
nichd_spline <- FALSE
bs_type <- "cs"
select <- FALSE

# parámetros para guardado según formula usada
suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

# clasificación de etapas con alto reporte según dinámica
stage_class <- rbind(
  data.table(nichd = niveles_nichd, dynamic = "uniform", class = 1),
  data.table(nichd = niveles_nichd, dynamic = "increase", class = c(0, 0, NA, NA, NA, 1, 1)),
  data.table(nichd = niveles_nichd, dynamic = "decrease", class = c(1, 1, NA, NA, NA, 0, 0)),
  data.table(nichd = niveles_nichd, dynamic = "plateau", class = c(0, NA, 1, 1, 1, NA, 0)),
  data.table(nichd = niveles_nichd, dynamic = "inverse_plateau", class = c(1, NA, 0, 0, 0, NA, 1))
)

# rutas
ruta_base_sensitivity <- paste0("./results/", suffix, "/augmentation_results/")
output_dir <- paste0("./results/", suffix, "/metrics_results/")

ruta_coadmin_pos <- paste0("./results/", suffix, "/augmentation_results/positive_coadmin_by_stage.csv")
ruta_coadmin_neg <- paste0("./results/", suffix, "/augmentation_results/negative_coadmin_by_stage.csv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# niveles de reducción
reduction_levels <- c(0, seq(10, 90, by = 10))

################################################################################
# Etiquetas para gráficos
################################################################################

# etiquetas para traducción 
dynamic_labels <- c(
  "uniform" = "Uniforme",
  "increase" = "Incremento",
  "decrease" = "Disminución",
  "plateau" = "Meseta",
  "inverse_plateau" = "Valle"
)

nichd_labels <- c(
  "term_neonatal" = "Neonato a término",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temprana",
  "late_adolescence" = "Adolescencia tardía"
)

metrics_to_plot <- c("sensitivity", "specificity", "PPV", "NPV", "F1")
metric_labels <- c(
  "sensitivity" = "Sensibilidad",
  "specificity" = "Especificidad",
  "PPV" = "Valor Predictivo Positivo",
  "NPV" = "Valor Predictivo Negativo",
  "F1" = "F1-Score"
)

method_pairs <- list(
  list(gam = "GAM-logIOR", classic = "Estratificado-IOR", label = "IOR"),
  list(gam = "GAM-RERI", classic = "Estratificado-RERI", label = "RERI"),
  list(gam = "GAM-Doble", classic = "Estratificado-Doble", label = "Doble")
)

################################################################################
# Carga de umbrales nulos
################################################################################
ruta_null_dist <- paste0("./results/", suffix, "/null_distribution_results/null_distribution.csv")

null_thresholds_ior <- fread(paste0("./results/", suffix, "/null_distribution_results/null_thresholds.csv"))
thresh_col <- paste0("threshold_", percentil)
null_thresholds_ior <- null_thresholds_ior[, .(stage, threshold_ior = get(thresh_col))]

null_thresholds_reri <- fread(paste0("./results/", suffix, "/null_distribution_results/null_thresholds_reri.csv"))
null_thresholds_reri <- null_thresholds_reri[, .(stage, threshold_reri = get(thresh_col))]

null_thresholds <- merge(null_thresholds_ior, null_thresholds_reri, by = "stage")

# datos de coadministración
coadmin_stage_pos <- fread(ruta_coadmin_pos)
coadmin_stage_neg <- fread(ruta_coadmin_neg)
setnames(coadmin_stage_pos, "nichd_num", "stage_num")
setnames(coadmin_stage_neg, "nichd_num", "stage_num")

################################################################################
# Carga de datos sin reducción 
################################################################################

data_baseline <- expandir_datos(0)  # carga datos sin reducción para cálculos de umbrales de poder
pos_high_base <- data_baseline$pos_high
neg_high_base <- data_baseline$neg_high

################################################################################
# Gráficos de comparación de distribuciones
################################################################################

null_distribution <- fread(ruta_null_dist)

# sampleo distribución nula para graficar
null_sample <- null_distribution[sample(.N, min(.N, 50000))]

comparison_data <- rbind(
  null_sample[, .(stage_num = stage, log_ior_lower90 = log_lower90, reri_lower90 = reri_lower90, source = "Distribución nula")],
  pos_high_base[, .(stage_num, log_ior_lower90 = gam_log_ior_lower90, reri_lower90 = gam_reri_lower90, source = "Positivos")],
  neg_high_base[, .(stage_num, log_ior_lower90 = gam_log_ior_lower90, reri_lower90 = gam_reri_lower90, source = "Negativos")]
)
comparison_data[, stage_name := factor(stage_num, levels = 1:7, labels = nichd_labels)] # ordeno para las facetas

# corto data para mostrar en gráfico (esto lo debería hacer limitando ejes en el gráfico igual)
comparison_data <- comparison_data[ is.finite(log_ior_lower90)]
comparison_data <- comparison_data[ is.finite(reri_lower90)]

# Paleta de colores consistente
color_palette <- c(
  "Distribución nula" = "gray60",
  "Positivos" = "#4DAF4A",
  "Negativos" = "#E41A1C"
)

# Gráfico combinado con facets
p_null_vs_obs_ior <- ggplot(
  comparison_data, 
  aes(x = log_ior_lower90, fill = source)
) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~ stage_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(limits = c(-5, 5)) + # limito para visualizar mejor las distribuciones
  labs(
    title = sprintf("Distribución nula vs señales detectadas (%s)", percentil),
    x = "Log(IOR) - Límite inferior IC 90%",
    y = "Densidad",
    fill = "Fuente"
  )

# Guardar gráfico combinado
ggsave(
  paste0(output_dir, "fig_null_vs_observed_ior.png"), 
  p_null_vs_obs_ior,
  width = 16, 
  height = 12, 
  dpi = 300
)

print(p_null_vs_obs_ior)

# Gráfico combinado con facets
p_null_vs_obs_reri <- ggplot(
  comparison_data, 
  aes(x = reri_lower90, fill = source)
) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  facet_wrap(~ stage_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(limits = c(-0.25, 0.25)) + # limito para visualizar mejor las distribuciones
  labs(
    title = sprintf("Distribución nula vs señales detectadas (%s)", percentil),
    x = "RERI - Límite inferior IC 90%",
    y = "Densidad",
    fill = "Fuente"
  ) 
# Guardar gráfico combinado
ggsave(
  paste0(output_dir, "fig_null_vs_observed_reri.png"), 
  p_null_vs_obs_reri,
  width = 16, 
  height = 12, 
  dpi = 300
)

print(p_null_vs_obs_reri)

message(sprintf("Dataset base: %d positivos, %d negativos", nrow(pos_high_base), nrow(neg_high_base)))

################################################################################
# Cálculo de subsets de poder para cada método
################################################################################

# GAM-logIOR
message("GAM-logIOR")
power_gam_ior <- calculate_power_gam(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,
  use_threshold_ior = TRUE,
  use_threshold_reri = FALSE,
  detection = "ior"
)

# GAM-RERI
message("GAM-RERI")
power_gam_reri <- calculate_power_gam(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,
  use_threshold_ior = FALSE,
  use_threshold_reri = TRUE,
  detection = "reri"
)

# Estratificado-IOR
message("Estratificado-IOR")
power_cls_ior <- calculate_power_classic(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "ior",
  na_remove = TRUE
)

# Estratificado-RERI
message("Estratificado-RERI")
power_cls_reri <- calculate_power_classic(
  data_pos = pos_high_base,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "reri",
  na_remove = TRUE
)

power_ior_list <- list(
  "GAM" = power_gam_ior,
  "Estratificado" = power_cls_ior
)

p_surface_ior <- plot_power_surface(
  power_results_list = power_ior_list,
  target_power = target_power, 
  detection = "IOR",
  grid_size = grid_res,
  t_range = c(0, tij_max), 
  n_range = c(0, n_max)
)

ggsave(
  paste0(output_dir, "fig_power_surface_ior_combined.png"), 
  p_surface_ior,
  width = 14,      
  height = 7,     
  dpi = 300
)

print(p_surface_ior)

power_reri_list <- list(
  "GAM" = power_gam_reri,
  "Estratificado" = power_cls_reri
)

p_surface_reri <- plot_power_surface(
  power_results_list = power_reri_list,
  target_power = target_power,
  detection = "RERI", 
  grid_size = grid_res,
  t_range = c(0, tij_max),
  n_range = c(0, n_max)
)

ggsave(
  paste0(output_dir, "fig_power_surface_reri_combined.png"),
  p_surface_reri, 
  width = 14,
  height = 7,
  dpi = 300
)

print(p_surface_reri)
  
# Guardo IDs de tripletes para cada subset de poder
power_ids <- list(
  "GAM-logIOR" = unique(power_gam_ior$superset_pos[class == 1]$triplet_id),
  "GAM-RERI" = unique(power_gam_reri$superset_pos[class == 1]$triplet_id),
  "Estratificado-IOR" = unique(power_cls_ior$superset_pos[class == 1]$triplet_id),
  "Estratificado-RERI" = unique(power_cls_reri$superset_pos[class == 1]$triplet_id)
)

# subsets de intersección
triplets_intersection_reri <- intersect(power_ids[["GAM-RERI"]], power_ids[["Estratificado-RERI"]])
message(sprintf("\nIntersección GAM-RERI & Estratificado-RERI: %d tripletes", length(triplets_intersection_reri)))

triplets_intersection_ior <- intersect(power_ids[["GAM-logIOR"]], power_ids[["Estratificado-IOR"]])
message(sprintf("\nIntersección GAM-logIOR & Estratificado-IOR: %d tripletes", length(triplets_intersection_ior)))

# Resumen de subsets
summary_power <- data.table(
  method = names(power_ids),
  n_triplets = sapply(power_ids, length),
  t_star = c(power_gam_ior$t_star, power_gam_reri$t_star,
             power_cls_ior$t_star, power_cls_reri$t_star),
  n_star = c(power_gam_ior$n_star, power_gam_reri$n_star,
             power_cls_ior$n_star, power_cls_reri$n_star)
)
print(summary_power)

################################################################################
# Carga de datos con reducción 
################################################################################

# expando datos según niveles de reducción
datos_por_reduccion <- lapply(reduction_levels, expandir_datos)
names(datos_por_reduccion) <- as.character(reduction_levels)

################################################################################
# Definición de métodos 
################################################################################

# métodos con tipo de score mapeado 
metodos <- list(
  list(nombre = "GAM-logIOR", tipo = "IOR", es_gam = TRUE, 
       score_type = "gam_log_ior_lower90", score_type_auc = "gam_log_ior"),
  list(nombre = "GAM-RERI", tipo = "RERI", es_gam = TRUE, 
       score_type = "gam_reri_lower90", score_type_auc = "gam_reri"),
  list(nombre = "Estratificado-IOR", tipo = "IOR", es_gam = FALSE, 
       score_type = "classic_log_ior_lower90", score_type_auc = "classic_log_ior"),
  list(nombre = "Estratificado-RERI", tipo = "RERI", es_gam = FALSE, 
       score_type = "classic_reri_lower90", score_type_auc = "classic_reri")
)

################################################################################
# Cálculo de métricas
################################################################################

# Métricas originales sin reducción
res_global_original <- list()
res_dinamica_original <- list()
res_etapa_original <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo de niveles de reducción
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in metodos) {
    
    message(sprintf("  %s", met$nombre))
    
    # Métricas globales
    # Filtrar negativos solo a etapas de alto reporte
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    dt_global <- rbind(datos$pos_high, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    # Agregar a nivel triplete para métricas globales
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE,  # acá TRUE para evitar doble conteo
    score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`(
      method = met$nombre,
      reduction_pct = red_pct,
      dataset = "original"
    )]
    res_global_original[[length(res_global_original) + 1]] <- metrics_global
    
    # Métricas por dinámica
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # remuevo uniform (creo que redundante)
     
    for (dyn in dinamicas) {
      # Etapas altas para esta dinámica
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd] 
      pos_dyn <- datos$pos_high[nichd %in% etapas_altas & dynamic == dyn]
  
      # Negativos en las MISMAS etapas altas
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]
  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE)
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc) 
      metrics_dyn[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        dynamic = dyn,
        dataset = "original"
      )]
      res_dinamica_original[[length(res_dinamica_original) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa (NO AGREGAR POR TRIPLETE)
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      
      # Para positivos: solo etapas clasificadas como altas (redundante)
      pos_stage <- datos$pos_high[stage_num == s & class == 1]
      
      # Negativos de la misma etapa
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      # NO agregar - cada fila es una etapa específica
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        stage_num = s,
        nichd = nichd_label,
        dataset = "original"
      )]
      res_etapa_original[[length(res_etapa_original) + 1]] <- metrics_stage
    }
  }
}

metrics_global_original <- rbindlist(res_global_original, fill = TRUE)
metrics_dynamic_original <- rbindlist(res_dinamica_original, fill = TRUE)
metrics_stage_original <- rbindlist(res_etapa_original, fill = TRUE)

fwrite(metrics_global_original, paste0(output_dir, "metrics_global_original.csv"))
fwrite(metrics_dynamic_original, paste0(output_dir, "metrics_dynamic_original.csv"))
fwrite(metrics_stage_original, paste0(output_dir, "metrics_stage_original.csv"))

################################################################################
# Cálculo de métricas en subset de poder estadístico
################################################################################

res_global_filtered <- list()
res_dinamica_filtered <- list()
res_etapa_filtered <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in metodos) {
    
    message(sprintf("  %s", met$nombre))
    
    ids_filtrar <- power_ids[[met$nombre]]
    
    # Métricas globales filtradas
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    pos_fil <- datos$pos_high[triplet_id %in% ids_filtrar]
    dt_global <- rbind(pos_fil, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE, 
      score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`(
      method = met$nombre,
      reduction_pct = red_pct,
      dataset = "filtered"
    )]
    res_global_filtered[[length(res_global_filtered) + 1]] <- metrics_global
    
    # Métricas por dinámica filtradas
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # (redundante)
    
    for (dyn in dinamicas) {
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd]
      pos_dyn <- datos$pos_high[triplet_id %in% ids_filtrar & nichd %in% etapas_altas & dynamic == dyn]
      
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE) 
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_dyn[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        dynamic = dyn,
        dataset = "filtered"
      )]
      res_dinamica_filtered[[length(res_dinamica_filtered) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa filtradas
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      pos_stage <- datos$pos_high[triplet_id %in% ids_filtrar & stage_num == s & class == 1]
      
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`(
        method = met$nombre,
        reduction_pct = red_pct,
        stage_num = s,
        nichd = nichd_label,
        dataset = "filtered"
      )]
      res_etapa_filtered[[length(res_etapa_filtered) + 1]] <- metrics_stage
    }
  }
}

metrics_global_filtered <- rbindlist(res_global_filtered, fill = TRUE)
metrics_dynamic_filtered <- rbindlist(res_dinamica_filtered, fill = TRUE)
metrics_stage_filtered <- rbindlist(res_etapa_filtered, fill = TRUE)

fwrite(metrics_global_filtered, paste0(output_dir, "metrics_global_filtered.csv"))
fwrite(metrics_dynamic_filtered, paste0(output_dir, "metrics_dynamic_filtered.csv"))
fwrite(metrics_stage_filtered, paste0(output_dir, "metrics_stage_filtered.csv"))

################################################################################
# Cálculo de métricas en subset de intersección
################################################################################

pares_interseccion <- list(
  list(
    metodos = metodos[sapply(metodos, function(m) m$nombre %in% c("GAM-RERI", "Estratificado-RERI"))],
    triplets = triplets_intersection_reri,
    label = "RERI"
  ),
  list(
    metodos = metodos[sapply(metodos, function(m) m$nombre %in% c("GAM-logIOR", "Estratificado-IOR"))],
    triplets = triplets_intersection_ior,
    label = "IOR"
  )
)

res_global_inter <- list()
res_dinamica_inter <- list()
res_etapa_inter <- list()

for (red_pct in reduction_levels) {
  
  message(sprintf("\nReducción %d%%", red_pct)) # trackeo
  datos <- datos_por_reduccion[[as.character(red_pct)]]
  
  for (met in par$metodos) {
    message(sprintf("%s [intersección %s]", met$nombre, par$label)))
    
    pos_int <- datos$pos_high[triplet_id %in% par$triplets]
    
    # Métricas a nivel global intersección
    etapas_alto_reporte <- stage_class[class == 1, unique(nichd)]
    neg_global <- datos$neg_high[nichd %in% etapas_alto_reporte]

    pos_int <- datos$pos_high[triplet_id %in% triplets_intersection]
    dt_global <- rbind(pos_int, neg_global, fill = TRUE)
    dt_global <- aplicar_deteccion(dt_global, met$nombre, met$tipo, use_null = met$es_gam)
    
    metrics_global <- calcular_metricas_simple(dt_global, n_boot, agregar_por_triplete = TRUE, score_type = met$score_type, score_type_auc = met$score_type_auc)
    metrics_global[, `:=`( method = met$nombre, reduction_pct = red_pct, dataset = "intersection" )]
    res_global_inter[[length(res_global_inter) + 1]] <- metrics_global
    
    # Métricas por dinámica en intersección 
    dinamicas <- setdiff(unique(datos$pos_high$dynamic), "uniform") # redundante
    
    for (dyn in dinamicas) {
      etapas_altas <- stage_class[dynamic == dyn & class == 1, nichd]
      pos_dyn <- pos_int[nichd %in% etapas_altas & dynamic == dyn]
      neg_dyn <- datos$neg_high[nichd %in% etapas_altas]  
  
      dt_dyn <- rbind(pos_dyn, neg_dyn, fill = TRUE)
      dt_dyn <- aplicar_deteccion(dt_dyn, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_dyn <- calcular_metricas_simple(dt_dyn, n_boot, agregar_por_triplete = TRUE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_dyn[, `:=`( method = met$nombre, reduction_pct = red_pct, dynamic = dyn, dataset = "intersection")]
      res_dinamica_inter[[length(res_dinamica_inter) + 1]] <- metrics_dyn
    }
    
    # Métricas por etapa en intersección 
    for (s in 1:7) {
      nichd_label <- niveles_nichd[s]
      pos_stage <- pos_int[stage_num == s & class == 1] # redundante class == 1
      
      neg_stage <- datos$neg_high[stage_num == s]
      
      dt_stage <- rbind(pos_stage, neg_stage, fill = TRUE)
      dt_stage <- aplicar_deteccion(dt_stage, met$nombre, met$tipo, use_null = met$es_gam)
      
      metrics_stage <- calcular_metricas_simple(dt_stage, n_boot, agregar_por_triplete = FALSE, 
        score_type = met$score_type, score_type_auc = met$score_type_auc)
      metrics_stage[, `:=`( method = met$nombre, reduction_pct = red_pct, stage_num = s, nichd = nichd_label, dataset = "intersection")]
      res_etapa_inter[[length(res_etapa_inter) + 1]] <- metrics_stage
    }
  }
}

metrics_global_intersection <- rbindlist(res_global_inter, fill = TRUE)
metrics_dynamic_intersection <- rbindlist(res_dinamica_inter, fill = TRUE)
metrics_stage_intersection <- rbindlist(res_etapa_inter, fill = TRUE)

fwrite(metrics_global_intersection, paste0(output_dir, "metrics_global_intersection.csv"))
fwrite(metrics_dynamic_intersection, paste0(output_dir, "metrics_dynamic_intersection.csv"))
fwrite(metrics_stage_intersection, paste0(output_dir, "metrics_stage_intersection.csv"))


