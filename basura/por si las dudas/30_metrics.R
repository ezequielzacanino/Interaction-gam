################################################################################
# Script de análisis de métricas y sensibilidad 
# Script: 30_metrics
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
# Configuraicón
################################################################################

percentil <- "p95"
n_boot <- 2000

# parámetros para calculo de subset de poder
target_power <- 0.80
tij_max <- 0.25  
n_max <- 300 
grid_res <- 30

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
suffix_sensitivity <- paste0(
  "sens_",
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

# parámetros de guardado según formula
suffix_base <- paste0(
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
ruta_base_sensitivity <- paste0("./results/", suffix_sensitivity, "/augmentation_results/")
output_dir <- paste0("./results/", suffix_sensitivity, "/metrics_results/")

ruta_coadmin_pos <- paste0("./results/", suffix_sensitivity, "/augmentation_results/positive_coadmin_by_stage.csv")
ruta_coadmin_neg <- paste0("./results/", suffix_sensitivity, "/augmentation_results/negative_coadmin_by_stage.csv")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# niveles de reducción
reduction_levels <- c(0, seq(10, 90, by = 10))

################################################################################
# CONFIGURACIÓN DE MÉTODOS
################################################################################

# clasificación de métodos
method_configs <- list(
  list(name = "GAM-logIOR", type = "IOR", is_gam = TRUE),
  list(name = "GAM-RERI", type = "RERI", is_gam = TRUE),
  list(name = "GAM-Doble", type = "Doble", is_gam = TRUE),
  list(name = "Estratificado-IOR", type = "IOR", is_gam = FALSE),
  list(name = "Estratificado-RERI", type = "RERI", is_gam = FALSE),
  list(name = "Estratificado-Doble", type = "Doble", is_gam = FALSE)
)

# helper para parámetros de threshold
get_threshold_params <- function(cfg, use_null_global) {
  list(
    use_t_ior = cfg$is_gam && use_null_global && (cfg$type %in% c("IOR", "Doble")),
    use_t_reri = cfg$is_gam && use_null_global && (cfg$type %in% c("RERI", "Doble"))
  )
}

################################################################################
# CARGA DE DATOS
################################################################################

null_thresholds_ior <- fread(paste0("./results/", suffix_base, "/null_distribution_results/null_thresholds.csv"))
thresh_col <- paste0("threshold_", percentil)
null_thresholds_ior <- null_thresholds_ior[, .(stage, threshold_ior = get(thresh_col))]

null_thresholds_reri <- fread(paste0("./results/", suffix_base, "/null_distribution_results/null_thresholds_reri.csv"))
null_thresholds_reri <- null_thresholds_reri[, .(stage, threshold_reri = get(thresh_col))]

null_thresholds <- merge(null_thresholds_ior, null_thresholds_reri, by = "stage")

# datos de coadministración
coadmin_stage_pos <- fread(ruta_coadmin_pos)
coadmin_stage_neg <- fread(ruta_coadmin_neg)

# renombro nichd_num a stage_num para consistencia con el resto del script
setnames(coadmin_stage_pos, "nichd_num", "stage_num")
setnames(coadmin_stage_neg, "nichd_num", "stage_num")

################################################################################
# CARGA DE DATOS BASE
################################################################################

# Carga dataset sin reducción
# filtra por injection_sucess = TRUE y 
# filtra uniform
data_baseline <- load_and_expand_data(red_pct = 0)

# pos_high y neg_high ya tienen todas las columnas necesarias:
# t_ij, n_coadmin (para GAM)
# n_coadmin_stage (ya viene de augmentation_results via expand_clean_all_metrics)
# gam_log_ior_lower90, gam_reri_lower90, classic_log_ior_lower90, classic_reri_lower90
# threshold_ior, threshold_reri (ya merged en expand_clean_all_metrics)
# triplet_id, stage_num, nichd, dynamic, class, label

pos_high_baseline <- data_baseline$pos_high
neg_high_baseline <- data_baseline$neg_high

pos_high_baseline <- merge(pos_high_baseline, 
  coadmin_stage_pos[, .(triplet_id, stage_num, n_coadmin_stage)], 
  by = c("triplet_id", "stage_num"), 
  all.x = TRUE
)

# innecesario pero por las dudas
neg_high_baseline <- merge(neg_high_baseline, 
  coadmin_stage_neg[, .(triplet_id, stage_num, n_coadmin_stage)],   
  by = c("triplet_id", "stage_num"), 
  all.x = TRUE
)

message(sprintf("Dataset baseline cargado: %d positivos, %d negativos", 
                nrow(pos_high_baseline), nrow(neg_high_baseline)))

###########
# subset de poder para GAM-logIOR
###########

message("\nsubset GAM-logIOR")
result_gam_power_logior <- calculate_power_gam(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,  
  use_threshold_ior = use_threshold_ior,    
  use_threshold_reri = FALSE,  # solo IOR
  detection = "ior"   
)

###########
# subset de poder para GAM-RERI
###########

message("\nsubset GAM-RERI")
result_gam_power_reri <- calculate_power_gam(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,  
  use_threshold_ior = FALSE,  # solo RERI
  use_threshold_reri = use_threshold_reri,
  detection = "reri"   
)

###########
# subset de poder para GAM-Doble
###########

message("\nsubset GAM-Doble")
result_gam_power <- calculate_power_gam(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = null_thresholds,
  metric_n = "n_coadmin",
  grid_resolution = grid_res,  
  use_threshold_ior = use_threshold_ior,    
  use_threshold_reri = use_threshold_reri,
  detection = "double"   
)

###########
# subset de poder para Estratificado-IOR
###########

message("\nsubset Estratificado-IOR")
result_classic_power_ior <- calculate_power_classic(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",  # clásico usa n_coadmin_stage
  grid_resolution = grid_res,
  detection = "ior",
  na_remove = TRUE
)

###########
# subset de poder para Estratificado-RERI
###########

message("\nEstratificado-RERI")
result_classic_power_reri <- calculate_power_classic(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "reri",
  na_remove = TRUE
)

###########
# subset de poder para Estratificado-Doble
###########

message("\nsubset Estratificado-Doble")
result_classic_power <- calculate_power_classic(
  data_pos = pos_high_baseline,
  target_power = target_power,
  null_thresholds = NULL,
  metric_n = "n_coadmin_stage",
  grid_resolution = grid_res,
  detection = "double",
  na_remove = TRUE
)

###########
# construcción de los datasets filtrados por poder
###########

# lista de datasets filtrados 
filtered_datasets <- list(
  "GAM-logIOR" = list(
    pos = result_gam_power_logior$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_gam_power_logior$t_star,
    n_star = result_gam_power_logior$n_star
  ),
  "GAM-RERI" = list(
    pos = result_gam_power_reri$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_gam_power_reri$t_star,
    n_star = result_gam_power_reri$n_star
  ),
  "GAM-Doble" = list(
    pos = result_gam_power$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_gam_power$t_star,
    n_star = result_gam_power$n_star
  ),
  "Estratificado-IOR" = list(
    pos = result_classic_power_ior$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_classic_power_ior$t_star,
    n_star = result_classic_power_ior$n_star
  ),
  "Estratificado-RERI" = list(
    pos = result_classic_power_reri$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_classic_power_reri$t_star,
    n_star = result_classic_power_reri$n_star
  ),
  "Estratificado-Doble" = list(
    pos = result_classic_power$superset_pos[class == 1],
    neg = neg_high_baseline,
    t_star = result_classic_power$t_star,
    n_star = result_classic_power$n_star
  )
)

###########
# triplet_ids para compatibilidad con análisis de intersección
###########

# Creo ids de tripletes que cumplen con la configuración mínima de poder
power_ids <- list()
for (method_name in names(filtered_datasets)) {
  ds <- filtered_datasets[[method_name]]
  pos_ids <- unique(ds$pos$triplet_id)
  neg_ids <- unique(ds$neg$triplet_id)
  
  power_ids[[method_name]] <- data.table(
    triplet_id = c(pos_ids, neg_ids),
    label = c(rep(1, length(pos_ids)), rep(0, length(neg_ids)))
  )
  message(sprintf("  %s: %d tripletes", method_name, nrow(power_ids[[method_name]])))
}

# calculo de intersección de tripletes 
triplets_intersection <- integer(0)
if (all(c("GAM-Doble", "Estratificado-Doble") %in% names(power_ids))) {
  triplets_gam_all <- power_ids[["GAM-Doble"]]$triplet_id
  triplets_cls_all <- power_ids[["Estratificado-Doble"]]$triplet_id
  triplets_intersection <- intersect(triplets_gam_all, triplets_cls_all)
  message(sprintf("  Intersección GAM-Doble ∩ Estratificado-Doble: %d tripletes", 
                  length(triplets_intersection)))
}

###########
# Resumen de subsets
###########

summary_power <- rbindlist(lapply(names(filtered_datasets), function(m) {
  ds <- filtered_datasets[[m]]
  data.table(
    method = m,
    t_star = ds$t_star,
    n_star = ds$n_star,
    n_pos_retained = nrow(ds$pos),
    n_neg_retained = nrow(ds$neg),
    n_triplets_pos = uniqueN(ds$pos$triplet_id)
  )
}))

print(summary_power)

################################################################################
# Carga de todos los datos con reducción 
################################################################################

# load_and_expand_data carga los datasets por nivel de reducción
# el primer argumento pasa los niveles de reducción 
# filtra por injection_sucess = TRUE y 
# filtra uniform
expanded_data_by_reduction <- lapply(reduction_levels, load_and_expand_data)
names(expanded_data_by_reduction) <- as.character(reduction_levels)

################################################################################
# Métricas originales (globales, por dinámica y por etapa)
################################################################################

# helper para aplicar a todos los métodos y niveles
# paso parametro de distribución nula en argumento use_null
# aggregation_level se determina según la función evaluada
run_evaluation <- function(eval_fn, use_null = use_threshold_ior, aggregation_level = "triplet") {
  rbindlist(lapply(expanded_data_by_reduction, function(data_exp) {
    rbindlist(lapply(method_configs, function(cfg) {
      eval_fn(data_exp, cfg, use_null, aggregation_level) 
    }))
  }))
}

# uso el helper para los wrappers que contienen evaluate_detection()
# evaluate_detection() contiene add_detection_status para clasificación de señal en detectada
# y compute_bootstrap_metrics para IC95 de métricas

# globales: agregar a nivel triplete (sin inflación por etapas)
metrics_global_original <- run_evaluation(evaluate_global_wrapper, 
                                            use_null = use_threshold_ior, 
                                            aggregation_level = "triplet")

# por dinámica: cada etapa es observación independiente  
metrics_dynamic_original <- run_evaluation(evaluate_by_dynamic,
                                           use_null = use_threshold_ior,
                                           aggregation_level = "stage")

# por etapa: cada etapa es observación independiente
metrics_stage_original <- run_evaluation(evaluate_by_stage_wrapper,
                                         use_null = use_threshold_ior,
                                         aggregation_level = "stage")

fwrite(metrics_global_original, paste0(output_dir, "metrics_global_original.csv"))
fwrite(metrics_dynamic_original, paste0(output_dir, "metrics_dynamic_original.csv"))
fwrite(metrics_stage_original, paste0(output_dir, "metrics_stage_original.csv"))

message(sprintf("  Global: %d filas", nrow(metrics_global_original)))
message(sprintf("  Dinámica: %d filas", nrow(metrics_dynamic_original)))
message(sprintf("  Etapa: %d filas", nrow(metrics_stage_original)))

################################################################################
# Métricas originales en dataset filtrado (global, por dinámica y por etapa)
################################################################################

# helper para crear datasets filtrados para cada método y nivel de reducción
create_filtered_data <- function(data_exp, method_name, ids_list) {
  if (!(method_name %in% names(ids_list))) return(NULL)
  
  ids <- ids_list[[method_name]]
  list(
    pos_high = data_exp$pos_high[triplet_id %in% ids[label == 1, triplet_id]],
    neg_high = data_exp$neg_high,  # negativos no se filtran
    reduction_pct = data_exp$reduction_pct
  )
}

# evaluao con datos filtrados
results_filtered <- lapply(expanded_data_by_reduction, function(data_exp) {
  
  # para cada método crea dataset filtrado y evalua

  lapply(method_configs, function(cfg) {

    ds_filtered <- create_filtered_data(data_exp, cfg$name, power_ids)

    if (is.null(ds_filtered)) return(NULL)

    list(
      global = evaluate_global_wrapper(ds_filtered, cfg, use_null = TRUE, aggregation_level = "triplet"),
      dynamic = evaluate_by_dynamic(ds_filtered, cfg, use_null = TRUE, aggregation_level = "stage"),
      stage = evaluate_by_stage_wrapper(ds_filtered, cfg, use_null = TRUE, aggregation_level = "stage")      
    ) 
  })
})

# combinación de resultados
metrics_global_filtered <- rbindlist(lapply(results_filtered, function(x) {
  rbindlist(lapply(x, function(y) y$global))
}))
metrics_dynamic_filtered <- rbindlist(lapply(results_filtered, function(x) {
  rbindlist(lapply(x, function(y) y$dynamic))
}))
metrics_stage_filtered <- rbindlist(lapply(results_filtered, function(x) {
  rbindlist(lapply(x, function(y) y$stage))
}))

# guardado 
fwrite(metrics_global_filtered, paste0(output_dir, "metrics_global_filtered.csv"))
fwrite(metrics_dynamic_filtered, paste0(output_dir, "metrics_dynamic_filtered.csv"))
fwrite(metrics_stage_filtered, paste0(output_dir, "metrics_stage_filtered.csv"))

################################################################################
# Métricas originales en dataset de intersección (global, por dinámica y por etapa)
################################################################################

# métodos a usar para dataset de intersección
configs_intersection <- list(
  list(name = "GAM-Doble", type = "Doble", is_gam = TRUE),
  list(name = "Estratificado-Doble", type = "Doble", is_gam = FALSE)
)
  
results_intersection <- lapply(expanded_data_by_reduction, function(data_exp) {
  # dataset de intersección
  ds_int <- list(
    pos_high = data_exp$pos_high[triplet_id %in% triplets_intersection],  # solo positivos en intersección
    neg_high = data_exp$neg_high,  # todos los negativos, sin filtrar
    reduction_pct = data_exp$reduction_pct
  )
  lapply(configs_intersection, function(cfg) {
    list(  # wrappers para detección de señal y cálculo de IC95 por bootstrap
      global = evaluate_global_wrapper(ds_int, cfg, use_null = cfg$is_gam, aggregation_level = "triplet"),
      dynamic = evaluate_by_dynamic(ds_int, cfg, use_null = cfg$is_gam, aggregation_level = "stage"),
      stage = evaluate_by_stage_wrapper(ds_int, cfg, use_null = cfg$is_gam, aggregation_level = "stage")
    )
  })
})

metrics_global_intersection <- rbindlist(lapply(results_intersection, function(x) {
  rbindlist(lapply(x, function(y) y$global))
}))
metrics_dynamic_intersection <- rbindlist(lapply(results_intersection, function(x) {
  rbindlist(lapply(x, function(y) y$dynamic))
}))
metrics_stage_intersection <- rbindlist(lapply(results_intersection, function(x) {
  rbindlist(lapply(x, function(y) y$stage))
}))

# guardado
fwrite(metrics_global_intersection, paste0(output_dir, "metrics_global_intersection.csv"))
fwrite(metrics_dynamic_intersection, paste0(output_dir, "metrics_dynamic_intersection.csv"))
fwrite(metrics_stage_intersection, paste0(output_dir, "metrics_stage_intersection.csv"))

################################################################################
# TRADUCCIONES DE ETIQUETAS
################################################################################

# traducción de dinámicas
dynamic_labels <- c(
  "uniform" = "Uniforme",
  "increase" = "Incremento",
  "decrease" = "Disminución",
  "plateau" = "Meseta",
  "inverse_plateau" = "Valle"
)

# traducción de etapas
nichd_labels <- c(
  "term_neonatal" = "Neonato a término",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temprana",
  "late_adolescence" = "Adolescencia tardía"
)

# pares de métodos a comparar
method_pairs <- list(
  list(gam = "GAM-logIOR", classic = "Estratificado-IOR", label = "IOR"),
  list(gam = "GAM-RERI", classic = "Estratificado-RERI", label = "RERI"),
  list(gam = "GAM-Doble", classic = "Estratificado-Doble", label = "Doble")
)

# métricas a graficar
metrics_to_plot <- c("sensitivity", "specificity", "PPV", "NPV", "F1")
metric_labels <- c(
  "sensitivity" = "Sensibilidad",
  "specificity" = "Especificidad", 
  "PPV" = "Valor Predictivo Positivo",
  "NPV" = "Valor Predictivo Negativo",
  "F1" = "F1-Score"
)

################################################################################
# Gráficos - Métricas globales
################################################################################

# loop para generar gráficos comparando métricas globales por método
for (pair in method_pairs) {
  
  data_pair <- metrics_global_original[method %in% c(pair$gam, pair$classic)]
  
  # datos en formato largo
  data_long <- rbindlist(lapply(metrics_to_plot, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    
    data.table(
      method = data_pair$method,
      reduction_pct = data_pair$reduction_pct,
      metric = m,
      value = data_pair[[m]],
      lower = data_pair[[lower_col]],
      upper = data_pair[[upper_col]]
    )
  }))
  
  data_long[, metric_label := metric_labels[metric]]
  data_long[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  
  p <- ggplot(data_long, aes(x = reduction_pct, y = value, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = reduction_levels) +
    labs(
      title = sprintf("Métricas Globales - %s", pair$label),
      subtitle = "Dataset original, IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Valor de la Métrica",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_global_metrics_", tolower(pair$label), "_original.png"),
    p, width = 12, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_global_metrics_", tolower(pair$label), "_original.svg"),
    p, width = 12, height = 10, device = svglite
  )
  print(p)
}

# Sensibilidad para todos los métodos
data_sens <- metrics_global_original[, .(
  method, reduction_pct, 
  sensitivity, sensitivity_lower, sensitivity_upper
)]

data_sens[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
data_sens[, Detection := fcase(
  grepl("logIOR|IOR", method), "IOR",
  grepl("RERI", method), "RERI",
  grepl("Doble", method), "Doble"
)]

p_sens_all <- ggplot(data_sens, 
                     aes(x = reduction_pct, y = sensitivity, 
                         color = method, linetype = Type)) +
  geom_point(size = 4, position = position_dodge(width = 6)) +
  geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
                width = 4, alpha = 0.5, position = position_dodge(width = 6)) +
  scale_color_manual(values = c(
    "GAM-logIOR" = "#27AE60", 
    "Estratificado-IOR" = "#E67E22",
    "GAM-RERI" = "#2ECC71", 
    "Estratificado-RERI" = "#E74C3C",
    "GAM-Doble" = "#16A085", 
    "Estratificado-Doble" = "#C0392B"
  )) +
  scale_linetype_manual(values = c("GAM" = "solid", "Estratificado" = "dashed")) +
  scale_x_continuous(breaks = reduction_levels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Sensibilidad - Todos los Métodos - Nivel Global",
    subtitle = "Dataset original - IC 95% Bootstrap",
    x = "Reducción del Dataset (%)",
    y = "Sensibilidad",
    color = "Método",
    linetype = "Tipo"
  )

print(p_sens_all)

ggsave(paste0(output_dir, "fig_sensitivity_all_methods_original.png"),
       p_sens_all, width = 12, height = 8, dpi = 300)

ggsave(paste0(output_dir, "fig_sensitivity_all_methods_original.svg"),
       p_sens_all, width = 12, height = 8, device = svglite)

################################################################################
# Gráficos - Métricas por dinámica
################################################################################

# loop para generar gráficos comparando métricas por dinámica, por método
for (pair in method_pairs) {
  
  data_pair <- metrics_dynamic_original[method %in% c(pair$gam, pair$classic)]
  
  if (nrow(data_pair) == 0) next
  
  data_pair[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  
  p_dyn_facet <- ggplot(data_pair, 
                        aes(x = reduction_pct, y = sensitivity, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
                  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ dynamic, scales = "free_y", ncol = 2, labeller = labeller(dynamic = dynamic_labels)) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = seq(0, 90, 30)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Sensibilidad por Dinámica - %s", pair$label),
      subtitle = "IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Sensibilidad",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_dynamic_", tolower(pair$label), "_original.png"),
    p_dyn_facet, width = 12, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_dynamic_", tolower(pair$label), "_original.svg"),
    p_dyn_facet, width = 12, height = 10, device = svglite
  )
  print(p_dyn_facet)
}

################################################################################
# Gráficos - Métricas por etapa
################################################################################

# loop para generar gráficos comparando métricas por etapa, por método
for (pair in method_pairs) {
  
  data_pair <- metrics_stage_original[method %in% c(pair$gam, pair$classic)]
  
  if (nrow(data_pair) == 0) next
  
  data_pair[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  data_pair[, nichd := factor(nichd, levels = niveles_nichd)]
  
  p_stage_facet <- ggplot(data_pair, 
                          aes(x = reduction_pct, y = sensitivity, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
                  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ nichd, scales = "free_y", ncol = 4, labeller = labeller(nichd = nichd_labels)) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = seq(0, 90, 30)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Sensibilidad por Etapa - %s", pair$label),
      subtitle = "IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Sensibilidad",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_stage_", tolower(pair$label), "_original.png"),
    p_stage_facet, width = 16, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_stage_", tolower(pair$label), "_original.svg"),
    p_stage_facet, width = 16, height = 10, device = svglite
  )
  print(p_stage_facet)
}

################################################################################
# Gráficos - Métricas globales filtradas
################################################################################

# loop para generar gráficos comparando métricas globales por método (filtrado)
for (pair in method_pairs) {
  
  data_pair <- metrics_global_filtered[method %in% c(pair$gam, pair$classic)]
  
  if (nrow(data_pair) == 0) next
  
  data_long <- rbindlist(lapply(metrics_to_plot, function(m) {
    data.table(
      method = data_pair$method,
      reduction_pct = data_pair$reduction_pct,
      metric = m,
      value = data_pair[[m]],
      lower = data_pair[[paste0(m, "_lower")]],
      upper = data_pair[[paste0(m, "_upper")]]
    )
  }))
  
  data_long[, metric_label := metric_labels[metric]]
  data_long[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  
  p <- ggplot(data_long, aes(x = reduction_pct, y = value, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = reduction_levels) +
    labs(
      title = sprintf("Métricas Globales - %s ", pair$label),
      subtitle = "Dataset filtrado por poder - IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Valor de la Métrica",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_global_metrics_", tolower(pair$label), "_filtered.png"),
    p, width = 12, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_global_metrics_", tolower(pair$label), "_filtered.svg"),
    p, width = 12, height = 10, device = svglite
  )
  print(p)
}

################################################################################
# Gráficos - Métricas por dinámica filtradas
################################################################################

# loop para generar gráficos comparando métricas por dinámica, por método (filtrado)
for (pair in method_pairs) {
  
  data_pair <- metrics_dynamic_filtered[method %in% c(pair$gam, pair$classic)]
  
  if (nrow(data_pair) == 0) next
  
  data_pair[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  
  p_dyn_facet <- ggplot(data_pair, 
                        aes(x = reduction_pct, y = sensitivity, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
                  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ dynamic, scales = "free_y", ncol = 2, labeller = labeller(dynamic = dynamic_labels)) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = seq(0, 90, 30)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Sensibilidad por Dinámica - %s", pair$label),
      subtitle = "Dataset filtrado por poder - IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Sensibilidad",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_dynamic_", tolower(pair$label), "_filtered.png"),
    p_dyn_facet, width = 12, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_dynamic_", tolower(pair$label), "_filtered.svg"),
    p_dyn_facet, width = 12, height = 10, device = svglite
  )
  print(p_dyn_facet)
}

################################################################################
# Gráficos - Métricas por etapa filtradas
################################################################################

# loop para generar gráficos comparando métricas por etapa, por método (filtrado)
for (pair in method_pairs) {
  
  data_pair <- metrics_stage_filtered[method %in% c(pair$gam, pair$classic)]
  
  if (nrow(data_pair) == 0) next
  
  data_pair[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
  data_pair[, nichd := factor(nichd, levels = niveles_nichd)]
  
  p_stage_facet <- ggplot(data_pair, 
                          aes(x = reduction_pct, y = sensitivity, color = Type)) +
    geom_point(size = 4, position = position_dodge(width = 6)) +
    geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
                  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
    facet_wrap(~ nichd, scales = "free_y", ncol = 4, labeller = labeller(nichd = nichd_labels)) +
    scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
    scale_x_continuous(breaks = seq(0, 90, 30)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Sensibilidad por Etapa - %s ", pair$label),
      subtitle = "Dataset filtrado por poder - IC 95%",
      x = "Reducción del Dataset (%)",
      y = "Sensibilidad",
      color = "Método"
    ) 
  
  ggsave(
    paste0(output_dir, "fig_stage_", tolower(pair$label), "_filtered.png"),
    p_stage_facet, width = 16, height = 10, dpi = 300
  )
  ggsave(
    paste0(output_dir, "fig_stage_", tolower(pair$label), "_filtered.svg"),
    p_stage_facet, width = 16, height = 10, device = svglite
  )
  print(p_stage_facet)
  message(sprintf("  Guardado: fig_stage_%s_filtered.png", tolower(pair$label)))
}

################################################################################
# Gráficos - Métricas en dataset de intersección
################################################################################
  
# Métricas globales
data_inter_global <- metrics_global_intersection
  
data_long <- rbindlist(lapply(metrics_to_plot, function(m) {
  data.table(
    method = data_inter_global$method,
    reduction_pct = data_inter_global$reduction_pct,
    metric = m,
    value = data_inter_global[[m]],
    lower = data_inter_global[[paste0(m, "_lower")]],
    upper = data_inter_global[[paste0(m, "_upper")]]
  )
}))
  
data_long[, metric_label := metric_labels[metric]]
data_long[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]

p_inter_global <- ggplot(data_long, 
  aes(x = reduction_pct, y = value, color = Type)) +
  geom_point(size = 4, position = position_dodge(width = 6)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
  scale_x_continuous(breaks = reduction_levels) +
  labs(
    title = "Métricas - Dataset de Intersección",
    subtitle = sprintf("Dataset filtrado por poder - IC 95 (N= %d tripletes comunes)", 
    length(triplets_intersection)),
    x = "Reducción del Dataset (%)",
    y = "Valor de la Métrica",
    color = "Método"
  ) 
print(p_inter_global)

ggsave(paste0(output_dir, "fig_intersection_global_metrics.png"),
p_inter_global, width = 12, height = 10, dpi = 300)
ggsave(paste0(output_dir, "fig_intersection_global_metrics.svg"),
p_inter_global, width = 12, height = 10, device = svglite)

# por dinámica
data_inter_dyn <- metrics_dynamic_intersection
data_inter_dyn[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]

p_inter_dyn <- ggplot(data_inter_dyn, 
  aes(x = reduction_pct, y = sensitivity, color = Type)) +
  geom_point(size = 4, position = position_dodge(width = 6)) +
  geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
  facet_wrap(~ dynamic, scales = "free_y", ncol = 2, labeller = labeller(dynamic = dynamic_labels)) +
  scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
  scale_x_continuous(breaks = seq(0, 90, 30)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Sensibilidad por Dinámica - Dataset de Intersección",
    subtitle = "IC 95%",
    x = "Reducción del Dataset (%)",
    y = "Sensibilidad",
    color = "Método"
  ) 
print(p_inter_dyn)

ggsave(paste0(output_dir, "fig_intersection_dynamic.png"),
p_inter_dyn, width = 12, height = 10, dpi = 300)
ggsave(paste0(output_dir, "fig_intersection_dynamic.svg"),
p_inter_dyn, width = 12, height = 10, device = svglite)

# por etapa 
data_inter_stage <- metrics_stage_intersection
data_inter_stage[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
data_inter_stage[, nichd := factor(nichd, levels = niveles_nichd)]

p_inter_stage <- ggplot(data_inter_stage, 
  aes(x = reduction_pct, y = sensitivity, color = Type)) +
  geom_point(size = 4, position = position_dodge(width = 6)) +
  geom_errorbar(aes(ymin = sensitivity_lower, ymax = sensitivity_upper), 
  width = 4, alpha = 0.6, position = position_dodge(width = 6)) +
  facet_wrap(~ nichd, scales = "free_y", ncol = 4, labeller = labeller(nichd = nichd_labels)) +
  scale_color_manual(values = c("GAM" = "#16A085", "Estratificado" = "#C0392B")) +
  scale_x_continuous(breaks = seq(0, 90, 30)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Sensibilidad por Etapa - Dataset de Intersección",
    subtitle = "IC 95%",
    x = "Reducción del Dataset (%)",
    y = "Sensibilidad",
    color = "Método"  )

print(p_inter_stage)

ggsave(paste0(output_dir, "fig_intersection_stage.png"),
p_inter_stage, width = 16, height = 10, dpi = 300)
ggsave(paste0(output_dir, "fig_intersection_stage.svg"),
p_inter_stage, width = 16, height = 10, device = svglite)

################################################################################
# Gráfico - Degradación de sensibilidad
################################################################################

# calculo degradación relativa (% de pérdida desde el basal)
data_degr <- metrics_global_original[, .(
  method, reduction_pct, sensitivity
)]

baseline_sens <- data_degr[reduction_pct == 0, .(method, baseline = sensitivity)]
data_degr <- merge(data_degr, baseline_sens, by = "method")
data_degr[, degradation_pct := 100 * (1 - sensitivity / baseline)]
data_degr[, Type := ifelse(grepl("GAM", method), "GAM", "Estratificado")]
data_degr[, Detection := fcase(
  grepl("logIOR|IOR", method), "IOR",
  grepl("RERI", method), "RERI",
  grepl("Doble", method), "Doble"
)]

p_degr <- ggplot(data_degr, 
                 aes(x = reduction_pct, y = degradation_pct, 
                     color = Detection, linetype = Type)) +
  geom_point(size = 4, position = position_dodge(width = 6)) +
  scale_color_manual(values = c(
    "IOR" = "#E74C3C", "RERI" = "#3498DB", "Doble" = "#2ECC71"
  )) +
  scale_linetype_manual(values = c("GAM" = "solid", "Estratificado" = "dashed")) +
  scale_x_continuous(breaks = reduction_levels) +
  labs(
    title = "Degradación de Sensibilidad por Reducción del Dataset",
    subtitle = "Pérdida porcentual respecto al basal (0% reducción)",
    x = "Reducción del Dataset (%)",
    y = "Degradación de Sensibilidad (%)",
    color = "Detección",
    linetype = "Método"
  ) 

print(p_degr)

ggsave(paste0(output_dir, "fig_sensitivity_degradation.png"),
       p_degr, width = 12, height = 8, dpi = 300)
ggsave(paste0(output_dir, "fig_sensitivity_degradation.svg"),
       p_degr, width = 12, height = 8, device = svglite)
################################################################################
# Archivos generados
################################################################################
summary_txt <- sprintf("

================================================================================
Archivos:
  - metrics_global_original.csv
  - metrics_dynamic_original.csv
  - metrics_stage_original.csv
  - metrics_global_filtered.csv
  - metrics_dynamic_filtered.csv
  - metrics_stage_filtered.csv
  - metrics_global_intersection.csv 
  - metrics_dynamic_intersection.csv 
  - metrics_stage_intersection.csv 

Directorio: %s
================================================================================
", output_dir)

cat(summary_txt)
writeLines(summary_txt, paste0(output_dir, "RESUMEN.txt"))
