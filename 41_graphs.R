################################################################################
# Faceted plot visualization script
# Script: 41_graphs.R
################################################################################

source("00_functions.R", local = TRUE)

output_dir <- paste0("./results/", suffix, "/metrics_results/")
fig_output_dir <- paste0(output_dir, "facet_figures/")

# Create output directory if it does not exist
dir.create(fig_output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Preprocessing
################################################################################

# Display labels for NICHD developmental stages
nichd_labels <- c(
  "term_neonatal" = "Neonato a term.",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temp.",
  "late_adolescence" = "Adolescencia tardía"
)

# Metric display labels
metric_labels <- c(
  "AUC" = "AUC",
  "sensitivity" = "Sensibilidad",
  "specificity" = "Especificidad",
  "PPV" = "VPP",
  "NPV" = "VPN",
  "F1" = "F1-Score"
)

# Signal dynamic display labels
dynamic_labels <- c(
  "uniform" = "Uniforme",
  "increase" = "Aumento",
  "plateau" = "Meseta",
  "decrease" = "Disminución",
  "inverse_plateau" = "Valle"
)

# Method pairs to compare (GAM vs. stratified, for IOR and RERI)
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

# Dataset versions to process
dataset_versions <- c("original", "filtered", "intersection")

# Helper: readable title for each dataset version
version_title_label <- function(version) {
  switch(version,
    "original" = "Dataset Original",
    "filtered" = "Dataset Filtrado por Poder",
    "intersection" = "Dataset de Intersección",
    version
  )
}

# Build color and linetype scales for a given method pair
make_scales <- function(pair) {
  list(
    color = setNames(c("#16A085", "#C0392B"), c(pair$gam, pair$classic)),
    linetype = setNames(c("solid", "dashed"), c(pair$gam, pair$classic)),
    labels = setNames(c(pair$gam_label, pair$classic_label), c(pair$gam, pair$classic))
  )
}

################################################################################
# Dynamic pattern plot
################################################################################

# Plot of tangential weight functions by signal dynamic type
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
# Data loading
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

# Add method_type column and convert nichd/dynamic to factors across all versions
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
# Data preparation for faceted plots
################################################################################

# Reshape stage-level metrics from wide to long format for ggplot2 faceting.
#
# Filters the data to the two methods in the given pair, then pivots
# all six performance metrics into a single long-format data.table.
#
# Arguments:
#  dt : data.table with per-stage metrics in wide format
#  pair : list with method identifiers and display labels (gam, classic, etc.)
#
# Returns: long-format data.table ready for ggplot2

prepare_facet_data <- function(dt, pair) {
  # Keep only the two methods belonging to this pair
  dt_filtered <- dt[method %in% c(pair$gam, pair$classic)]
  # Metrics to include in the plot
  metrics <- c("AUC", "sensitivity", "specificity", "PPV", "NPV", "F1")
  # Pivot to long format, one row per method × stage × metric
  dt_long <- rbindlist(lapply(metrics, function(m) {
    lower_col <- paste0(m, "_lower")
    upper_col <- paste0(m, "_upper")
    
    # Skip metric if confidence interval columns are missing
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
  # Assign display label for each method
  dt_long[, method_label := ifelse(
    method == pair$gam,
    pair$gam_label,
    pair$classic_label
  )]
  # Set metric factor order to match the desired panel sequence
  dt_long[, metric_label := factor(
    metric_label,
    levels = metric_labels
  )]
  return(dt_long)
}

# Reshape dynamic-level metrics from wide to long format for faceting
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

# Reshape global (unstratified) metrics from wide to long format for faceting
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
# Main plot generation function
################################################################################

# Generate all faceted stage-level plots across all dataset versions and method pairs.
# This definition overrides the one above, consolidating logic via make_scales().

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
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Valor de la Métrica",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs( title = sprintf("Métricas por Etapa - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version))
}

# Faceted plot: metrics × signal dynamics
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
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Valor de la Métrica",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs( title = sprintf("Métricas por Dinámica - %s vs %s", pair$gam_label, pair$classic_label), subtitle = version_title_label(version)) 
}

# Faceted plot: global (unstratified) metrics via facet_wrap per metric
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
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 4),
      name = "Valor de la Métrica",
      # small multiplicative padding so errorbars don't clip at panel edge
      expand = expansion(mult = 0.08)
    ) +
    labs(
      title = sprintf("Métricas Globales - %s vs %s", pair$gam_label, pair$classic_label),
      subtitle = version_title_label(version)
    ) 
}

################################################################################
# Plot saving
################################################################################

save_plot <- function(p, file_suffix, width, height) {
  png_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".png")
  svg_path <- paste0(fig_output_dir, "fig_facet_", file_suffix, ".svg")
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(svg_path, p, width = width, height = height, device = svglite)
}

# Stage-stratified faceted plots — iterate over all versions and method pairs
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

# Dynamic-stratified faceted plots — iterate over all versions and method pairs
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

# Global (unstratified) faceted plots — iterate over all versions and method pairs
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

# Execute all faceted plot generators
generate_all_facet_plots()
generate_all_facet_plots_dynamic()
generate_all_facet_plots_global()

################################################################################
# Positive triplet plot generation
################################################################################

# File paths for positive triplet data
ruta_positive_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_info <- "./drug.csv"
output_dir_positive <- paste0(fig_output_dir, "positive_triplets/")

# Positives detected in the original (non-simulated) dataset, from 40_network.R
ruta_network_graphs <- paste0("./results/", suffix, "/network/network_triplets_for_graphs.rds")

# Verify that all required input files exist before proceeding
check_positive_files <- function() {
  archivos_requeridos <- c(ruta_positive_results, ruta_ade_raw, ruta_drug_info)
  nombres <- c("positive_triplets_results.rds", "ade_raw.csv", "drug.csv")
  return(TRUE)
}

################################################################################
# Loading and preprocessing functions
################################################################################

# Load and preprocess positive triplet data.
#
# Arguments:
#   ruta_results : path to the RDS file with augmentation results
#   ruta_ade: path to the CSV file with raw ADE reports
#
# Returns: list with two elements:
#   $results: data.table filtered to baseline cases (reduction_pct == 0, successful)
#   $ade_raw: preprocessed ADE data.table
#
# Only baseline, fully successful triplets are retained for plotting

carga_datos_positive <- function(ruta_results, ruta_ade) {
  # Load augmentation results
  dt <- readRDS(ruta_results)
  
  # Keep only baseline (no reduction) and fully successful cases
  dt <- dt[reduction_pct == 0 & model_success == TRUE & injection_success == TRUE]
  
  # Load raw ADE reports
  ade_raw_dt <- fread(ruta_ade)
  
  # Process sex variable if sex-stratified analysis is enabled
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

# Build a canonical drug ID translation table.
#
# Arguments:
#  ruta_drug_info : path to the drug CSV file
#
# Returns: data.table with columns atc_concept_id -> canonical_id.
#
# Drug variants sharing the same base name (text before the first ; or ,)
# are collapsed to a single canonical ID (the minimum numeric ID in the group)

translation_table <- function(ruta_drug_info) {
  drug_info_original <- fread(ruta_drug_info)
  drug_info_original[, atc_concept_id := as.character(atc_concept_id)]
  drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]
  
  # Build canonical mapping: one representative ID per base drug name
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

# Apply drug ID canonicalization to the ADE dataset and set keys for fast lookup.
#
# Arguments:
#  ade_raw_dt: raw ADE data.table
#  translation_table: data.table mapping atc_concept_id -> canonical_id
#
# Returns: processed data.table with key set on (atc_concept_id, meddra_concept_id, nichd_num)

translate_ade <- function(ade_raw_dt, translation_table) {
  # Cast ID to character and merge with the translation table
  ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]
  
  ade_raw_dt <- merge(
    ade_raw_dt,
    translation_table[, .(atc_concept_id, canonical_id)],
    by = "atc_concept_id",
    all.x = TRUE
  )
  
  # Replace original IDs with their canonical counterparts
  ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
  ade_raw_dt[, canonical_id := NULL]
  
  # Remove duplicate report–drug–event combinations
  nrow_antes <- nrow(ade_raw_dt)
  cols_unicos <- c("safetyreportid", "atc_concept_id", "meddra_concept_id")
  if (include_sex) cols_unicos <- c(cols_unicos, "sex")
  
  ade_raw_dt <- unique(ade_raw_dt, by = cols_unicos)
  message(sprintf("  Duplicados eliminados: %d", nrow_antes - nrow(ade_raw_dt)))
  
  # Create ordered NICHD factor and a numeric version for indexing
  ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
  ade_raw_dt[, nichd_num := as.integer(nichd)]
  
  # Set key for efficient subset lookups downstream
  setkey(ade_raw_dt, atc_concept_id, meddra_concept_id, nichd_num)
  
  return(ade_raw_dt)
}

################################################################################
# Report count calculation
################################################################################

# Count spontaneous reports by NICHD stage for a specific drug-drug-event triplet.
#
# Arguments:
#  drug_a: ATC concept ID of drug A
#  drug_b: ATC concept ID of drug B
#  meddra_event: MedDRA concept ID of the adverse event
#  ade_dt: preprocessed ADE data.table (output of translate_ade)
#
# Returns: data.table with 7 rows (one per stage) and columns
#   n_a, n_b, n_ab, n_evento, n_evento_ab
#
# Note: computing counts on the fly here is highly inefficient for the full
# pipeline; ideally this should have been precomputed in 10_augmentation.R.
# Re-running that script would take ~2.5 days, so we calculate counts here instead

calculate_triplet_counts <- function(drug_a, drug_b, meddra_event, ade_dt) {
  # Unique report IDs for each drug and the event
  ids_a <- unique(ade_dt[atc_concept_id == drug_a, safetyreportid])
  ids_b <- unique(ade_dt[atc_concept_id == drug_b, safetyreportid])
  ids_event <- unique(ade_dt[meddra_concept_id == meddra_event, safetyreportid])
  
  # Co-exposure and exclusive-exposure report sets
  ids_ab <- intersect(ids_a, ids_b)
  ids_a_only <- setdiff(ids_a, ids_b)
  ids_b_only <- setdiff(ids_b, ids_a)
  
  # Helper: count unique reports per NICHD stage for a given ID subset
  count_per_stage <- function(ids_subset) {
    if (length(ids_subset) == 0) {
      return(data.table(nichd_num = 1:7, n = 0L))
    }
    counts <- ade_dt[safetyreportid %in% ids_subset, .(n = uniqueN(safetyreportid)), by = nichd_num]
    # Fill in any missing stages with zero counts
    etapas_completas <- data.table(nichd_num = 1:7)
    counts <- merge(etapas_completas, counts, by = "nichd_num", all.x = TRUE)
    counts[is.na(n), n := 0L]
    return(counts[order(nichd_num)])
  }
  
  # Compute per-stage counts for each exposure group
  n_a_dt <- count_per_stage(ids_a_only)
  n_b_dt <- count_per_stage(ids_b_only)
  n_ab_dt <- count_per_stage(ids_ab)
  n_event_dt <- count_per_stage(ids_event)
  n_event_ab_dt <- count_per_stage(intersect(ids_event, ids_ab))
  
  # Combine into a single result table
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
# Triplet data expansion
################################################################################

# Expand a single triplet result row into a stage-level data.table with
# metrics and report counts merged together.
#
# Arguments:
#  row: single-row data.table from the results object
#  ade_dt: preprocessed ADE data.table (can be NULL if counts are embedded)
#  precomputed_counts: optional data.table of precomputed counts keyed by triplet_id
#
# Returns: long data.table with one row per NICHD stage, containing GAM and
#  classic metric estimates, confidence bounds, and report counts.
#
# Metrics are extracted from list columns in the RDS results.
# Injected counts (from the augmentation diagnostics) are added to the
# observed event counts before plotting

expand_triplets_counts <- function(row, ade_dt, precomputed_counts = NULL) {
  # Extract NICHD stage indices
  stages <- unlist(row$stage)
  
  # GAM metrics
  gam_log_ior <- unlist(row$log_ior)
  gam_log_ior_lower <- unlist(row$log_ior_lower90)
  gam_log_ior_upper <- gam_log_ior + (gam_log_ior - gam_log_ior_lower)
  
  gam_reri <- unlist(row$reri_values)
  gam_reri_lower <- unlist(row$reri_lower90)
  gam_reri_upper <- unlist(row$reri_upper90)
  
  # Classic (stratified) metrics
  cls_log_ior <- unlist(row$log_ior_classic)
  cls_log_ior_lower <- unlist(row$log_ior_classic_lower90)
  cls_log_ior_upper <- cls_log_ior + (cls_log_ior - cls_log_ior_lower)
  
  cls_reri <- unlist(row$RERI_classic)
  cls_reri_lower <- unlist(row$RERI_classic_lower90)
  cls_reri_upper <- unlist(row$RERI_classic_upper90)
  
  # Resolve report counts: supports precomputed, embedded (from 40_network.R),
  # or on-the-fly calculation from ade_dt
  if (!is.null(precomputed_counts)) {
    counts_dt <- precomputed_counts[triplet_id == row$triplet_id]
    if (nrow(counts_dt) == 0)
      counts_dt <- calculate_triplet_counts(row$drugA, row$drugB, row$meddra, ade_dt)
  } else if (!is.null(row$counts_by_stage[[1]])) {
    # Counts are already embedded in the network results row
    counts_dt <- row$counts_by_stage[[1]]
  } else {
    counts_dt <- calculate_triplet_counts(row$drugA, row$drugB, row$meddra, ade_dt)
  }
  
  # Extract injected counts per stage from diagnostics$injection_by_stage
  # (columns: nichd_num, N); fill missing stages with 0
    inj_by_stage <- tryCatch({
    diag <- row$diagnostics[[1]]
    if (!is.null(diag) && !is.null(diag$injection_by_stage)) {
      ibs <- merge(data.table(nichd_num = 1:7), diag$injection_by_stage, by = "nichd_num", all.x = TRUE)
      ibs[is.na(N), N := 0L]
      ibs[order(nichd_num)]$N
    } else { rep(0L, 7) }
  }, error = function(e) rep(0L, 7))
  
  # Add injected reports to observed event counts for plotting purposes
  counts_dt[, n_evento_ab := n_evento_ab + inj_by_stage]
  counts_dt[, n_evento := n_evento + inj_by_stage]

  # Guard against mismatched vector lengths across sources
  n <- min(length(stages), length(gam_log_ior), nrow(counts_dt))
  if (n == 0) return(NULL)
  
  # Assemble final stage-level data.table
  result <- data.table(
    triplet_id = row$triplet_id,
    dynamic = row$dynamic,
    stage_num = stages[1:n],
    # GAM metrics
    gam_log_ior = gam_log_ior[1:n],
    gam_log_ior_lower = gam_log_ior_lower[1:n],
    gam_log_ior_upper = gam_log_ior_upper[1:n],
    gam_reri = gam_reri[1:n],
    gam_reri_lower = gam_reri_lower[1:n],
    gam_reri_upper = gam_reri_upper[1:n],
    # Classic metrics
    cls_log_ior = cls_log_ior[1:n],
    cls_log_ior_lower = cls_log_ior_lower[1:n],
    cls_log_ior_upper = cls_log_ior_upper[1:n],
    cls_reri = cls_reri[1:n],
    cls_reri_lower = cls_reri_lower[1:n],
    cls_reri_upper = cls_reri_upper[1:n],
    # Report counts
    n_a = counts_dt$n_a[1:n],
    n_b = counts_dt$n_b[1:n],
    n_ab = counts_dt$n_ab[1:n],
    n_evento = counts_dt$n_evento[1:n],
    n_evento_ab = counts_dt$n_evento_ab[1:n]
  )
  
  return(result)
}

################################################################################
# Triplet-level visualization
################################################################################

# Generate a dual-axis plot for a single triplet: metric trajectory (primary Y)
# overlaid on report count bars (secondary log Y)
#
# Arguments:
#  plot_dt: expanded stage-level data.table for this triplet
#  metric_col: column name of the main metric values
#  lower_col: column name of the lower confidence bound
#  upper_col: column name of the upper confidence bound
#  y_label: label for the primary Y axis
#  file_suffix: identifier string used in the output filename
#  y_limit: symmetric Y range for the primary axis (default ±10)
#  max_count: expected maximum report count, used to calibrate the log scale (default 5000)
#  plot_title: optional plot title (defaults to dynamic + triplet ID)
#  plot_subtitle: optional subtitle (defaults to "A + B")

graph_metrics_counts <- function(plot_dt, metric_col, lower_col, upper_col,
                                  y_label, file_suffix, y_limit = 10, max_count = 5000,
                                  plot_title = NULL, plot_subtitle = NULL) {
  # Pivot count columns to long format for bar rendering
  counts_long <- melt(
    plot_dt,
    id.vars = "stage_num",
    measure.vars = c("n_a", "n_b", "n_ab", "n_evento", "n_evento_ab"),
    variable.name = "metric",
    value.name = "count"
  )
  
  # Legend labels for each count group
  metric_labels <- c(
    n_a = "A", n_b = "B", n_ab = "A-B",
    n_evento = "Evento", n_evento_ab = "A-B-Evento"
  )
  counts_long[, metric_label := factor(metric_labels[as.character(metric)], levels = unname(metric_labels))]
  
  # Compute dodged bar X positions manually (5 groups, evenly spaced within each stage)
  bar_w <- 0.12
  counts_long[, metric_idx := as.integer(metric)]
  counts_long[, x_center := stage_num + (metric_idx - 3) * bar_w]
  counts_long[, xmin := x_center - bar_w/2]
  counts_long[, xmax := x_center + bar_w/2]
  
  # Log-linear transformation for the secondary axis:
  # maps count = 1 -> -y_limit and count = max_count -> +y_limit
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
  
  # Build the composite plot
  p <- ggplot() +
    # Report count bars (secondary axis, log-scaled)
    geom_rect(
      data = counts_long,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = metric_label),
      alpha = 0.6, color = NA
    ) +
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # Confidence interval ribbon for the metric
    geom_ribbon(
      data = plot_dt,
      aes(x = stage_num, ymin = .data[[lower_col]], ymax = .data[[upper_col]]),
      alpha = 0.3, fill = "#2c3e50"
    ) +
    # Metric trajectory line
    geom_line(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", linewidth = 1.2
    ) +
    # Metric value points
    geom_point(
      data = plot_dt,
      aes(x = stage_num, y = .data[[metric_col]]),
      color = "#2c3e50", size = 4, fill = "white", shape = 21, stroke = 1.5
    ) +
    # Axis scales
    scale_x_continuous(
      breaks = 1:7,
      labels = nichd_labels,
      name = NULL
    ) +
    scale_y_continuous(
      name = y_label,
      limits = c(-y_limit, y_limit),
      breaks = {
        # Compute readable breaks regardless of the axis range
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
    # Plot labels
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
# Positive triplet plot generation
################################################################################

# Save a single positive triplet plot to disk.
#
# Arguments:
#  p: ggplot object
#  triplet_id: numeric triplet identifier
#  dynamic: dynamic label string (used in the filename)
#  suffix: method/metric identifier string (e.g. "GAM_LogIOR")
#  dir: output directory (defaults to output_dir_positive)

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

# Generate plots for all positive triplets from the augmentation pipeline.
#
# Source: positive_triplets_results.rds (baseline only: reduction_pct == 0)
# Counts: computed on the fly from ade_raw.csv via calculate_triplet_counts()
# Metrics: GAM-LogIOR, GAM-RERI, Classic-LogIOR, Classic-RERI
# (each metric is only plotted if its CI bounds are finite)

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

# Generate plots for triplets from the ORIGINAL dataset detected as positive in 40_network.R.
#
# Source : network_triplets_for_graphs.rds
#  already filtered to positives (triplet_id %in% positives$triplet_id in 40_network.R)
#  drugA_name, drugB_name, meddra_name resolved by 40_network.R
#  counts_by_stage embedded per row (not recomputed here)
# Title: "<drugA_name> + <drugB_name> | Triplete: <id>" / "<meddra_name>"
# Metric: GAM-LogIOR only (RERI and classic metrics are NULL in this pipeline; skipped)
# Output: subfolder "network/" inside output_dir_positive

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

      # drugA_name, drugB_name, and meddra_name are guaranteed to be present (set by 40_network.R)
      ptitle <- paste0(row$drugA_name, " + ", row$drugB_name,
                          " | Triplete: ", row$triplet_id)
      psubtitle <- row$meddra_name

      # Counts come from counts_by_stage embedded in the row; ade_dt is not needed here
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

# Run both positive triplet plot generators
generate_positive_graphs_from_results()
generate_positive_graphs_from_network()

