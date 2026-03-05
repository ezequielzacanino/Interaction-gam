################################################################################
# Análisis de redes 
# Script: 40_network.R
################################################################################

pacman::p_load(
  data.table, tidyverse, mgcv, parallel, doParallel, foreach, doRNG,
  igraph, ggraph, tidygraph, scales, RColorBrewer, patchwork, graphlayouts,
  ggrepel, networkD3, htmlwidgets
)

set.seed(9427)
setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)   
source("01_theme.R")                   
theme_set(theme_base())               

################################################################################
# Configuración
################################################################################

percentil <- "p95"          
use_null_threshold <- TRUE  

min_reports_triplet <- 5   
min_stages_positive <- 1    
n_random_networks <- 500    # n redes aleatorias modelo nulo 
jaccard_threshold <- 0.1    # umbral Jaccard para considerar soporte biológico

# Parámetros grafo bipartito 
n_drugs_bipartite <- 50    # drogas (top por degree en g_S) en grafo bipartito
n_genes_bipartite <- 10     # genes (top por degree en g_B) en grafo bipartito

spline_individuales <- TRUE  
include_sex <- FALSE         
include_stage_sex <- FALSE
k_spline <- 7             
include_nichd <- FALSE       
nichd_spline <- FALSE         
bs_type <- "cs"            
select <- FALSE               
method <- "fREML"            

suffix <- paste0(
  if (spline_individuales) "si" else "",
  if (include_sex) "s" else "",
  if (include_stage_sex) "ss" else "",
  if (include_nichd) "n" else "",
  if (nichd_spline) "ns" else "",
  bs_type
)

ruta_ade_raw <- "./ade_raw.csv"
ruta_drug_gene <- "./drug_gene.csv"
ruta_drug_info <- "./drug.csv"
ruta_concept <- "./vocabulary/concept.csv"
ruta_rel <- "./vocabulary/concept_relationship.csv"
ruta_twosides <- "./twosides/TWOSIDES.csv.gz"

output_dir <- paste0("./results/", suffix, "/network/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, floor(detectCores() * 0.75))

ruta_null_dist <- paste0("./results/", suffix, "/null_distribution_results/")

null_thresholds_ior <- fread(paste0(ruta_null_dist, "null_thresholds.csv"))
thresh_col <- paste0("threshold_", percentil)  # e.g. "threshold_p95"
null_thresholds_ior <- null_thresholds_ior[, .(stage, threshold_ior = get(thresh_col))]

# Umbrales para RERI (por las dudas)
null_thresholds_reri <- fread(paste0(ruta_null_dist, "null_thresholds_reri.csv"))
null_thresholds_reri <- null_thresholds_reri[, .(stage, threshold_reri = get(thresh_col))]

# Combina ambos umbrales
null_thresholds <- merge(null_thresholds_ior, null_thresholds_reri, by = "stage")

################################################################################
# Función helper
################################################################################

# Calcula cuántas etapas NICHD son positivas para un triplete dado.
#
# Parámetros:
# lower90_vec: vector numérico con el límite inferior del IC 90% de log_IOR (una entrada por etapa)
# null_thresh: data.table con columnas stage (1:7) y threshold_ior (umbral nulo por etapa)
# use_null: lógico; TRUE = exige superar el umbral nulo además de lower90 > 0
#
# Retorna: entero con el número de etapas donde el triplete es positivo

count_positive_stages_triplet <- function(lower90_vec, null_thresh, use_null) {
  if (is.null(lower90_vec) || length(lower90_vec) == 0) return(0L)
  stages <- seq_along(lower90_vec)  
  if (use_null) {
    # busco umbral nulo por etapa con indice
    thresh <- null_thresh$threshold_ior[match(stages, null_thresh$stage)]
    as.integer(sum(lower90_vec > 0 & lower90_vec > thresh, na.rm = TRUE))
  } else {
    # si no nulo, solo nominal
    as.integer(sum(lower90_vec > 0, na.rm = TRUE))
  }
}

################################################################################
# Unificación de ids y carga de vocabulario
################################################################################

# carga vocabulario completo de OMOP (concept.csv)
concept <- fread(ruta_concept, quote = "")
# minusculas para evitar problemas de joint
concept[, `:=`(
  vocabulary_id = tolower(vocabulary_id),
  concept_class_id = tolower(concept_class_id)
)]

# Extrae ATC del vocabulario (para mapear a RxNorm)
atc_concepts <- concept[vocabulary_id == "atc", .(
  concept_id = as.character(concept_id),
  atc_code   = concept_code
)]
# extraigo grupo anatómico (letra de código ATC, para colorear nodos)
atc_concepts[, atc_group := substr(atc_code, 1, 1)]

# tabla de drogas con sus atc ids y nombres
drug_info <- fread(ruta_drug_info)
drug_info[, `:=`(
  atc_concept_id = as.character(atc_concept_id),
  # Limpio nombres
  base_name = tolower(trimws(sub("[;,].*", "", atc_concept_name)))
)]

# Para drogas con el mismo nombre (base_name) pero distinto atc id,
# elige el id mínimo como id canónico para evitar duplicados
canonical_map <- drug_info[, .(canonical_id = min(atc_concept_id)), by = base_name]

# tabla de traducción con id canónico de cada atc concept
translation_table <- merge(drug_info[, .(atc_concept_id, base_name)], canonical_map, by = "base_name")

# mapa de id canónico a nombre en chr
drug_names_map <- unique(translation_table[, .(
  atc_concept_id = canonical_id,
  drug_name = base_name
)])
setkey(drug_names_map, atc_concept_id)

# Agrego grupo ATC (letra) al mapa de nombres para colorear los nodos
drug_names_map <- merge(
  drug_names_map,
  atc_concepts[, .(concept_id, atc_group)],
  by.x = "atc_concept_id", by.y = "concept_id",
  all.x = TRUE
)

################################################################################
# Carga de datos y preprocesamiento
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

# cambio atc_concept_id por id canónica 
ade_raw_dt <- merge(ade_raw_dt, translation_table[, .(atc_concept_id, canonical_id)],
                    by = "atc_concept_id", all.x = TRUE)
ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]  # aplica el reemplazo
ade_raw_dt[, canonical_id := NULL]    # limpia columna auxiliar

# saco duplicados
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

###########
# Carga de datos droga-gen
###########

drug_gene <- fread(ruta_drug_gene)
drug_gene[, atc_concept_id := as.character(atc_concept_id)]

# unifico ids canónicos en la tabla droga-gen
drug_gene <- merge(drug_gene, translation_table[, .(atc_concept_id, canonical_id)],
                   by = "atc_concept_id", all.x = TRUE)
drug_gene[!is.na(canonical_id), atc_concept_id := canonical_id]
drug_gene[, canonical_id := NULL]
# elimino duplicados
drug_gene <- unique(drug_gene, by = c("atc_concept_id", "gene_symbol"))

message(sprintf("  Drogas únicas: %d  Asociaciones droga-gen: %d  Genes únicos: %d",
                uniqueN(drug_names_map$atc_concept_id),
                nrow(drug_gene),
                uniqueN(drug_gene$gene_symbol)))

################################################################################
# Construcción de tripletes
################################################################################

# agrego por reporte: lista de drogas, lista de eventos, y etapa 
reports_meta <- ade_raw_dt[, .(
  drugs     = list(unique(atc_concept_id[!is.na(atc_concept_id)])),
  events    = list(unique(meddra_concept_id[!is.na(meddra_concept_id)])),
  nichd_num = unique(nichd_num[!is.na(nichd_num)])[1]   # toma la primera etapa si hay varias
), by = safetyreportid]

# Genera todos los pares A-B por eventos para cada reporte
# make_triplets() de 00_functions.R
triplets_list <- lapply(seq_len(nrow(reports_meta)), function(i) {
  # definición de triplete
  if (length(reports_meta$drugs[[i]]) >= 2 && length(reports_meta$events[[i]]) >= 1)
    make_triplets(reports_meta$drugs[[i]], reports_meta$events[[i]],
                  reports_meta$safetyreportid[i], reports_meta$nichd_num[i])
})
triplets_dt <- rbindlist(triplets_list, use.names = TRUE, fill = TRUE)
rm(triplets_list, reports_meta); gc()   # libera memoria 

# filtro por requisitos mínimos
candidatos_triplets <- triplets_dt[, .N, by = .(drugA, drugB, meddra)][N >= min_reports_triplet]
candidatos_triplets[, triplet_id := .I]   # id de triplete
message(sprintf(" Tripletes candidatos: %d", nrow(candidatos_triplets)))

fwrite(candidatos_triplets, paste0(output_dir, "triplets_metadata.csv"))

################################################################################
# Ajuste GAM
################################################################################

# carga si ya se hizo
cache_file <- paste0(output_dir, "all_triplets_results.rds")

if (file.exists(cache_file)) {
  all_triplets_results <- readRDS(cache_file)
} else {
  all_results <- list()
  batch_size <- 500
  batch_num <- 0L
  n_total <- nrow(candidatos_triplets)

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("fit_gam", "ade_raw_dt",
                      "spline_individuales", "include_sex", "include_stage_sex",
                      "k_spline", "nichd_spline", "bs_type", "select"),
                envir = environment())
  clusterEvalQ(cl, { library(data.table); library(mgcv); library(MASS) })

  for (start in seq(1, n_total, by = batch_size)) {
    batch_num <- batch_num + 1L
    batch_rows <- candidatos_triplets[start:min(start + batch_size - 1, n_total)]

    message(sprintf(" Lote %d tripletes %d-%d de %d",
                    batch_num, start, start + nrow(batch_rows) - 1, n_total))

    batch_res <- foreach(idx = seq_len(nrow(batch_rows)),
                         .packages = c("data.table", "mgcv", "MASS"),
                         .errorhandling = "pass") %dopar% {
      rowt <- batch_rows[idx]
      tryCatch(
        fit_gam(
          drugA_id = rowt$drugA, drugB_id = rowt$drugB, event_id = rowt$meddra,
          ade_data = ade_raw_dt, spline_individuales = spline_individuales,
          include_sex = include_sex, include_stage_sex = include_stage_sex,
          k_spline = k_spline, nichd_spline = nichd_spline, include_nichd = include_nichd,
          bs_type = bs_type, select = select
        ),
        error = function(e) list(success = FALSE)
      )
    }

    batch_dt <- rbindlist(lapply(seq_along(batch_res), function(i) {
      res <- batch_res[[i]]
      rowt <- batch_rows[i]
      base <- data.table(triplet_id = rowt$triplet_id, drugA = rowt$drugA,
                         drugB = rowt$drugB, meddra = rowt$meddra,
                         n_reports = rowt$N)
      if (is.null(res$success) || !res$success)
        return(cbind(base, model_success = FALSE))

      cbind(base, data.table(
        model_success = TRUE,
        n_events = res$n_events,
        n_coadmin = res$n_coadmin,
        n_stages_significant  = res$n_stages_significant,
        model_aic = res$model_aic,
        stage = list(1:7),
        log_ior = list(res$log_ior),
        gam_log_ior_lower90 = list(res$log_ior_lower90)
      ))
    }), fill = TRUE)

    all_results[[batch_num]] <- batch_dt
    saveRDS(batch_dt, paste0(output_dir, "checkpoint_batch_", sprintf("%04d", batch_num), ".rds"))
    rm(batch_res, batch_dt); gc()
  }

  stopCluster(cl)
  all_triplets_results <- rbindlist(all_results, fill = TRUE)
  saveRDS(all_triplets_results, cache_file)
}

message(sprintf("Exitosos: %d/%d (%.1f%%)",
                sum(all_triplets_results$model_success, na.rm = TRUE),
                nrow(all_triplets_results),
                100 * mean(all_triplets_results$model_success, na.rm = TRUE)))

################################################################################
# Funciones para exportar objeto 
################################################################################

# crea objeto exportable para la función de graficado de 31_metrics_graphs

# Función de conteos por etapa (equivale a calculate_triplet_counts)
count_by_stage <- function(ids, ade_dt) {
  if (length(ids) == 0) return(data.table(nichd_num = 1:7, n = 0L)) 
  dt <- ade_dt[safetyreportid %in% ids, .(n = uniqueN(safetyreportid)), by = nichd_num] # numero de reportes por etapa para droga
  dt <- merge(data.table(nichd_num = 1:7), dt, by = "nichd_num", all.x = TRUE)
  dt[is.na(n), n := 0L]
  dt[order(nichd_num)]
}

# creo objetos con conteos y les aplico función anterior para tenerlos por etapa
calculate_counts_for_triplet <- function(drug_a, drug_b, meddra_event, ade_dt) {
  ids_a <- unique(ade_dt[atc_concept_id == drug_a, safetyreportid])
  ids_b <- unique(ade_dt[atc_concept_id == drug_b, safetyreportid])
  ids_event <- unique(ade_dt[meddra_concept_id == meddra_event, safetyreportid])
  ids_ab <- intersect(ids_a, ids_b)

  data.table(
    nichd_num = 1:7,
    n_a = count_by_stage(setdiff(ids_a, ids_b), ade_dt)$n,
    n_b = count_by_stage(setdiff(ids_b, ids_a), ade_dt)$n,
    n_ab = count_by_stage(ids_ab, ade_dt)$n,
    n_evento = count_by_stage(ids_event, ade_dt)$n,
    n_evento_ab = count_by_stage(intersect(ids_event, ids_ab), ade_dt)$n
  )
}

# Solo tripletes con modelo exitoso y con IOR calculado
dt_graph_base <- all_triplets_results[
  model_success == TRUE & !sapply(gam_log_ior_lower90, is.null)
]

# calculo de conteos por etapa para los tripletes con las funciones anteriores
counts_list <- vector("list", nrow(dt_graph_base))
for (i in seq_len(nrow(dt_graph_base))) {
  counts_list[[i]] <- calculate_counts_for_triplet(
    dt_graph_base$drugA[i],
    dt_graph_base$drugB[i],
    dt_graph_base$meddra[i],
    ade_raw_dt
  )
}

# Arma el objeto final con los nombres que espera expand_triplets_counts()
# tambien columnas dummy para las métricas que no existen en este pipeline
network_triplets_for_graphs <- dt_graph_base[, .(
  triplet_id,
  drugA,
  drugB,
  meddra,
  n_reports,
  model_success,
  reduction_pct = 0,       # requerido por el filtro de carga_datos_positive()
  injection_success = TRUE,     # ídem
  dynamic = "network",   # sin dinámica asignada; label neutro
  stage,
  # renombro para que coincida con lo que lee expand_triplets_counts()
  log_ior = log_ior,
  log_ior_lower90 = gam_log_ior_lower90,
  # métricas que no existen en este pipeline → NULL, los if() del script las ignoran
  reri_values = list(NULL),
  reri_lower90 = list(NULL),
  reri_upper90 = list(NULL),
  log_ior_classic = list(NULL),
  log_ior_classic_lower90 = list(NULL),
  RERI_classic = list(NULL),
  RERI_classic_lower90 = list(NULL),
  RERI_classic_upper90 = list(NULL),
  diagnostics = list(NULL)   # sin inyección. tryCatch devuelve rep(0L,7)
)]

# Agrega conteos como columna-lista
network_triplets_for_graphs[, counts_by_stage := counts_list]
meddra_names <- concept[
  vocabulary_id == "meddra",
  .(meddra = as.character(concept_id), meddra_name = concept_name)
]

# nombres de drogas y meddra event
network_triplets_for_graphs[, drugA_name := drug_names_map[as.character(drugA), drug_name]]
network_triplets_for_graphs[, drugB_name := drug_names_map[as.character(drugB), drug_name]]
network_triplets_for_graphs[, meddra := as.character(meddra)]

network_triplets_for_graphs <- merge(
  network_triplets_for_graphs,
  meddra_names,
  by = "meddra", all.x = TRUE
)

rm(dt_graph_base, counts_list); gc()

################################################################################
# Clasificación de señales
################################################################################

# positivo cada triplete con min_stages_positive etapas positivas (1 pero por si quiero cambiar)
positives <- if (
  "gam_log_ior_lower90" %in% names(all_triplets_results) &&
  any(all_triplets_results$model_success == TRUE, na.rm = TRUE)
) {
  all_triplets_results[
    model_success == TRUE & !sapply(gam_log_ior_lower90, is.null),
    .(n_stages_positive = count_positive_stages_triplet(
        gam_log_ior_lower90[[1]], null_thresholds, use_null_threshold)),
    by = .(triplet_id, drugA, drugB, meddra)
  ][n_stages_positive >= min_stages_positive]
} else {
  data.table(triplet_id = integer(0), drugA = character(0),
             drugB = character(0), meddra = character(0),
             n_stages_positive = integer(0))
}

# negativos
negatives <- all_triplets_results[
  model_success == TRUE & !(triplet_id %in% positives$triplet_id),
  .(triplet_id, drugA, drugB, meddra)
]

message(sprintf("Positivos: %d | Negativos: %d", nrow(positives), nrow(negatives)))
network_triplets_for_graphs <- network_triplets_for_graphs[triplet_id %in% positives$triplet_id]

saveRDS(network_triplets_for_graphs,
        paste0(output_dir, "network_triplets_for_graphs.rds"))

################################################################################
# Capa S (señal/statistical)
################################################################################

# Grafo no dirigido de interacciones farmacológicas positivas
dt_edges_stat  <- unique(positives[, .(from = pmin(drugA, drugB), to = pmax(drugA, drugB))])
vec_nodes_stat <- unique(c(dt_edges_stat$from, dt_edges_stat$to))

graph_stat_ddi <- graph_from_data_frame(dt_edges_stat, directed = FALSE)

################################################################################
# Capa B (biológica) grafo bipartito droga-Gen
################################################################################

# reduzco las conexiones a las drogas presentes en la red estadística 
dt_edges_bio <- unique(drug_gene[atc_concept_id %in% vec_nodes_stat, .(from = atc_concept_id, to = gene_symbol)])
dt_nodes_bio <- data.table(
  name = c(vec_nodes_stat, unique(dt_edges_bio$to)),
  type = c(rep("drug", length(vec_nodes_stat)), rep("gene", uniqueN(dt_edges_bio$to)))
)

graph_bio_gene <- graph_from_data_frame(dt_edges_bio, directed = FALSE, vertices = dt_nodes_bio)

################################################################################
# métricas intercapas
################################################################################

get_drug_genes <- function(id_drug) dt_edges_bio[from == id_drug, to]

dt_edge_metrics <- rbindlist(lapply(seq_len(ecount(graph_stat_ddi)), function(i) {
  nodes_edge <- ends(graph_stat_ddi, i)
  drug_a <- nodes_edge[1]; drug_b <- nodes_edge[2]
  
  genes_a <- get_drug_genes(drug_a)
  genes_b <- get_drug_genes(drug_b)
  
  n_intersect <- length(intersect(genes_a, genes_b))
  n_union     <- length(union(genes_a, genes_b))
  val_jaccard <- if (n_union > 0) n_intersect / n_union else 0
  
  val_distance <- tryCatch(distances(graph_bio_gene, v = drug_a, to = drug_b, mode = "all")[1, 1], 
                           error = function(e) Inf)
  
  data.table(drugA = drug_a, drugB = drug_b, jaccard = val_jaccard, 
             n_shared = n_intersect, bio_distance = val_distance, 
             closes_triangle = n_intersect > 0)
}))

# atributos intercapa a las aristas del grafo estadístico
E(graph_stat_ddi)$jaccard <- dt_edge_metrics$jaccard
E(graph_stat_ddi)$n_shared <- dt_edge_metrics$n_shared
E(graph_stat_ddi)$bio_distance <- dt_edge_metrics$bio_distance
E(graph_stat_ddi)$closes_triangle <- dt_edge_metrics$closes_triangle
E(graph_stat_ddi)$has_bio_support <- dt_edge_metrics$jaccard > 0

vec_degree_stat <- degree(graph_stat_ddi)
fwrite(dt_edge_metrics, paste0(output_dir, "edge_metrics_interlayer.csv"))
print(dt_edge_metrics)

################################################################################
# Modelo nulo con rewiring paralelo preservando grados
################################################################################

file_null_cache <- paste0(output_dir, "null_model_cache.rds")

# verifico si ya se corrio el rewiring
if (file.exists(file_null_cache)) {
  dt_null_results <- readRDS(file_null_cache)
} else {
  # paso data.table de aristas biológicas a data.frame para paralelización
  df_edges_bio <- as.data.frame(dt_edges_bio) 
  
  # cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("graph_stat_ddi", "graph_bio_gene", "df_edges_bio"), envir = environment()) # objetos
  clusterEvalQ(cl, library(igraph))
  
  # bucle paralelizado
  dt_null_results <- foreach(i = seq_len(n_random_networks), .combine = rbind, 
                             .packages = "igraph", .options.RNG = 9427) %dorng% {
    graph_permuted <- rewire(graph_stat_ddi, keeping_degseq(niter = ecount(graph_stat_ddi) * 100)) # rewiring
    
    metrics <- sapply(seq_len(ecount(graph_permuted)), function(j) {
      # nodos contectados por arista j
      nodes_edge <- ends(graph_permuted, j)
      # genes asociados a nodo drogaA
      genes_a <- df_edges_bio$to[df_edges_bio$from == nodes_edge[1]]
      # genes asociados a nodo drogaB
      genes_b <- df_edges_bio$to[df_edges_bio$from == nodes_edge[2]]
      
      # Jaccard entre conjuntos de genes
      n_intersect <- length(intersect(genes_a, genes_b))
      n_union     <- length(union(genes_a, genes_b))
      val_jaccard <- if (n_union > 0) n_intersect / n_union else 0
      # distancia en el grafo biológico entre los dos fármacos 
      val_dist    <- tryCatch(distances(graph_bio_gene, v = nodes_edge[1], to = nodes_edge[2], mode = "all")[1, 1], 
                              error = function(e) Inf)
      # vector con jaccard, indicador de intersección y distancia
      c(val_jaccard, n_intersect > 0, val_dist)
    })
    
    data.frame(mean_jaccard = mean(metrics[1, ]), prop_closure = mean(metrics[2, ]), 
               mean_bio_dist = mean(metrics[3, is.finite(metrics[3, ])]))
  }
  stopCluster(cl)
  dt_null_results <- as.data.table(dt_null_results)
  saveRDS(dt_null_results, file_null_cache)
}

################################################################################
# Test de fisher 
################################################################################

# pares positivos únicos
dt_pairs_pos <- unique(positives[, .(drugA = pmin(drugA, drugB), drugB = pmax(drugA, drugB))])
dt_pairs_pos[, has_signal := TRUE]

# pares negativos únicos
dt_pairs_neg <- unique(negatives[, .(drugA = pmin(drugA, drugB), drugB = pmax(drugA, drugB))])

# saco de negativos los que son positivos para otro evento
dt_pairs_neg <- dt_pairs_neg[!dt_pairs_pos, on = c("drugA", "drugB")]
dt_pairs_neg[, has_signal := FALSE]

# uno pares
dt_all_pairs_unique <- rbind(dt_pairs_pos, dt_pairs_neg)

# lista de genes por droga usando drug_gene original
genes_por_droga <- split(drug_gene$gene_symbol, drug_gene$atc_concept_id)

# calculo el jaccard y si comparten gen para todos los pares
dt_all_pairs_unique[, c("jaccard", "has_shared") := {
  
  resultados <- lapply(seq_len(.N), function(i) {
    genes_a <- genes_por_droga[[ drugA[i] ]]
    genes_b <- genes_por_droga[[ drugB[i] ]]
    
    # si alguna de las drogas no tiene genes en la base de datos, jaccard 0
    if (is.null(genes_a) || is.null(genes_b)) {
      return(list(0, FALSE))
    }
    
    n_intersect <- length(intersect(genes_a, genes_b))
    n_union <- length(union(genes_a, genes_b))
    
    val_jaccard <- if (n_union > 0) n_intersect / n_union else 0
    
    list(val_jaccard, n_intersect > 0)
  })
  
  # devuelve como columnas
  list(sapply(resultados, `[[`, 1), sapply(resultados, `[[`, 2))
}]

# proporción de negativos que comparten genes
print(paste("Porcentaje de Negativos que comparten genes:", 
            round(dt_all_pairs_unique[has_signal == FALSE, mean(has_shared) * 100], 2), "%"))

print(paste("Porcentaje de Positivos que comparten genes:", 
            round(dt_all_pairs_unique[has_signal == TRUE, mean(has_shared) * 100], 2), "%"))

# test de fisher sobre pares unicos
tabla_fisher <- table(Signal = dt_all_pairs_unique$has_signal, 
                      SharedGenes = dt_all_pairs_unique$has_shared)
print("Tabla de contingencia")
print(tabla_fisher)

res_fisher <- fisher.test(tabla_fisher)
print(res_fisher)

################################################################################
# resumen de métricas
################################################################################

# métricas del moodelo
n_triplets_evaluados <- nrow(all_triplets_results[model_success == TRUE])
n_positivos <- nrow(positives)
n_negativos <- nrow(negatives)

# métricas de la capa S
n_nodos_S   <- vcount(graph_stat_ddi)
n_aristas_S <- ecount(graph_stat_ddi)

# capa B
n_con_soporte_bio <- sum(E(graph_stat_ddi)$has_bio_support)
pct_soporte_bio <- (n_con_soporte_bio / max(1, n_aristas_S)) * 100
mean_jaccard_obs  <- mean(E(graph_stat_ddi)$jaccard)

# modelo nulo
mean_jaccard_null <- mean(dt_null_results$mean_jaccard)
sd_jaccard_null <- sd(dt_null_results$mean_jaccard)
z_score_jaccard <- ifelse(sd_jaccard_null > 0, (mean_jaccard_obs - mean_jaccard_null) / sd_jaccard_null, NA)
p_val_empirico  <- sum(dt_null_results$mean_jaccard >= mean_jaccard_obs) / nrow(dt_null_results)

# Fisher
or_fisher <- res_fisher$estimate
pval_fisher <- res_fisher$p.value

# Tabla resumen
dt_resumen_metricas <- data.table(
  Metrica = c(
    "Tripletes GAM ajustados exitosamente",
    "Señales Positivas Detectadas (DDI-Evento)",
    "Nodos en Red S (Fármacos)",
    "Nodos-Gen en la Red B (genes)",
    "Aristas en Red S (Pares DDI únicos)",
    "Pares DDI con Soporte Biológico (Genes compartidos)",
    "Porcentaje de la Red con Soporte Biológico (%)",
    "Índice Jaccard Medio (Observado)",
    "Índice Jaccard Medio (Modelo Nulo azar)",
    "Z-Score (Jaccard Observado vs Nulo)",
    "P-valor empírico (Soporte biológico de la red)",
    "Odds Ratio (Test Fisher: Señal vs Gen compartido)",
    "P-valor (Test Fisher)"
  ),
  Valor = c(
    as.character(n_triplets_evaluados),
    as.character(n_positivos),
    as.character(n_nodos_S),
    uniqueN(drug_gene$gene_symbol),
    as.character(n_aristas_S),
    as.character(n_con_soporte_bio),
    sprintf("%.2f %%", pct_soporte_bio),
    sprintf("%.4f", mean_jaccard_obs),
    sprintf("%.4f", mean_jaccard_null),
    sprintf("%.2f", z_score_jaccard),
    sprintf("%.5f", p_val_empirico),
    sprintf("%.2f", or_fisher),
    sprintf("%.2e", pval_fisher)
  )
)
# 
print(dt_resumen_metricas, row.names = FALSE, justify = "left")

fwrite(dt_resumen_metricas, paste0(output_dir, "00_resumen_metricas_globales.csv"))

################################################################################
# atributos
################################################################################

# nodos / vertices
# asigno nombre de fármacos a nodos mapeando desde drug_names_map
V(graph_stat_ddi)$drug_name <- drug_names_map[V(graph_stat_ddi)$name, drug_name]
# asigno atc group con mismo mapeo
V(graph_stat_ddi)$atc_group <- drug_names_map[V(graph_stat_ddi)$name, atc_group]
# asigno conectividad entre nodos
V(graph_stat_ddi)$degree_s  <- vec_degree_stat

# ordeno grupos atc
vec_atc_levels <- sort(unique(V(graph_stat_ddi)$atc_group))
# paleta de colores para atc
pal_atc <- setNames(
  colorRampPalette(brewer.pal(min(8, length(vec_atc_levels)), "Set2"))(length(vec_atc_levels)),
  vec_atc_levels
)
# asigno colores
V(graph_stat_ddi)$color_node <- pal_atc[V(graph_stat_ddi)$atc_group]

################################################################################
# visualización en circulo
################################################################################

# filtro red estadística para visualizar solo pares con soporte genético (Jaccard > 0)
vec_nodes_bio <- unique(c(dt_edge_metrics[jaccard > 0, drugA], dt_edge_metrics[jaccard > 0, drugB]))
graph_bio_support <- induced_subgraph(graph_stat_ddi, vids = vec_nodes_bio)

# selección de los 100 nodos principales por grado de conectividad para reducir densidad visual
vec_degree_bio <- degree(graph_bio_support)
vec_top100_bio <- names(sort(vec_degree_bio, decreasing = TRUE))[1:min(100, length(vec_degree_bio))]
graph_bio_support <- induced_subgraph(graph_bio_support, vids = vec_top100_bio)
vec_degree_bio <- degree(graph_bio_support)

# parametrización de etiquetas mostrando solo fármacos en el cuartil superior de conectividad
val_label_thr <- quantile(vec_degree_bio, 0.75)
vec_labels_bio <- ifelse(vec_degree_bio >= val_label_thr, V(graph_bio_support)$drug_name, NA_character_)
# escalado dinámico del tamaño de fuente según el grado 
vec_text_sizes <- rescale(pmax(vec_degree_bio - val_label_thr + 1, 0), to = c(2.8, 5.5))

set.seed(9427)
mat_coords <- layout_in_circle(graph_bio_support) * 0.60

plot_circle <- ggraph(as_tbl_graph(graph_bio_support), layout = "manual", x = mat_coords[, 1], y = mat_coords[, 2]) +
  geom_edge_arc(aes(alpha = jaccard, color = jaccard), strength = 0.15, width = 0.6) +
  geom_node_point(aes(size = vec_degree_bio, color = atc_group), alpha = 0.92) +
  geom_node_text(aes(label = vec_labels_bio, angle = -((-node_angle(x, y) + 90) %% 180) + 90, size = vec_text_sizes),
                 hjust = "outward", fontface = "bold", color = "#2C3E50", na.rm = TRUE) +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  scale_edge_alpha(range = c(0.10, 0.75), guide = "none") +
  scale_edge_color_gradient(low = "#AED6F1", high = "#C0392B", name = "Jaccard") +
  scale_color_manual(values = pal_atc, name = "ATC") +
  coord_cartesian(clip = "off") + theme_graph(base_family = "sans") + theme(legend.position = "right")

ggsave(paste0(output_dir, "fig_01_circle_bio.png"), plot_circle, width = 15, height = 15, dpi = 300)

################################################################################
# barras con hubs
################################################################################

dt_hub_drugs <- data.table(
  drug_name = V(graph_stat_ddi)$drug_name,
  degree = vec_degree_stat,
  atc_group = V(graph_stat_ddi)$atc_group
)[order(-degree)][1:min(30, .N)]

plot_hubs <- ggplot(dt_hub_drugs, aes(x = reorder(drug_name, degree), y = degree, fill = atc_group)) +
  geom_col(alpha = 0.85, width = 0.75) +
  geom_text(aes(label = degree), hjust = -0.3, size = 3, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pal_atc, name = "ATC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Top 30 Drugs by DDI Degree", x = NULL, y = "Degree") + theme(legend.position = "right")

ggsave(paste0(output_dir, "fig_04_hubs_bar.png"), plot_hubs, width = 12, height = 10, dpi = 300)

################################################################################
# visualización droga-gen
################################################################################

v_n_drugs <- 20
v_n_genes <- 10

vec_top_genes <- names(head(sort(degree(graph_bio_gene)[V(graph_bio_gene)$type == "gene"], decreasing = TRUE), v_n_genes))
vec_drugs_with_genes <- unique(dt_edges_bio$from)
vec_deg_drugs <- vec_degree_stat[names(vec_degree_stat) %in% vec_drugs_with_genes]
vec_top_drugs <- names(head(sort(vec_deg_drugs, decreasing = TRUE), v_n_drugs))

# Unir aristas Mixtas
dt_edges_dg <- dt_edges_bio[from %in% vec_top_drugs & to %in% vec_top_genes, .(from, to, type = "Drug-Gene")]
dt_edges_dd <- dt_edges_stat[from %in% vec_top_drugs & to %in% vec_top_drugs, .(from, to, type = "DDI")]
dt_edges_mix <- rbind(dt_edges_dg, dt_edges_dd)

# Grafo Mixto
graph_mixed <- graph_from_data_frame(dt_edges_mix, directed = FALSE)
V(graph_mixed)$type <- ifelse(V(graph_mixed)$name %in% vec_top_genes, "gene", "drug")
V(graph_mixed)$degree <- degree(graph_mixed)

vec_colors_mix <- ifelse(V(graph_mixed)$type == "gene", "#95A5A6", pal_atc[V(graph_mixed)$name])
vec_labels_mix <- ifelse(V(graph_mixed)$type == "gene", V(graph_mixed)$name, 
                         drug_names_map[match(V(graph_mixed)$name, atc_concept_id), drug_name])

###########
# función base para distintas visualizaciones
###########

plot_mixed_base <- function(layout_type) {
  set.seed(9427)
  ggraph(as_tbl_graph(graph_mixed), layout = layout_type) +
    geom_edge_link(aes(color = type, width = type, alpha = type)) +
    geom_node_point(aes(size = degree), color = vec_colors_mix, alpha = 0.9) +
    geom_node_text(aes(label = vec_labels_mix), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c("Drug-Gene" = "#E67E22", "DDI" = "#BDC3C7")) +
    scale_edge_width_manual(values = c("Drug-Gene" = 0.4, "DDI" = 0.4), guide = "none") +
    scale_edge_alpha_manual(values = c("Drug-Gene" = 0.35, "DDI" = 0.75), guide = "none") +
    scale_size(range = c(2, 10)) + theme_graph()
}

# kamada kawai
plot_mix_kk <- plot_mixed_base("kk") 
ggsave(paste0(output_dir, "fig_07a_mixed_kk.png"), plot_mix_kk, width = 16, height = 12, dpi = 300)

# fruchterman reingold
plot_mix_fr <- plot_mixed_base("fr") 
ggsave(paste0(output_dir, "fig_07b_mixed_fr.png"), plot_mix_fr, width = 16, height = 12, dpi = 300)

# Circular 
plot_mix_circle <- plot_mixed_base("circle")
ggsave(paste0(output_dir, "fig_07c_mixed_circle.png"), plot_mix_circle, width = 16, height = 12, dpi = 300)

################################################################################
# Sankey
################################################################################

g_S <- graph_stat_ddi
top_ids <- vec_top_drugs

positive_edges <- as.data.table(as_edgelist(g_S, names = TRUE))
names(positive_edges) <- c("drugA", "drugB")
positive_edges[, edge_id := .I]

gene_bridge_dt <- rbindlist(lapply(positive_edges$edge_id, function(idx) {
  dA <- positive_edges[edge_id == idx, drugA]
  dB <- positive_edges[edge_id == idx, drugB]
  
  genesA <- drug_gene[atc_concept_id == dA, gene_symbol]
  genesB <- drug_gene[atc_concept_id == dB, gene_symbol]
  
  common <- intersect(genesA, genesB)
  
  if (length(common) > 0) {
    data.table(drugA = dA, drugB = dB, gene_symbol = common)
  } else {
    NULL
  }
}))

# paquete para gráfico sankey
pacman::p_load(ggalluvial)

gene_counts <- gene_bridge_dt[, .(n_ddis = .N), by = gene_symbol][order(-n_ddis)]

# Filtramos las conexiones droga-gen que involucran a estas drogas top
sankey_links_dt <- rbind(
  gene_bridge_dt[drugA %in% top_ids, .(drug = drugA, gene = gene_symbol)],
  gene_bridge_dt[drugB %in% top_ids, .(drug = drugB, gene = gene_symbol)]
)

# pesos: cantidad de señales DDI que esta droga tiene mediadas por este gen
sankey_links_dt <- sankey_links_dt[, .(weight = .N), by = .(drug, gene)]
sankey_links_dt[, drug_name := drug_names_map[drug, drug_name]]

# Limpiar casos donde no hay gen (si las hay)
sankey_links_dt <- sankey_links_dt[!is.na(drug_name) & !is.na(gene)]

# Nodos del Sankey: drogas (izquierda) y genes (derecha)
node_labels <- unique(c(sankey_links_dt$drug_name, sankey_links_dt$gene))
node_df <- data.frame(name = node_labels)

# índices
sankey_links_dt$source <- match(sankey_links_dt$drug_name, node_labels) - 1
sankey_links_dt$target <- match(sankey_links_dt$gene, node_labels) - 1

# filtro top 15
top_n_static <- 15

# top 15 fármacos y genes basados en el volumen de conexiones (pesos)
top_drugs_static <- sankey_links_dt[, .(total = sum(weight)), by = drug_name][order(-total)][1:min(top_n_static, .N), drug_name]
top_genes_static <- sankey_links_dt[, .(total = sum(weight)), by = gene][order(-total)][1:min(top_n_static, .N), gene]

# Subconjunto de datos para el plot estático
static_sankey_dt <- sankey_links_dt[drug_name %in% top_drugs_static & gene %in% top_genes_static]

# sankey estático
p_static_sankey <- ggplot(static_sankey_dt, 
                          aes(axis1 = drug_name, axis2 = gene, y = weight)) +
  # Las cintas (flujos)
  geom_alluvium(aes(fill = drug_name), alpha = 0.65, curve_type = "cubic", width = 0.2) +
  # Las cajas (nodos)
  geom_stratum(fill = "#ECF0F1", color = "#2C3E50", width = 0.2, linewidth = 0.8) +
  # Los textos de los nodos
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.8, fontface = "bold", color = "#2C3E50") +
  # Escalas y ejes
  scale_x_discrete(limits = c("Fármacos (Top 15)", "Farmacogenes (Top 15)"), 
                   expand = c(0.15, 0.15)) +
  # Paleta de colores vibrante para los flujos
  scale_fill_viridis_d(option = "turbo", guide = "none") + 
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "#7F8C8D", margin = margin(b = 15))
  )

# p_static_sankey 
ggsave(paste0(output_dir, "fig_03b_sankey_static.png"), 
       p_static_sankey, width = 12, height = 10, dpi = 300)

################################################################################
# Comparación TWOSIDES
################################################################################

# cargo dataset twosides 
dt_twosides <- fread(ruta_twosides)
setnames(dt_twosides, tolower(names(dt_twosides)))
setnames(dt_twosides, "drug_1_rxnorn_id", "drug_1_rxnorm_id", skip_absent = TRUE)

dt_twosides[, rxnorm1 := as.character(drug_1_rxnorm_id)]
dt_twosides[, rxnorm2 := as.character(drug_2_rxnorm_id)]

dt_rel <- fread(ruta_rel, quote = "", colClasses = list(character = c("concept_id_1", "concept_id_2")))
dt_rel[, relationship_id := tolower(relationship_id)]

# extraigo ingredientes RxNorm
dt_rxnorm <- concept[vocabulary_id == "rxnorm" & concept_class_id == "ingredient",
                     .(concept_id = as.character(concept_id), rxnorm_code = as.character(concept_code))]

# mapeo usando el vocabulario OMOP. de ATC a RxNorm ingredient
vec_rel_ids <- c("atc - rxnorm", "atc - rxnorm pr lat", "atc - rxnorm eq", "atc - rxnorm sec", "maps to", "subsumes")
dt_atc_to_rxnorm <- dt_rel[relationship_id %in% vec_rel_ids]

# Mapeo por OMOP vocabulary y coincidencia de nombres
dt_raw_map <- dt_atc_to_rxnorm[atc_concepts, on = .(concept_id_1 = concept_id), nomatch = 0][
  dt_rxnorm, on = .(concept_id_2 = concept_id), nomatch = 0,
  .(atc_concept_id = as.character(concept_id_1), rxnorm_id = rxnorm_code)
]

dt_map_omop <- unique(merge(dt_raw_map, translation_table[, .(atc_concept_id, canonical_id)], by = "atc_concept_id", all.x = TRUE)[
  !is.na(canonical_id), .(atc_concept_id = canonical_id, rxnorm_id)])

dt_rxnorm_names <- concept[vocabulary_id == "rxnorm" & concept_class_id == "ingredient",
                           .(rxnorm_id = as.character(concept_code), drug_name = tolower(concept_name))]
dt_map_name <- merge(drug_names_map, dt_rxnorm_names, by = "drug_name")[, .(atc_concept_id, rxnorm_id)]

dt_mapping <- unique(rbind(dt_map_omop, dt_map_name))

# Filtro interacciones de twosides con punto de corte significativo (PRR > 2)
dt_tw_sig <- dt_twosides[prr > 2 & !is.na(rxnorm1) & !is.na(rxnorm2)]
dt_tw_pairs <- unique(dt_tw_sig[, .(pair_key = paste(pmin(rxnorm1, rxnorm2), pmax(rxnorm1, rxnorm2), sep = "__"))])
dt_tw_pairs[, in_twosides := TRUE]
setkey(dt_tw_pairs, pair_key)
rm(dt_tw_sig, dt_twosides); gc()

################################################################################
# Comparación Univariada de Pares Mapeados
################################################################################

vec_drugs_pos <- unique(c(as.character(positives$drugA), as.character(positives$drugB)))
dt_map_sub <- dt_mapping[atc_concept_id %in% vec_drugs_pos]

dt_pos_strong <- unique(positives[, .(atc_A = as.character(pmin(drugA, drugB)), atc_B = as.character(pmax(drugA, drugB)))])
dt_map_A <- dt_map_sub[, .(atc_A = as.character(atc_concept_id), rxnormA = rxnorm_id)]
dt_map_B <- dt_map_sub[, .(atc_B = as.character(atc_concept_id), rxnormB = rxnorm_id)]

# Cruzar mapeos de cada fármaco y definir key
dt_pos_mapped <- merge(dt_pos_strong, dt_map_A, by = "atc_A", all.x = TRUE, allow.cartesian = TRUE)
dt_pos_mapped <- merge(dt_pos_mapped, dt_map_B, by = "atc_B", all.x = TRUE, allow.cartesian = TRUE)
dt_pos_mapped <- dt_pos_mapped[!is.na(rxnormA) & !is.na(rxnormB)]
dt_pos_mapped[, pair_key := paste(pmin(rxnormA, rxnormB), pmax(rxnormA, rxnormB), sep = "__")]
setkey(dt_pos_mapped, pair_key)

n_mapped_unique <- uniqueN(dt_pos_mapped[, .(atc_A, atc_B)])

# limito twosides a las interacciones que se cruzan
vec_my_rxnorm <- unique(c(dt_pos_mapped$rxnormA, dt_pos_mapped$rxnormB))
dt_tw_restricted <- dt_tw_pairs[, c("r1", "r2") := tstrsplit(pair_key, "__", fixed = TRUE)][r1 %in% vec_my_rxnorm & r2 %in% vec_my_rxnorm]
dt_tw_restricted[, c("r1", "r2") := NULL]

dt_rxnorm_pairs <- unique(dt_pos_mapped[, .(pair_key)])
dt_rxnorm_pairs[, in_my_model := TRUE]
setkey(dt_rxnorm_pairs, pair_key)

# Evaluación Comparativa
dt_comp <- merge(dt_rxnorm_pairs, dt_tw_restricted, by = "pair_key", all = TRUE)
dt_comp[, status := fcase(
  in_my_model == TRUE & in_twosides == TRUE, "Concordante",
  in_my_model == TRUE & is.na(in_twosides), "Nuevo",
  is.na(in_my_model) & in_twosides == TRUE, "Solo TWOSIDES"
)]

# clasificación de pares para matriz de confusión vs TWOSIDES
val_tp <- sum(dt_comp$status == "Concordante", na.rm = TRUE)
val_fp <- sum(dt_comp$status == "Nuevo", na.rm = TRUE)
val_fn <- sum(dt_comp$status == "Solo TWOSIDES", na.rm = TRUE)
# cálculo de indicadores estadísticos de la red
val_prec <- val_tp / max(val_tp + val_fp, 1)
val_rec <- val_tp / max(val_tp + val_fn, 1)
val_f1 <- ifelse(val_prec + val_rec > 0, 2 * val_prec * val_rec / (val_prec + val_rec), 0)

# Exportar tabla unificada detallada
dt_comp_details <- merge(dt_comp, unique(dt_pos_mapped[, .(pair_key, atc_A, atc_B)]), by = "pair_key", all.x = TRUE)
dt_comp_details <- merge(dt_comp_details, drug_names_map[, .(atc_A = as.character(atc_concept_id), name_A = drug_name)], by = "atc_A", all.x = TRUE)
dt_comp_details <- merge(dt_comp_details, drug_names_map[, .(atc_B = as.character(atc_concept_id), name_B = drug_name)], by = "atc_B", all.x = TRUE)

fwrite(dt_comp_details, paste0(output_dir, "twosides_comparison_positives.csv"))
fwrite(unique(dt_comp_details[status == "Nuevo", .(atc_A, atc_B, name_A, name_B)]), paste0(output_dir, "seniales_novedosas_pares.csv"))

################################################################################
# Resumen de métricas
################################################################################

dt_metrics_summary <- data.table(
  grupo = "positives_strong_all",
  n_mapeados = n_mapped_unique,
  tp = val_tp,
  fp = val_fp,
  fn = val_fn,
  precision = round(val_prec, 3),
  recall = round(val_rec,  3),
  f1 = round(val_f1,   3)
)

print(dt_metrics_summary)
fwrite(dt_metrics_summary, paste0(output_dir, "metrics_versus_twosides.csv"))
