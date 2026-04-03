################################################################################
# Análisis de redes para plausibilidad biológica
# Script 40_network
################################################################################

pacman::p_load("data.table", "tidyverse", "mgcv", "parallel", 
               "doParallel", "pbapply", "igraph", "ggraph", 
               "visNetwork", "htmlwidgets", "scales") 

set.seed(9427)
setwd("D:/Bioestadística/gam-farmacovigilancia")
source("00_functions.R", local = TRUE)
source("01_theme.R")
theme_set(base_theme())

################################################################################
# Configuración
################################################################################

min_reports_triplet <- 5
min_stages_positive <- 1
n_random_networks <- 500
jaccard_threshold <- 0.1

spline_individuales <- TRUE 
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7    
include_nichd = FALSE
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
output_dir <- paste0("./results/", suffix, "/network_multiplex/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, floor(detectCores() * 0.75))

################################################################################
# Unificación de IDs canónicos
################################################################################

drug_info <- fread(ruta_drug_info)
drug_info[, `:=`(
  atc_concept_id = as.character(atc_concept_id),
  base_name = tolower(trimws(sub("[;,].*", "", atc_concept_name)))
)]

# Mapeo de IDs canónicos
canonical_map <- drug_info[, .(canonical_id = min(atc_concept_id)), by = base_name]
translation_table <- merge(drug_info[, .(atc_concept_id, base_name)], 
                           canonical_map, by = "base_name")

# Aplicar a dataset principal
ade_raw_dt <- fread(ruta_ade_raw)
ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]
ade_raw_dt <- merge(ade_raw_dt, 
                    translation_table[, .(atc_concept_id, canonical_id)], 
                    by = "atc_concept_id", all.x = TRUE)
ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

if (include_sex) {
  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
}

# Aplicar a drug_gene
drug_gene <- fread(ruta_drug_gene)
drug_gene[, atc_concept_id := as.character(atc_concept_id)]
drug_gene <- merge(drug_gene, translation_table[, .(atc_concept_id, canonical_id)],
                   by = "atc_concept_id", all.x = TRUE)
drug_gene[!is.na(canonical_id), atc_concept_id := canonical_id]
drug_gene[, canonical_id := NULL]
drug_gene <- unique(drug_gene, by = c("atc_concept_id", "gene_symbol"))

# Clasificación ATC (nivel 1 = grupo terapéutico)
if (file.exists(ruta_drug_info)) {
  atc_class <- fread(ruta_drug_info)
  atc_class[, atc_concept_id := as.character(atc_concept_id)]
  atc_class <- merge(atc_class, translation_table[, .(atc_concept_id, canonical_id)],
                     by = "atc_concept_id", all.x = TRUE)
  atc_class[!is.na(canonical_id), atc_concept_id := canonical_id]
  atc_class[, canonical_id := NULL]
  atc_class <- unique(atc_class[, .(atc_concept_id, atc1_concept_name)])
} else {
  atc_class <- unique(drug_info[, .(
    atc_concept_id = canonical_id,
    atc1_concept_name = substr(atc_concept_id, 1, 1)
  )])
}

# Diccionario de nombres para visualización
drug_names_map <- unique(translation_table[, .(
  atc_concept_id = canonical_id, 
  drug_name = base_name
)])
setkey(drug_names_map, atc_concept_id)

message(sprintf("  Drogas únicas: %d", uniqueN(drug_info$canonical_id)))
message(sprintf("  Asociaciones droga-gen: %d", nrow(drug_gene)))
message(sprintf("  Genes únicos: %d", uniqueN(drug_gene$gene_symbol)))

################################################################################
# Construcción de tripletes y ajuste de modelos
################################################################################

reports_meta <- ade_raw_dt[, .(
  drugs = list(unique(atc_concept_id[!is.na(atc_concept_id)])),
  events = list(unique(meddra_concept_id[!is.na(meddra_concept_id)])),
  nichd_num = unique(nichd_num[!is.na(nichd_num)])[1]
), by = safetyreportid]

triplets_list <- pblapply(seq_len(nrow(reports_meta)), function(i) {
  drug_list <- reports_meta$drugs[[i]]
  event_list <- reports_meta$events[[i]]
  stage <- reports_meta$nichd_num[i]
  
  if(length(drug_list) >= 2 && length(event_list) >= 1) {
    make_triplets(drug = drug_list, event = event_list, report_id = reports_meta$safetyreportid[i], nichd_stage = stage)
  }
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE, fill = TRUE)
rm(triplets_list, reports_meta); gc()

trip_counts <- triplets_dt[, .N, by = .(drugA, drugB, meddra)]
candidatos_triplets <- trip_counts[N >= min_reports_triplet]
candidatos_triplets[, triplet_id := .I]

message(sprintf("  Tripletes candidatos: %d", nrow(candidatos_triplets)))

fwrite(candidatos_triplets, paste0(output_dir, "triplets_metadata.csv"))

################################################################################
# Ajuste de modelos GAM (con caché)
################################################################################

cache_file <- paste0(output_dir, "all_triplets_results.rds")

if (file.exists(cache_file)) {
  all_triplets_results <- readRDS(cache_file)
} else {  
  batch_size <- 500
  n_batches <- ceiling(nrow(candidatos_triplets) / batch_size)
  all_results <- list()
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(candidatos_triplets))
    batch_indices <- start_idx:end_idx
    
    message(sprintf("  Procesando lote %d/%d", batch, n_batches))
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    clusterExport(cl, c("fit_differential_gam", "ade_raw_dt", "candidatos_triplets", 
                        "spline_individuales", "include_sex", "include_stage_sex", 
                        "k_spline", "nichd_spline", "bs_type", "select"), 
                  envir = environment())
    clusterEvalQ(cl, {library(data.table); library(mgcv)})
    
    batch_results <- foreach(idx = batch_indices, 
                             .packages = c("data.table", "mgcv"), 
                             .errorhandling = "pass") %dopar% {
      rowt <- candidatos_triplets[idx]
      tryCatch({
        fit_differential_gam(drugA_id = rowt$drugA, drugB_id = rowt$drugB,
                            event_id = rowt$meddra, ade_data = ade_raw_dt,
                            spline_individuales = spline_individuales,
                            include_sex = include_sex,
                            include_stage_sex = include_stage_sex,
                            k_spline = k_spline, nichd_spline = nichd_spline,
                            bs_type = bs_type, select = select)
      }, error = function(e) list(success = FALSE))
    }
    stopCluster(cl)
    
    # Procesar resultados del lote
    batch_dt_list <- lapply(seq_along(batch_results), function(i) {
      res <- batch_results[[i]]
      rowt <- candidatos_triplets[batch_indices[i]]
      base_res <- data.table(triplet_id = rowt$triplet_id, drugA = rowt$drugA, 
                             drugB = rowt$drugB, meddra = rowt$meddra, n_reports = rowt$N)
      
      if (is.null(res$success) || !res$success) {
        return(cbind(base_res, model_success = FALSE))
      } else {
        cols_metrics <- data.table(
          model_success = TRUE, n_events = res$n_events, n_coadmin = res$n_coadmin,
          n_stages_significant = res$n_stages_significant, model_aic = res$model_aic,
          stage = list(1:7), log_ior = list(res$log_ior), 
          log_ior_lower90 = list(res$log_ior_lower90)
        )
        return(cbind(base_res, cols_metrics))
      }
    })
    
    batch_dt <- rbindlist(batch_dt_list, fill = TRUE)
    all_results[[batch]] <- batch_dt
    
    # Checkpoint intermedio
    saveRDS(batch_dt, paste0(output_dir, "checkpoint_batch_", sprintf("%04d", batch), ".rds"))
    rm(batch_results, batch_dt, batch_dt_list); gc()
  }
  
  all_triplets_results <- rbindlist(all_results, fill = TRUE)
  saveRDS(all_triplets_results, cache_file)
}

message(sprintf("  Modelos exitosos: %d/%d (%.1f%%)",
                sum(all_triplets_results$model_success, na.rm = TRUE),
                nrow(all_triplets_results),
                100 * mean(all_triplets_results$model_success, na.rm = TRUE)))

################################################################################
# Clasificación de señales positivas y negativas
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CLASIFICACIÓN DE SEÑALES")
message(paste(rep("=", 80), collapse = ""))

# Positivos: Al menos min_stages_positive etapas con señal
positives <- all_triplets_results[model_success == TRUE, {
  lowers <- unlist(log_ior_lower90)
  n_pos <- if(length(lowers) == 0) 0 else sum(lowers > 0, na.rm=TRUE)
  list(n_stages_positive = n_pos)
}, by = .(triplet_id, drugA, drugB, meddra)]

positives <- positives[n_stages_positive >= min_stages_positive]

# Negativos: Resto de tripletes exitosos sin señal
negatives <- all_triplets_results[
  model_success == TRUE & 
  !(triplet_id %in% positives$triplet_id),
  .(triplet_id, drugA, drugB, meddra)
]

message(sprintf("  Tripletes positivos (señal detectada): %d", nrow(positives)))
message(sprintf("  Tripletes negativos (sin señal): %d", nrow(negatives)))

################################################################################
# CAPA S (Statistical - Farmacovigilancia)
################################################################################

# 1- Aristas de capa S: Solo DDIs con señal positiva
stat_edges <- unique(positives[, .(
  from = pmin(drugA, drugB),
  to = pmax(drugA, drugB),
  has_signal = TRUE
)])

# 2- Nodos de capa S: Drogas involucradas
drugs_in_S <- unique(c(stat_edges$from, stat_edges$to))

# 3- Construcción del grafo GS
g_S <- graph_from_data_frame(stat_edges[, .(from, to)], directed = FALSE)
E(g_S)$has_signal <- TRUE

message(sprintf("  Nodos (drogas): %d", vcount(g_S)))
message(sprintf("  Aristas (DDIs positivas): %d", ecount(g_S)))

################################################################################
# CAPA B (Biological - Molecular)
################################################################################

# 1- Filtrar solo drogas presentes en capa S
bio_edges <- unique(drug_gene[atc_concept_id %in% drugs_in_S,
                               .(from = atc_concept_id, to = gene_symbol)])

# 2- Nodos de capa B: Drogas + Genes
all_nodes_B <- data.table(
  name = c(drugs_in_S, unique(bio_edges$to)),
  type = c(rep("drug", length(drugs_in_S)), 
           rep("gene", uniqueN(bio_edges$to)))
)

# 3- Construcción del grafo GB (bipartito)
g_B <- graph_from_data_frame(bio_edges, directed = FALSE, vertices = all_nodes_B)

message(sprintf("  Nodos totales: %d (drogas: %d, genes: %d)",
                vcount(g_B), sum(V(g_B)$type == "drug"), sum(V(g_B)$type == "gene")))
message(sprintf("  Aristas droga-gen: %d", ecount(g_B)))

################################################################################
# MÉTRICAS INTRACAPA - Capa S
################################################################################
message("MÉTRICAS INTRACAPA (CAPA S)")

###########
# 1- Degree distribution
###########

degree_S <- degree(g_S)
degree_summary_S <- data.table(
  metric = c("Media", "Mediana", "Max", "Min"),
  value = c(mean(degree_S), median(degree_S), max(degree_S), min(degree_S))
)

message("\nDistribución de grado (Capa S):")
print(degree_summary_S)

# Identificar hubs (percentil 90)
hub_threshold <- quantile(degree_S, 0.90)
hubs_S <- names(degree_S)[degree_S >= hub_threshold]

message(sprintf("  Hubs identificados (>p90): %d", length(hubs_S)))

###########
# 2- Assortativity por clase ATC
###########

# Asignar clase ATC a nodos de GS
drug_atc_map <- merge(
  data.table(atc_concept_id = V(g_S)$name),
  atc_class,
  by = "atc_concept_id",
  all.x = TRUE
)

# Convertir a factor numérico para assortativity
atc_numeric <- as.numeric(as.factor(drug_atc_map$atc1_concept_name))
names(atc_numeric) <- drug_atc_map$atc_concept_id

# Calcular assortativity
assortativity_atc <- assortativity_nominal(
  g_S,
  types = atc_numeric[V(g_S)$name],
  directed = FALSE
)

message(sprintf("\nAssortativity por clase ATC (Capa S): %.4f", assortativity_atc))

# Guardar métricas intracapa
metrics_intralayer <- data.table(
  metric = c("n_nodes", "n_edges", "mean_degree", "median_degree", 
             "assortativity_atc", "n_hubs"),
  value = c(vcount(g_S), ecount(g_S), mean(degree_S), median(degree_S),
            assortativity_atc, length(hubs_S))
)

fwrite(metrics_intralayer, paste0(output_dir, "metrics_intralayer_S.csv"))

################################################################################
# MÉTRICAS INTERCAPA - S condicionada a B
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("MÉTRICAS INTERCAPA (S CONDICIONADA A B)")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Índice de Jaccard biológico
###########

# Para cada arista en GS, calcular Jaccard de genes compartidos
jaccard_results <- rbindlist(pblapply(seq_len(ecount(g_S)), function(i) {
  
  edge_ends <- ends(g_S, i)
  drugA <- edge_ends[1]
  drugB <- edge_ends[2]
  
  # Genes de cada droga en capa B
  genes_A <- bio_edges[from == drugA, to]
  genes_B <- bio_edges[from == drugB, to]
  
  # Jaccard
  n_shared <- length(intersect(genes_A, genes_B))
  n_union <- length(union(genes_A, genes_B))
  
  jaccard <- if (n_union > 0) n_shared / n_union else 0
  
  data.table(
    drugA = drugA,
    drugB = drugB,
    n_genes_A = length(genes_A),
    n_genes_B = length(genes_B),
    n_shared_genes = n_shared,
    jaccard_index = jaccard
  )
}))

# Agregar a aristas de GS
E(g_S)$jaccard <- jaccard_results$jaccard_index
E(g_S)$n_shared_genes <- jaccard_results$n_shared_genes

message(sprintf("  Media Jaccard: %.4f", mean(jaccard_results$jaccard_index)))
message(sprintf("  Mediana Jaccard: %.4f", median(jaccard_results$jaccard_index)))
message(sprintf("  DDIs con soporte biológico (Jaccard > 0): %d (%.1f%%)",
                sum(jaccard_results$jaccard_index > 0),
                100 * mean(jaccard_results$jaccard_index > 0)))

fwrite(jaccard_results, paste0(output_dir, "jaccard_biological_index.csv"))

###########
# 2- Shortest path biológico
###########
shortest_paths_results <- rbindlist(pblapply(seq_len(ecount(g_S)), function(i) {
  
  edge_ends <- ends(g_S, i)
  drugA <- edge_ends[1]
  drugB <- edge_ends[2]
  
  # Calcular distancia en GB
  path_length <- tryCatch({
    distances(g_B, v = drugA, to = drugB, mode = "all")[1,1]
  }, error = function(e) Inf)
  
  data.table(
    drugA = drugA,
    drugB = drugB,
    bio_distance = path_length
  )
}))

E(g_S)$bio_distance <- shortest_paths_results$bio_distance

message(sprintf("  Media distancia biológica: %.2f", 
                mean(shortest_paths_results$bio_distance[is.finite(shortest_paths_results$bio_distance)])))
message(sprintf("  DDIs con distancia = 2 (gen compartido): %d (%.1f%%)",
                sum(shortest_paths_results$bio_distance == 2),
                100 * mean(shortest_paths_results$bio_distance == 2)))
message(sprintf("  DDIs desconectadas (distancia = Inf): %d",
                sum(is.infinite(shortest_paths_results$bio_distance))))

fwrite(shortest_paths_results, paste0(output_dir, "biological_distances.csv"))

###########
# 3- Ratio de triadic closure
###########

# Para cada arista en GS, verificar si existe un gen que cierra el triángulo
triadic_closure_results <- rbindlist(pblapply(seq_len(ecount(g_S)), function(i) {
  
  edge_ends <- ends(g_S, i)
  drugA <- edge_ends[1]
  drugB <- edge_ends[2]
  
  # Genes conectados a ambas drogas en GB
  neighbors_A <- neighbors(g_B, drugA)
  genes_A <- V(g_B)[neighbors_A]$name[V(g_B)[neighbors_A]$type == "gene"]
  
  neighbors_B <- neighbors(g_B, drugB)
  genes_B <- V(g_B)[neighbors_B]$name[V(g_B)[neighbors_B]$type == "gene"]
  
  # ¿Existe al menos un gen compartido?
  closes_triangle <- length(intersect(genes_A, genes_B)) > 0
  
  data.table(
    drugA = drugA,
    drugB = drugB,
    closes_triangle = closes_triangle
  )
}))

E(g_S)$closes_triangle <- triadic_closure_results$closes_triangle

ratio_closure <- mean(triadic_closure_results$closes_triangle)

message(sprintf("  Ratio de triadic closure: %.4f", ratio_closure))
message(sprintf("  DDIs que cierran triángulo: %d/%d",
                sum(triadic_closure_results$closes_triangle),
                nrow(triadic_closure_results)))

fwrite(triadic_closure_results, paste0(output_dir, "triadic_closure.csv"))

###########
# Resumen de métricas intercapa
###########

metrics_interlayer <- data.table(
  metric = c("mean_jaccard", "median_jaccard", "prop_bio_support",
             "mean_bio_distance", "prop_distance_2", "prop_disconnected",
             "triadic_closure_ratio"),
  value = c(mean(jaccard_results$jaccard_index),
            median(jaccard_results$jaccard_index),
            mean(jaccard_results$jaccard_index > 0),
            mean(shortest_paths_results$bio_distance[is.finite(shortest_paths_results$bio_distance)]),
            mean(shortest_paths_results$bio_distance == 2),
            mean(is.infinite(shortest_paths_results$bio_distance)),
            ratio_closure)
)

fwrite(metrics_interlayer, paste0(output_dir, "metrics_interlayer.csv"))

################################################################################
# MODELO NULO - Degree-preserving rewiring (CON CACHÉ)
################################################################################


# Archivo de caché para modelo nulo
null_cache_file <- paste0(output_dir, "null_model_cache.rds")

# Métricas observadas
obs_mean_jaccard <- mean(jaccard_results$jaccard_index)
obs_median_jaccard <- median(jaccard_results$jaccard_index)
obs_prop_closure <- ratio_closure
obs_mean_bio_distance <- mean(shortest_paths_results$bio_distance[is.finite(shortest_paths_results$bio_distance)])

if (file.exists(null_cache_file)) {
  message(sprintf("  Cargando %d redes permutadas ", n_random_networks))
  null_results <- readRDS(null_cache_file)
  
} else {
  message(sprintf(" Permutación de %d redes", n_random_networks))
  
  # Null distributions
  null_results <- data.table(
    mean_jaccard = numeric(n_random_networks),
    median_jaccard = numeric(n_random_networks),
    prop_closure = numeric(n_random_networks),
    mean_bio_distance = numeric(n_random_networks)
  )
  
  pb <- txtProgressBar(max = n_random_networks, style = 3)
  
  for (i in 1:n_random_networks) {
    
    # Rewiring de GS preservando grados (keeping_degseq)
    g_S_perm <- rewire(g_S, keeping_degseq(niter = ecount(g_S) * 10))
    
    # Recalcular métricas en red permutada
    perm_jaccard <- sapply(seq_len(ecount(g_S_perm)), function(j) {
      edge_ends <- ends(g_S_perm, j)
      drugA <- edge_ends[1]
      drugB <- edge_ends[2]
      
      genes_A <- bio_edges[from == drugA, to]
      genes_B <- bio_edges[from == drugB, to]
      
      n_shared <- length(intersect(genes_A, genes_B))
      n_union <- length(union(genes_A, genes_B))
      
      if (n_union > 0) n_shared / n_union else 0
    })
    
    perm_closure <- sapply(seq_len(ecount(g_S_perm)), function(j) {
      edge_ends <- ends(g_S_perm, j)
      drugA <- edge_ends[1]
      drugB <- edge_ends[2]
      
      neighbors_A <- neighbors(g_B, drugA)
      genes_A <- V(g_B)[neighbors_A]$name[V(g_B)[neighbors_A]$type == "gene"]
      
      neighbors_B <- neighbors(g_B, drugB)
      genes_B <- V(g_B)[neighbors_B]$name[V(g_B)[neighbors_B]$type == "gene"]
      
      length(intersect(genes_A, genes_B)) > 0
    })
    
    perm_distances <- sapply(seq_len(ecount(g_S_perm)), function(j) {
      edge_ends <- ends(g_S_perm, j)
      drugA <- edge_ends[1]
      drugB <- edge_ends[2]
      
      tryCatch({
        distances(g_B, v = drugA, to = drugB, mode = "all")[1,1]
      }, error = function(e) Inf)
    })
    
    # Guardar métricas del null
    null_results$mean_jaccard[i] <- mean(perm_jaccard)
    null_results$median_jaccard[i] <- median(perm_jaccard)
    null_results$prop_closure[i] <- mean(perm_closure)
    null_results$mean_bio_distance[i] <- mean(perm_distances[is.finite(perm_distances)])
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Guardar caché
  saveRDS(null_results, null_cache_file)
}

# Calcular p-valores (one-tailed, observado > nulo)
p_mean_jaccard <- mean(null_results$mean_jaccard >= obs_mean_jaccard)
p_median_jaccard <- mean(null_results$median_jaccard >= obs_median_jaccard)
p_closure <- mean(null_results$prop_closure >= obs_prop_closure)
p_bio_distance <- mean(null_results$mean_bio_distance <= obs_mean_bio_distance)

message("\n", paste(rep("=", 80), collapse = ""))
message("RESULTADOS DEL MODELO NULO")
message(paste(rep("=", 80), collapse = ""))

cat(sprintf("
Jaccard medio:
  Observado: %.4f
  Null: %.4f ± %.4f
  p-valor: %.4f

Jaccard mediano:
  Observado: %.4f
  Null: %.4f ± %.4f
  p-valor: %.4f

Triadic closure:
  Observado: %.4f
  Null: %.4f ± %.4f
  p-valor: %.4f

Distancia biológica media:
  Observado: %.2f
  Null: %.2f ± %.2f
  p-valor: %.4f
",
obs_mean_jaccard, mean(null_results$mean_jaccard), sd(null_results$mean_jaccard), p_mean_jaccard,
obs_median_jaccard, mean(null_results$median_jaccard), sd(null_results$median_jaccard), p_median_jaccard,
obs_prop_closure, mean(null_results$prop_closure), sd(null_results$prop_closure), p_closure,
obs_mean_bio_distance, mean(null_results$mean_bio_distance, na.rm = TRUE), 
sd(null_results$mean_bio_distance, na.rm = TRUE), p_bio_distance
))

# Guardar distribuciones nulas
fwrite(null_results, paste0(output_dir, "null_distributions.csv"))

null_summary <- data.table(
  metric = c("mean_jaccard", "median_jaccard", "prop_closure", "mean_bio_distance"),
  observed = c(obs_mean_jaccard, obs_median_jaccard, obs_prop_closure, obs_mean_bio_distance),
  null_mean = c(mean(null_results$mean_jaccard), mean(null_results$median_jaccard),
                mean(null_results$prop_closure), mean(null_results$mean_bio_distance, na.rm = TRUE)),
  null_sd = c(sd(null_results$mean_jaccard), sd(null_results$median_jaccard),
              sd(null_results$prop_closure), sd(null_results$mean_bio_distance, na.rm = TRUE)),
  p_value = c(p_mean_jaccard, p_median_jaccard, p_closure, p_bio_distance)
)

fwrite(null_summary, paste0(output_dir, "null_model_summary.csv"))

################################################################################
# RED DE VISUALIZACIÓN - Proyección unimodal enriquecida
################################################################################

###########
# 1- Proyección: Solo drogas, genes como atributos
###########

# Usar aristas de GS con atributos biológicos ya calculados
vis_edges <- data.table(
  from = ends(g_S, E(g_S))[,1],
  to = ends(g_S, E(g_S))[,2],
  jaccard = E(g_S)$jaccard,
  n_shared_genes = E(g_S)$n_shared_genes,
  bio_distance = E(g_S)$bio_distance,
  closes_triangle = E(g_S)$closes_triangle
)

# Clasificar aristas: con/sin soporte biológico
vis_edges[, has_bio_support := jaccard > 0]
vis_edges[, orphan := jaccard == 0]

# Preparar nodos
vis_nodes <- data.table(
  id = V(g_S)$name,
  degree = degree(g_S)
)

# Añadir información de clasificación ATC y nombres
vis_nodes <- merge(vis_nodes, drug_names_map, 
                   by.x = "id", by.y = "atc_concept_id", all.x = TRUE)
vis_nodes <- merge(vis_nodes, atc_class,
                   by.x = "id", by.y = "atc_concept_id", all.x = TRUE)

# Identificar hubs para etiquetado
vis_nodes[, is_hub := degree >= hub_threshold]

message(sprintf("  Nodos (drogas): %d", nrow(vis_nodes)))
message(sprintf("  Aristas: %d", nrow(vis_edges)))
message(sprintf("  Aristas con soporte biológico: %d (%.1f%%)",
                sum(vis_edges$has_bio_support),
                100 * mean(vis_edges$has_bio_support)))
message(sprintf("  Aristas huérfanas (sin genes compartidos): %d (%.1f%%)",
                sum(vis_edges$orphan),
                100 * mean(vis_edges$orphan)))

###########
# 2- Construcción del grafo de visualización
###########
vis_nodes <- unique(vis_nodes, by = names(vis_nodes)[1])
g_vis <- graph_from_data_frame(vis_edges[, .(from, to)], 
                                directed = FALSE, 
                                vertices = vis_nodes)

# Transferir atributos de aristas
E(g_vis)$jaccard <- vis_edges$jaccard
E(g_vis)$n_shared_genes <- vis_edges$n_shared_genes
E(g_vis)$bio_distance <- vis_edges$bio_distance
E(g_vis)$has_bio_support <- vis_edges$has_bio_support
E(g_vis)$orphan <- vis_edges$orphan

###########
# 3- Codificación visual (Gramática de Gráficos de Red)
###########

# NODOS:
# - Color: Por grupo terapéutico ATC
atc_colors <- setNames(
  rainbow(uniqueN(vis_nodes$atc1_concept_name, na.rm = TRUE)),
  unique(vis_nodes$atc1_concept_name)
)
V(g_vis)$color <- atc_colors[V(g_vis)$atc1_concept_name]

# - Tamaño: Degree (importancia clínica)
V(g_vis)$size <- rescale(V(g_vis)$degree, to = c(3, 20))

# - Label: Nombre común (solo hubs)
V(g_vis)$label <- ifelse(V(g_vis)$is_hub, V(g_vis)$drug_name, NA)
V(g_vis)$label.cex <- 0.7
V(g_vis)$label.color <- "black"

# ARISTAS:
# - Grosor: Proporcional a Jaccard
E(g_vis)$width <- rescale(E(g_vis)$jaccard, to = c(0.5, 5))

# - Estilo: Sólida (con soporte) vs Punteada (huérfana)
E(g_vis)$lty <- ifelse(E(g_vis)$has_bio_support, 1, 3)

# - Color: Escala de grises según distancia biológica
# Normalizar distancia: más cercano = más oscuro
max_finite_dist <- max(E(g_vis)$bio_distance[is.finite(E(g_vis)$bio_distance)])
normalized_dist <- ifelse(
  is.finite(E(g_vis)$bio_distance),
  E(g_vis)$bio_distance / max_finite_dist,
  1  # Desconectadas = gris claro
)
E(g_vis)$color <- gray(normalized_dist)

###########
# 4- Layout optimizado
###########

set.seed(9427)
layout_vis <- layout_with_fr(g_vis, niter = 500)
layout_vis <- norm_coords(layout_vis, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

################################################################################
# VISUALIZACIONES
################################################################################

###########
# 1- Red de visualización principal
###########

viz_params <- list(
  # Layout
  layout_iterations = 100,
  layout_repulsion = 25.0,
  
  # Nodos
  node_size_range = c(2, 4),
  node_border_width = 0.5,
  node_border_color = "#2C3E50",
  node_alpha = 0.85,
  
  # Etiquetas
  label_size = 0.2,
  label_color = "#2C3E50",
  label_dist = 0.5,
  label_degree_only = TRUE,
  
  # Aristas
  edge_width_range = c(1, 5),
  edge_alpha_bio = 0.5,
  edge_alpha_orphan = 0.2,
  edge_curvature = 0.15,
  
  # Colores
  color_palette = "Set2",
  color_bio_support = "#27AE60",
  color_orphan = "#95A5A6",
  
  # Imagen
  img_width = 8000,
  img_height = 8000,
  img_res = 400,
  bg_color = "#FFFFFF"
)

###########
# 1- Layout optimizado con mayor spacing
###########

set.seed(9427)

# Fruchterman-Reingold con parámetros ajustados
layout_vis_improved <- layout_with_fr(
  g_vis, 
  niter = viz_params$layout_iterations,
  grid = "nogrid",
  area = vcount(g_vis)^2 * viz_params$layout_repulsion
)

# Normalizar coordenadas
layout_vis_improved <- norm_coords(
  layout_vis_improved, 
  ymin = -1, ymax = 1, 
  xmin = -1, xmax = 1
)

###########
# 2- Sistema de colores profesional
###########

library(RColorBrewer)

n_atc_groups <- uniqueN(V(g_vis)$atc1_concept_name, na.rm = TRUE)

# Seleccionar paleta según número de grupos
if (n_atc_groups <= 8) {
  atc_colors_improved <- brewer.pal(
    max(3, n_atc_groups), 
    viz_params$color_palette
  )[1:n_atc_groups]
} else {
  base_palette <- brewer.pal(8, viz_params$color_palette)
  atc_colors_improved <- colorRampPalette(base_palette)(n_atc_groups)
}

# Mapear colores a grupos ATC
atc_groups <- sort(unique(V(g_vis)$atc1_concept_name))
names(atc_colors_improved) <- atc_groups

# Aplicar colores con transparencia
V(g_vis)$color_improved <- alpha(
  atc_colors_improved[V(g_vis)$atc1_concept_name],
  viz_params$node_alpha
)

###########
# 3- Tamaños de nodos escalados
###########

V(g_vis)$size_improved <- rescale(
  V(g_vis)$degree,
  to = viz_params$node_size_range
)

###########
# 4- Etiquetas inteligentes
###########

if (viz_params$label_degree_only) {
  degree_threshold <- quantile(V(g_vis)$degree, 0.90)
  V(g_vis)$label_improved <- ifelse(
    V(g_vis)$degree >= degree_threshold,
    V(g_vis)$drug_name,
    NA
  )
} else {
  V(g_vis)$label_improved <- V(g_vis)$drug_name
}

###########
# 5- Sistema de colores de aristas mejorado
###########

E(g_vis)$color_improved <- ifelse(
  E(g_vis)$has_bio_support,
  alpha(viz_params$color_bio_support, viz_params$edge_alpha_bio),
  alpha(viz_params$color_orphan, viz_params$edge_alpha_orphan)
)

# Grosor escalado
E(g_vis)$width_improved <- rescale(
  E(g_vis)$jaccard,
  to = viz_params$edge_width_range
)

###########
# 6- Generación de gráfico principal
###########

png(
  paste0(output_dir, "fig_network_visualization.png"),
  width = viz_params$img_width,
  height = viz_params$img_height,
  res = viz_params$img_res,
  bg = viz_params$bg_color
)

par(mar = c(1, 1, 4, 1), bg = viz_params$bg_color)

plot(
  g_vis,
  layout = layout_vis_improved,
  
  # Nodos
  vertex.color = V(g_vis)$color_improved,
  vertex.size = V(g_vis)$size_improved,
  vertex.frame.color = viz_params$node_border_color,
  vertex.frame.width = viz_params$node_border_width,
  
  # Etiquetas
  vertex.label = V(g_vis)$label_improved,
  vertex.label.cex = viz_params$label_size,
  vertex.label.color = viz_params$label_color,
  vertex.label.dist = viz_params$label_dist,
  vertex.label.family = "sans",
  vertex.label.font = 2,
  
  # Aristas
  edge.width = E(g_vis)$width_improved,
  edge.lty = E(g_vis)$lty,
  edge.color = E(g_vis)$color_improved,
  edge.curved = viz_params$edge_curvature,
  edge.arrow.mode = 0,
  
  # Layout
  rescale = FALSE,
  asp = 0,
  xlim = c(-1.2, 1.2),
  ylim = c(-1.2, 1.2)
)

# Título mejorado
title(
  main = "Red de Interacciones Droga-Droga (DDIs)",
  sub = sprintf(
    "Proyección con Soporte Biológico | n=%d drogas, %d DDIs detectadas",
    vcount(g_vis), ecount(g_vis)
  ),
  cex.main = 2.0,
  cex.sub = 1.2,
  font.main = 2,
  col.main = "#2C3E50",
  col.sub = "#7F8C8D"
)

# Leyenda mejorada
legend(
  "bottomright",
  legend = c(
    "Droga (nodo)",
    "DDI con soporte biológico",
    "DDI huérfana (sin genes)",
    paste0("Tamaño ∝ Degree (", 
           min(V(g_vis)$degree), "-", 
           max(V(g_vis)$degree), ")"),
    paste0("Grosor ∝ Jaccard (0-", 
           sprintf("%.2f", max(E(g_vis)$jaccard)), ")")
  ),
  lty = c(NA, 1, 3, NA, NA),
  lwd = c(NA, 3, 2, NA, NA),
  pch = c(21, NA, NA, 21, NA),
  pt.cex = c(2, NA, NA, c(1, 2), NA),
  pt.bg = c(alpha("#3498DB", 0.7), NA, NA, NA, NA),
  col = c(
    viz_params$node_border_color,
    viz_params$color_bio_support,
    viz_params$color_orphan,
    viz_params$node_border_color,
    NA
  ),
  cex = 1.3,
  bty = "n",
  bg = alpha("white", 0.9),
  box.col = "#BDC3C7"
)

dev.off()
message(sprintf("    - Paleta: %s", viz_params$color_palette))
message(sprintf("    - Layout iterations: %d", viz_params$layout_iterations))
message(sprintf("    - Nodos etiquetados: %d", sum(!is.na(V(g_vis)$label_improved))))

###########
# 2- Distribución de Jaccard
###########

p_jaccard_dist <- ggplot(jaccard_results, aes(x = jaccard_index)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, fill = "#3498DB", alpha = 0.7, color = "black") +
  geom_density(color = "#E74C3C", linewidth = 1.5) +
  geom_vline(xintercept = obs_mean_jaccard, 
             linetype = "dashed", color = "#2C3E50", linewidth = 1) +
  annotate("text", x = obs_mean_jaccard + 0.05, y = Inf, vjust = 1.5,
           label = sprintf("Media: %.3f", obs_mean_jaccard),
           fontface = "bold", size = 4) +
  labs(
    title = "Distribución del Índice de Jaccard Biológico",
    subtitle = sprintf("n = %d pares de drogas con DDI detectada", nrow(jaccard_results)),
    x = "Índice de Jaccard (genes compartidos)",
    y = "Densidad"
  )

ggsave(paste0(output_dir, "fig_jaccard_distribution.png"),
       p_jaccard_dist, width = 10, height = 7, dpi = 300)

###########
# 3- Validación null model
###########

null_long <- rbindlist(list(
  data.table(metric = "Jaccard medio", 
             value = null_results$mean_jaccard,
             observed = obs_mean_jaccard),
  data.table(metric = "Triadic closure",
             value = null_results$prop_closure,
             observed = obs_prop_closure)
))

p_null_validation <- ggplot(null_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 40, fill = "#95A5A6", alpha = 0.7, color = "black") +
  geom_density(color = "#3498DB", linewidth = 1.2) +
  geom_vline(aes(xintercept = observed), 
             color = "#E74C3C", linewidth = 1.5, linetype = "dashed") +
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  labs(
    title = "Validación vs Modelo Nulo",
    subtitle = sprintf("n = %d permutaciones con degree-preserving rewiring", n_random_networks),
    x = "Valor",
    y = "Densidad"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))

ggsave(paste0(output_dir, "fig_null_model_validation.png"),
       p_null_validation, width = 10, height = 10, dpi = 300)

message("  ✓ Validación null model guardada")

###########
# 4- Top genes por participación en DDIs
###########

# Contar participación de genes en DDIs
gene_participation <- bio_edges[from %in% V(g_S)$name, .N, by = to]
setorder(gene_participation, -N)
setnames(gene_participation, c("to", "N"), c("gene_symbol", "n_ddis"))

p_gene_participation <- ggplot(gene_participation[1:min(30, .N)],
                                aes(x = reorder(gene_symbol, n_ddis), 
                                    y = n_ddis)) +
  geom_col(aes(fill = n_ddis), alpha = 0.8, color = "black") +
  scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
  coord_flip() +
  labs(
    title = "Top Genes por Participación en DDIs",
    subtitle = "Genes conectados a drogas involucradas en interacciones detectadas",
    x = NULL,
    y = "Número de DDIs",
    fill = NULL
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

ggsave(paste0(output_dir, "fig_top_genes_participation.png"),
       p_gene_participation, width = 10, height = 12, dpi = 300)

message("  ✓ Participación de genes guardada")

fwrite(gene_participation, paste0(output_dir, "gene_participation_ranking.csv"))

################################################################################
# ANÁLISIS COMPLEMENTARIOS
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("ANÁLISIS COMPLEMENTARIOS")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Enriquecimiento de genes compartidos (Fisher test)
###########

# Tabla de contingencia: Positivos vs Negativos × Con genes vs Sin genes
all_pairs <- rbind(
  positives[, .(drugA, drugB, has_signal = TRUE)],
  negatives[, .(drugA, drugB, has_signal = FALSE)]
)

# Añadir información de genes compartidos
all_pairs_genes <- merge(
  all_pairs,
  jaccard_results[, .(drugA = pmin(drugA, drugB),
                      drugB = pmax(drugA, drugB),
                      jaccard_index)],
  by = c("drugA", "drugB"),
  all.x = TRUE
)

all_pairs_genes[is.na(jaccard_index), jaccard_index := 0]
all_pairs_genes[, has_shared_genes := jaccard_index > 0]

# Tabla de contingencia
contingency_table <- table(
  Signal = all_pairs_genes$has_signal,
  SharedGenes = all_pairs_genes$has_shared_genes
)

fisher_enrichment <- fisher.test(contingency_table)

message("\nEnriquecimiento de genes compartidos:")
message(sprintf("  Fisher OR: %.3f", fisher_enrichment$estimate))
message(sprintf("  p-valor: %.4e", fisher_enrichment$p.value))

enrichment_summary <- as.data.table(contingency_table)
fwrite(enrichment_summary, paste0(output_dir, "enrichment_contingency.csv"))

# Resumen de proporciones
proportion_summary <- all_pairs_genes[, .(
  prop_shared_genes = mean(has_shared_genes),
  n_pairs = .N
), by = has_signal]

message("\nProporciones de genes compartidos:")
print(proportion_summary)

fwrite(proportion_summary, paste0(output_dir, "proportion_shared_genes.csv"))

################################################################################
# Gráfico de grafo global
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CONSTRUCCIÓN RED GLOBAL COLAPSADA")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Construcción de nodos globales
###########

# Nodos de drogas con información completa
global_drugs <- data.table(
  id = V(g_S)$name,
  label = drug_names_map[V(g_S)$name, drug_name],
  group = atc_class[match(V(g_S)$name, atc_concept_id), atc1_concept_name],
  degree_ddi = degree(g_S),  # Degree solo en capa DDI
  type = "drug"
)

# Calcular degree de genes en capa biológica
genes_in_network <- unique(bio_edges[from %in% V(g_S)$name, to])
gene_degree <- bio_edges[from %in% V(g_S)$name, .N, by = to]

# Nodos de genes
global_genes <- data.table(
  id = genes_in_network,
  label = genes_in_network,
  group = "Gene",
  degree_bio = gene_degree[match(genes_in_network, to), N],
  type = "gene"
)

# Combinar nodos
global_nodes <- rbindlist(list(
  global_drugs[, .(id, label, group, degree = degree_ddi, type)],
  global_genes[, .(id, label, group, degree = degree_bio, type)]
), fill = TRUE)

message(sprintf("  Nodos totales: %d (drogas: %d, genes: %d)",
                nrow(global_nodes),
                nrow(global_drugs),
                nrow(global_genes)))

###########
# 2- Construcción de aristas globales
###########

# Aristas DDI (droga-droga) con atributos
global_edges_ddi <- data.table(
  from = ends(g_S, E(g_S))[,1],
  to = ends(g_S, E(g_S))[,2],
  weight = E(g_S)$jaccard,
  edge_type = "DDI",
  has_bio_support = E(g_S)$has_bio_support
)

# Aristas biológicas (droga-gen)
global_edges_bio <- bio_edges[from %in% V(g_S)$name, .(
  from = from,
  to = to,
  weight = 1,
  edge_type = "biological",
  has_bio_support = TRUE
)]

# Combinar aristas
global_edges <- rbind(global_edges_ddi, global_edges_bio, fill = TRUE)

message(sprintf("  Aristas totales: %d (DDI: %d, droga-gen: %d)",
                nrow(global_edges),
                nrow(global_edges_ddi),
                nrow(global_edges_bio)))

###########
# 3- Construcción del grafo global
###########

g_global <- graph_from_data_frame(
  global_edges[, .(from, to)],
  directed = FALSE,
  vertices = global_nodes
)

# Transferir atributos de aristas
E(g_global)$weight <- global_edges$weight
E(g_global)$edge_type <- global_edges$edge_type
E(g_global)$has_bio_support <- global_edges$has_bio_support

message(sprintf("  Grafo global: %d nodos, %d aristas",
                vcount(g_global), ecount(g_global)))

###########
# 4- Codificación visual para red global
###########

# Sistema de colores para grupos ATC + genes
n_atc_groups <- uniqueN(global_drugs$group, na.rm = TRUE)

if (n_atc_groups <= 8) {
  atc_colors_global <- brewer.pal(max(3, n_atc_groups), "Set2")
} else {
  atc_colors_global <- colorRampPalette(brewer.pal(8, "Set2"))(n_atc_groups)
}

atc_groups_global <- sort(unique(global_drugs$group))
names(atc_colors_global) <- atc_groups_global

# Color para genes (gris claro)
node_colors_global <- ifelse(
  V(g_global)$type == "drug",
  atc_colors_global[V(g_global)$group],
  "#ECF0F1"
)
V(g_global)$color <- alpha(node_colors_global, 0.85)

# Forma diferente para genes vs drogas
V(g_global)$shape <- ifelse(V(g_global)$type == "drug", "circle", "square")

# Tamaño proporcional al degree
V(g_global)$size <- rescale(V(g_global)$degree, to = c(2, 8))

# Etiquetas solo para hubs
degree_threshold_global <- quantile(V(g_global)$degree, 0.95)
V(g_global)$label <- ifelse(
  V(g_global)$degree >= degree_threshold_global,
  V(g_global)$label,
  NA
)
V(g_global)$label.cex <- 0.6
V(g_global)$label.color <- "#2C3E50"

# Colores de aristas
E(g_global)$color <- ifelse(
  E(g_global)$edge_type == "DDI",
  ifelse(E(g_global)$has_bio_support, 
         alpha("#27AE60", 0.5),  # Verde para DDI con soporte
         alpha("#95A5A6", 0.3)), # Gris para DDI huérfana
  alpha("#3498DB", 0.2)          # Azul para droga-gen
)

# Grosor de aristas
E(g_global)$width <- ifelse(
  E(g_global)$edge_type == "DDI",
  rescale(E(g_global)$weight, to = c(0.5, 3)),
  0.3
)

###########
# 5- Layout optimizado para red global
###########

set.seed(9427)

# Fruchterman-Reingold con más iteraciones para red compleja
# CORRECCIÓN: layout_with_fr no acepta weights = 0
# Opción 1: Usar layout SIN pesos (recomendado para redes heterogéneas)
layout_global <- layout_with_fr(
  g_global,
  niter = 200,
  grid = "nogrid",
  weights = NA,  # <--- CAMBIO: Desactivar uso de pesos
  area = vcount(g_global)^2 * 50
)

layout_global <- norm_coords(
    layout_global,
      ymin = -1, ymax = 1,
        xmin = -1, xmax = 1
        )


###########
# 6- Generación de gráfico global
###########

png(
  paste0(output_dir, "fig_network_global_collapsed.png"),
  width = 10000,
  height = 10000,
  res = 400,
  bg = "#FFFFFF"
)

par(mar = c(1, 1, 4, 1), bg = "#FFFFFF")

plot(
  g_global,
  layout = layout_global,
  
  # Nodos
  vertex.color = V(g_global)$color,
  vertex.size = V(g_global)$size,
  vertex.shape = V(g_global)$shape,
  vertex.frame.color = "#2C3E50",
  vertex.frame.width = 0.5,
  
  # Etiquetas
  vertex.label = V(g_global)$label,
  vertex.label.cex = V(g_global)$label.cex,
  vertex.label.color = V(g_global)$label.color,
  vertex.label.dist = 0.5,
  vertex.label.family = "sans",
  vertex.label.font = 2,
  
  # Aristas
  edge.width = E(g_global)$width,
  edge.color = E(g_global)$color,
  edge.curved = 0.1,
  edge.arrow.mode = 0,
  
  # Layout
  rescale = FALSE,
  asp = 0,
  xlim = c(-1.2, 1.2),
  ylim = c(-1.2, 1.2)
)

# Título
title(
  main = "Red Global Multiplex: Drogas + Genes",
  sub = sprintf(
    "Capa DDI + Capa Biológica | %d drogas (círculos), %d genes (cuadrados)",
    nrow(global_drugs), nrow(global_genes)
  ),
  cex.main = 2.5,
  cex.sub = 1.5,
  font.main = 2,
  col.main = "#2C3E50",
  col.sub = "#7F8C8D"
)

# Leyenda
legend(
  "bottomright",
  legend = c(
    "Droga",
    "Gen",
    "DDI con soporte biológico",
    "DDI huérfana",
    "Conexión droga-gen"
  ),
  pch = c(21, 22, NA, NA, NA),
  lty = c(NA, NA, 1, 1, 1),
  lwd = c(NA, NA, 3, 2, 1),
  pt.cex = c(2, 2, NA, NA, NA),
  pt.bg = c(
    alpha("#3498DB", 0.7),
    alpha("#ECF0F1", 0.7),
    NA, NA, NA
  ),
  col = c(
    "#2C3E50",
    "#2C3E50",
    alpha("#27AE60", 0.8),
    alpha("#95A5A6", 0.6),
    alpha("#3498DB", 0.5)
  ),
  cex = 1.5,
  bty = "n",
  bg = alpha("white", 0.9),
  box.col = "#BDC3C7"
)

dev.off()

message("  ✓ Red global colapsada guardada")

# Guardar métricas de red global
global_metrics <- data.table(
  metric = c(
    "n_nodes_total", "n_drugs", "n_genes",
    "n_edges_total", "n_ddi_edges", "n_bio_edges",
    "density", "mean_degree",
    "transitivity", "diameter"
  ),
  value = c(
    vcount(g_global),
    nrow(global_drugs),
    nrow(global_genes),
    ecount(g_global),
    nrow(global_edges_ddi),
    nrow(global_edges_bio),
    edge_density(g_global),
    mean(degree(g_global)),
    transitivity(g_global, type = "global"),
    diameter(g_global, directed = FALSE)
  )
)

fwrite(global_metrics, paste0(output_dir, "metrics_global_network.csv"))

message(sprintf("    - Densidad: %.4f", global_metrics[metric == "density", value]))
message(sprintf("    - Degree medio: %.2f", global_metrics[metric == "mean_degree", value]))
message(sprintf("    - Transitividad: %.4f", global_metrics[metric == "transitivity", value]))

################################################################################
# VISUALIZACIÓN INTERACTIVA CON visNetwork (CORREGIDA)
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("CONSTRUCCIÓN DE RED INTERACTIVA")
message(paste(rep("=", 80), collapse = ""))

###########
# 1- Preparación de nodos para visNetwork
###########

# Nodos de drogas
vis_drugs <- data.table(
  id = V(g_S)$name,
  label = drug_names_map[V(g_S)$name, drug_name],
  group = atc_class[match(V(g_S)$name, atc_concept_id), atc1_concept_name],
  value = degree(g_S),
  type = "drug",
  shape = "dot",
  title = sprintf(
    "<b>%s</b><br>ID: %s<br>Grupo ATC: %s<br>DDIs detectadas: %d<br>Hub: %s",
    drug_names_map[V(g_S)$name, drug_name],
    V(g_S)$name,
    atc_class[match(V(g_S)$name, atc_concept_id), atc1_concept_name],
    degree(g_S),
    ifelse(degree(g_S) >= hub_threshold, "Sí", "No")
  )
)

# Nodos de genes
vis_genes <- data.table(
  id = genes_in_network,
  label = genes_in_network,
  group = "Gene",
  value = gene_degree[match(genes_in_network, to), N],
  type = "gene",
  shape = "triangle",
  title = sprintf(
    "<b>%s</b><br>Tipo: Gen<br>Conectado a: %d drogas",
    genes_in_network,
    gene_degree[match(genes_in_network, to), N]
  )
)

# Combinar nodos
vis_nodes_full <- rbind(vis_drugs, vis_genes, fill = TRUE)

# Escalar tamaños
vis_nodes_full[type == "drug", value := rescale(value, to = c(10, 40))]
vis_nodes_full[type == "gene", value := rescale(value, to = c(5, 20))]

###########
# 2- Preparación de aristas para visNetwork
###########

# Aristas DDI
vis_edges_ddi <- data.table(
  from = ends(g_S, E(g_S))[,1],
  to = ends(g_S, E(g_S))[,2],
  value = rescale(E(g_S)$jaccard, to = c(1, 8)),
  dashes = !as.logical(E(g_S)$has_bio_support),  # <--- CORRECCIÓN: Convertir explícitamente a lógico
  title = sprintf(
    "<b>DDI</b><br>Jaccard: %.3f<br>Genes compartidos: %d<br>Distancia biológica: %.1f<br>Cierra triángulo: %s",
    E(g_S)$jaccard,
    E(g_S)$n_shared_genes,
    ifelse(is.finite(E(g_S)$bio_distance), E(g_S)$bio_distance, Inf),
    ifelse(E(g_S)$closes_triangle, "Sí", "No")
  ),
  layer = "DDI"
)

vis_edges_bio <- bio_edges[from %in% V(g_S)$name, .(
  from = from,
  to = to,
  value = 1,
  dashes = FALSE,
  title = sprintf("<b>Droga-Gen</b><br>%s → %s", from, to),
  layer = "biological"
)]

# Combinar aristas
vis_edges_full <- rbind(
  vis_edges_ddi[, .(from, to, value, dashes, title, layer)],
  vis_edges_bio[, .(from, to, value, dashes, title, layer)],
  fill = TRUE
)

###########
# 3- Construcción de red interactiva con visNetwork
###########

colores_mates <- c(
  "#B8C5D6",  # Azul grisáceo
  "#D4C5B8",  # Beige
  "#C5D6C5",  # Verde grisáceo
  "#D6C5D4",  # Rosa grisáceo
  "#D6D4C5",  # Amarillo grisáceo
  "#C5D6D4",  # Cyan grisáceo
  "#D6C5C5",  # Rojo grisáceo
  "#C5C5D6"   # Violeta grisáceo
)

# Extender paleta si hay más grupos
if (length(atc_groups_global) > length(colores_mates)) {
  colores_mates <- rep(colores_mates, length.out = length(atc_groups_global))
}

# Asignar colores mates a grupos
names(colores_mates) <- atc_groups_global[1:length(colores_mates)]

# Aplicar colores a nodos
vis_nodes_full_matte <- copy(vis_nodes_full)
vis_nodes_full_matte[type == "drug", color := colores_mates[group]]
vis_nodes_full_matte[type == "gene", color := "#D5DBDB"]  # Gris claro para genes

# Crear visNetwork con mayor repulsión
vis_net <- visNetwork(
  nodes = vis_nodes_full_matte[, .(id, label, group, value, shape, title, color)],
  edges = vis_edges_full[, .(from, to, value, title)],
  width = "100%",
  height = "900px"
) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(
      gravitationalConstant = -150,      # Mucha más repulsión
      centralGravity = 0.005,            # Menos gravedad central
      springLength = 200,                # Resortes más largos
      springConstant = 0.04,             # Resortes más débiles
      damping = 0.5,
      avoidOverlap = 1.0                 # Máximo evitar solapamiento
    ),
    stabilization = list(
      enabled = TRUE,
      iterations = 2000                  # Más iteraciones para estabilizar
    )
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    hover = TRUE,
    zoomView = TRUE,
    dragView = TRUE
  ) %>%
  visLegend()


# Guardar red interactiva
visSave(
  vis_net,
  file = paste0(output_dir, "interactive_network.html"),
  selfcontained = TRUE
)

message("  ✓ Red interactiva guardada: interactive_network.html")
message(sprintf("    - Nodos totales: %d (drogas: %d, genes: %d)",
                nrow(vis_nodes_full),
                nrow(vis_drugs),
                nrow(vis_genes)))
message(sprintf("    - Aristas totales: %d (DDIs: %d, droga-gen: %d)",
                nrow(vis_edges_full),
                nrow(vis_edges_ddi),
                nrow(vis_edges_bio)))

################################################################################
# RESUMEN EJECUTIVO
################################################################################

message("\n", paste(rep("=", 80), collapse = ""))
message("GENERANDO RESUMEN EJECUTIVO")
message(paste(rep("=", 80), collapse = ""))

summary_text <- sprintf("
================================================================================
ANÁLISIS DE PLAUSIBILIDAD BIOLÓGICA - ARQUITECTURA MULTIPLEX
CONFIGURACIÓN:
Modelo: %s
Drogas únicas: %d
Genes únicos: %d
Tripletes evaluados: %d
CAPA S (STATISTICAL - FARMACOVIGILANCIA):
Nodos (drogas): %d
Aristas (DDIs positivas): %d
Degree medio: %.2f
Assortativity ATC: %.4f
CAPA B (BIOLOGICAL - MOLECULAR):
Nodos totales: %d (drogas: %d, genes: %d)
Aristas droga-gen: %d
MÉTRICAS INTERCAPA (S CONDICIONADA A B):
Índice de Jaccard biológico:
  Media: %.4f
  Mediana: %.4f
  DDIs con soporte (Jaccard > 0): %.1f%%
Shortest path biológico:
  Distancia media: %.2f
  DDIs con gen compartido (dist=2): %.1f%%
  DDIs desconectadas: %.1f%%
Triadic closure:
  Ratio: %.4f
  DDIs que cierran triángulo: %d/%d (%.1f%%)
VALIDACIÓN VS MODELO NULO:
Jaccard medio:
  Observado: %.4f
  Null: %.4f ± %.4f
  p-valor: %.4f %s
Triadic closure:
  Observado: %.4f
  Null: %.4f ± %.4f
  p-valor: %.4f %s
Distancia biológica:
  Observado: %.2f
  Null: %.2f ± %.2f
  p-valor: %.4f %s
ENRIQUECIMIENTO DE GENES COMPARTIDOS:
Fisher OR: %.3f
p-valor: %.4e %s
Señales con genes compartidos: %.1f%%
No-señales con genes compartidos: %.1f%%
RED GLOBAL COLAPSADA:
Nodos totales: %d (drogas: %d, genes: %d)
Aristas totales: %d (DDI: %d, droga-gen: %d)
Densidad: %.4f
Transitividad: %.4f
Fecha: %s
",
suffix,
uniqueN(drug_info$canonical_id),
uniqueN(drug_gene$gene_symbol),
nrow(candidatos_triplets),
vcount(g_S),
ecount(g_S),
mean(degree_S),
assortativity_atc,
vcount(g_B),
sum(V(g_B)$type == "drug"),
sum(V(g_B)$type == "gene"),
ecount(g_B),
obs_mean_jaccard,
obs_median_jaccard,
100 * mean(jaccard_results$jaccard_index > 0),
obs_mean_bio_distance,
100 * mean(shortest_paths_results$bio_distance == 2),
100 * mean(is.infinite(shortest_paths_results$bio_distance)),
obs_prop_closure,
sum(triadic_closure_results$closes_triangle),
nrow(triadic_closure_results),
100 * obs_prop_closure,
obs_mean_jaccard,
mean(null_results$mean_jaccard),
sd(null_results$mean_jaccard),
p_mean_jaccard,
ifelse(p_mean_jaccard < 0.05, "***", "ns"),
obs_prop_closure,
mean(null_results$prop_closure),
sd(null_results$prop_closure),
p_closure,
ifelse(p_closure < 0.05, "***", "ns"),
obs_mean_bio_distance,
mean(null_results$mean_bio_distance, na.rm = TRUE),
sd(null_results$mean_bio_distance, na.rm = TRUE),
p_bio_distance,
ifelse(p_bio_distance < 0.05, "***", "ns"),
fisher_enrichment$estimate,
fisher_enrichment$p.value,
ifelse(fisher_enrichment$p.value < 0.05, "***", "ns"),
100 * proportion_summary[has_signal == TRUE, prop_shared_genes],
100 * proportion_summary[has_signal == FALSE, prop_shared_genes],
vcount(g_global),
nrow(global_drugs),
nrow(global_genes),
ecount(g_global),
nrow(global_edges_ddi),
nrow(global_edges_bio),
global_metrics[metric == "density", value],
global_metrics[metric == "transitivity", value],
format(Sys.time(), "%Y-%m-%d %H:%M")
)

cat(summary_text)
writeLines(summary_text, paste0(output_dir, "resumen_multiplex_analysis.txt"))

################################################################################
# GUARDADO DE OBJETOS FINALES
################################################################################
saveRDS(list(
  g_S = g_S,
  g_B = g_B,
  g_vis = g_vis,
  g_global = g_global,
  jaccard_results = jaccard_results,
  shortest_paths_results = shortest_paths_results,
  triadic_closure_results = triadic_closure_results,
  null_results = null_results,
  metrics_intralayer = metrics_intralayer,
  metrics_interlayer = metrics_interlayer,
  global_metrics = global_metrics
), 
paste0(output_dir, "network_multiplex_objects.rds"))

message("\n", paste(rep("=", 80), collapse = ""))
message("ANÁLISIS COMPLETADO EXITOSAMENTE")
message(paste(rep("=", 80), collapse = ""))
