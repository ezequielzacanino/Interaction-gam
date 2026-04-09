################################################################################
# Semi-synthetic data generation script with sensitivity analysis
# Script 10_augmentation
################################################################################

source("00_functions.R", local = TRUE)

################################################################################
# Configuration
################################################################################

# Control set parameters
n_pos <- 500
n_neg <- 10000
lambda_fc <- 0.75
dinamicas <- c("uniform","increase","decrease","plateau","inverse_plateau")
z90 <- qnorm(0.95)

# Filtering parameters for triplet construction
min_reports_triplet <- 2         
min_nichd_with_rep <- 2          
all_nichd_rep <- FALSE           
max_events_per_pair <- 10000

# Sensitivity analysis parameters
reduction_levels <- seq(10, 90, by = 10)  # 10%, 20% ... 90%

# Batch processing parameters
batch_size_pos <- 50    # Batch size for positive triplets
save_interval <- 1   # Save every X completed batches

n_null_reports <- 100000  

output_dir <- paste0("./results/", suffix, "/augmentation_results/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Data loading
################################################################################

ade_raw_dt <- fread(ruta_ade_raw)

cols_req <- c("safetyreportid", "atc_concept_id", "meddra_concept_id", "nichd")

# In case sex is included in the GAM formula
if (include_sex) {
  cols_req <- c(cols_req, "sex")
}
if (include_sex) {
  ade_raw_dt[, sex := toupper(trimws(sex))]
  ade_raw_dt[sex == "M", sex := "MALE"]
  ade_raw_dt[sex == "F", sex := "FEMALE"]
  ade_raw_dt[, sex := factor(sex, levels = c("MALE", "FEMALE"))]
  
  sex_summary <- ade_raw_dt[, .(n = .N), by = sex]
  message("\nDistribución de sexo:")
  print(sex_summary)
}

message(sprintf("Dataset %s filas", format(nrow(ade_raw_dt), big.mark = ",")))

################################################################################
# Preprocessing 
################################################################################

drug_info_original <- fread(ruta_drug_info)
drug_info_original[, atc_concept_id := as.character(atc_concept_id)]
drug_info_original[, base_name := tolower(trimws(sub("[;,].*", "", atc_concept_name)))]

canonical_map <- drug_info_original[, .(
  canonical_id = min(atc_concept_id)
), by = base_name]

translation_table <- merge(
  drug_info_original[, .(atc_concept_id, base_name)],
  canonical_map,
  by = "base_name"
)

cat(sprintf("IDs originales: %d\n", uniqueN(translation_table$atc_concept_id)))
cat(sprintf("IDs únicos: %d\n", uniqueN(translation_table$canonical_id)))

ade_raw_dt[, atc_concept_id := as.character(atc_concept_id)]

# Unify drugs sharing the same active compound but different formulations
ade_raw_dt <- merge(
  ade_raw_dt, 
  translation_table[, .(atc_concept_id, canonical_id)], 
  by = "atc_concept_id", 
  all.x = TRUE
)

ade_raw_dt[!is.na(canonical_id), atc_concept_id := canonical_id]
ade_raw_dt[, canonical_id := NULL]

nrow_before <- nrow(ade_raw_dt)
ade_raw_dt <- unique(ade_raw_dt, by = c("safetyreportid", "atc_concept_id", "meddra_concept_id"))

ade_raw_dt[, nichd := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
ade_raw_dt[, nichd_num := as.integer(nichd)]

################################################################################
# Candidate triplet construction
################################################################################

# Extract unique drugs and events per report
drugs_by_report <- unique(ade_raw_dt[, .(safetyreportid, atc_concept_id)])
events_by_report <- unique(ade_raw_dt[, .(safetyreportid, meddra_concept_id)])

reports <- unique(ade_raw_dt[, .(safetyreportid, nichd, nichd_num)])
drugs_list <- drugs_by_report[, .(drugs = list(unique(atc_concept_id))), by = safetyreportid]
events_list <- events_by_report[, .(events = list(unique(meddra_concept_id))), by = safetyreportid]

reports_meta <- merge(reports, drugs_list, by = "safetyreportid", all.x = TRUE)
reports_meta <- merge(reports_meta, events_list, by = "safetyreportid", all.x = TRUE)

report_combo <- copy(reports_meta)  

triplets_list <- pblapply(seq_len(nrow(report_combo)), function(i) {
  rowi <- report_combo[i]
  make_triplets(
    drug = rowi$drugs[[1]], 
    event = rowi$events[[1]], 
    report_id = rowi$safetyreportid, 
    nichd_stage = rowi$nichd
  )
})

triplets_dt <- rbindlist(triplets_list, use.names = TRUE)
rm(triplets_list); gc()

trip_counts_by_stage <- unique(triplets_dt[, .(drugA, drugB, meddra, nichd_num, safetyreportid)])[
  , .N, by = .(drugA, drugB, meddra, nichd_num)
]

trip_summary <- trip_counts_by_stage[, .(
  N = sum(N),
  n_stages = uniqueN(nichd_num),
  stages_with_data = list(nichd_num)
), by = .(drugA, drugB, meddra)]

if (all_nichd_rep) {
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet & 
    n_stages == 7
  ]
  message(sprintf("Filtro estricto: todas las etapas deben tener reportes"))
} else {
  candidatos_pos <- trip_summary[
    N >= min_reports_triplet & 
    n_stages >= min_nichd_with_rep
  ]
  message(sprintf("Filtro flexible: mínimo %d etapas con reportes", min_nichd_with_rep)) # option to apply strict filter requiring all stages
}

# Filtering outcome summary
message("\nEstadísticas de filtrado:")
message(sprintf("  Total candidatos iniciales: %s", 
                format(nrow(trip_summary), big.mark = ",")))
message(sprintf("  Candidatos que cumplen filtro de reportes (N >= %d): %s",
                min_reports_triplet,
                format(sum(trip_summary$N >= min_reports_triplet), big.mark = ",")))
message(sprintf("  Candidatos que cumplen filtro de etapas: %s",
                format(nrow(candidatos_pos), big.mark = ",")))

stage_distribution <- candidatos_pos[, .N, by = n_stages][order(n_stages)]
message("\nDistribución de candidatos por número de etapas:")
print(stage_distribution)

candidatos_pos <- candidatos_pos[, 
  .SD[sample(.N, min(.N, max_events_per_pair))], 
  by = .(drugA, drugB)
]
message(sprintf("Tripletes post diversificación: %d", nrow(candidatos_pos)))

n_pos_final <- min(n_pos, nrow(candidatos_pos))
positivos_sel <- candidatos_pos[sample(.N, n_pos_final)]

# Summary of selected positive candidates
message(sprintf("\nTripletes positivos seleccionados: %d", nrow(positivos_sel)))
message(sprintf("  Reportes promedio: %.1f", mean(positivos_sel$N)))
message(sprintf("  Etapas promedio: %.1f", mean(positivos_sel$n_stages)))

################################################################################
# Positive triplet selection
################################################################################

pos_meta_base <- positivos_sel[, .(drugA, drugB, meddra, N)]
pos_meta_base[, base_triplet_id := 1:.N]

pos_meta <- pos_meta_base[, {
  data.table(
    drugA = drugA,
    drugB = drugB,
    meddra = meddra,
    N = N,
    base_triplet_id = base_triplet_id,
    dynamic = dinamicas
  )
}, by = base_triplet_id]

pos_meta[, fold_change := fold_change(.N, lambda = lambda_fc)]
pos_meta[, triplet_id := 1:.N]

top30_ids <- pos_meta[order(-N)][1:min(.N, 30), triplet_id]
pos_meta[, is_top30 := triplet_id %in% top30_ids]

fwrite(pos_meta, paste0(output_dir, "positive_triplets_metadata.csv"))

################################################################################
# Co-administration counts by NICHD stage for positive triplets
################################################################################

coadmin_stage_pos_list <- lapply(1:nrow(pos_meta), function(i) {
  row <- pos_meta[i]
  result <- coadmin_by_stage(
    drugA = row$drugA, 
    drugB = row$drugB, 
    meddra = row$meddra, 
    ade_data = ade_raw_dt
  )
  result[, `:=`(
    triplet_id = row$triplet_id,
    drugA = row$drugA,
    drugB = row$drugB,
    meddra = row$meddra
  )]
  return(result)
})

coadmin_stage_pos <- rbindlist(coadmin_stage_pos_list)
setcolorder(coadmin_stage_pos, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "nichd", "n_coadmin_stage"))

fwrite(coadmin_stage_pos, paste0(output_dir, "positive_coadmin_by_stage.csv"))

rm(coadmin_stage_pos_list, coadmin_stage_pos)
gc()

################################################################################
# Batch processing of positive triplets with sensitivity analysis
################################################################################

# Set up parallel cluster for batch processing
cl <- makeCluster(n_cores)
registerDoParallel(cl)
  
# Export required functions and data to cluster workers
clusterExport(cl, c("process_single_positive", "fit_reduced_model", 
  "reduce_dataset_by_stage", "inject_signal", "fit_gam",
  "calculate_classic_ior", "calculate_classic_reri", "calc_basic_counts",
  "pos_meta", "ade_raw_dt", "reduction_levels", "niveles_nichd",
  "spline_individuales", "include_sex", "include_stage_sex",
  "k_spline", "bs_type", "select", "nichd_spline", "z90",
  "generate_dynamic", "fold_change"), 
  envir = environment())

clusterEvalQ(cl, {
  library(data.table)
  library(mgcv)
  library(MASS)
  library(doRNG)
})

n_batches <- ceiling(nrow(pos_meta) / batch_size_pos)

# Detect existing checkpoints from previous runs
existing_cp_pos <- list.files(
  output_dir,
  pattern = "checkpoint_positives_batch_\\d+\\.rds",
  full.names = TRUE
)

# If checkpoints exist, load all but the last one and re-run the last batch in case it was corrupted
if (length(existing_cp_pos) > 0) {
  cp_batches_pos <- sort(as.integer(
    gsub(".*checkpoint_positives_batch_(\\d+)\\.rds", "\\1", existing_cp_pos)
  ))
  last_cp_pos <- max(cp_batches_pos)
  redo_batch <- last_cp_pos           # rerun the last batch in case it was corrupted
  load_batches <- cp_batches_pos[cp_batches_pos < redo_batch]
  
  message(sprintf(
    "Último checkpoint detectado: lote %d",
    last_cp_pos
  ))
  
  # Load batches prior to the one being rerun
  if (length(load_batches) > 0) {
    positives_scores_list <- setNames(
      lapply(load_batches, function(b)
        readRDS(paste0(output_dir, "checkpoint_positives_batch_", b, ".rds"))
      ),
      as.character(load_batches)
    )
  } else {
    positives_scores_list <- list()
  }
  
  start_batch_pos <- redo_batch
} else {
  # No checkpoints found: start from scratch
  positives_scores_list <- list()
  start_batch_pos <- 1
}


for (batch in start_batch_pos:n_batches) {
  
  start_idx <- (batch - 1) * batch_size_pos + 1
  end_idx <- min(batch * batch_size_pos, nrow(pos_meta))
  batch_indices <- start_idx:end_idx

  message(sprintf("Lote %d / %d (tripletes %d-%d)", batch, n_batches, start_idx, end_idx))  
  
  # Process batch in parallel
  batch_results <- foreach(
    idx = batch_indices,
    .packages = c("data.table", "mgcv"),
    .errorhandling = "pass",
    .verbose = FALSE,
    .options.RNG = 7113
  ) %dorng% {    # switched from %dopar% to %dorng% for reproducible parallel RNG
    process_single_positive(
      idx, pos_meta, ade_raw_dt, reduction_levels,
      spline_individuales, include_sex, include_stage_sex,
      k_spline, bs_type, select, nichd_spline, z90,
      base_seed = 7113)
  }
  
  # Clean and normalize batch results
  batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)
  
  if (length(batch_results_clean) > 0) {
    batch_dt <- rbindlist(batch_results_clean, fill = TRUE)
    positives_scores_list[[batch]] <- batch_dt
    # Progress info for each processed batch
    message(sprintf("Lote %d completado: %d tripletes exitosos (base)", batch, sum(batch_dt$model_success & batch_dt$reduction_pct == 0, na.rm = TRUE)))
    
    # Save checkpoint every save_interval batches or at the last one
    if (batch %% save_interval == 0 || batch == n_batches) {
      checkpoint_file <- paste0(output_dir, "checkpoint_positives_batch_", batch, ".rds")
      saveRDS(batch_dt, checkpoint_file)
      message(sprintf("Checkpoint guardado: %s", checkpoint_file))
    }
  }
  # Free memory after each batch
  rm(batch_results, batch_results_clean, batch_dt)
  gc(verbose = FALSE)
}
stopCluster(cl)

# Combine all batch results into a single data.table
positives_scores <- rbindlist(positives_scores_list, fill = TRUE)
rm(positives_scores_list)
gc()

# Summary of injection results
message(sprintf("Total tripletes procesados: %d", nrow(positives_scores[reduction_pct == 0])))
message(sprintf("Exitosos (base): %d", sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE)))

# Save results separately for each reduction level
for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if(red_pct == 0) "" else paste0("_", red_pct)
  
  subset_data <- positives_scores[reduction_pct == red_pct]
  
  if (nrow(subset_data) > 0) {
    # Save as RDS
    saveRDS(subset_data, paste0(output_dir, "positive_triplets_results", suffix_file, ".rds"))
    
    # CSV version
    subset_csv <- copy(subset_data)
    
    if ("diagnostics" %in% names(subset_csv)) {
      subset_csv[, diagnostics := NULL]
    }
    
    # Flatten list columns to comma-separated strings
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]
    
    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }
    fwrite(subset_csv, paste0(output_dir, "positive_triplets_results", suffix_file, ".csv"))
  }
}

gc()

pos_success <- positives_scores[reduction_pct == 0 & injection_success == TRUE]

# Summary of successful injections
message("\nPositivos exitosos (base): ", sum(positives_scores$model_success & positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE))

################################################################################
# Negative triplet selection
################################################################################

# Build negative pool from same drugs and events that positive pool
drugs_from_pos <- unique(c(pos_meta$drugA, pos_meta$drugB))
events_from_pos <- unique(pos_meta$meddra)

# Initial pool summary
message(sprintf("Pool de drogas: %d", length(drugs_from_pos)))
message(sprintf("Pool de eventos: %d", length(events_from_pos)))

# Build a lookup set of positive triplet identifiers to exclude from the
# negative pool — prevents contamination between the two sets
pos_triplet_ids <- paste(
  pmin(pos_meta$drugA, pos_meta$drugB),
  pmax(pos_meta$drugA, pos_meta$drugB),
  pos_meta$meddra,
  sep = "_"
)
pos_triplet_set <- unique(pos_triplet_ids)

# Generate all possible drug-pair x event combinations using chunked CJ to
# avoid a single massive cross-join that would exhaust available memory
chunk_size <- 50
n_drugs <- length(drugs_from_pos)
n_chunks <- ceiling(n_drugs / chunk_size)

candidatos_neg_list <- list()
pb <- txtProgressBar(max = n_chunks, style = 3)

for (chunk_idx in 1:n_chunks) {
  start_idx <- (chunk_idx - 1) * chunk_size + 1
  end_idx <- min(chunk_idx * chunk_size, n_drugs)
  
  drugs_chunk <- drugs_from_pos[start_idx:end_idx]
  
  # Cross-join a drug chunk against the full drug and event pools
  chunk_combinations <- CJ(
    drugA = drugs_chunk,
    drugB = drugs_from_pos,
    meddra = events_from_pos
  )
  # Enforce canonical ordering (drugA < drugB) to eliminate mirror-image duplicates
  chunk_combinations[, `:=`(
    drugA_ord = pmin(drugA, drugB),
    drugB_ord = pmax(drugA, drugB)
  )]
  chunk_combinations <- chunk_combinations[
    drugA == drugA_ord
  ]
  chunk_combinations[, `:=`(
    drugA = drugA_ord,
    drugB = drugB_ord,
    drugA_ord = NULL,
    drugB_ord = NULL
  )]

  # Remove self-pairs and any triplet already used as a positive
  chunk_combinations <- chunk_combinations[drugA != drugB]
  chunk_combinations[, triplet_id := paste(drugA, drugB, meddra, sep = "_")]
  chunk_combinations <- chunk_combinations[!triplet_id %in% pos_triplet_set]
  # Keep only triplets that appear in the observed data with enough reports
  chunk_candidates <- merge(
    chunk_combinations[, .(drugA, drugB, meddra, triplet_id)],
    trip_summary,
    by = c("drugA", "drugB", "meddra"),
    all.x = FALSE
  )

  chunk_candidates <- chunk_candidates[N >= min_reports_triplet]
  
  if (nrow(chunk_candidates) > 0) {
    candidatos_neg_list[[chunk_idx]] <- chunk_candidates
  }
  
  rm(chunk_combinations, chunk_candidates)
  gc(verbose = FALSE)
  
  setTxtProgressBar(pb, chunk_idx)
}
close(pb)

candidatos_neg <- rbindlist(candidatos_neg_list, use.names = TRUE)
rm(candidatos_neg_list)
gc()

candidatos_neg <- unique(candidatos_neg, by = "triplet_id")

# Count observed reports per NICHD stage for each negative candidate to apply
# the minimum-stages filter, mirroring the filter used for positive candidates
neg_ids <- paste(candidatos_neg$drugA, candidatos_neg$drugB, candidatos_neg$meddra, sep = "_")

neg_counts_by_stage <- unique(
  triplets_dt[paste(drugA, drugB, meddra, sep = "_") %in% neg_ids, 
              .(drugA, drugB, meddra, nichd_num, safetyreportid)]
)[, .(n_reports = .N), by = .(drugA, drugB, meddra, nichd_num)]

# Rename column
setnames(neg_counts_by_stage, "nichd_num", "nichd")

# Map NICHD stage names to numeric indices (niveles_nichd is defined in 00_functions.R)
nichd_to_num <- setNames(1:7, niveles_nichd)
neg_counts_by_stage[, nichd_num := nichd_to_num[nichd]]

# 4. Compute n_stages_with_data directly (no dcast needed)
stage_counts_summary <- neg_counts_by_stage[, .(
  n_stages_with_data = uniqueN(nichd_num[n_reports > 0])
), by = .(drugA, drugB, meddra)]

# Merge stage counts back into negative candidate table
setkey(candidatos_neg, drugA, drugB, meddra)
setkey(stage_counts_summary, drugA, drugB, meddra)

candidatos_neg_full <- merge(
  candidatos_neg,
  stage_counts_summary,
  by = c("drugA", "drugB", "meddra"),
  all.x = TRUE
)
# Apply the same minimum-stages criterion used for positive candidates
candidatos_neg_filtered <- candidatos_neg_full[n_stages_with_data >= min_nichd_with_rep]

# Summary of filtered negative candidates
message(sprintf("\nCandidatos negativos filtrados: %d (de %d)", nrow(candidatos_neg_filtered), nrow(candidatos_neg)))

# Random sample from the filtered pool up to the configured n_neg ceiling
n_neg_final <- min(n_neg, nrow(candidatos_neg_filtered))
selected_negatives <- candidatos_neg_filtered[sample(.N, n_neg_final)]

selected_negatives[, triplet_id := 1:.N]

# Summary of selected negative pool
message(sprintf("\nNegativos seleccionados: %d", nrow(selected_negatives)))
message(sprintf("  Reportes promedio: %.1f", mean(selected_negatives$N)))
message(sprintf("  Etapas promedio: %.1f", mean(selected_negatives$n_stages_with_data)))

################################################################################
# Co-administration counts by NICHD stage for negative triplets
################################################################################

# Computes how many co-administration reports exist per NICHD stage for each
# selected negative triplet; saved for use in downstream analysis
coadmin_stage_neg_list <- lapply(1:nrow(selected_negatives), function(i) {
  row <- selected_negatives[i]
  result <- coadmin_by_stage(
    drugA = row$drugA, 
    drugB = row$drugB, 
    meddra = row$meddra, 
    ade_data = ade_raw_dt
  )
  result[, `:=`(
    triplet_id = row$triplet_id,
    drugA = row$drugA,
    drugB = row$drugB,
    meddra = row$meddra
  )]
  return(result)
})

coadmin_stage_neg <- rbindlist(coadmin_stage_neg_list)
setcolorder(coadmin_stage_neg, c("triplet_id", "drugA", "drugB", "meddra", "nichd_num", "n_coadmin_stage"))

fwrite(coadmin_stage_neg, paste0(output_dir, "negative_coadmin_by_stage.csv"))

rm(coadmin_stage_neg_list, coadmin_stage_neg)
gc()

# Pre-compute reduced datasets for each sensitivity level; passed to workers later
reduced_datasets <- setNames(
  lapply(reduction_levels, function(red_pct) {
    reduce_dataset_by_stage(ade_raw_dt, red_pct, seed = 7113)
  }),
  as.character(reduction_levels)
)

################################################################################
# Batch processing of negative triplets
################################################################################

batch_size_neg <- 50
n_batches      <- ceiling(nrow(selected_negatives) / batch_size_neg)

# Helper: detect existing checkpoints and return the batch to start from 
detect_neg_checkpoints <- function(out_dir, pass_tag) {
  pattern  <- sprintf("checkpoint_neg_%s_batch_\\d+\\.rds", pass_tag)
  cp_files <- list.files(out_dir, pattern = pattern, full.names = TRUE)

  if (length(cp_files) == 0) return(1L)

  cp_nums <- sort(as.integer(
    gsub(sprintf(".*checkpoint_neg_%s_batch_(\\d+)\\.rds", pass_tag), "\\1", cp_files)
  ))
  last_cp <- max(cp_nums)
  message(sprintf("  [%s] Checkpoint found — resuming from batch %d", pass_tag, last_cp))
  return(last_cp)
}

# Helper: run all batches for a single reduction-level pass 
# pass_tag : string used in checkpoint filenames ("base", "red10", "red20", ...)
# ade_pass : the dataset to use for this pass (full or reduced)
# red_pct  : numeric value stored in the reduction_pct output column
#
# The parallel backend must be registered (registerDoParallel) before calling
# this function. GAM parameters are read from the enclosing environment via
# lexical scoping (spline_individuales, k_spline, etc.).

run_neg_batch_pass <- function(pass_tag, ade_pass, red_pct, cl) {

  # Export ade_pass from this function's local environment to all workers
  # clusterExport by default looks in .GlobalEnv, so envir must be set explicitly here since ade_pass is a function argument, not a global variable.
  clusterExport(cl, "ade_pass", envir = environment())

  start_batch <- detect_neg_checkpoints(output_dir, pass_tag)

  for (batch in start_batch:n_batches) {

    start_idx <- (batch - 1) * batch_size_neg + 1
    end_idx <- min(batch * batch_size_neg, nrow(selected_negatives))
    batch_indices <- start_idx:end_idx

    message(sprintf("\n[%s] Batch %d / %d  (triplets %d-%d)",
                    pass_tag, batch, n_batches, start_idx, end_idx))
    batch_results <- foreach(
      idx = batch_indices,
      .packages = c("data.table", "mgcv"),
      .errorhandling = "pass",
      .verbose = FALSE,
      .options.RNG = 7113
    ) %dorng% {

      set.seed(7113 + idx)

      rowt <- selected_negatives[idx]
      rowt$type <- "negative"

      # Basic co-occurrence counts used as fallback when the model fails
      counts <- tryCatch(
        calc_basic_counts(ade_pass, rowt$drugA, rowt$drugB, rowt$meddra),
        error = function(e) list(n_events_coadmin = NA, n_events_total = NA, n_coadmin = NA)
      )

      tryCatch({

        model_res <- fit_gam(
          drugA_id = rowt$drugA,
          drugB_id = rowt$drugB,
          event_id = rowt$meddra,
          ade_data = ade_pass,
          spline_individuales = spline_individuales,
          include_sex = include_sex,
          include_stage_sex = include_stage_sex,
          k_spline = k_spline,
          bs_type = bs_type,
          select = select,
          nichd_spline = nichd_spline
        )

        classic_res <- calculate_classic_ior(
          drugA_id = rowt$drugA,
          drugB_id = rowt$drugB,
          event_id = rowt$meddra,
          ade_data = ade_pass
        )

        classic_reri <- calculate_classic_reri(
          drugA_id = rowt$drugA,
          drugB_id = rowt$drugB,
          event_id = rowt$meddra,
          ade_data = ade_pass
        )

        if (!model_res$success) {
          # GAM failed .return NA placeholders for all model-derived columns
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = counts$n_events_coadmin,
            model_success = FALSE,
            n_events = counts$n_events_total,
            n_coadmin = counts$n_coadmin,
            n_stages_significant = NA_integer_,
            max_ior = NA_real_,
            mean_ior = NA_real_,
            model_aic = NA_real_,
            stage = list(1:7),
            log_ior = list(rep(NA_real_, 7)),
            log_ior_lower90 = list(rep(NA_real_, 7)),
            ior_values = list(rep(NA_real_, 7)),
            formula_used = if (!is.null(model_res$formula_attempted)) model_res$formula_attempted else NA_character_,
            message = if (!is.null(model_res$error_msg))         model_res$error_msg         else NA_character_,
            classic_success = classic_res$success,
            log_ior_classic = list(rep(NA_real_, 7)),
            log_ior_classic_lower90 = list(rep(NA_real_, 7)),
            ior_classic = list(rep(NA_real_, 7)),
            reri_classic_success = FALSE,
            reri_values = list(rep(NA_real_, 7)),
            reri_lower90 = list(rep(NA_real_, 7)),
            reri_upper90 = list(rep(NA_real_, 7)),
            n_stages_reri_significant = NA_integer_,
            RERI_classic = list(rep(NA_real_, 7)),
            RERI_classic_lower90 = list(rep(NA_real_, 7)),
            RERI_classic_upper90 = list(rep(NA_real_, 7)),
            RERI_classic_se = list(rep(NA_real_, 7))
          )
        } else {
          # GAM succeeded. store all per-stage estimates
          data.table(
            triplet_id = rowt$triplet_id,
            drugA = rowt$drugA,
            drugB = rowt$drugB,
            meddra = rowt$meddra,
            type = "negative",
            reduction_pct = red_pct,
            N = model_res$n_events_coadmin,
            model_success = TRUE,
            n_coadmin  = model_res$n_coadmin,
            n_events = model_res$n_events_total,
            n_stages_significant = model_res$n_stages_significant,
            max_ior = model_res$max_ior,
            mean_ior = model_res$mean_ior,
            model_aic = model_res$model_aic,
            stage = list(1:7),
            log_ior = list(model_res$log_ior),
            log_ior_lower90 = list(model_res$log_ior_lower90),
            ior_values = list(model_res$ior_values),
            classic_success = classic_res$success,
            log_ior_classic = if (classic_res$success) list(classic_res$results_by_stage$log_ior_classic) else list(rep(NA_real_, 7)),
            log_ior_classic_lower90 = if (classic_res$success) list(classic_res$results_by_stage$log_ior_classic_lower90) else list(rep(NA_real_, 7)),
            ior_classic = if (classic_res$success) list(classic_res$results_by_stage$ior_classic) else list(rep(NA_real_, 7)),
            reri_classic_success = classic_reri$success,
            reri_values = list(model_res$reri_values),
            reri_lower90 = list(model_res$reri_lower90),
            reri_upper90 = list(model_res$reri_upper90),
            n_stages_reri_significant = model_res$n_stages_reri_significant,
            RERI_classic = if (classic_reri$success) list(classic_reri$results_by_stage$RERI_classic) else list(rep(NA_real_, 7)),
            RERI_classic_lower90 = if (classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_lower90) else list(rep(NA_real_, 7)),
            RERI_classic_upper90 = if (classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_upper90) else list(rep(NA_real_, 7)),
            RERI_classic_se = if (classic_reri$success) list(classic_reri$results_by_stage$RERI_classic_se) else list(rep(NA_real_, 7))
          )
        }

      }, error = function(e) {
        data.table(
          triplet_id = rowt$triplet_id,
          drugA = rowt$drugA,
          drugB = rowt$drugB,
          meddra = rowt$meddra,
          type = "negative",
          reduction_pct = red_pct,
          N = counts$n_events_coadmin,
          model_success = FALSE,
          n_events = counts$n_events_total,
          n_coadmin = counts$n_coadmin,
          error_msg = paste("Unhandled error:", e$message)
        )
      })

    } # end foreach

    # Filter out foreach-level condition objects (distinct from per-triplet
    # errors already handled by the inner tryCatch above)
    batch_results_clean <- Filter(function(x) !inherits(x, "error"), batch_results)

    if (length(batch_results_clean) > 0) {
      batch_dt <- rbindlist(batch_results_clean, fill = TRUE)

      # Write to disk immediately — this is the only stored copy of this batch;
      # freed right after to prevent accumulation across batches
      checkpoint_file <- sprintf(
        "%scheckpoint_neg_%s_batch_%d.rds", output_dir, pass_tag, batch
      )
      saveRDS(batch_dt, checkpoint_file)

      message(sprintf("  Checkpoint saved: %s  (%d/%d successful)",
                      basename(checkpoint_file),
                      sum(batch_dt$model_success, na.rm = TRUE),
                      nrow(batch_dt)))
      rm(batch_dt)
    }

    rm(batch_results, batch_results_clean)
    gc(verbose = FALSE)

  } # end batch loop
}

# Workers receive only ade_raw_dt. no reduced copies exported here
message(sprintf("\nPass 1 / %d: base dataset (reduction = 0%%) ",
                length(reduction_levels) + 1))

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c(
  "fit_gam", "ade_raw_dt", "z90",
  "selected_negatives", "niveles_nichd",
  "spline_individuales", "include_sex", "include_stage_sex",
  "k_spline", "nichd_spline", "include_nichd",
  "bs_type", "select", "method",
  "calculate_classic_ior", "calculate_classic_reri",
  "calc_basic_counts"
), envir = environment())

clusterEvalQ(cl, {
  library(data.table)
  library(mgcv)
  library(MASS)
  library(doRNG)
})

run_neg_batch_pass(pass_tag = "base", ade_pass = ade_raw_dt, red_pct = 0, cl = cl)

stopCluster(cl)
gc()

# Pass 2
# Each iteration:
#  Builds a single reduced copy of the dataset for this level only
#  Starts a fresh cluster and exports only that one copy
#  Processes all batches for this level
#  Stops the cluster and frees the reduced dataset before the next iteration

for (red_pct in reduction_levels) {

  pass_tag <- sprintf("red%d", red_pct)
  n_pass <- which(reduction_levels == red_pct) + 1

  message(sprintf("\nPass %d / %d: reduced dataset (%d%% reduction)",
                  n_pass, length(reduction_levels) + 1, red_pct))

  # Build this level's reduced dataset. never co-exists with another reduced copy
  ade_reduced <- reduce_dataset_by_stage(ade_raw_dt, red_pct, seed = 7113)

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  clusterExport(cl, c(
    "fit_gam", "z90", "niveles_nichd",
    "selected_negatives",
    "spline_individuales", "include_sex", "include_stage_sex",
    "k_spline", "nichd_spline", "include_nichd",
    "bs_type", "select", "method",
    "calculate_classic_ior", "calculate_classic_reri",
    "calc_basic_counts",
    "ade_reduced"          # single reduced copy, not a list of all levels
  ), envir = environment())

  clusterEvalQ(cl, {
    library(data.table)
    library(mgcv)
    library(MASS)
    library(doRNG)
  })

  run_neg_batch_pass(pass_tag = pass_tag, ade_pass = ade_reduced, red_pct = red_pct, cl = cl)

  stopCluster(cl)

  # Release the reduced dataset before the next iteration creates a new one
  rm(ade_reduced)
  gc()
}

################################################################################
# Assemble final results from checkpoints and save per reduction level
################################################################################
# All pass results live on disk as per-batch checkpoint files. This is the
# first time the complete dataset is materialised in memory, and only after
# every cluster has been stopped and its memory released.

all_pass_tags <- c("base", sprintf("red%d", reduction_levels))

negatives_scores <- rbindlist(
  lapply(all_pass_tags, function(tag) {
    cp_files <- sort(list.files(
      output_dir,
      pattern = sprintf("checkpoint_neg_%s_batch_\\d+\\.rds", tag),
      full.names = TRUE
    ))
    if (length(cp_files) == 0) {
      warning(sprintf("No checkpoint files found for pass '%s'", tag))
      return(NULL)
    }
    rbindlist(lapply(cp_files, readRDS), fill = TRUE)
  }),
  fill = TRUE
)

message(sprintf("Total negative rows assembled: %s",
                format(nrow(negatives_scores), big.mark = ",")))
message(sprintf("Successful (base): %d",
                sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0,
                    na.rm = TRUE)))

# Save one RDS and one flat CSV per reduction level for downstream scripts
for (red_pct in c(0, reduction_levels)) {
  suffix_file <- if (red_pct == 0) "" else paste0("_", red_pct)

  subset_data <- negatives_scores[reduction_pct == red_pct]

  if (nrow(subset_data) > 0) {

    saveRDS(subset_data,
            paste0(output_dir, "negative_triplets_results", suffix_file, ".rds"))

    # Flatten list columns to comma-separated strings for CSV compatibility
    subset_csv <- copy(subset_data)
    list_cols <- names(subset_csv)[sapply(subset_csv, is.list)]

    for (col in list_cols) {
      subset_csv[, (col) := sapply(get(col), function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ",")
      })]
    }

    fwrite(subset_csv,
           paste0(output_dir, "negative_triplets_results", suffix_file, ".csv"))
    rm(subset_csv)
  }
}

################################################################################
# Null pool creation
################################################################################
# Simple random sampling; drug/event labels are permuted later to destroy real associations

# Sampling parameters
n_null_reports <- 100000             

# Get all unique report IDs available in the dataset
all_reports <- unique(ade_raw_dt$safetyreportid)

# Reports potentially overlapping with the positive set are not excluded because:
# Drug/event permutation destroys any real associations in this pool
# This maximizes the available sample size for the null distribution
selected_null_reports <- sample(all_reports, n_null_reports)

# Null pool metadata: report ID and NICHD stage
null_pool_meta <- unique(ade_raw_dt[
  safetyreportid %in% selected_null_reports,
  .(safetyreportid, nichd, nichd_num)
])

# Verify report distribution across NICHD stages
stage_distribution <- null_pool_meta[, .N, by = nichd][order(nichd)]
message("\nDistribución de reportes por etapa NICHD:")
print(stage_distribution)

# Additional diagnostics could be added to verify that selected reports are representative (especially by stage)
# Given the large dataset size, representativeness is probabilistically very likely, but worth noting

# Save null pool metadata
fwrite(null_pool_meta, paste0(output_dir, "null_pool_reports_metadata.csv"))

message("Total de reportes: ", nrow(null_pool_meta))

################################################################################
# Positive vs. Negative comparison
################################################################################

pos_success <- positives_scores[model_success == TRUE & injection_success == TRUE & reduction_pct == 0]
neg_success <- negatives_scores[model_success == TRUE & reduction_pct == 0]

message("\nTriplets exitosos (base):")
message(" Positivos: ", nrow(pos_success))
message(" Negativos: ", nrow(neg_success))

if (nrow(pos_success) > 0 && nrow(neg_success) > 0) {
  
  pos_expanded <- pos_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type, dynamic)]
  
  neg_expanded <- neg_success[, {
    stages <- unlist(stage)
    log_iors <- unlist(log_ior)
    log_ior_l90s <- unlist(log_ior_lower90)
    
    n <- min(length(stages), length(log_iors), length(log_ior_l90s))
    if (n == 0) {
      data.table()
    } else {
      data.table(
        stage = stages[1:n],
        log_ior = log_iors[1:n],
        log_ior_lower90 = log_ior_l90s[1:n]
      )
    }
  }, by = .(triplet_id, type)]
  
  message("\nEstadísticas de log-IOR (base):")
  message("  Positivos. Media: ", round(mean(pos_expanded$log_ior, na.rm = TRUE), 4))
  message("  Negativos.  Media: ", round(mean(neg_expanded$log_ior, na.rm = TRUE), 4))
  
  if (nrow(pos_expanded) > 30 && nrow(neg_expanded) > 30) {
    test_result <- wilcox.test(
      pos_expanded$log_ior,
      neg_expanded$log_ior,
      alternative = "greater"
    )
    message("\nWilcoxon Test (Positivos > Negativos):")
    message("  p-value: ", format.pval(test_result$p.value, digits = 3))
  }
  
  all_expanded <- rbind(
    pos_expanded[, .(triplet_id, type, log_ior, stage)],
    neg_expanded[, .(triplet_id, type, log_ior, stage)]
  )
  
  p1 <- ggplot(all_expanded, aes(x = type, y = log_ior, fill = type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Distribución de Log-IOR por tipo de Triplete (Base)",
      x = "Tipo",
      y = "Log-IOR"
    )
  
  ggsave(paste0(output_dir, "comparison_log_ior_distribution.png"),
         p1, width = 8, height = 6, dpi = 300)
  
}

################################################################################
# Final summary
################################################################################

summary_stats <- data.table(
  metric = c(
    "config_min_reports_triplet",
    "config_min_nichd_with_rep",
    "config_all_nichd_rep",
    "reduction_levels_tested",
    "n_positive_total",
    "n_positive_injected",
    "n_positive_modeled_base",
    "n_negative_total",
    "n_negative_modeled_base",
    "mean_ior_positive_base",
    "mean_ior_negative_base",
    "mean_stages_sig_positive_base",
    "mean_stages_sig_negative_base"
  ),
  value = c(
    min_reports_triplet,
    min_nichd_with_rep,
    as.numeric(all_nichd_rep),
    paste(reduction_levels, collapse = ","),
    nrow(pos_meta),
    sum(positives_scores$injection_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    sum(positives_scores$model_success & positives_scores$reduction_pct == 0, na.rm = TRUE),
    nrow(selected_negatives),
    sum(negatives_scores$model_success & negatives_scores$reduction_pct == 0, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$mean_ior, na.rm = TRUE),
    mean(positives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE),
    mean(negatives_scores[reduction_pct == 0]$n_stages_significant, na.rm = TRUE)
  )
)

print(summary_stats)
fwrite(summary_stats, paste0(output_dir, "summary_statistics.csv"))

################################################################################
# Aggregated sensitivity analysis
################################################################################

sensitivity_summary <- rbind(
  positives_scores[, .(
    type = "positive",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct],
  negatives_scores[, .(
    type = "negative",
    reduction_pct = reduction_pct[1],
    n_total = .N,
    n_success = sum(model_success, na.rm = TRUE),
    mean_ior = mean(mean_ior, na.rm = TRUE),
    mean_stages_sig = mean(n_stages_significant, na.rm = TRUE)
  ), by = reduction_pct]
)

print(sensitivity_summary)
fwrite(sensitivity_summary, paste0(output_dir, "sensitivity_analysis_summary.csv"))

################################################################################
# Dynamics detection
################################################################################

# Checks whether the model correctly recovers the shape of the injected dynamic

# Skip rerunning the full pipeline from scratch
ruta_pos_results <- paste0("./results/", suffix, "/augmentation_results/positive_triplets_results.rds")
positives_scores <- readRDS(ruta_pos_results)

# Filter to successfully processed positive triplets
pos_for_dynamics <- positives_scores[
  model_success == TRUE & 
  injection_success == TRUE &
  reduction_pct == 0 &   
  !is.na(dynamic)
]

# Expand per-stage data from list columns
pos_dynamics_expanded <- pos_for_dynamics[, {
  stages <- unlist(stage)
  log_iors <- unlist(log_ior)
  
  # Validate list lengths before indexing
  n <- min(length(stages), length(log_iors))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      log_ior = log_iors[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

# Add stage name labels after expansion
pos_dynamics_expanded[, stage_name := niveles_nichd[stage]]

# Compute mean log-IOR by injected dynamic and developmental stage
dynamics_summary <- pos_dynamics_expanded[, .(
  mean_log_ior = mean(log_ior, na.rm = TRUE),
  sd_log_ior = sd(log_ior, na.rm = TRUE),
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

# Add stage name labels to summary table
dynamics_summary[, stage_name := niveles_nichd[stage]]

# Use the uniform dynamic as the reference baseline
uniform_baseline <- dynamics_summary[dynamic == "uniform", .(
  stage, 
  baseline_log_ior = mean_log_ior
)]

# Compute differences relative to the uniform baseline
dynamics_diff <- merge(
  dynamics_summary[dynamic != "uniform"],
  uniform_baseline,
  by = "stage",
  all.x = TRUE
)

dynamics_diff[, log_ior_diff := mean_log_ior - baseline_log_ior]

print(dynamics_diff[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_log_ior, 
  baseline = baseline_log_ior,
  difference = log_ior_diff
)])

################################################################################
# Bootstrap confidence intervals
################################################################################

n_boot <- 100

# Apply bootstrap to all dynamic-stage combinations
dynamics_nonuniform <- unique(pos_dynamics_expanded[dynamic != "uniform", dynamic])
stages <- 1:7

bootstrap_results <- rbindlist(pblapply(dynamics_nonuniform, function(dyn) {
  rbindlist(lapply(stages, function(s) {
    boot_res <- bootstrap_dynamic_diff(pos_dynamics_expanded, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

# Merge bootstrap CIs back into the main differences table
dynamics_with_ci <- merge(
  dynamics_diff,
  bootstrap_results,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci, paste0(output_dir, "dynamics_recovery_analysis.csv"))

recovery_stats <- dynamics_with_ci[, .(
  mean_difference = mean(log_ior_diff, na.rm = TRUE),
  max_difference = max(abs(log_ior_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats)

################################################################################
# Visualization
################################################################################

# Color palette for injected dynamics
dynamic_colors <- c(
  "increase" = "#E41A1C",
  "decrease" = "#377EB8", 
  "plateau" = "#4DAF4A",
  "inverse_plateau" = "#984EA3"
)

nichd_labels <- c(
  "term_neonatal" = "Neonato a term.",
  "infancy" = "Lactante",
  "toddler" = "Deambulador",
  "early_childhood" = "Preescolar",
  "middle_childhood" = "Escolar",
  "early_adolescence" = "Adolescencia temp.",
  "late_adolescence" = "Adolescencia tardía"
)

# Delta log-IOR vs. uniform dynamic with bootstrap confidence intervals
p_dynamics_diff <- ggplot(
  dynamics_with_ci[!is.na(mean_diff)],
  aes(x = stage_name, y = log_ior_diff, color = dynamic, group = dynamic)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  labs(
    x = "Etapa",
    y = "Δ Log(IOR) (vs uniforme)",
    color = "Dinámica inyectada",
    fill = "Dinámica inyectada"
  ) 

ggsave(
  paste0(output_dir, "fig_dynamics_recovery.png"),
  p_dynamics_diff,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_dynamics_diff)

################################################################################
# Dynamics detection using RERI
################################################################################

# Checks whether the model correctly recovers the shape of the injected dynamic
# using GAM-derived RERI (absolute risk scale) rather than IOR (odds ratio scale)
 
# Expand per-stage data including GAM-derived RERI values
pos_dynamics_expanded_reri <- pos_for_dynamics[, {
  stages <- unlist(stage)
  reri_vals <- unlist(reri_values)  # GAM-derived RERI
  
  # Validate list lengths
  n <- min(length(stages), length(reri_vals))
  
  if (n > 0) {
    data.table(
      stage = stages[1:n],
      reri = reri_vals[1:n]
    )
  } else {
    data.table()
  }
}, by = .(triplet_id, dynamic, fold_change)]

# Add stage name labels after expansion
pos_dynamics_expanded_reri[, stage_name := niveles_nichd[stage]]

# Compute median RERI by injected dynamic and developmental stage
dynamics_summary_reri <- pos_dynamics_expanded_reri[, .(
  mean_reri = median(reri, na.rm = TRUE),   # median, not mean
  sd_reri = mad(reri, na.rm = TRUE),        # MAD instead of SD
  n_triplets = uniqueN(triplet_id)
), by = .(dynamic, stage)]

# Add stage name labels to RERI summary
dynamics_summary_reri[, stage_name := niveles_nichd[stage]]

# Compute uniform baseline for RERI (mirrors the log-IOR baseline block above)
uniform_baseline_reri <- dynamics_summary_reri[dynamic == "uniform", .(
  stage,
  baseline_reri = mean_reri   # already median; keeping name for consistency
)]

# Differences vs. uniform baseline (delta RERI)
dynamics_diff_reri <- merge(
  dynamics_summary_reri[dynamic != "uniform"],
  uniform_baseline_reri,
  by = "stage",
  all.x = TRUE
)

dynamics_diff_reri[, reri_diff := mean_reri - baseline_reri]

print(dynamics_diff_reri[order(dynamic, stage), .(
  dynamic, 
  stage_name, 
  mean_reri, 
  baseline = baseline_reri,
  difference = reri_diff
)])

################################################################################
# Bootstrap confidence intervals for RERI
################################################################################

# Apply bootstrap to all dynamic-stage combinations (RERI)
dynamics_nonuniform_reri <- unique(pos_dynamics_expanded_reri[dynamic != "uniform", dynamic])
stages_reri <- 1:7

bootstrap_results_reri <- rbindlist(pblapply(dynamics_nonuniform_reri, function(dyn) {
  rbindlist(lapply(stages_reri, function(s) {
    boot_res <- bootstrap_dynamic_diff_reri(pos_dynamics_expanded_reri, dyn, s, n_boot)
    cbind(data.table(dynamic = dyn, stage = s), boot_res)
  }))
}))

# Merge bootstrap CIs into the RERI differences table
dynamics_with_ci_reri <- merge(
  dynamics_diff_reri,
  bootstrap_results_reri,
  by = c("dynamic", "stage"),
  all.x = TRUE
)

dynamics_with_ci_reri[, stage_name := factor(stage_name, levels = niveles_nichd)]

fwrite(dynamics_with_ci_reri, paste0(output_dir, "dynamics_recovery_analysis_reri.csv"))

recovery_stats_reri <- dynamics_with_ci_reri[, .(
  mean_difference = mean(reri_diff, na.rm = TRUE),
  max_difference = max(abs(reri_diff), na.rm = TRUE),
  stages_significant = sum(ci_lower > 0 | ci_upper < 0, na.rm = TRUE)
), by = dynamic]

print(recovery_stats_reri)

################################################################################
# Delta RERI plot (vs. uniform dynamic)
################################################################################F

# Delta RERI vs. uniform dynamic with bootstrap confidence intervals

p_dynamics_diff_reri <- ggplot(
  dynamics_with_ci_reri[!is.na(mean_diff)],
  aes(x = stage_name, y = reri_diff, color = dynamic, group = dynamic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = dynamic),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  scale_x_discrete(labels = nichd_labels) +
  scale_fill_manual(
    values = dynamic_colors,
    labels = c(
      "increase" = "Creciente",
      "decrease" = "Decreciente",
      "plateau" = "Meseta",
      "inverse_plateau" = "Valle"
    )
  ) +
  labs(
    x = "Etapa",
    y = "Δ RERI (vs uniforme)",
    color = "Dinámica inyectada",
    fill = "Dinámica inyectada"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
  
ggsave(
  paste0(output_dir, "fig_dynamics_recovery_reri.png"),
  p_dynamics_diff_reri,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_dynamics_diff_reri)
