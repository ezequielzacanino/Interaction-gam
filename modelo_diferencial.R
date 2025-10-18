# script 2
# Librerías ---------------------------------------------------------------
library(data.table)
library(mgcv)
library(tidyverse)
library(gtsummary)
library(broom)
library(gt)
library(rmarkdown)
library(knitr)
library(kableExtra)

# -------------------------------------------------------------------------
# Carga de datos
# -------------------------------------------------------------------------

# directorio 
dir_datos <- "./"

# carga de datos 
ade_raw_dt <- fread(paste0(dir_datos, "ade_raw.csv"))
drug_dt <- fread(paste0(dir_datos, "drug.csv"))
event_dt <- fread(paste0(dir_datos, "event.csv"))

# Merge 
# mapeo meddra_concept_4
# más eficiente con data.table que con dyplr, porque dyplr crea nuevo objeto
ade_raw_dt[event_dt, on = .(meddra_concept_id), meddra_concept_4 := i.meddra_concept_code_4]
ade_raw_dt[event_dt, on = .(meddra_concept_id), meddra_concept_2 := i.meddra_concept_code_2]


# =========================================================================
# SCRIPT DE ESCANEO: Identificación de fármacos con reportes conjuntos
# =========================================================================

# Definir el fármaco y evento adverso de interés
id_droga_referencia <- 21604422  # Cambiar por el atc_concept_id de interés
id_evento_referencia <- 36718321  # Cambiar por el meddra_concept_4 de interés

# Obtener nombre del fármaco y evento de referencia
nombre_referencia <- drug_dt[atc_concept_id == id_droga_referencia, atc_concept_name][1]
nombre_evento_ref <- event_dt[meddra_concept_id == id_evento_referencia, meddra_concept_name_1][1]

message("\n", paste(rep("=", 70), collapse = ""))
message("ESCANEO DE FÁRMACOS CON REPORTES CONJUNTOS")
message(paste(rep("=", 70), collapse = ""))
message(paste("Fármaco de referencia:", nombre_referencia, "(ID:", id_droga_referencia, ")"))
message(paste("Evento adverso:", nombre_evento_ref, "(ID:", id_evento_referencia, ")"))
message(paste(rep("=", 70), collapse = ""))

# Identificar reportes con el fármaco de referencia
reportes_ref <- unique(ade_raw_dt[atc_concept_id == id_droga_referencia, safetyreportid])
message(paste("\nReportes totales con", nombre_referencia, ":", length(reportes_ref)))

# Identificar reportes con el evento adverso
reportes_evento <- unique(ade_raw_dt[meddra_concept_id == id_evento_referencia, safetyreportid])
message(paste("Reportes totales con", nombre_evento_ref, ":", length(reportes_evento)))

# Reportes que tienen AMBOS (fármaco + evento)
reportes_ref_con_evento <- intersect(reportes_ref, reportes_evento)
message(paste("Reportes con", nombre_referencia, "Y", nombre_evento_ref, ":", 
              length(reportes_ref_con_evento)))

# Encontrar todos los otros fármacos que aparecen con el fármaco de referencia
# en reportes que también tienen el evento adverso
message("\nBuscando fármacos coadministrados en reportes con el evento...")

# Obtener todos los fármacos presentes en reportes con fármaco_ref + evento
otros_farmacos <- ade_raw_dt[
  safetyreportid %in% reportes_ref_con_evento & 
  atc_concept_id != id_droga_referencia,
  .(N_reportes_conjuntos = uniqueN(safetyreportid)),
  by = .(atc_concept_id)
][order(-N_reportes_conjuntos)]

# Agregar nombres de fármacos
otros_farmacos[drug_dt, on = .(atc_concept_id), 
               atc_concept_name := i.atc_concept_name]

# Calcular métricas adicionales
otros_farmacos[, N_reportes_totales := uniqueN(
  ade_raw_dt[atc_concept_id == .BY[[1]], safetyreportid]
), by = atc_concept_id]

otros_farmacos[, porcentaje_conjunto := 
                 round(100 * N_reportes_conjuntos / N_reportes_totales, 2)]

# Reordenar columnas
otros_farmacos <- otros_farmacos[, .(
  atc_concept_id, 
  atc_concept_name,
  N_reportes_conjuntos,
  N_reportes_totales,
  porcentaje_conjunto
)]

# Mostrar top 20
message("\n", paste(rep("-", 70), collapse = ""))
message("TOP 20 FÁRMACOS CON MÁS REPORTES CONJUNTOS")
message(paste(rep("-", 70), collapse = ""))
print(head(otros_farmacos, 20))

# -------------------------------------------------------------------------

# triplete candidato
id_droga_a <- 21604415  #clonazepam
id_droga_b <- 21604422  #acido valproico
id_evento_adverso <- 36718321  #somnoliencia

# -------------------------------------------------------------------------
# Preprocesado
# -------------------------------------------------------------------------


# mapeo de nombres desde el dataset "drug"
nombre_a <- drug_dt[atc_concept_id == id_droga_a, atc_concept_name][1]
nombre_b <- drug_dt[atc_concept_id == id_droga_b, atc_concept_name][1]

# clasificación de reportes
reportes_droga_a <- unique(ade_raw_dt[atc_concept_id == id_droga_a, safetyreportid]) # junto id de reportes con cada fármaco
reportes_droga_b <- unique(ade_raw_dt[atc_concept_id == id_droga_b, safetyreportid])
reportes_conjuntos <- intersect(reportes_droga_a, reportes_droga_b) # intersección para ver administración conjunta

message(paste("Reportes conjuntos:", length(reportes_conjuntos)))


# preparo dataset para GAM 
# modificar cuando se agreguen más covariables al modelo
datos_modelo <- unique(ade_raw_dt[, .(safetyreportid, nichd)]) # no uso id duplicadas porque ya se cuales corresponden a administración conjunta

# variable respuesta: EA ocurrió/no
reportes_ea <- unique(ade_raw_dt[meddra_concept_id == id_evento_adverso, safetyreportid])
datos_modelo[, ea_ocurrio := ifelse(safetyreportid %in% reportes_ea, 1, 0)]  # ":=" computo columna según expresión. %in% elementos izq en lado der

# variables de exposición 
datos_modelo[, droga_a := ifelse(safetyreportid %in% reportes_droga_a, 1, 0)]
datos_modelo[, droga_b := ifelse(safetyreportid %in% reportes_droga_b, 1, 0)]

# variable exposición factorizada
datos_modelo[, exposicion := "Ninguno"]
datos_modelo[droga_a == 1 & droga_b == 0, exposicion := paste("Solo", nombre_a)]
datos_modelo[droga_a == 0 & droga_b == 1, exposicion := paste("Solo", nombre_b)]
datos_modelo[droga_a == 1 & droga_b == 1, exposicion := paste(nombre_a, "+", nombre_b)]
niveles_exposicion <- c("Ninguno",
                        paste("Solo", nombre_a),
                        paste("Solo", nombre_b),
                        paste(nombre_a, "+", nombre_b))
datos_modelo[, exposicion := factor(exposicion, levels = niveles_exposicion)] # paso a factor e indico niveles

# nichd como factor ordenado
niveles_nichd <- c("term_neonatal", "infancy", "toddler", "early_childhood",
                   "middle_childhood", "early_adolescence", "late_adolescence")
datos_modelo[, nichd_ord := factor(nichd, levels = niveles_nichd, ordered = TRUE)]
datos_modelo[, nichd_num := as.numeric(nichd_ord)]

# resumen de candidato
print(datos_modelo[, .N, by = .(exposicion, ea_ocurrio)])


# ========================================================================
# MODELO DIFERENCIAL (EA ~ DrogaA + DrogaB + s(nichd) + s(nichd, by=AB))
# ========================================================================

# indicador administración conjunta
datos_modelo[, droga_ab := as.integer(droga_a == 1 & droga_b == 1)]

# ajuste modelo diferencial
modelo_dif <- gam(
  ea_ocurrio ~ droga_a + droga_b +
    s(nichd_num, k = 7) +
    s(nichd_num, bs = "cs", by = droga_ab, k = 7),
  data   = datos_modelo,
  family = binomial(link = "logit"),
  method = "REML"
)

plot(modelo_dif)
summary(modelo_dif)

# ------------------------------------------------------------------------
# Predicciones por etapa y estado de exposición (None, A, B, A+B)
# ------------------------------------------------------------------------
grid_dif <- CJ(
  nichd_ord = factor(niveles_nichd, levels = niveles_nichd, ordered = TRUE),
  exposicion = factor(c("Ninguno",
                        paste("Solo", nombre_a),
                        paste("Solo", nombre_b),
                        paste(nombre_a, "+", nombre_b)),
                      levels = c("Ninguno",
                                 paste("Solo", nombre_a),
                                 paste("Solo", nombre_b),
                                 paste(nombre_a, "+", nombre_b))),
  sorted = FALSE
)
grid_dif[, nichd_num := as.numeric(nichd_ord)]
grid_dif[, `:=`(
  droga_a  = as.integer(exposicion %in% c(paste("Solo", nombre_a), paste(nombre_a, "+", nombre_b))),
  droga_b  = as.integer(exposicion %in% c(paste("Solo", nombre_b), paste(nombre_a, "+", nombre_b))),
  droga_ab = as.integer(exposicion == paste(nombre_a, "+", nombre_b))
)]

pred_dif <- predict(modelo_dif, newdata = grid_dif, type = "link", se.fit = TRUE)
z90 <- qnorm(0.95)

grid_dif[, `:=`(
  lp = pred_dif$fit,
  se = pred_dif$se.fit
)]
grid_dif[, `:=`(
  li   = lp - z90 * se,
  ls   = lp + z90 * se,
  prob = plogis(lp)
)]

# ------------------------------------------------------------------------
# GRÁFICO: log-odds por NICHD y exposición (modelo diferencial)
# ------------------------------------------------------------------------
nombre_evento_elegido <- event_dt[meddra_concept_id == id_evento_adverso, meddra_concept_name_1][1]


p_dif <- ggplot(grid_dif, aes(x = nichd_ord, y = lp, group = exposicion, color = exposicion)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = li, ymax = ls, fill = exposicion), alpha = 0.2, linetype = "dashed") +
  labs(
    title = paste("Modelo diferencial - Riesgo de", nombre_evento_elegido),
    subtitle = paste("Efectos por etapa NICHD para", nombre_a, "y", nombre_b),
    x = "Etapa de Desarrollo (NICHD)", y = "Log-Odds", color = "Exposición", fill = "Exposición"
  ) +
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
print(p_dif)

# ------------------------------------------------------------------------
# OR por exposición (vs 'Ninguno') en escala logit -> OR = exp(delta LP)
# ------------------------------------------------------------------------
ref_dif <- grid_dif[exposicion == "Ninguno", .(nichd_ord, ref_lp = lp, ref_li = li, ref_ls = ls)]
or_dif  <- merge(grid_dif[exposicion != "Ninguno"], ref_dif, by = "nichd_ord")
or_dif[, `:=`(
  OR    = exp(lp - ref_lp),
  OR_li = exp(li - ref_li),
  OR_ls = exp(ls - ref_ls)
)]
y_breaks <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.6, 2, 3, 4, 5, 10, 25, 50, 100)

p_or_dif <- ggplot(or_dif, aes(x = nichd_ord, y = OR, group = exposicion, color = exposicion)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = OR_li, ymax = OR_ls, fill = exposicion), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans = "log10", breaks = y_breaks, labels = y_breaks) +
  labs(
    title = paste("Modelo diferencial - OR de", nombre_evento_elegido),
    subtitle = "Comparado con 'Ninguno' (IC 90%)",
    x = "Etapa NICHD", y = "Odds Ratio (OR)", color = "Exposición", fill = "Exposición"
  ) +
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
print(p_or_dif)

# ------------------------------------------------------------------------
# IOR (del modelo diferencial)
# Nota: en este especificación, log(IOR) = LP_AB - LP_A - LP_B + LP_0,
# que coincide con el término suave diferencial s(nichd_num, by=droga_ab)
# ------------------------------------------------------------------------
# Pasamos a wide para combinar por etapa
w_lp <- dcast(grid_dif[, .(nichd_ord, exposicion, lp, se)], nichd_ord ~ exposicion, value.var = c("lp","se"))
col_0  <- "Ninguno"
col_A  <- paste("Solo", nombre_a)
col_B  <- paste("Solo", nombre_b)
col_AB <- paste(nombre_a, "+", nombre_b)

# log-IOR y su SE (asumiendo independencia de pred varianza en escala link).
# Para SE correcto usamos lpmatrix y covarianza más abajo (como en tu IOR previo).
logior_dif <- w_lp[[paste0("lp_", col_AB)]] - w_lp[[paste0("lp_", col_A)]] -
              w_lp[[paste0("lp_", col_B)]] + w_lp[[paste0("lp_", col_0)]]

# SE correcto vía delta multivariado usando lpmatrix
Xp_dif <- predict(modelo_dif, newdata = grid_dif, type = "lpmatrix")
Vb_dif <- vcov(modelo_dif)
cov_link_dif <- Xp_dif %*% Vb_dif %*% t(Xp_dif)  # cov entre LPs
# índices por etapa/exp para armar varianza de (AB - A - B + 0)
idx_fun <- function(et, expo) which(grid_dif$nichd_ord == et & grid_dif$exposicion == expo)
logior_se <- sapply(levels(grid_dif$nichd_ord), function(et) {
  iAB <- idx_fun(et, col_AB); iA <- idx_fun(et, col_A); iB <- idx_fun(et, col_B); i0 <- idx_fun(et, col_0)
  cvec <- rep(0, nrow(grid_dif)); cvec[c(iAB, iA, iB, i0)] <- c(1, -1, -1, 1)
  sqrt(max(as.numeric(t(cvec) %*% cov_link_dif %*% cvec), 0))
})
datos_ior_dif <- data.table(
  nichd_ord = factor(levels(grid_dif$nichd_ord), levels = niveles_nichd, ordered = TRUE),
  log_ior = as.numeric(logior_dif),
  se_log_ior = as.numeric(logior_se)
)
datos_ior_dif[, `:=`(
  log_ior_li = log_ior - z90 * se_log_ior,
  log_ior_ls = log_ior + z90 * se_log_ior)]
datos_ior_dif[, `:=`(
  ior    = exp(log_ior),
  ior_li = exp(log_ior_li),
  ior_ls = exp(log_ior_ls))]
datos_ior_dif[, `:=`(
  significativo = !(ior_li <= 1 & ior_ls >= 1)
)]

p_ior_dif <- ggplot(datos_ior_dif, aes(x = nichd_ord, y = ior, group = 1)) +
  geom_ribbon(aes(ymin = ior_li, ymax = ior_ls), alpha = 0.2) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans = "log10", breaks = y_breaks, labels = y_breaks) +
  coord_cartesian(ylim = c(0.2, 100)) +
  labs(
    title = paste("Modelo diferencial - IOR para", nombre_evento_elegido),
    subtitle = "IC 90% (línea punteada = sin interacción)",
    x = "Etapa NICHD", y = "Interaction Odds Ratio (IOR)"
  ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_ior_dif)

message("Modelo diferencial - Etapas con interacción significativa (IC90):")
print(datos_ior_dif[significativo == TRUE, .(nichd_ord, ior, ior_li, ior_ls)])

# ------------------------------------------------------------------------
# RERI (modelo diferencial) con covarianza correcta en escala prob
# ------------------------------------------------------------------------
# Probabilidades y covarianza en prob usando delta multivariante
grid_dif[, prob := plogis(lp)]
D_dif <- diag(as.numeric(grid_dif$prob * (1 - grid_dif$prob)), nrow = nrow(grid_dif))
cov_prob_dif <- D_dif %*% cov_link_dif %*% D_dif

idx_row_dif <- function(nichd, expo) which(grid_dif$nichd_ord == nichd & grid_dif$exposicion == expo)

res_reri_dif <- data.table(
  nichd_ord = character(), p0 = numeric(), pA = numeric(), pB = numeric(), pAB = numeric(),
  RERI = numeric(), RERI_li = numeric(), RERI_ls = numeric()
)
for (et in niveles_nichd) {
  i0  <- idx_row_dif(et, col_0); iA <- idx_row_dif(et, col_A)
  iB  <- idx_row_dif(et, col_B); iAB <- idx_row_dif(et, col_AB)
  if (length(i0)!=1 || length(iA)!=1 || length(iB)!=1 || length(iAB)!=1) next

  p0  <- grid_dif$prob[i0]; pA <- grid_dif$prob[iA]; pB <- grid_dif$prob[iB]; pAB <- grid_dif$prob[iAB]
  idxs <- c(iAB, iA, iB, i0)
  cov_sub <- cov_prob_dif[idxs, idxs]
  cvec <- c(1, -1, -1, 1)
  se_RERI <- sqrt(max(as.numeric(t(cvec) %*% cov_sub %*% cvec), 0))
  RERI_val <- pAB - pA - pB + p0
  res_reri_dif <- rbind(res_reri_dif, list(et, p0, pA, pB, pAB,
                                           RERI_val, RERI_val - z90 * se_RERI, RERI_val + z90 * se_RERI))
}
setnames(res_reri_dif, c("nichd_ord","p0","pA","pB","pAB","RERI","RERI_li","RERI_ls"))
res_reri_dif[, nichd_ord := factor(nichd_ord, levels = niveles_nichd, ordered = TRUE)]

p_reri_dif <- ggplot(res_reri_dif, aes(x = nichd_ord, y = RERI, group = 1)) +
  geom_ribbon(aes(ymin = RERI_li, ymax = RERI_ls), alpha = 0.2) +
  geom_line(size = 1) + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = paste("Modelo diferencial - RERI para", nombre_evento_elegido),
    x = "Etapa NICHD", y = "RERI (IC 90%)",
    caption = "RERI > 0: sinergia aditiva; < 0: antagonismo aditivo"
  ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_reri_dif)

# ------------------------------------------------------------------------
# TABLAS GT (modelo diferencial)
# ------------------------------------------------------------------------
s_dif <- summary(modelo_dif)

df_metricas_dif <- tibble::tibble(
  Métrica = c("AIC", "Deviance", "Dev. explicado", "EDF total"),
  Valor   = c(round(AIC(modelo_dif), 3),
              round(deviance(modelo_dif), 3),
              round(s_dif$dev.expl, 3),
              round(sum(s_dif$edf), 3))
)

P_dif <- s_dif$p.table
df_parametricos_dif <- tibble::tibble(
  Tipo         = "Paramétrico",
  Término      = rownames(P_dif),
  Estimación   = round(P_dif[,"Estimate"], 3),
  `Std. Error` = round(P_dif[,"Std. Error"], 3),
  `p-valor`    = round(P_dif[,"Pr(>|z|)"], 3)
)

S_dif <- s_dif$s.table
df_suavizados_dif <- tibble::tibble(
  Tipo     = "Suavizado",
  Término  = rownames(S_dif),
  EDF      = round(S_dif[,"edf"], 3),
  `Ref.df` = round(S_dif[,"Ref.df"], 3),
  `p-valor`= round(S_dif[,"p-value"], 3)
)

df_coefs_dif <- bind_rows(df_parametricos_dif, df_suavizados_dif)

df_ior_dif <- datos_ior_dif %>%
  transmute(
    Etapa        = nichd_ord,
    IOR          = round(ior, 3),
    `IC 90% inf` = round(ior_li, 3),
    `IC 90% sup` = round(ior_ls, 3),
    Signif       = ifelse(significativo, "Sí", "No")
  )

df_reri_dif <- res_reri_dif %>%
  transmute(
    Etapa        = nichd_ord,
    RERI         = round(RERI, 3),
    `IC 90% inf` = round(RERI_li, 3),
    `IC 90% sup` = round(RERI_ls, 3)
  )

tbl_metricas_dif <- gt(df_metricas_dif) %>% tab_header("Modelo diferencial - Métricas")
tbl_coefs_dif    <- gt(df_coefs_dif)    %>% tab_header("Modelo diferencial - Coeficientes GAM") %>%
  cols_label(
    Tipo="Tipo", Término="Término", Estimación="Estimación",
    `Std. Error`="Error Estándar", EDF="EDF", `Ref.df`="Ref.df", `p-valor`="p-valor"
  )
tbl_ior_dif  <- gt(df_ior_dif)  %>% tab_header("Modelo diferencial - IOR por NICHD (IC 90%)")
tbl_reri_dif <- gt(df_reri_dif) %>% tab_header("Modelo diferencial - RERI por NICHD (IC 90%)")

print(tbl_metricas_dif)
print(tbl_coefs_dif)
print(tbl_ior_dif)
print(tbl_reri_dif)
# ========================================================================

# 1. Crear un nombre de archivo seguro y descriptivo
nombre_archivo_base <- paste(
    "Reporte_Diferencial",
    gsub("[^A-Za-z0-9]", "", nombre_a),
    gsub("[^A-Za-z0-9]", "", nombre_b),
    gsub("[^A-Za-z0-9]", "_", strtrim(nombre_evento_elegido, 40)),
    format(Sys.Date(), "%Y%m%d"),
    sep = "_"
)
nombre_archivo_html <- paste0(nombre_archivo_base, ".html")

# 2. Preparar datos para el reporte (que faltaban en script 2)
tabla_resumen_datos <- datos_modelo[, .N, by = .(exposicion, ea_ocurrio)]
tabla_resumen_completa <- datos_modelo[, .N, by = .(exposicion)]
datos_ior_signif_dif <- datos_ior_dif[significativo == TRUE]

# 3. Convertir tablas gt a data.frames para kableExtra
df_metricas_dif_tabla <- data.frame(
    Métrica = c("AIC", "Deviance", "Dev. explicado", "EDF total"),
    Valor = c(
        round(AIC(modelo_dif), 3),
        round(deviance(modelo_dif), 3),
        round(s_dif$dev.expl, 3),
        round(sum(s_dif$edf), 3)
    )
)

# Tabla de coeficientes (reutilizando el código del script 1)
P_dif <- s_dif$p.table
df_parametricos_dif_tabla <- data.frame(
    Tipo = "Paramétrico",
    Término = rownames(P_dif),
    Estimación = round(P_dif[, "Estimate"], 3),
    `Std. Error` = round(P_dif[, "Std. Error"], 3),
    `p-valor` = round(P_dif[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE, check.names = FALSE
)

S_dif <- s_dif$s.table
df_suavizados_dif_tabla <- data.frame(
    Tipo = "Suavizado",
    Término = rownames(S_dif),
    EDF = round(S_dif[, "edf"], 3),
    `Ref.df` = round(S_dif[, "Ref.df"], 3),
    `p-valor` = round(S_dif[, "p-value"], 3),
    stringsAsFactors = FALSE, check.names = FALSE
)

# Usar dplyr::bind_rows para combinar data.frames con columnas diferentes
df_coefs_dif_completo <- dplyr::bind_rows(df_parametricos_dif_tabla, df_suavizados_dif_tabla)

df_ior_tabla_dif <- datos_ior_dif %>%
    transmute(
        Etapa = nichd_ord,
        IOR = round(ior, 3),
        `IC 90% inf` = round(ior_li, 3),
        `IC 90% sup` = round(ior_ls, 3),
        Signif = ifelse(significativo, "Sí", "No")
    )

df_reri_tabla_dif <- res_reri_dif %>%
    transmute(
        Etapa = factor(nichd_ord, levels = niveles_nichd),
        RERI = round(RERI, 3),
        `IC 90% inf` = round(RERI_li, 3),
        `IC 90% sup` = round(RERI_ls, 3)
    ) %>%
    arrange(Etapa)

# 4. Definir la plantilla R Markdown COMPLETA para HTML
plantilla_vector <- c(
'---',
'title: "Análisis de Interacción Farmacológica mediante Modelo GAM Diferencial"',
'subtitle: "`r params$nombre_a` + `r params$nombre_b` - `r params$nombre_evento`"',
'author: "Autor: Zacañino Ezequiel"',
'date: "`r format(Sys.Date(), \'%d de %B de %Y\')`"',
'output:',
'  html_document:',
'    toc: true',
'    toc_float: true',
'    number_sections: true',
'    code_folding: "hide"',
'    theme: journal',
'    fig_width: 8',
'    fig_height: 6',
'params:',
'  nombre_a: ""',
'  nombre_b: ""',
'  nombre_evento: ""',
'  id_droga_a: 0',
'  id_droga_b: 0',
'  id_evento_adverso: 0',
'  reportes_conjuntos_n: 0',
'  tabla_resumen_completa: NULL',
'  tabla_resumen_datos: NULL',
'  datos_ior_signif: NULL',
'  p_log_odds: NULL',
'  p_or: NULL',
'  p_ior: NULL',
'  p_reri: NULL',
'  df_metricas: NULL',
'  df_coefs: NULL',
'  df_ior_tabla: NULL',
'  df_reri_tabla: NULL',
'---',
'',
'```{r setup, include=FALSE}',
'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align = "center")',
'library(knitr)',
'library(kableExtra)',
'library(dplyr)',
'library(tidyr)',
'```',
'',
'# Resumen de análisis (Modelo Diferencial):',
'',
'Este reporte muestra los resultados del análisis de interacción farmacológica entre **`r params$nombre_a`** y **`r params$nombre_b`** en relación al evento adverso **`r params$nombre_evento`**, utilizando un modelo GAM **diferencial** que permite evaluar efectos específicos de la interacción por etapa ontogénica.',
'',
'## Especificación del Modelo Diferencial',
'',
'El modelo diferencial utilizado tiene la forma:',
'```',
'EA ~ DrogaA + DrogaB + s(nichd) + s(nichd, by = AB)',
'```',
'',
'donde:',
'- `s(nichd)` captura el efecto basal por edad',
'- `s(nichd, by = AB)` captura el efecto **diferencial** específico de la administración conjunta',
'',
'## Hallazgos Principales',
'',
'- **Fármacos analizados**: `r params$nombre_a` (ID: `r params$id_droga_a`) y `r params$nombre_b` (ID: `r params$id_droga_b`)',
'- **Evento adverso**: `r params$nombre_evento` (ID: `r params$id_evento_adverso`)',
'- **Reportes con administración conjunta**: `r params$reportes_conjuntos_n`',
'- **Total de reportes analizados**: `r sum(params$tabla_resumen_completa$N)`',
'',
'```{r significancia-summary, results="asis"}',
'if(nrow(params$datos_ior_signif) > 0) {',
'  cat("- **Etapas con interacción estadísticamente significativa (IC 90%)**: ")',
'  cat(paste(params$datos_ior_signif$nichd_ord, collapse = ", "))',
'} else {',
'  cat("- **No se encontraron interacciones estadísticamente significativas** en ninguna etapa NICHD.")',
'}',
'```',
'',
'# Descripción de los Datos',
'',
'## Distribución de Reportes por Grupo de Exposición',
'```{r tabla-resumen-completa}',
'kable(params$tabla_resumen_completa, caption = "Distribución total de reportes por grupo de exposición.", col.names = c("Grupo de Exposición", "N° de Reportes")) %>% kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)',
'```',
'',
'## Ocurrencia del Evento Adverso por Grupo',
'```{r tabla-resumen-evento}',
'params$tabla_resumen_datos %>%',
'  mutate(evento_label = ifelse(ea_ocurrio == 1, "Evento Ocurrió", "Sin Evento")) %>%',
'  select(-ea_ocurrio) %>%',
'  tidyr::pivot_wider(names_from = evento_label, values_from = N, values_fill = 0) %>%',
'  mutate(Total = `Sin Evento` + `Evento Ocurrió`, `Tasa (%)` = round((`Evento Ocurrió` / Total) * 100, 2)) %>%',
'  kable(caption = "Incidencia del evento adverso por grupo de exposición.") %>%',
'  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%',
'  add_header_above(c(" " = 1, "Reportes" = 2, "Resumen" = 2))',
'```',
'',
'# Resultados del Modelo GAM Diferencial',
'',
'## Métricas de Ajuste del Modelo',
'```{r metricas-gam}',
'kable(params$df_metricas, caption = "Métricas globales del modelo GAM diferencial ajustado.") %>% kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)',
'```',
'',
'## Coeficientes del Modelo',
'```{r coeficientes-gam}',
'kable(params$df_coefs, caption = "Coeficientes paramétricos y términos suavizados del modelo GAM diferencial.") %>%',
'  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE) %>%',
'  pack_rows("Términos Paramétricos", 1, sum(params$df_coefs$Tipo == "Paramétrico")) %>%',
'  pack_rows("Términos Suavizados", sum(params$df_coefs$Tipo == "Paramétrico") + 1, nrow(params$df_coefs))',
'```',
'',
'# Análisis Visual de la Interacción',
'## Evolución del Riesgo por Etapa de Desarrollo (Log-Odds)',
'```{r grafico-logodds, fig.cap="Log-odds del evento adverso por etapa NICHD y grupo de exposición (modelo diferencial). Las bandas representan intervalos de confianza al 90%. Se observa cómo el riesgo basal varía entre grupos y a lo largo del desarrollo."}',
'print(params$p_log_odds)',
'```',
'Este gráfico muestra la evolución del riesgo (en escala log-odds) del evento adverso a lo largo de las diferentes etapas del desarrollo según NICHD usando el **modelo diferencial**. Las diferencias entre las curvas sugieren efectos diferenciales de los tratamientos según la edad.',
'',
'## Odds Ratios por Exposición',
'```{r grafico-or, fig.cap="Odds Ratios por etapa NICHD y exposición comparado con el grupo \'Ninguno\'. El modelo diferencial permite identificar cambios específicos en el riesgo asociado a la administración conjunta."}',
'print(params$p_or)',
'```',
'',
'## Interaction Odds Ratios (IOR)',
'```{r grafico-ior, fig.cap="Evolución del Interaction Odds Ratio (IOR) según modelo diferencial. La escala logarítmica permite visualizar tanto efectos sinérgicos (> 1) como antagónicos (< 1). La línea punteada horizontal marca la ausencia de interacción (IOR = 1)."}',
'print(params$p_ior)',
'```',
'',
'```{r tabla-ior}',
'kable(params$df_ior_tabla, caption = "Interaction Odds Ratio (IOR) por etapa NICHD con IC 90% - Modelo Diferencial.") %>%',
'  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%',
'  column_spec(which(colnames(params$df_ior_tabla) == "Signif"), bold = TRUE) %>%',
'  row_spec(which(params$df_ior_tabla$Signif == "Sí"), background = "#fcf8e3")',
'```',
'',
'## Exceso de Riesgo Relativo por Interacción (RERI)',
'```{r grafico-reri, fig.cap="Exceso de Riesgo Relativo por Interacción (RERI) según modelo diferencial. Valores positivos indican sinergia aditiva, valores negativos antagonismo aditivo. La línea en cero marca la ausencia de interacción aditiva."}',
'print(params$p_reri)',
'```',
'',
'```{r tabla-reri}',
'kable(params$df_reri_tabla, caption = "Exceso de Riesgo Relativo por Interacción (RERI) por etapa NICHD con IC 90% - Modelo Diferencial.") %>%',
'  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)',
'```',
'',
'# Interpretación de Resultados',
'',
'## Ventajas del Modelo Diferencial',
'',
'El modelo diferencial presenta ventajas específicas:',
'',
'1. **Eficiencia paramétrica**: Al modelar directamente el efecto diferencial de la interacción, puede ser más eficiente estadísticamente.',
'2. **Interpretación directa**: El término `s(nichd, by = AB)` representa directamente la variación del efecto de interacción por edad.',
'3. **Menor número de parámetros**: Comparado con modelos factorizados, puede requerir menos parámetros para capturar la interacción.',
'',
'## Comparación con Enfoques Alternativos',
'',
'Este modelo diferencial permite una interpretación complementaria a los modelos factorizados tradicionales, especialmente útil cuando se sospecha que el efecto de interacción tiene una estructura suave específica a lo largo del desarrollo ontogénico.',
''
)

# 5. Escribir plantilla a archivo temporal y renderizar
plantilla_rmd <- paste(plantilla_vector, collapse = "\n")
nombre_plantilla <- "reporte_diferencial_completo_template.Rmd"
writeLines(plantilla_rmd, nombre_plantilla)

message("Generando reporte HTML completo (modelo diferencial)...")

tryCatch({
    rmarkdown::render(
        input = nombre_plantilla,
        output_format = "html_document",
        output_file = nombre_archivo_html,
        params = list(
            nombre_a = nombre_a,
            nombre_b = nombre_b,
            nombre_evento = nombre_evento_elegido,
            id_droga_a = id_droga_a,
            id_droga_b = id_droga_b,
            id_evento_adverso = id_evento_adverso,
            reportes_conjuntos_n = length(reportes_conjuntos),
            tabla_resumen_completa = tabla_resumen_completa,
            tabla_resumen_datos = tabla_resumen_datos,
            datos_ior_signif = datos_ior_signif_dif,
            p_log_odds = p_dif,
            p_or = p_or_dif,
            p_ior = p_ior_dif,
            p_reri = p_reri_dif,
            df_metricas = df_metricas_dif_tabla,
            df_coefs = df_coefs_dif_completo,
            df_ior_tabla = df_ior_tabla_dif,
            df_reri_tabla = df_reri_tabla_dif
        ),
        quiet = TRUE # Poner en FALSE para ver el log de pandoc
    )
    
    message(paste("\n\u2713 Reporte HTML completo generado exitosamente:", nombre_archivo_html))
    message(paste("   Ubicación:", file.path(getwd(), nombre_archivo_html)))
    if (file.exists(nombre_archivo_html)) {
        size_mb <- round(file.info(nombre_archivo_html)$size / 1024^2, 2)
        message(paste("   Tamaño:", size_mb, "MB"))
        # Opcional: abrir el archivo automáticamente
        # utils::browseURL(nombre_archivo_html)
    }

}, error = function(e) {
    message("\n\u2717 ERROR AL GENERAR EL REPORTE HTML COMPLETO")
    message(paste(rep("=", 50), collapse = ""))
    message("Ocurrió un error durante la renderización del archivo R Markdown.")
    message("Revisa la consola en busca de mensajes de error de 'knitr' o 'pandoc'.")
    message(paste(rep("=", 50), collapse = ""))
    message("\nDetalle del error:")
    message(e$message)
    message("\nPuedes revisar los gráficos y tablas generados en R mientras tanto.")
})

# 6. Limpiar archivos temporales
if (file.exists(nombre_plantilla)) {
    file.remove(nombre_plantilla)
    message("\n\u2713 Archivo temporal de plantilla eliminado")
}



p_dif <- ggplot(grid_dif, aes(x = factor(nichd_ord,
                                         levels = c("term_neonatal", "infancy", "toddler",
                                                    "early_childhood", "middle_childhood",
                                                    "early_adolescence", "late_adolescence"),
                                         labels = c("Neonato", "Lactante", "Niñez temprana",
                                                    "Primera infancia", "Segunda infancia",
                                                    "Adolescencia temprana", "Adolescencia tardía")),
                               y = lp, group = exposicion, color = exposicion)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = li, ymax = ls, fill = exposicion), alpha = 0.2, linetype = "dashed") +
  labs(
    title = paste("Modelo diferencial - Riesgo de", nombre_evento_elegido),
    subtitle = paste("Efectos por etapa NICHD para", nombre_a, "y", nombre_b),
    x = "Etapa de desarrollo (NICHD)", 
    y = "Log-Odds", 
    color = "Exposición", 
    fill = "Exposición"
  ) +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 18) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p_dif)
