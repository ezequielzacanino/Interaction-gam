library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

################################################################################
# Parámetros 
################################################################################
set.seed(9427)
n_reports <- 150000
n_stages  <- 7

# Riesgos basales DESADOS (verdadero mecanismo)
# Configuración: A y B tienen efecto moderado, juntos tienen efecto ALTO (sinergia)
pA_real  <- 0.05      # Riesgo droga A sola
pB_real  <- 0.05     # Riesgo droga B sola  
pAB_real <- 0.01      # Riesgo conjunto
pC_real  <- 0.04      # Riesgo baseline población general

cat(sprintf("Configuración real:\n"))
cat(sprintf("p_A = %.3f, p_B = %.3f, p_AB = %.3f\n", pA_real, pB_real, pAB_real))
cat(sprintf("Suma individual: %.3f | Exceso de sinergia: %.3f\n", 
            pA_real + pB_real, pAB_real - (pA_real + pB_real)))

################################################################################
# 1- Simulación dataset base
################################################################################
dt <- data.table(
  id    = 1:n_reports,
  stage = sample(1:n_stages, n_reports, replace = TRUE),
  A     = rbinom(n_reports, 1, 0.12),
  B     = rbinom(n_reports, 1, 0.10)
)

# Asignar p_base verdadero por combinación A/B
dt[, p_base := fifelse(A==1 & B==1, pAB_real,
                fifelse(A==1 & B==0, pA_real,
                fifelse(A==0 & B==1, pB_real, pC_real)))]
dt[, E := rbinom(.N, 1, p_base)]

# Estimados empíricos (con muestreo)
pA  <- dt[A==1 & B==0, mean(E)]
pB  <- dt[A==0 & B==1, mean(E)]
pAB <- dt[A==1 & B==1, mean(E)]
p0  <- dt[A==0 & B==0, mean(E)]
p_global <- dt[, mean(E)]

cat(sprintf("\nEstimados empíricos del dataset:\n"))
cat(sprintf("p_A=%.4f, p_B=%.4f, p_AB=%.4f, p_0=%.4f\n", pA, pB, pAB, p0))

################################################################################
# 2- dinámica y fold-change
################################################################################
f <- tanh(seq(-pi, pi, length.out = n_stages))
dyn <- (f - min(f)) / (max(f) - min(f))
names(dyn) <- as.character(1:n_stages)

sample_fold_change <- function(n, lambda = 0.75) 1 + rexp(n, rate = lambda)
fc_global <- sample_fold_change(1, lambda = 0.75)[1]

methods <- c("Global", "Multiplicativo", "Aditivo", "Max_Relativo")

################################################################################
# 3- inyección con métricas de interacción
################################################################################

inject_and_analyze <- function(name, dt_input, fold_change) {
  dt2 <- copy(dt_input)
  dt2[, p_new := NA_real_]
  
  # Calcular Pbase según método
  if (name == "Global") {
    # Original: usa p_global (riesgo poblacional)
    Pbase <- p_global
    formula_desc <- "p_global"
  } else if (name == "Multiplicativo") {
    # Independencia multiplicativa: OR_AB = OR_A * OR_B
    # En escala de probabilidad: p = (pA*pB/p0) cappeado
    Pbase <- (pA * pB) / p0
    formula_desc <- "(pA*pB)/p0"
  } else if (name == "Aditivo") {
    # Independencia aditiva: p_AB = p_A + p_B - p_A*p_B
    Pbase <- pA + pB - (pA * pB)
    formula_desc <- "pA + pB - pA*pB"
  } else if (name == "Max_Relativo") {
    # Dominancia: mejor efecto individual
    Pbase <- max(pA, pB)
    formula_desc <- "max(pA, pB)"
  }
  
  # Inyección en coadministración (A=1, B=1)
  sel <- which(dt2$A == 1 & dt2$B == 1)
  if (length(sel) > 0) {
    stages_sel <- as.character(dt2$stage[sel])
    dyn_sel <- dyn[stages_sel]
    
    # Probabilidad dinámica: Pbase * FC * (1 + dinámica normalizada)
    rhs <- (Pbase * fold_change) * (1 + dyn_sel)
    rhs <- pmin(pmax(rhs, 0.001), 0.999)
    dt2$p_new[sel] <- rhs
    
    # Simular eventos
    dt2$E[sel] <- rbinom(length(sel), 1, rhs)
  }
  
  ###########
  # Métricas de interacción por etapa
  ###########
  results_by_stage <- dt2[, {
    # Proporciones observadas post-inyección
    p_A_obs  <- mean(E[A==1 & B==0])
    p_B_obs  <- mean(E[A==0 & B==1])
    p_AB_obs <- mean(E[A==1 & B==1])
    p_0_obs  <- mean(E[A==0 & B==0])
    
    # RERI (Relative Excess Risk due to Interaction)
    # RERI = p_AB - p_A - p_B + p_0
    # > 0: Sinergia aditiva, = 0: Aditivo, < 0: Sub-aditivo/antagonismo
    reri <- p_AB_obs - p_A_obs - p_B_obs + p_0_obs
    
    # IOR (Interaction Odds Ratio) - escala multiplicativa
    # IOR = (OR_AB * OR_0) / (OR_A * OR_B)
    # log(IOR) = logit(p_AB) - logit(p_A) - logit(p_B) + logit(p_0)
    safe_logit <- function(p) {
      p <- pmin(pmax(p, 0.001), 0.999)
      log(p/(1-p))
    }
    
    logit_AB <- safe_logit(p_AB_obs)
    logit_A  <- safe_logit(p_A_obs)
    logit_B  <- safe_logit(p_B_obs)
    logit_0  <- safe_logit(p_0_obs)
    
    log_ior <- logit_AB - logit_A - logit_B + logit_0
    ior <- exp(log_ior)
    
    # Clasificación
    interaccion_tipo <- ifelse(
      reri > 0.01 & ior > 1.1, "Sinérgica (multiplicativa)",
      ifelse(reri > 0.01, "Sinérgica Aditiva",
             ifelse(abs(reri) <= 0.01 & abs(ior-1) <= 0.1, "Aditiva/Independiente",
                    ifelse(reri < -0.01, "Sub-aditiva/Antagónica", "Indeterminada")))
    )
    
    .(
      p_A = p_A_obs, p_B = p_B_obs, p_AB = p_AB_obs, p_0 = p_0_obs,
      Pbase_used = Pbase,
      reri = reri,
      log_ior = log_ior,
      ior = ior,
      interaccion_tipo = interaccion_tipo[1],
      n_coadmin = sum(A==1 & B==1)
    )
  }, by = stage][order(stage)]
  
  results_by_stage[, method := name]
  results_by_stage[, formula := formula_desc]
  return(results_by_stage)
}

# -----------------------
# 4) Ejecutar para todos los métodos
# -----------------------
res_list <- lapply(methods, inject_and_analyze, dt_input = dt, fold_change = fc_global)
res <- rbindlist(res_list)

# -----------------------
# 5) Resumen comparativo
# -----------------------
cat(sprintf("\n\n=== COMPARACIÓN DE MÉTODOS (Fold-change = %.2f) ===\n", fc_global))

summary_table <- res[, .(
  Pbase = mean(Pbase_used),
  p_AB_mean = mean(p_AB),
  RERI_mean = mean(reri),
  IOR_mean = mean(ior),
  Tipo_Interaccion = interaccion_tipo[1],
  Clasificacion = case_when(
    mean(reri) > 0.005 & mean(ior) > 1.05 ~ "✓ SINÉRGICO",
    abs(mean(reri)) <= 0.005 & abs(mean(ior)-1) <= 0.05 ~ "○ ADITIVO (Independiente)",
    mean(reri) < -0.005 ~ "✗ SUB-ADITIVO",
    TRUE ~ "~ MIXTO"
  )
), by = .(method, formula)]

print(summary_table)

################################################################################
# 6- Visualizaciones
################################################################################

# Paleta de colores por tipo de interacción
colores_tipo <- c(
  "Sinérgica (multiplicativa)" = "#E74C3C",
  "Sinérgica Aditiva" = "#E67E22", 
  "Aditiva/Independiente" = "#3498DB",
  "Sub-aditiva/Antagónica" = "#9B59B6",
  "Indeterminada" = "#95A5A6"
)

# Gráfico 1: Probabilidades por método y etapa
gg_probs <- ggplot(res, aes(x = factor(stage), y = p_AB, color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = pA + pB, linetype = "dashed", color = "gray50", 
             alpha = 0.7, size = 1) +
  geom_hline(yintercept = pAB_real, linetype = "solid", color = "darkgreen", 
             alpha = 0.7, size = 1) +
  annotate("text", x = 1, y = pA + pB + 0.01, label = "Umbral Aditivo (pA+pB)", 
           hjust = 0, color = "gray50", size = 3) +
  annotate("text", x = 1, y = pAB_real + 0.01, label = "Verdadero p_AB (sinérgico)", 
           hjust = 0, color = "darkgreen", size = 3) +
  facet_wrap(~method, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(title = "Probabilidad Observada P(E|A,B) por Método",
       subtitle = sprintf("Fold-change = %.2f | Línea verde = efecto real sinérgico | Línea gris = aditividad pura",
                         fc_global),
       x = "Etapa NICHD", y = "Probabilidad p(E|A,B)",
       color = "Método") +
  theme(legend.position = "none")

print(gg_probs)

# Gráfico 2: RERI (escala aditiva)
gg_reri <- ggplot(res, aes(x = factor(stage), y = reri, fill = interaccion_tipo)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_hline(yintercept = pAB_real - (pA + pB), color = "darkgreen", 
             linetype = "dashed", size = 1) +
  facet_wrap(~method, ncol = 2) +
  scale_fill_manual(values = colores_tipo) +
  theme_bw(base_size = 12) +
  labs(title = "RERI: Exceso de Riesgo por Interacción (Escala Aditiva)",
       subtitle = "RERI = p_AB - p_A - p_B + p_0 | >0: Sinergia | =0: Aditivo | <0: Antagonismo",
       y = "RERI", x = "Etapa", fill = "Tipo") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, 0.1))

print(gg_reri)

# Gráfico 3: Log-IOR (escala multiplicativa)
gg_ior <- ggplot(res, aes(x = factor(stage), y = log_ior, fill = interaccion_tipo)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  # Calcular log-IOR verdadero
  geom_hline(yintercept = log((pAB_real/(1-pAB_real)) / 
                               ((pA/(1-pA)) * (pB/(1-pB)) / (p0/(1-p0)))),
             color = "darkgreen", linetype = "dashed", size = 1) +
  facet_wrap(~method, ncol = 2) +
  scale_fill_manual(values = colores_tipo) +
  theme_bw(base_size = 12) +
  labs(title = "Log-IOR: Medida de Interacción (Escala Multiplicativa)",
       subtitle = "Log(IOR) = logit(p_AB) - logit(p_A) - logit(p_B) + logit(p_0) | >0: Sinergia",
       y = "Log(IOR)", x = "Etapa", fill = "Tipo") +
  theme(legend.position = "bottom")

print(gg_ior)

# Gráfico 4: Comparación directa IOR vs RERI (scatter por método)
# Promediar por método para simplificar
res_avg <- res[, .(
  reri = mean(reri),
  ior = mean(ior),
  p_AB = mean(p_AB)
), by = method]

gg_compare <- ggplot(res_avg, aes(x = reri, y = ior-1, color = method, size = p_AB)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           fill = "green", alpha = 0.1) + # Cuadrante sinérgico
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           fill = "red", alpha = 0.1) +   # Cuadrante sub-aditivo
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(title = "Clasificación de Interacción: RERI vs IOR",
       subtitle = "Cuadrante superior-derecho (verde): Sinergia | Inferior-izquierdo (rojo): Sub-aditividad",
       x = "RERI (escala aditiva)", y = "IOR - 1 (escala multiplicativa)",
       color = "Método", size = "p(E|A,B)") +
  theme(legend.position = "right")

print(gg_compare)

################################################################################
# 7- Análisis de sensibilidad a fold-change
################################################################################
cat("\n\n=== ANÁLISIS DE SENSIBILIDAD A FOLD-CHANGE ===\n")

fc_values <- seq(1.0, 3.0, by = 0.5)
sensitivity_results <- list()

for(fc in fc_values) {
  temp_res <- rbindlist(lapply(methods, inject_and_analyze, dt_input = dt, fold_change = fc))
  temp_summary <- temp_res[, .(
    reri = mean(reri),
    ior = mean(ior),
    sinergia_aditiva = mean(reri > 0.01),
    sinergia_mult = mean(ior > 1.1)
  ), by = method]
  temp_summary[, fc := fc]
  sensitivity_results[[as.character(fc)]] <- temp_summary
}

sens_df <- rbindlist(sensitivity_results)

# Gráfico de sensibilidad
gg_sens <- ggplot(sens_df, aes(x = fc, y = reri, color = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = pAB_real - (pA + pB), color = "darkgreen", 
             linetype = "dotted", size = 1) +
  scale_x_continuous(breaks = fc_values) +
  theme_bw(base_size = 12) +
  labs(title = "Sensibilidad del RERI al Fold-Change",
       subtitle = "Línea punteada verde: RERI objetivo (sinergia real)",
       x = "Fold-change", y = "RERI promedio", color = "Método")

print(gg_sens)

# Tabla final resumen
final_table <- sens_df[, .(
  `RERI@FC=1.0` = reri[fc==1.0],
  `RERI@FC=2.0` = reri[fc==2.0],
  `Tipo@FC=1.0` = ifelse(reri[fc==1.0] > 0.01, "Sinérgico", 
                         ifelse(abs(reri[fc==1.0])<0.01, "Aditivo", "Sub-aditivo")),
  `Tipo@FC=2.0` = ifelse(reri[fc==2.0] > 0.01, "Sinérgico", 
                         ifelse(abs(reri[fc==2.0])<0.01, "Aditivo", "Sub-aditivo"))
), by = method]

cat("\nTabla resumen por método y fold-change:\n")
print(final_table)

# Combinar plots
composite_plot <- (gg_probs / gg_reri) | (gg_ior / gg_compare)
composite_plot + plot_annotation(
  title = 'Análisis Comparativo de Métodos de Inyección: Aditivo vs Sinérgico',
  subtitle = sprintf('Configuración: p_A=%.2f, p_B=%.2f, p_AB_real=%.2f (Sinérgico)', 
                     pA_real, pB_real, pAB_real),
  theme = theme(plot.title = element_text(size = 16, face = 'bold'))
)