# Modelado de la Dinámica Ontogénica de Interacciones Farmacológicas en Pediatría

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)


## Descripción

Este repositorio contiene los scripts utilizados para la validación de modelos aditivos generalizados (GAM) para la detección de señales de interacción farmacológica con modelado de dinámicas ontogénicas, en sistemas de reporte espontáneo de efectos adversos en pediatría.

## Conceptos

- **Triplete**: Combinación de DrogaA + DrogaB + Evento adverso
- **Modulación ontogénica**: Variaciones en el odds de efecto adverso según etapa del desarrollo (NICHD stages)
- **Dinámica**: Patrón temporal de la modulación (uniforme, aumento, disminución, meseta, meseta inversa)
- **IOR (Interaction Odds Ratio)**: Métrica que cuantifica la interacción sinérgica o antagónica entre drogas
- **RERI (Relative Excess Risk Due to Interaction)**: Métrica que cuantifica la interacción sinérgica o antagónica entre drogas
- **Métodos estratificados**: Se refiere al cálculo aritmético de las medidas de interacción IOR/RERI estratificando dataset por etapas

## Workflow

```
00_functions.R               # Funciones del pipeline
01_theme.R                   # Configuración para gráficos
02_descriptive.R             # Análisis descriptivo del dataset
10_augmentation.R            # Generación de datos y resultados
11_augmentation_base.R       # Generación de datos y resultados sin reducción (ligero)
20_null.R                    # Construcción de distribución nula
30_metrics.R                 # Análisis de resultados
31_metrics_base.R            # Análisis de resultados sin reducción (ligero)
32_metrics_facet.R           # Generación de gráficos facetados
40_network.R                 # Análisis de redes
```

## Metodología

### 1. Generación de Datos Semisintéticos (`10_augmentation.R`)

**Controles Positivos:**
- Armado de tripletes con requisitos mínimos (trabajo final requiere al menos 2 reportes A-B con al menos 2 etapas con presencia de reporte A-B)
- Selección aleatoria de 500 tripletes candidatos
- Inyección de señal con dinámicas ontogénicas simuladas (increase, decrease, plateau, inverse-plateau, uniform)
- Asignación de fold-changes muestreados de distribución exponencial (λ = 0.75)
- Ajustado de positivos inyectados en copia independiente de dataset original 
- 5 sets de 500 tripletes positivos (se excluye uniform de analisis final) 

**Controles Negativos:**
- Selección aleatoria de 5000 tripletes no inyectados (mutuamente exclusivos con positivos) con características de reporte similares a positivos
- Ajustado de negativos 

**Fórmula de inyección:**
```
P(evento|etapa) = fold_change × e_j + dinámica(etapa)

siendo e_j = P(A) + P(B) - P(A)×P(B)  
```

### 2. Distribución Nula (`20_null.R`)

- Pool aleatorio de 100,000 reportes
- Permutación estratificada por etapa de desarrollo y droga (rompe asociación droga-evento)
- Generación de 10000 tripletes permutados
- Cálculo de percentiles (P90, P95, P99) del IC90 inferior por etapa

### 3. Criterio de Detección 

Una señal es **positiva** si al menos una etapa cumple:
1. **Criterio nominal**: IC90 inferior del log-IOR > 0 (único criterio utilizado en métodos estratificados)
2. **Criterio null model**: IC90 inferior del log-IOR > Percentil de distribución nula (seleccionable, utilizado en análisis final)

### 4. Modelo GAM

**Fórmula de base:**
```r
evento ~ s(nichd_num, by=drugA, bs="cs", k=7) +
         s(nichd_num, by=drugB, bs="cs", k=7) +
         s(nichd_num, by=drugA_drugB, bs="cs", k=7) 
         
```

**Métricas de interacción:**
```
log(IOR) = lp₁₁ - lp₁₀ - lp₀₁ + lp₀₀           # predicciones obtenidas a través de predict()
RERI = p₁₁ - p₁₀ - p₀₁ + p₀₀                   # probabilidades obtenidas a través de plogis()
```

Donde:
- ₁₁: DrogaA + DrogaB
- ₁₀: DrogaA sola
- ₀₁: DrogaB sola
- ₀₀: Ninguna


### Parámetros para ajuste de fórmula

Configuración utilizada en presentación final:
```r
spline_individuales <- TRUE    # Splines para efectos individuales
include_sex <- FALSE           # Incluir sexo como covariable
include_stage_sex <- FALSE     # Interacción stage-sex
nichd_spline <- FALSE          # Spline para efecto base de stage
bs_type <- "cs"                # Tipo de spline: "cs", "tp", "cr"
k_spline <- 7                  # Número de knots (= etapas NICHD)
select <- FALSE                # Penalización hasta cero
```

## Estructura de Datos

Columnas requeridas:
- `safetyreportid`: ID único del reporte
- `atc_concept_id`: Código ATC de fármaco
- `meddra_concept_id`: Código MedDRA de evento adverso
- `nichd`: Etapa NICHD (factor ordenado)
- `nichd_num`: Etapa NICHD numérica (1-7)
- `sex`: Sexo (opcional, si `include_sex = TRUE`)


