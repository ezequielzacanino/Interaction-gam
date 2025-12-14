# Modelado de la Dinámica Ontogénica de Interacciones Farmacológicas en Pediatría

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)


## Descripción

Este repositorio contiene los scripts utilizados para la validación de modelos aditivos generalizados (GAM) para la detección de señales de interacción farmacológica con modelado de dinámicas ontogénicas, en sistemas de reporte espontáneo de efectos adversos en pediatría.

## Conceptos

- **Triplete**: Combinación de DrogaA + DrogaB + Evento adverso
- **Modulación ontogénica**: Variaciones en el odds de efecto adverso según etapa del desarrollo (NICHD stages)
- **Dinámica**: Patrón temporal de la modulación (uniforme, aumento, disminución, meseta, meseta inversa)
- **IOR (Interaction Odds Ratio)**: Métrica principal que cuantifica la interacción sinérgica o antagónica entre drogas

## Workflow

```
00_functions.R               # Funciones del pipeline
01_augmentation.R            # Generación de datos semisintéticos
02_nulldistribution.R        # Construcción de distribución nula
032_validation_param.R       # Validación y métricas de rendimiento
033_validation_comparison.R  # Comparación entre IOR clásico y GAM
```

## Metodología

### 1. Generación de Datos Semisintéticos (`01_augmentation.R`)

**Controles Positivos:**
- Armado de tripletes con requisitos mínimos (10 reportes conjuntos)
- Selección aleatoria de 750 tripletes como candidatos
- Inyección de señal con dinámica ontogénica simulada (increase, decrease, plateau, inverse-plateau, uniform)
- Fold-changes muestreados de distribución exponencial (λ = 0.75)
- Ajustado de positivos inyectados

**Controles Negativos:**
- Selección aleatoria de 2500 tripletes no inyectados (mutuamente exclusivos con positivos) con características de reporte similares a positivos
- Ajustado de negativos 

**Fórmula de inyección:**
```
P(evento|etapa) = fold_change × e_j + dinámica(etapa)

siendo e_j = P(A) + P(B) - P(A)×P(B)  
```

### 2. Distribución Nula (`02_nulldistribution.R`)

- Pool aleatorio de 100,000 reportes
- Permutación estratificada por etapa de desarrollo (rompe asociación droga-evento)
- Generación de 10000 tripletes permutados
- Cálculo de percentiles (P90, P95, P99) del IC90 inferior por etapa

### 3. Criterio de Detección 

Una señal es **positiva** si al menos una etapa cumple:
1. **Criterio nominal**: IC90 inferior del log-IOR > 0
2. **Criterio null model**: IC90 inferior del log-IOR > Percentil de distribución nula (seleccionable)

### 4. Modelo GAM

**Fórmula de base:**
```r
evento ~ drugA + drugB +
         s(nichd_num, bs="cs", k=7) +
         s(nichd_num, by=drugA_drugB, bs="cs", k=7) 
         
```

**Métrica principal:**
```
log(IOR) = logit(P₁₁) - logit(P₁₀) - logit(P₀₁) + logit(P₀₀)
```

Donde:
- P₁₁: Prob(evento | DrogaA + DrogaB)
- P₁₀: Prob(evento | DrogaA sola)
- P₀₁: Prob(evento | DrogaB sola)
- P₀₀: Prob(evento | ninguna)


### Parámetros para ajuste de fórmula

```r
spline_individuales <- FALSE    # Splines para efectos individuales
include_sex <- FALSE            # Incluir sexo como covariable
include_stage_sex <- FALSE      # Interacción stage-sex
nichd_spline <- TRUE            # Spline para efecto base de stage
bs_type <- "cs"                 # Tipo de spline: "cs", "tp", "cr"
k_spline <- 7                   # Número de knots (= etapas NICHD)
select <- FALSE                 # Penalización hasta cero
```

## Estructura de Datos

Columnas requeridas:
- `safetyreportid`: ID único del reporte
- `atc_concept_id`: Código ATC de fármaco
- `meddra_concept_id`: Código MedDRA de evento adverso
- `nichd`: Etapa NICHD (factor ordenado)
- `nichd_num`: Etapa NICHD numérica (1-7)
- `sex`: Sexo (opcional, si `include_sex = TRUE`)


