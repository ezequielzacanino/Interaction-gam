# Modelado de la Dinámica Ontogénica de las Interacciones Farmacológicas en Pediatría

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)

## Descripción

Este repositorio implementa un flujo de validación metodológica para detectar señales de interacciones farmacológicas (drug-drug interactions, DDI) con dinámica ontogénica en población pediátrica usando modelos aditivos generalizados (GAM) sobre sistemas de reporte espontáneo.

El pipeline reproduce la lógica del manuscrito adjunto:

- genera datasets semisintéticos a partir de FAERS pediátrico curado;
- compara un enfoque GAM contra métodos estratificados clásicos basados en `IOR` y `RERI`;
- construye una distribución nula empírica por permutación;
- evalúa sensibilidad, especificidad, AUC, F1 y métricas relacionadas bajo escasez progresiva de datos;
- analiza plausibilidad biológica mediante redes droga-gen;
- compara señales detectadas contra `TWOSIDES`.

La unidad de análisis es el **triplete** `drugA - drugB - event`, definido a partir de reportes con coadministración y un evento adverso reportado.

## Objetivo analítico

El objetivo es evaluar si un GAM puede detectar mejor que los métodos estratificados patrones no lineales de desproporcionalidad a lo largo de las 7 etapas NICHD del desarrollo pediátrico.

En la implementación actual:

- los controles positivos se construyen por inyección semisintética de dinámicas ontogénicas;
- los controles negativos se seleccionan a partir de combinaciones alternativas de los mismos fármacos y eventos;
- la clasificación final del GAM usa un criterio doble: señal nominal positiva y superación de un umbral derivado de la distribución nula;
- la comparación metodológica se realiza en escala multiplicativa (`log(IOR)`) y aditiva (`RERI`).

## Conceptos clave

- **Triplete**: combinación `drugA`, `drugB`, `meddra`.
- **Modulación ontogénica**: cambio del patrón de desproporcionalidad según etapa NICHD.
- **Dinámica**: forma del patrón inyectado a través de las 7 etapas (`uniform`, `increase`, `decrease`, `plateau`, `inverse_plateau`).
- **IOR**: medida de interacción en escala multiplicativa.
- **RERI**: medida de interacción en escala aditiva.
- **Distribución nula**: distribución empírica construida por permutación estratificada por etapa.
- **Subset calibrado por poder**: subconjunto de tripletes retenidos para asegurar comparaciones más justas entre métodos.

## Estructura del repositorio

```text
00_functions.R               # Funciones base del pipeline, configuración global y utilidades
01_theme.R                   # Tema gráfico común
02_descriptive.R             # Análisis descriptivo exploratorio
03_descriptive.R             # Variante/adaptación del análisis descriptivo
10_augmentation.R            # Generación semisintética, positivos, negativos y sensibilidad
20_null.R                    # Construcción de distribución nula empírica
30_metrics.R                 # Cálculo de métricas, poder y comparaciones entre métodos
40_network.R                 # Análisis de redes, soporte biológico y comparación con TWOSIDES
41_graphs.R                  # Generación de figuras facetadas finales
50_supuestos.R               # Análisis adicionales de supuestos/diagnóstico
dataset_drugbank.R           # Preparación de insumos de DrugBank
simulacion_inyeccion.R       # Simulación auxiliar para validar la estrategia de inyección
drug.csv                     # Información de drogas
drug_gene.csv                # Relaciones droga-gen/proteína
ade_raw.csv                  # Base curada de reportes espontáneos
twosides/                    # Datos de referencia para comparación externa
vocabulary/                  # Vocabulario OMOP/MedDRA/ATC/RxNorm
results/                     # Resultados generados por distintas configuraciones
```

## Datos de entrada

### Archivo principal

`ade_raw.csv` debe contener, como mínimo:

- `safetyreportid`: identificador único del reporte.
- `atc_concept_id`: identificador ATC de la droga.
- `meddra_concept_id`: identificador MedDRA del evento.
- `nichd`: etapa del desarrollo pediátrico según NICHD.

Si se habilita `include_sex`, además:

- `sex`: sexo del paciente.

### Otros insumos requeridos

- `drug.csv`: nombres y metadatos de drogas.
- `drug_gene.csv`: pares droga-gen/proteína para soporte biológico.
- `vocabulary/concept.csv`: mapeos OMOP, MedDRA, ATC y RxNorm.
- `vocabulary/concept_relationship.csv`: relaciones de vocabulario.
- `twosides/TWOSIDES.csv.gz`: base de comparación externa.

## Configuración global

La configuración central vive en `00_functions.R`.

Parámetros activos en la versión actual del pipeline:

```r
spline_individuales <- TRUE
include_sex <- FALSE
include_stage_sex <- FALSE
include_nichd <- FALSE
nichd_spline <- FALSE
k_spline <- 7
bs_type <- "cs"
select <- FALSE
method <- "fREML"
percentil <- "p95"
```

Con esta combinación, el sufijo de resultados por defecto es:

```r
suffix <- "sics"
```

Por lo tanto, los outputs de la configuración actual se escriben en `results/sics/`.

## Resumen metodológico

### 1. Generación semisintética (`10_augmentation.R`)

- unifica IDs de drogas con mismo compuesto activo;
- construye tripletes candidatos desde reportes con al menos dos drogas y un evento;
- filtra positivos con al menos `2` reportes y presencia en al menos `2` etapas;
- selecciona `500` tripletes base positivos;
- inyecta dinámicas ontogénicas usando cinco perfiles (`uniform` incluida como comparador);
- asigna tamaños de efecto mediante una exponencial negativa (`lambda_fc = 0.75`);
- ajusta el GAM y calcula `log(IOR)` y `RERI` para cada triplete;
- genera un set negativo de hasta `10000` tripletes con fármacos y eventos del mismo universo;
- ejecuta reducción iterativa del dataset (`10%` a `90%`) para el análisis de sensibilidad;
- selecciona `100000` reportes para construir el pool usado por la distribución nula.

### 2. Distribución nula (`20_null.R`)

- toma el pool guardado por `10_augmentation.R`;
- permuta drogas y eventos dentro de cada etapa del desarrollo, preservando la estructura por etapa;
- genera tripletes permutados y reintroduce esas permutaciones en el dataset original;
- ajusta el GAM para los tripletes del universo nulo;
- calcula, por etapa, percentiles empíricos (`p90`, `p95`, `p99`) del límite inferior del IC90 para `log(IOR)` y `RERI`.

### 3. Métricas y calibración por poder (`30_metrics.R`)

- expande los resultados positivos y negativos por etapa;
- compara las distribuciones observadas frente a la distribución nula;
- calcula sensibilidad, especificidad, `PPV`, `NPV`, `F1` y `AUC`;
- estima intervalos por bootstrap no paramétrico (`n_boot = 2000`);
- define subsets calibrados a un poder objetivo del `80%`;
- evalúa resultados a nivel global, por dinámica y por etapa;
- analiza tres escenarios: `original`, `filtered` e `intersection`.

### 4. Redes y validación externa (`40_network.R`)

- ajusta candidatos del dataset original con al menos `5` reportes por triplete;
- detecta señales positivas con criterio GAM-IOR y umbral nulo por etapa;
- integra relaciones droga-gen/proteína derivadas de DrugBank;
- evalúa soporte biológico mediante genes compartidos e índice de Jaccard;
- contrasta la red observada frente a redes nulas generadas por rewiring;
- compara pares y tripletes positivos contra `TWOSIDES`;
- separa señales concordantes, novedosas y perdidas.

### 5. Figuras finales (`41_graphs.R`)

- genera figuras de dinámicas simuladas;
- produce figuras facetadas por etapa, métrica y versión del dataset;
- exporta gráficos en `png` y `svg` para resultados principales.

## Criterios de detección

### GAM

Un triplete se clasifica como positivo si existe al menos una etapa donde:

- el límite inferior del IC90 de `log(IOR)` o `RERI` es mayor a `0`; y
- además supera el umbral empírico de la distribución nula para esa etapa.

Por defecto, el percentil usado es `p95`.

### Métodos estratificados

Los métodos clásicos usan criterio nominal por etapa:

- `log(IOR)_lower90 > 0`
- `RERI_lower90 > 0`

## Orden recomendado de ejecución

```r
source("00_functions.R", local = TRUE)
source("10_augmentation.R", local = TRUE)
source("20_null.R", local = TRUE)
source("30_metrics.R", local = TRUE)
source("40_network.R", local = TRUE)
source("41_graphs.R", local = TRUE)
```

## Principales salidas

Bajo la configuración actual, los resultados se organizan en `results/sics/`:

```text
results/sics/
  augmentation_results/
  null_distribution_results/
  metrics_results/
```

Archivos relevantes generados por el pipeline:

- `augmentation_results/positive_triplets_metadata.csv`
- `augmentation_results/positive_triplets_results*.csv`
- `augmentation_results/negative_triplets_results*.csv`
- `augmentation_results/null_pool_reports_metadata.csv`
- `null_distribution_results/null_distribution.csv`
- `null_distribution_results/null_thresholds.csv`
- `null_distribution_results/null_thresholds_reri.csv`
- `metrics_results/metrics_global_original.csv`
- `metrics_results/metrics_dynamic_original.csv`
- `metrics_results/metrics_stage_original.csv`
- `metrics_results/fig_power_surface_ior_combined.png`
- `metrics_results/fig_power_surface_reri_combined.png`

En ejecuciones orientadas al manuscrito también aparecen resultados bajo `results/sics_manuscrito/`, incluyendo:

- `network/edge_metrics_interlayer.csv`
- `network/concordant_pairs_full_summary.csv`
- `network/novel_pairs_full_summary.csv`
- `network/metrics_versus_twosides.csv`
- `network/twosides_comparison_triplets.csv`

## Dependencias

El pipeline carga librerías mediante `pacman::p_load()`. Entre las dependencias principales:

- `data.table`
- `tidyverse`
- `mgcv`
- `MASS`
- `parallel`
- `doParallel`
- `foreach`
- `doRNG`
- `pROC`
- `akima`
- `igraph`
- `ggraph`
- `tidygraph`
- `patchwork`
- `networkD3`
- `htmlwidgets`
- `DHARMa`
- `svglite`

## Consideraciones prácticas

- `ade_raw.csv` es grande, por lo que se recomienda trabajar con suficiente memoria RAM.
- La paralelización usa por defecto el `75%` de los núcleos disponibles.
- Los scripts están pensados para ejecutarse desde la raíz del proyecto.
- `00_functions.R` fija el directorio de trabajo a `D:/Bioestadística/gam-farmacovigilancia`.
- El pipeline guarda checkpoints intermedios para evitar perder progreso en corridas largas.

## Alcance del README

Este README prioriza describir la implementación real del repositorio actual. Cuando existen diferencias menores entre el borrador del manuscrito y los parámetros del código, se documenta el comportamiento observado en los scripts vigentes.
