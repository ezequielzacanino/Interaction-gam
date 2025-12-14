# Análisis de Poder Estadístico - Documentación

## Descripción

El script `04_power_analysis.R` implementa un análisis completo de poder estadístico para la detección de interacciones farmacológicas usando modelos GAM (Generalized Additive Models) y IOR (Information Odds Ratio) clásico.

## Funcionalidades Principales

### 1. Cálculo de Potencia
- **Objetivo**: Calcular poder estadístico para diferentes combinaciones de fold-change y número de reportes
- **Parámetros**: 
  - Fold-change: 1.1 a 5.0 (incrementos de 0.2)
  - Número de reportes: 10 a 200 (incrementos de 10)
  - Simulaciones por punto: 100

### 2. Identificación de Umbrales
- **Objetivo**: Determinar umbrales mínimos para 80% de poder estadístico
- **Métodos**: GAM e IOR clásico
- **Output**: Tabla con fold-change mínimo y número mínimo de reportes

### 3. Visualización con Heatmaps
- Heatmaps de poder por fold-change y número de reportes
- Contornos que marcan umbral del 80% de poder
- Gráficos separados para GAM e IOR clásico

### 4. Filtrado del Pool Positivo
- Aplica filtros basados en umbrales de 80% poder
- Genera tres pools: GAM, IOR clásico, y combinado
- Compara tamaños de pools filtrados vs original

### 5. Análisis de Validación
- Ejecuta validación completa en pools filtrados
- Calcula métricas de performance (AUROC, sensibilidad, especificidad, PPV, NPV, accuracy, F1-score)
- Análisis por etapa NICHD y por número de reportes

## Uso

### Requisitos Previos
```r
# Asegurar que existen los archivos:
- ./ade_raw.csv
- ./augmentation_results/positive_triplets_results.rds
- ./augmentation_results/negative_triplets_results.rds
- ./null_distribution_results/null_thresholds.csv
```

### Configuración
```r
# Parámetros principales (líneas 36-56 del script)
POWER_TARGET <- 0.80          # Objetivo de poder estadístico
fold_change_range <- seq(1.1, 5, by = 0.2)  # Rango fold-change
n_reports_range <- seq(10, 200, by = 10)    # Rango número de reportes
n_simulations <- 100          # Simulaciones por punto de grid

# Configuración de percentil
PERCENTILE_LEVEL <- "p95"     # Usado para umbrales nulos

# Configuración GAM (debe coincidir con otros scripts)
spline_individuales <- FALSE  
include_sex <- FALSE          
include_stage_sex <- FALSE    
k_spline <- 7                 
nichd_spline <- TRUE
bs_type <- "cs"
select <- TRUE
method <- "fREML" 
```

### Ejecución
```r
# Ejecutar script completo
source("04_power_analysis.R")

# El script creará automáticamente:
# - Directorio ./power_analysis_results/
# - Archivos CSV con resultados
# - Gráficos PNG con heatmaps y curvas ROC
```

## Archivos de Salida

### Datos
- `power_analysis_results.csv` - Resultados completos del análisis de poder
- `power_thresholds_80.csv` - Umbrales para 80% de poder
- `positive_pool_filtered_gam.csv` - Pool filtrado (GAM)
- `positive_pool_filtered_classic.csv` - Pool filtrado (IOR clásico)
- `positive_pool_filtered_combined.csv` - Pool filtrado (combinado)

### Métricas
- `metrics_validation_pool.csv` - Métricas globales pool filtrado
- `metrics_by_stage_validation.csv` - Métricas por etapa NICHD
- `metrics_by_reports_validation.csv` - Métricas por número de reportes

### Gráficos
- `heatmap_power_gam.png` - Heatmap poder para GAM
- `heatmap_power_classic.png` - Heatmap poder para IOR clásico
- `fig_roc_filtered_pool.png` - Curva ROC para pool filtrado

## Interpretación de Resultados

### Umbrales de Poder
Los umbrales identificados indican las características mínimas necesarias para detectar el 80% de las interacciones reales:

- **Fold-change**: Factor mínimo de aumento del efecto
- **Número de reportes**: Mínimo de coadministraciones requeridas

### Heatmaps
- **Colores**: Rojo intenso = alto poder, blanco = bajo poder
- **Contorno amarillo**: Marca el umbral del 80% de poder
- **Interpretación**: Valores en la esquina superior derecha tienen mayor poder

### Métricas de Validación
- **AUROC**: Capacidad discriminativa (0.5=aleatorio, 1.0=perfecto)
- **Sensibilidad**: Proporción de positivos detectados
- **Especificidad**: Proporción de negativos correctamente identificados
- **PPV/NPV**: Valores predictivos positivo y negativo

## Configuración Avanzada

### Modificar Rango de Análisis
```r
# Para análisis más granular
fold_change_range <- seq(1.1, 3, by = 0.1)      # Rango más estrecho
n_reports_range <- seq(20, 100, by = 5)        # Incrementos más pequeños

# Para análisis más amplio
fold_change_range <- seq(1.0, 10, by = 0.5)    # Rango más amplio
n_reports_range <- seq(5, 500, by = 25)        # Mayor número de reportes
```

### Ajustar Número de Simulaciones
```r
# Para mayor precisión (tiempo de ejecución mayor)
n_simulations <- 200

# Para ejecución más rápida (menor precisión)
n_simulations <- 50
```

### Cambiar Objetivo de Poder
```r
# Para 90% de poder
POWER_TARGET <- 0.90

# Para 70% de poder
POWER_TARGET <- 0.70
```

## Consideraciones Computacionales

### Tiempo de Ejecución
- **Estimación**: 1-3 horas para grid completo
- **Factores**: Número de simulaciones, tamaño del grid, complejidad de datos
- **Optimización**: El script incluye progress bars para monitorear progreso

### Uso de Memoria
- El script procesa datos en chunks para minimizar uso de memoria
- Se limpián objetos intermedios durante la ejecución
- Recomendado: 8GB+ RAM para datasets grandes

### Paralelización
- Utiliza procesamiento en paralelo cuando es posible
- Detecta automáticamente número de cores disponibles
- Configurable para entornos con restricciones de recursos

## Integración con Otros Scripts

### Dependencias
- `00_functions.R` - Funciones base (GAM, IOR, etc.)
- `01_augmentation.R` - Datos de augmentation
- `032_validation_param.R` - Funciones de validación

### Compatibilidad
- Mantiene mismos parámetros de configuración GAM
- Utiliza misma estructura de datos y formatos
- Genera archivos compatibles con pipeline existente

## Troubleshooting

### Errores Comunes
1. **Archivos no encontrados**: Verificar rutas en configuración
2. **Memoria insuficiente**: Reducir número de simulaciones o tamaño del grid
3. **Convergencia de modelos**: El script incluye manejo de errores robusto
4. **Datos insuficientes**: Verificar que existan datos de augmentation previos

### Validación de Resultados
- Los resultados se guardan automáticamente en múltiples formatos
- Verificar que los heatmaps se generen correctamente
- Confirmar que las métricas estén en rangos esperados (AUROC 0.5-1.0)

## Ejemplo de Uso Completo

```r
# 1. Verificar archivos necesarios
list.files("./augmentation_results/")
list.files("./null_distribution_results/")

# 2. Configurar parámetros
POWER_TARGET <- 0.80
fold_change_range <- seq(1.2, 4, by = 0.3)
n_reports_range <- seq(20, 150, by = 20)
n_simulations <- 50  # Reducido para ejemplo

# 3. Ejecutar análisis
source("04_power_analysis.R")

# 4. Revisar resultados
results <- fread("./power_analysis_results/power_analysis_results.csv")
thresholds <- fread("./power_analysis_results/power_thresholds_80.csv")

print("Umbrales para 80% poder:")
print(thresholds)
```

## Contacto y Soporte

Para dudas sobre el script de análisis de poder, revisar:
1. Comentarios en el código fuente
2. Logs de ejecución generados automáticamente
3. Archivos de resultados para diagnóstico

---
**Versión**: 1.0  
**Fecha**: Diciembre 2024  
**Autor**: Sistema de Análisis Farmacovigilancia GAM