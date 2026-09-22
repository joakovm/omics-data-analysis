# Análisis de expresión diferencial en la respuesta a SARS-CoV-2 mediante RNA-seq

Análisis de expresión génica diferencial a partir de datos de RNA-seq de sangre periférica, con el pipeline **limma-voom**, para caracterizar el perfil transcriptómico de la respuesta frente a distintas infecciones.

## Objetivo

Identificar diferencias en el perfil de expresión génica entre pacientes con COVID-19, neumonía bacteriana e individuos sanos, y estudiar en qué medida la raza actúa como variable de confusión en dicha respuesta.

## Datos

- **Fuente:** [GEO (Gene Expression Omnibus)](https://www.ncbi.nlm.nih.gov/geo/), identificador **GSE161731**.
- **Estudio original:** McClain, M.T. et al. (2021). *Dysregulated transcriptional responses to SARS-CoV-2 in the periphery*. Nature Communications 12(1). [https://doi.org/10.1038/s41467-021-21289-y](https://doi.org/10.1038/s41467-021-21289-y)
- Matriz de recuentos y metadatos de sangre periférica de sujetos con COVID-19, coronavirus estacional, gripe, neumonía bacteriana y controles sanos (`Datos/GSE161731_counts.csv`, `Datos/GSE161731_counts_key.csv`, `Datos/GSE161731_metadata.txt`). Para este análisis se seleccionó un subconjunto de 75 muestras (semilla fija para reproducibilidad) correspondientes a los cohortes COVID-19, neumonía bacteriana y sanos.

## Métodos

1. Construcción de un objeto `SummarizedExperiment` (con `rowRanges` vía `EnsDb.Hsapiens.v86`).
2. Limpieza de metadatos y selección aleatoria de 75 muestras.
3. Preprocesado: normalización por CPM, filtrado de genes poco expresados, normalización TMM.
4. Exploración: distribuciones de densidad, heatmap de distancias, escalado multidimensional (MDS) y PCA para detectar outliers y variables de confusión (raza, cohorte, lote).
5. Análisis de expresión diferencial con **limma-voom** (siete contrastes por raza y cohorte).
6. Análisis de significación biológica (sobrerrepresentación de términos GO) sobre los genes subexpresados con `clusterProfiler`.

## Cómo reproducirlo

1. Abrir `Analisis_expresion_diferencial.Rproj` en RStudio.
2. Renderizar `Informe_Analisis_Expresion_Diferencial.qmd` con Quarto (`quarto render` o el botón *Render*).

El informe ya renderizado está disponible en `Informe_Analisis_Expresion_Diferencial.pdf`. 

## Contenido de la carpeta

```
Informe_Analisis_Expresion_Diferencial.qmd   # Informe (código + texto)
Informe_Analisis_Expresion_Diferencial.pdf   # Informe renderizado
Datos/                                       # Matriz de recuentos y metadatos (GEO)
Resultados/                                  # Figuras exportadas (heatmap, MDS, PCA, volcanoplots)
# Referencias bibliográficas (no vinculadas al YAML)
```

