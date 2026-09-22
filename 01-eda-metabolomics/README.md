# EDA: de Metabolomics Workbench a SummarizedExperiment

Análisis exploratorio de datos (EDA) de metabolómica dirigido, entre otros objetivos, a familiarizarse con la clase `SummarizedExperiment` de Bioconductor y con las técnicas de análisis multivariante en R.

## Objetivo

Explorar si existen patrones diferenciados en la expresión de metabolitos entre ratones *wild type* (WT) y ratones con pérdida de función del gen **ZMYM2** (genotipo HET), un gen asociado a un mayor riesgo de esquizofrenia y trastornos del neurodesarrollo cuyo papel metabólico no está esclarecido.

## Datos

- **Fuente:** [Metabolomics Workbench](https://www.metabolomicsworkbench.org/), estudio **ST003805**.
- **Estudio original:** Huang, Wei-Chao, Kira Perzel Mandell, Sameer Aryal, et al. 2025. *Epigenetic Changes, Neuronal Dysregulation and Metabolomic Abnormalities in Zmym2 Mutant Mice, a Genetic Model of Schizophrenia and Neurodevelopmental Disorders*. bioRxiv. [https://doi.org/10.1101/2025.02.18.638656](https://doi.org/10.1101/2025.02.18.638656) 
- **Dataset utilizado:** resultado del análisis por espectrometría de masas de la cromatografía líquida de interacción hidrofílica en modo de ion positivo (HILIC positive ion mode), `ST003805_AN006254_Results.txt`, junto con sus metadatos (`coldata.txt`, `rowdata.txt`, `metadata.txt`).
- `data_pim.txt` y `data_untpim.txt` corresponden a otros datasets del mismo estudio combinado (otras cromatografías) que no se usan en el análisis final, conservados por trazabilidad de la descarga original.

## Métodos

1. Preprocesado de los datos crudos: eliminación de columnas con valores faltantes y transformación logarítmica (log2).
2. Construcción de un objeto `SummarizedExperiment` (assays, colData, rowData, metadata).
3. Análisis univariante de comparación de grupos (boxplot) para validar el preprocesado.
4. Análisis multivariante: análisis de componentes principales (PCA), agrupación jerárquica y heatmap de distancias.

## Cómo reproducirlo

1. Abrir `Analisis_Exploratorio.Rproj` en RStudio.
2. Renderizar `Informe_analisis_exploratorio.qmd` con Quarto (`quarto render` o el botón *Render*).

El informe ya renderizado está disponible en `Informe_analisis_exploratorio.pdf`.

## Contenido de la carpeta

```
Informe_analisis_exploratorio.qmd   # Informe (código + texto)
Informe_analisis_exploratorio.pdf   # Informe renderizado
ST003805_AN006254_Results.txt       # Matriz de expresión de metabolitos
coldata.txt / rowdata.txt / metadata.txt  # Metadatos del SummarizedExperiment
references.bib                      # Bibliografía citada en el informe
```
