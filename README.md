# Portfolio de análisis de datos ómicos

Colección de análisis estadísticos aplicados a datos ómicos públicos, desarrollados en R con [Quarto](https://quarto.org/) y presentados en formato de artículo científico. El objetivo de este repositorio es desarrollar el pipeline de un análisis ómico: desde la obtención e importación de datos públicos hasta el análisis exploratorio, el análisis multivariante y, en el segundo caso, la identificación de expresión diferencial.

Los dos análisis parten de datasets y preguntas biológicas distintas, pero comparten metodología (uso de `SummarizedExperiment` de Bioconductor, análisis exploratorio multivariante).

## Análisis incluidos

| Carpeta | Análisis | Dataset | Fuente | Técnicas |
|---|---|---|---|---|
| [`01-eda-metabolomics/`](01-eda-metabolomics) | Análisis exploratorio de datos (EDA) de metabolómica | ST003805 — metabolitos de ratones WT vs. HET (pérdida de función de ZMYM2) | [Metabolomics Workbench](https://www.metabolomicsworkbench.org/) | Estadística descriptiva, PCA, clustering jerárquico, heatmap |
| [`02-differential-expression/`](02-differential-expression) | Análisis de expresión diferencial (RNA-seq) | GSE161731 — sangre periférica de pacientes con COVID-19, neumonía bacteriana y controles sanos ([McClain et al., 2021](https://doi.org/10.1038/s41467-021-21289-y)) | [GEO (Gene Expression Omnibus)](https://www.ncbi.nlm.nih.gov/geo/) | EDA, MDS, PCA, limma-voom, análisis de sobrerrepresentación (GO) |

Cada carpeta es un análisis autocontenido, con su propio `.Rproj`, datos, informe en Quarto (`.qmd`) y PDF ya renderizado. Consulta el README de cada carpeta para el detalle de objetivos, métodos y cómo reproducirlo.

## Procedencia y uso de los datos

Todos los datos provienen de repositorios públicos financiados por el NIH y pensados explícitamente para su reutilización (GEO y Metabolomics Workbench). En cada informe y README se cita el identificador del dataset (accession) y el estudio original correspondiente.

## Requisitos generales

- R (≥ 4.x) y RStudio
- [Quarto](https://quarto.org/docs/get-started/)
- Paquetes de Bioconductor: `SummarizedExperiment`, `edgeR`, `limma`, `clusterProfiler`, `EnsDb.Hsapiens.v86`, `org.Hs.eg.db`, entre otros (ver el bloque de instalación al inicio de cada `.qmd`)

## Autor

Joaquín Villegas
