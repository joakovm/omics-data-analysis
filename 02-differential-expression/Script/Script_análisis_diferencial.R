#Librerias
if (!require(Rsubread, quietly = TRUE)) { BiocManager::install("Rsubread"); library(Rsubread) }
if (!require(edgeR, quietly = TRUE)) { BiocManager::install("edgeR"); library(edgeR) }
if (!require(limma, quietly = TRUE)) { BiocManager::install("limma"); library(limma) }
if (!require(pheatmap, quietly = TRUE)) { BiocManager::install("pheatmap"); library(pheatmap) }
if (!require(clusterProfiler, quietly = TRUE)) { BiocManager::install("clusterProfiler"); library(clusterProfiler) }
if (!require(factoextra, quietly = TRUE)) { install.packages("factoextra", dep=TRUE); library(factoextra) }
if (!require(org.Hs.eg.db, quietly = TRUE)) { BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db) }
if (!require(EnsDb.Hsapiens.v86, quietly = TRUE)) { BiocManager::install("EnsDb.Hsapiens.v86"); library(EnsDb.Hsapiens.v86) }
library(AnnotationDbi)

#Creación de directorios
workingDir<-getwd()
dataDir<-file.path(workingDir, "Datos")
resultadosDir<-file.path(workingDir, "Resultados")

#Lectura de datos
#Matriz de recuento
counts<-read.csv(file=file.path(dataDir, "GSE161731_counts.csv"), header = TRUE, row.names = 1,  check.names = FALSE)

#Metadatos
colData<-read.csv(file=file.path(dataDir,"GSE161731_counts_key.csv"), header = TRUE, row.names = 1)
metadata<-read.table(file=file.path(dataDir,"GSE161731_metadata.txt"), header = FALSE, sep = "\t", quote="", fill=TRUE)

#Identificación y limpieza de muestras diferentes
setdiff(colnames(counts), rownames(colData))
#Selección por nombre de columna (evita depender de la posición de las columnas extra)
counts1<- counts[, colnames(counts) %in% rownames(colData)]
#Comprobación muestras en común
stopifnot(all(colnames(counts1) == rownames(colData)))

#Creación del objeto rowRange
library(EnsDb.Hsapiens.v86)
rowranges <- genes(EnsDb.Hsapiens.v86, filter = GeneIdFilter(rownames(counts1)))
length(setdiff(rownames(counts1), names(rowranges)))
countsconid<-counts1[rownames(counts1) %in% names(rowranges),]
rowranges_ordenados <- rowranges[rownames(countsconid)]
stopifnot(all(rownames(countsconid) == names(rowranges_ordenados)))

#Creación del SummerizedExperiment
library(SummarizedExperiment)
sumex<-SummarizedExperiment(assays=list(count = as.matrix(countsconid)),
                            colData = colData, metadata = metadata,  
                            rowRanges = rowranges_ordenados)
sumex

#Preprocesado de los datos
#Separador comun "_"
library(stringr)
colData(sumex)$race<-str_replace_all(colData(sumex)$race, "[ /]", "_")
colData(sumex)$cohort<-str_replace_all(colData(sumex)$cohort, "[- ]", "_")

#Transformación de variables
colData(sumex)$age<- as.numeric(colData(sumex)$age)
colData(sumex)$gender<-as.factor(colData(sumex)$gender)
colData(sumex)$race<-as.factor(colData(sumex)$race)
colData(sumex)$cohort<-as.factor(colData(sumex)$cohort)
colData(sumex)$time_since_onset<-as.factor(colData(sumex)$time_since_onset)
colData(sumex)$hospitalized<-as.factor(colData(sumex)$hospitalized)
colData(sumex)$batch<-as.factor(colData(sumex)$batch)
str(colData(sumex))

#Selección de variables por cohortes
sumex_cohortes <- subset(sumex, select = cohort %in% c("COVID_19","Bacterial","healthy"))
#Eliminación de factores eliminados para facilitar graficar
colData(sumex_cohortes)$cohort<-droplevels(colData(sumex_cohortes)$cohort)
datossumex_coh<-assay(sumex_cohortes)
sumex_coh_nodup<- sumex_cohortes[!duplicated(datossumex_coh)]

#Selección de las 75 muestras al azar
#| include: true
#| echo: true
myseed <- sum(utf8ToInt("joaquinmartinezvillegas"))
set.seed(myseed)
seleccion<-sample(colnames(sumex_coh_nodup), 75)
sumex_coh_nodup_75<-sumex_coh_nodup[,seleccion]

# Normalización de los contajes para las 
#diferentes profundidades de secuenciación (CPM)
library(edgeR)
datos_75<-assay(sumex_coh_nodup_75)
counts.cpm<-cpm(datos_75)
umbral<-counts.cpm>0.5
filtro<-rowSums(umbral)>=5
datos_75_filtrado<-datos_75[filtro,]

#Transformación logarítmica de los datos
logcounts<-cpm(datos_75_filtrado, log = TRUE)
#Creación del objeto DGE
dgeobj<-DGEList(datos_75_filtrado)
#Normalización por composición de librerías (factor de normalización)
datos_75_filtrado_norm<-calcNormFactors(dgeobj)
#Matriz normalizada
logcounts_norm<-cpm(datos_75_filtrado_norm, log = TRUE)

#Cambio de formato ancho a vertical para graficar los datos
library(reshape2)
library(ggplot2)
plotdatos1<-melt(logcounts_norm)
plotdatos<-melt(logcounts)
#Gráficos de densidades
ggplot(plotdatos, aes(value, group = Var2)) + geom_density(alpha=0.05) + 
  labs(x="logCPM", y="Densidad", title = "Distribuciones de logCPM (Datos crudos)")
plotdatos1<-melt(logcounts_norm)
ggplot(plotdatos1, aes(value, group = Var2)) + geom_density(alpha=0.05) +
  labs(x="logCPM", y="Densidad", title = "Distribuciones de logCPM (Datos normalizados)")

#Cálculo de distancias euclídeas
library(factoextra)
distancias<-dist(t(logcounts_norm))
#Heatmap
fviz_dist((distancias), lab_size = 4)
#Resultado en pdf para ver mejor el gráfico (muchas muestras)
pdf(file = file.path(workingDir, "Resultados", "heatmap.pdf"))
fviz_dist((distancias))
dev.off()

#Creación de colores para PCA y MDS
#Se mapea por nombre de nivel (setNames), no por posición, para no depender
#del orden en que R ordene los niveles del factor
colorPorNivel<-function(factorVar, paleta){
  setNames(paleta, levels(factorVar))[as.character(factorVar)]
}
col.cohort<-colorPorNivel(colData(sumex_coh_nodup_75)$cohort, c("blue", "green", "red"))
head(data.frame(colData(sumex_coh_nodup_75)$cohort, col.cohort),10)
col.race<-colorPorNivel(colData(sumex_coh_nodup_75)$race, c("blue", "orange", "green", "red", "yellow", "darkgreen", "black"))
head(data.frame(colData(sumex_coh_nodup_75)$race, col.race),10)
col.batch<-colorPorNivel(colData(sumex_coh_nodup_75)$batch, c("blue", "red"))
head(data.frame(colData(sumex_coh_nodup_75)$batch, col.batch),10)
col.gender<-colorPorNivel(colData(sumex_coh_nodup_75)$gender, c("blue", "red"))
head(data.frame(colData(sumex_coh_nodup_75)$gender, col.gender),10)
col.hospitalized<-colorPorNivel(colData(sumex_coh_nodup_75)$hospitalized, c("blue", "red"))
head(data.frame(colData(sumex_coh_nodup_75)$hospitalized, col.hospitalized),10)
col.time_since_onset<-colorPorNivel(colData(sumex_coh_nodup_75)$time_since_onset, c("blue", "red", "green"))
head(data.frame(colData(sumex_coh_nodup_75)$time_since_onset, col.time_since_onset),10)

#MDS
limma::plotMDS(logcounts_norm,col=col.cohort, main="Cohort", cex=0.5)

#PCAs
#Selección de los 5 últimos digitos para poner nombres al graficar
names2plot<-substr(colnames(logcounts_norm),nchar(colnames(logcounts_norm)) - 4, nchar(colnames(logcounts_norm)))
#Función para graficar los dos primeros PCA
plotPCA<- function(X, labels=NULL, colors=NULL, dataDesc="",
                   scale=FALSE, formapunts=NULL, myCex=0.8){
  pcX<-prcomp(t(X), scale=FALSE)
  loads<-round(pcX$sdev^2/sum(pcX$sdev^2)*100, 1)
  xlab<-c(paste("PC1", loads[1], "%"))
  ylab<-c(paste("PC2", loads[2], "%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2], xlab=xlab, ylab=ylab, col=colors, pch=formapunts, main=dataDesc)
  text(pcX$x[,1], pcX$x[,2], names2plot, pos=1, cex=myCex)
}

#PCAs
plotPCA(logcounts_norm, colors=col.cohort, myCex=0.6, dataDesc="Plot of first 2 PCs_Cohorts")
plotPCA(logcounts_norm, colors=col.race, myCex=0.6, dataDesc="Plot of first 2 PCs_Race")
plotPCA(logcounts_norm, colors=col.batch, myCex=0.6, dataDesc="Plot of first 2 PCs_Batch")
plotPCA(logcounts_norm, colors=col.gender, myCex=0.6, dataDesc="Plot of first 2 PCs_Gender")
plotPCA(logcounts_norm, colors=col.hospitalized, myCex=0.6, dataDesc="Plot of first 2 PCs_Hospitalized")
plotPCA(logcounts_norm, colors=col.time_since_onset, myCex=0.6, dataDesc="Plot of first 2 PCs_Time since onset")

#Eliminación de outliers
#Identificados por inspección visual de los gráficos de PCA/MDS anteriores
#(se separan claramente del resto de muestras de su mismo grupo)
outliers <- c("95967", "DU18-02S0011625", "DU18-02S0011639", "DU09-02S0000151")
sumex_coh_nodup_75_so<-sumex_coh_nodup_75[, !(colnames(sumex_coh_nodup_75)
                                              %in% outliers)]

#Creación de la matriz de diseño
grupo<-paste(colData(sumex_coh_nodup_75_so)$race,
             colData(sumex_coh_nodup_75_so)$cohort, sep=".")
design=model.matrix(~0+grupo)
colnames(design)<-gsub("grupo", "", colnames(design))
rownames(design)<-rownames(colData(sumex_coh_nodup_75_so))

#Creación de la matriz de contraste
matrizcontraste<- makeContrasts(
  BlackBac.vs.WhiteBac = Black_African_American.Bacterial - White.Bacterial,
  BlackCov.vs.WhiteCov = Black_African_American.COVID_19 - White.COVID_19,
  WhiteCov.vs.Whiteheal = White.COVID_19 - White.healthy,
  AsianCov.vs.Asianheal = Asian.COVID_19 - Asian.healthy,
  BlackBac.vs.BlackCov = Black_African_American.Bacterial - Black_African_American.COVID_19,
  WhiteBac.vs.WhiteCov = White.Bacterial - White.COVID_19, 
  Bac.vs.healthy = Black_African_American.Bacterial - White.healthy,
  levels=design)

#Creación de otro objeto DGE sin outliers para ajustar los datos con VOOM
datos_sin_outliers<-assay(sumex_coh_nodup_75_so)
dgeobj1<-DGEList(datos_sin_outliers)
datos_sin_outliers_norm<-calcNormFactors(dgeobj1)

#Transformación de los datos con boom y ajuste al modelo
voomObj<-voom(datos_sin_outliers_norm, design)
fit<-lmFit(voomObj)
fit.contraste<-contrasts.fit(fit, matrizcontraste)
fit.contraste<-eBayes(fit.contraste)

#Resultados en topTable
toptab_BWbac<-topTable(fit.contraste, coef=1, sort.by = "p", number=nrow(fit.contraste),
                       adjust="fdr", lfc=1.5, p.value=0.05)
toptab_BWcov<-topTable(fit.contraste, coef=2, sort.by = "p", number=nrow(fit.contraste),
                       adjust="fdr", lfc=1.5, p.value=0.05)
toptab_Wcovheal<-topTable(fit.contraste, coef=3, sort.by = "p", number=nrow(fit.contraste),
                          adjust="fdr", lfc=1.5, p.value=0.05)
toptab_Acovheal<-topTable(fit.contraste, coef=4, sort.by = "p", number=nrow(fit.contraste), 
                          adjust="fdr", lfc=1.5, p.value=0.05)
toptab_Bbaccov<-topTable(fit.contraste, coef=5, sort.by = "p", number=nrow(fit.contraste), 
                         adjust="fdr", lfc=1.5, p.value=0.05)
toptab_Wbaccov<-topTable(fit.contraste, coef=6, sort.by = "p", number=nrow(fit.contraste), 
                         adjust="fdr", lfc=1.5, p.value=0.05)
toptab_bacheal<-topTable(fit.contraste, coef=7, sort.by = "p", number=nrow(fit.contraste), adjust="fdr", lfc=1.5, p.value=0.05)

#Guardado de resultados en disco para no perderlos al cerrar la sesión
listaToptabs<-list(BWbac=toptab_BWbac, BWcov=toptab_BWcov, Wcovheal=toptab_Wcovheal,
                    Acovheal=toptab_Acovheal, Bbaccov=toptab_Bbaccov,
                    Wbaccov=toptab_Wbaccov, bacheal=toptab_bacheal)
for (nombre in names(listaToptabs)) {
  write.csv(listaToptabs[[nombre]],
            file.path(resultadosDir, paste0("topTable_", nombre, ".csv")))
}
saveRDS(fit.contraste, file.path(resultadosDir, "fit_contraste.rds"))

#Volcanoplots
par(mar=c(2,2,2,2), mfrow=c(3,2))
volcanoplot(fit.contraste,coef=1,highlight=10,
            names=fit.contraste$genes$SYMBOL, main="BW.Bacterial")
volcanoplot(fit.contraste,coef=2,highlight=10,
            names=fit.contraste$genes$SYMBOL, main="BW.COVID_19_19")
volcanoplot(fit.contraste,coef=3,highlight=10,
            names=fit.contraste$genes$SYMBOL, main="WCOVID_19.vs.healthy")
volcanoplot(fit.contraste,coef=5,highlight=10,
            names=fit.contraste$genes$SYMBOL, main="Bbacterial.vs.COVID_19")
volcanoplot(fit.contraste,coef=6,highlight=10,
            names=fit.contraste$genes$SYMBOL, main="WBacterial.vs.COVID_19")

#Diagramas de Venn-Comparación Múltiple
diagramasVenn <- decideTests(fit.contraste, p.value = 0.05, lfc = 2)
summary(diagramasVenn)
vc<- vennCounts(diagramasVenn[,c(3, 7)])
vennDiagram(vc, include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue"), cex=c(1,1,1))

#Análisis significación biológica
#Creación de topTab sin restricciones para el universo
toptab_Wcovheal_universe<-topTable(fit.contraste, coef=3, sort.by = "p", number=nrow(fit.contraste))
topTabGO<-toptab_Wcovheal_universe
allEnsembl<-rownames(topTabGO)
selectedEnsemblDown<-rownames(subset(topTabGO, (logFC< -1.5) & (adj.P.Val < 0.05)))
#Análisis de sobreexpresión
library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(gene = selectedEnsemblDown, universe = allEnsembl, 
                keyType = "ENSEMBL", OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
#Comprobación de genes encontrados en la base de datos
sum(selectedEnsemblDown %in% keys(org.Hs.eg.db, keytype="ENSEMBL"))
#Gráfico resultado análisis de enriquecimiento
dotplot(ego, showCategory=12, font.size=4)