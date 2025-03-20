# ---------- Instalación de paquetes y carga de datos -----------
# Instalación de paquetes necesarios para el preprocesamiento RNA-Seq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("org.Mm.eg.db")
install.packages('gplots')
install.packages("RColorBrewer")

# Cargar paquetes
library(limma)
library(edgeR)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# -------------------- Cargar el dataset --------------------
# Cambiar el directorio de trabajo
setwd("/Users/jonathanquishpe/JoniDev/BIO/EPD4/Problema2")

# Cargar el dataset de RNA-Seq
seqdata <- read.delim("GSE290268_RNAseq_FPKM_data.txt", stringsAsFactors = FALSE)
head(seqdata)
dim(seqdata)

# Formatear los datos
# Eliminar la primera columna (identificadores de genes)
countdata <- seqdata[,-1]

# Almacenar los identificadores de genes como nombres de fila
rownames(countdata) <- seqdata[,1]
View(countdata)

##############################################################
# ---------------- Preprocesamiento CPM ---------------------#
##############################################################
# Calcular CPM (Counts Per Million)
myCPM <- cpm(countdata)
View(myCPM)

# -------- Detectar y Filtrar genes de baja expresión --------
# Definir un umbral de expresión (CPM > 0.5)
threshold <- myCPM > 0.5
head(threshold)

# Resumen de cuántos genes cumplen con el umbral en cada muestra
table(rowSums(threshold))

# Mantener genes con al menos 8 muestras que cumplen el umbral
keep <- rowSums(threshold) >= 8
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)

# -------- Normalizar y generar nuevo conjunto de datos --------
# Normalizar los datos usando CPM
counts.normalised <- myCPM[keep,]
dim(counts.normalised)

# Convertir a un objeto edgeR
dgeObj <- DGEList(counts.keep)
dgeObjNormalised <- DGEList(counts.normalised)


##############################################################
# ---------------- Preprocesamiento RPKM ---------------------#
##############################################################

# Primero necesitamos la longitud de los genes
# Asumiendo que tenemos un vector gene.length con las longitudes de los genes

# Calcular RPKM
myRPKM <- rpkm(countdata, gene.length = gene.length)
View(myRPKM)

# -------- Detectar y Filtrar genes de baja expresión --------
# Definir un umbral de expresión (RPKM > 0.5)
threshold <- myRPKM > 0.5
head(threshold)

# Resumen de cuántos genes cumplen con el umbral en cada muestra
table(rowSums(threshold))

# Mantener genes con al menos 8 muestras que cumplen el umbral
keep <- rowSums(threshold) >= 8
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)

# -------- Normalizar y generar nuevo conjunto de datos --------
# Normalizar los datos usando RPKM
counts.normalised <- myRPKM[keep,]
dim(counts.normalised)

# Convertir a un objeto edgeR
dgeObj <- DGEList(counts.keep)
dgeObjNormalised <- DGEList(counts.normalised)