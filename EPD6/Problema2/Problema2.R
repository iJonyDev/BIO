####################
# INSTALL PACKAGES #
####################
install.packages(c("cluster", "factoextra", "caret", "ggplot2", "edgeR"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "Biobase", "GEOquery"))

###########################
# LOAD DATA AND LIBRARIES #
###########################
library(edgeR)
library(cluster)
library(factoextra)
library(ggplot2)
library(Biobase)

####################################
# 1. CARGAR Y PREPROCESAR DATOS RNA-SEQ
####################################
setwd("/Users/jonathanquishpe/JoniDev/BIO/EPD6/Problema2")
# Cargar el dataset de RNA-Seq
seqdata <- read.delim("/Users/jonathanquishpe/JoniDev/BIO/EPD6/Problema2/GSE290268_RNAseq_FPKM_data.txt", stringsAsFactors = FALSE)

# Formatear los datos
countdata <- seqdata[,-1]
rownames(countdata) <- seqdata[,1]

#########################
##  PREPROCESAMIENTO  ###
#########################

# Calcular CPM (Counts Per Million)
myCPM <- cpm(countdata)

# Filtrar genes de baja expresión (CPM > 0.5 en al menos 8 muestras)
threshold <- myCPM > 0.5
keep <- rowSums(threshold) >= 8
counts.keep <- countdata[keep,]
counts.normalised <- cpm(counts.keep, log = TRUE) # Log2-CPM para mejor distribución

# Crear objeto DGEList para análisis posteriores
dgeObj <- DGEList(counts = counts.keep)
dgeObj <- calcNormFactors(dgeObj)

####################################
# 4. TRANSFORMACIÓN PARA CLUSTERING
####################################

# Escalar datos para clustering (centrar y escalar)
preprocessParams <- preProcess(t(counts.normalised), method = c("center", "scale"))
dataForClustering <- predict(preprocessParams, t(counts.normalised))

####################################
# CLUSTERING DE GENES (EXPRESIÓN)
####################################

# Seleccionar genes más variables para clustering (top 1000 por varianza)
gene_vars <- apply(counts.normalised, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]
top_gene_data <- counts.normalised[top_genes,]

# Escalar datos de genes
scaled_gene_data <- t(scale(t(top_gene_data)))

# Clustering jerárquico de genes
gene_dist <- dist(scaled_gene_data)
gene_clusters <- hclust(gene_dist, method = "ward.D2")

# Cortar árbol para obtener clusters (ejemplo con 5 clusters)
gene_cluster_cut <- cutree(gene_clusters, k = 5)

# Visualización de heatmap
heatmap(as.matrix(scaled_gene_data), 
        Rowv = as.dendrogram(gene_clusters),
        Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none")

####################################
# PREPARACIÓN PARA ENRIQUECIMIENTO #
####################################

# Guardar listas de genes por cluster para enriquecimiento
for (i in 1:max(gene_cluster_cut)) {
  genes_in_cluster <- names(gene_cluster_cut[gene_cluster_cut == i])
  write.table(genes_in_cluster, 
              file = paste0("gene_cluster_", i, ".txt"),
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE,
              quote = FALSE)
}

# Guardar todos los genes analizados
write.table(rownames(counts.normalised), 
            file = "all_genes.txt",
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

# Guardar datos normalizados para análisis posteriores
write.csv(counts.normalised, file = "normalized_expression_data.csv")

