####################
# INSTALL PACKAGES #
####################
install.packages(c("cluster","factoextra","biclust","caret","NbClust"))

# BioConductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("Biobase")
BiocManager::install("GEOquery")

###########################
# LOAD DATA AND LIBRARIES #
###########################
library(Biobase)
library(GEOquery)
library(caret)
library(cluster)
library(factoextra)

setwd("/Users/jonathanquishpe/JoniDev/BIO/EPD6/Problema1")
# Obtener los datos GSE62301
dataSetComplete <- getGEO("GSE62301", GSEMatrix=TRUE, AnnotGPL=TRUE)

# Para series GSE, los datos vienen en una lista, normalmente usamos el primer elemento
expressionSet <- dataSetComplete[[1]]

# Obtener la matriz de expresión
dataGSE <- exprs(expressionSet)

# Eliminar NAs si existen
dataGSE <- na.omit(dataGSE)

# Obtener información sobre genes y muestras
genesGSE <- rownames(dataGSE)
samplesGSE <- colnames(dataGSE)

# Crear copia para preprocesamiento
dataGSEPreProcess <- dataGSE

# Obtener número de muestras
samplesCount = ncol(dataGSEPreProcess)

# Obtener resumen de los datos
summary(dataGSEPreProcess)

############
# OUTLIERS #
############
for (i in 1:samplesCount) { # Por cada columna del dataset
  
  # Detección de outliers
  Q1 = quantile(dataGSEPreProcess[,i], 0.25, na.rm = TRUE)
  Q3 = quantile(dataGSEPreProcess[,i], 0.75, na.rm = TRUE)
  IQR = Q3-Q1
  leveDown = Q1-1.5*IQR
  leveUp = Q3 + 1.5*IQR
  
  media = mean(dataGSEPreProcess[,i], na.rm = TRUE)
  
  # Detectar tanto outliers superiores como inferiores
  outliers <- dataGSEPreProcess[,i] >= leveUp | dataGSEPreProcess[,i] <= leveDown
  dataGSEPreProcess[outliers,i] <- media
}

summary(dataGSEPreProcess)

#################
# PREPROCESSING #
#################
preprocessParams <- preProcess(dataGSEPreProcess, method=c("center", "scale")) 
dataGSEPreProcess <- predict(preprocessParams, dataGSEPreProcess)
summary(dataGSEPreProcess)

######################
# CLUSTERING: KMEANS #
######################
# Numero optimo de custers clusters KMeans
# Calcula el WSS inicial
wss <- (nrow(dataGSEPreProcess)-1)*sum(apply(dataGSEPreProcess,2,var))

# Calcula WSS para diferentes números de clusters
for (i in 2:15) 
  wss[i] <- sum(kmeans(dataGSEPreProcess, centers=i, iter.max=80)$withinss)

# Grafica el método del codo
plot(1:15, wss, type="b", xlab="Number of Clusters", 
     ylab="Suma de cuadrados dentro de los grupos",
     main="Método del Codo para K-means")

numClusters = 3 # Numero optimo de clusters para KMeans.
set.seed(100)
clusters <- kmeans(dataGSEPreProcess, centers = numClusters, iter.max=80)
clusters # 95.8 % de datos bien agrupados con 3 clusters.

# Visualización de los clusters
fviz_cluster(clusters, data = dataGSEPreProcess, geom = "point", ellipse.type = "convex")

# Seleccionamos aquellos clusters que tiene una suma de cuadrados mas cercana a 0.
# Almacenamos el cluster 1(1212.552), cluster 3(1347.270) que son los mejores clusters para enriquecerlos.

data_clus_1 <- dataGSEPreProcess[clusters$cluster == 1,]
data_clus_3 <- dataGSEPreProcess[clusters$cluster == 3,]

# Guardar los identificadores de los genes de los clusters seleccionados.
write.table(rownames(data_clus_1), file="cluster1.txt",sep = "\t", row.names=FALSE, col.names=FALSE,
            quote = FALSE)
write.table(rownames(data_clus_3), file="cluster3.txt",sep = "\t", row.names=FALSE, col.names=FALSE,
            quote = FALSE)
# Guardarmos todos los genes.
all_genes <- rownames(dataGSEPreProcess)
write.table(all_genes, file="all_genes.txt",sep = "\t", row.names=FALSE, col.names=FALSE,
            quote = FALSE)
