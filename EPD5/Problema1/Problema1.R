###########################
# INSTALACION Y LIBRERIAS #
###########################
install.packages(c("mlbench", "cluster", "factoextra", "biclust", "caret", "pheatmap", "ggplot2", "tidyr", "reshape2"))
# Load all required libraries
libraries <- c("mlbench", "factoextra", "pheatmap", "cluster", 
               "biclust", "caret", "ggplot2", "tidyr", "reshape2")
lapply(libraries, library, character.only = TRUE)

dataSet <- read.table("/Users/jonathanquishpe/JoniDev/BIO/EPD5/Problema1/GSE227018_norm_rma.txt", header = TRUE, sep = "\t", row.names = 1)

# --- Verificamos la distribucion de los datos antes del tratamiento de outliers y muestreo ---
str(dataSet)  # Estructura del dataset
summary(dataSet[,1:9]) # Resumen de las primeras 9 columnas
samplesCount = 9 # Cantidad de columnas

# Convertir el dataset a formato largo para visualizar los histogramas
data_long <- melt(dataSet)

# Crear los histogramas
ggplot(data_long, aes(x = value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_bw() +
  labs(title = "Distribución antes de tratamiento de outliers",
       x = "Valor",
       y = "Frecuencia")

############
# OUTLIERS #
############
for (i in 1:9){ # Por cada columna del dataset
  print(i)
  dataColumn <- dataSet[,i]
  # Deteccion de outliers
  Q1 = quantile(dataColumn, 0.25, na.rm = TRUE)
  Q3 = quantile(dataColumn, 0.75, na.rm = TRUE)
  IQR=Q3-Q1
  leveDown = Q1-1.5*IQR
  leveUp = Q3 + 1.5 * IQR
  
  # Recodificacion
  is.na(dataColumn) <- dataColumn >= leveUp
  is.na(dataColumn) <- dataColumn <= leveDown
  # Reemplazar los datos en el conjunto de datos preprocesado
  dataColumn[is.na(dataColumn)] <- mean(dataColumn, na.rm=TRUE)
  print(summary(dataColumn))
  dataSet[,i] <- dataColumn
}

# Resultado final despues del tratamiento de outliers
summary(dataSet[,1:9])

# --- Verificamos la distribucion de los datos despues de tratamiento de outliers y antes del muestreo ---

# Convertir el dataset a formato largo para visualizar los histogramas
data_long <- melt(dataSet)

# Crear los histogramas
ggplot(data_long, aes(x = value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_bw() +
  labs(title = "Distribución después de tratamiento de outliers y antes del muestreo",
       x = "Valor",
       y = "Frecuencia")


############
# MUESTREO #
############
# Tomamos una muestra porque el tamanio del dataset original es muy grande
# y requiere mucho tiempo y recursos de procesamiento.

# Establecer una semilla para reproducibilidad
set.seed(123)  # Puedes usar cualquier número entero como semilla

# Establecemos un tamanio de muestra
sample_size <- 1000

# Generar muestra aleatoria
dataSample <- dataSet[sample(nrow(dataSet), sample_size), ]
str(dataSample)

# --- Verificamos la distribucion de los datos despues de muestreo ---

# Convertir el dataset a formato largo para visualizar los histogramas
data_long <- melt(dataSample)

# Crear los histogramas de la muestra
ggplot(data_long, aes(x = value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_bw() +
  labs(title = "Distribución de todas las columnas de la muestra",
       x = "Valor",
       y = "Frecuencia")

# Resumen de la muestra
summary(dataSample)


####################
# PREPROCESAMIENTO #
####################
# Se aplica el metodo de estandarizacion ya que es el que mejor funciona para los algoritmos de clustering y biclustering.
preprocessParams <- preProcess(dataSample[,1:samplesCount], method=c("center", "scale"))
dataSample[,1:samplesCount] <- predict(preprocessParams, dataSample[,1:samplesCount])
summary(dataSample)

##########################
# CLUSTERING: DISTANCIAS #
##########################
distanciasDataSet <- get_dist(dataSample[,1:samplesCount], method = "pearson")               
head(round(as.matrix(distanciasDataSet), 2))[, 1:samplesCount]
distanciasDataSet
fviz_dist(distanciasDataSet, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
#######################
# CLUSTERING: K-MEANS #
#######################

# 1) Determinar el numero optimo de clusters.
fviz_nbclust(dataSample[,1:samplesCount], kmeans, method = "gap_stat")
fviz_nbclust(dataSample[,1:samplesCount], kmeans, method = "wss")
fviz_nbclust(dataSample[,1:samplesCount], kmeans, method = "silhouette")

# 2) Visualizar pheatmap con 2 clusters
pheatmap(t(dataSample[,1:samplesCount]), cutree_cols = 2)

# 3) Ejecutar metodo de clustering
set.seed(123)
clusters.km <- kmeans(dataSample[,1:samplesCount],2,nstart = 25)
clusters.km # datos agrupados en 2 clusters

data_clus_1 <- dataSample[clusters.km$cluster == 1,]
data_clus_2 <- dataSample[clusters.km$cluster == 2,]

# 4) Visualizar resultados
fviz_cluster(clusters.km, data = dataSample[,1:samplesCount], palette = "jco", ggtheme = theme_minimal())


##########################
# CLUSTERING: JERARQUICO #
##########################
# 1) Ejecutar metodo de clustering
clusters.hc <- hclust(dist(dataSample[,1:samplesCount]),  method = "ward.D2")

# 2) Visualizar resultados
clusters.hc
fviz_dend(clusters.hc, cex = 0.5, k = 4, palette = "jco")



#######################
# BICLUSTERING: BIMAX #
#######################

# Convertimos el dataset a matriz
dataMatrix <- as.matrix(dataSample) 
summary(dataMatrix) 

# 1) Ejecutar metodo de Biclustering
bicBimax <- biclust(dataMatrix, method=BCBimax(), minc = 2, minr = 2)

# 2) Visualizar resultados
bicBimax
bicluster(dataMatrix,bicBimax) # Visualizar los biclusters
bicluster(dataMatrix,bicBimax)[1] # Visualizar el primer bicluster
parallelCoordinates(dataMatrix, bicResult=bicBimax, number = 1) 
drawHeatmap(dataMatrix,bicResult=bicBimax, number = 1) 
drawHeatmap2(dataMatrix,bicResult=bicBimax, number = 1) 

# 3) Resumen
summary(bicBimax)

###########################
# BICLUSTERING: CC (CCCB) #
###########################

# 1) Ejecutar metodo de Biclustering
bicCC <- biclust(dataMatrix, method=BCCC(), number=10, alpha=1.5, delta=0.1)

# 2) Visualizar resultados
bicluster(dataMatrix, bicCC)
drawHeatmap(dataMatrix, bicResult=bicCC, number=1)
drawHeatmap2(dataMatrix, bicResult=bicCC, number=1)

# 3) Resumen
summary(bicCC)

