# ---------- Instalación de paquetes y carga de datos -----------
# Instalación de paquetes necesarios para el preprocesamiento RNA-Seq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

# Install all Bioconductor packages
BiocManager::install(c("edgeR", "limma", "Glimma", "org.Mm.eg.db","cluster","biclust","pheatmap"))

# Install CRAN packages
install.packages(c("gplots", "RColorBrewer", "factoextra", "caret"))

# Load all required libraries
libraries <- c("limma", "edgeR", "Glimma", "gplots", "caret", 
               "org.Mm.eg.db", "RColorBrewer", "factoextra","mlbench","cluster","biclust","pheatmap")
lapply(libraries, library, character.only = TRUE)

# Cargar el dataset de RNA-Seq
seqdata <- read.delim("/Users/jonathanquishpe/JoniDev/BIO/EPD5/Problema2/GSE261149_gene_count_matrix.txt", stringsAsFactors = FALSE)
head(seqdata)
dim(seqdata)
samplesCount = 9 

# Formatear los datos
# Eliminar la primera columna (identificadores de genes)
countdata <- seqdata[,-1]

# Almacenar los identificadores de genes como nombres de fila
rownames(countdata) <- seqdata[,1]
View(countdata)

####################
# PREPROCESAMIENTO #
####################
# Se aplica el metodo de estandarizacion ya que es el que mejor funciona para los algoritmos de clustering y biclustering.
preprocessParams <- preProcess(countdata[,1:samplesCount], method=c("center", "scale"))
countdata[,1:samplesCount] <- predict(preprocessParams, countdata[,1:samplesCount])
summary(countdata)


##########################
# CLUSTERING: DISTANCIAS #
##########################
distanciasDataSet <- get_dist(t(countdata[,1:samplesCount]), method = "pearson")             
head(round(as.matrix(distanciasDataSet), 2))[, 1:samplesCount]
distanciasDataSet
fviz_dist(distanciasDataSet, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Clustering jerárquico
hc <- hclust(distanciasDataSet, method = "ward.D2")

# Visualizar dendrograma
fviz_dend(hc, k = 3, cex = 0.5, palette = "jco")

#######################
# BICLUSTERING: BIMAX #
#######################

# Convertimos el dataset a matriz
dataMatrix <- as.matrix(countdata) 
summary(dataMatrix)

# 1) Ejecutar metodo de Biclustering
bicBimax <- biclust(dataMatrix, method=BCBimax(), minc = 2, minr = 2)

# 2) Visualizar resultados
bicBimax
bicluster(dataMatrix,bicBimax) #Obtendremos 100
bicluster(dataMatrix,bicBimax)[1] #Observamos el bicluster 1
parallelCoordinates(dataMatrix, bicResult=bicBimax, number = 1) 
drawHeatmap(dataMatrix,bicResult=bicBimax, number = 1) 
drawHeatmap2(dataMatrix,bicResult=bicBimax, number = 1) 

###########################
# BICLUSTERING: CC (CCCB) #
###########################

# Biclustering CC
bicCC <- biclust(dataMatrix, method=BCCC(), number=10, alpha=1.5, delta=0.1)

# Visualización
bicluster(dataMatrix, bicCC)
drawHeatmap(dataMatrix, bicResult=bicCC, number=1)
drawHeatmap2(dataMatrix, bicResult=bicCC, number=1)


#################################
# BICLUSTERING: Plaid (BCPlaid) #
#################################

# Biclustering Plaid
bicPlaid <- biclust(dataMatrix, method=BCPlaid(), cluster="b", fit.model = y ~ m + a + b,
                    background = TRUE, row.release = 0.7, col.release = 0.7,
                    shuffle = 3, back.fit = 5, max.layers = 10, iter.startup = 5,
                    iter.layer = 10)

# Visualización
bicluster(dataMatrix, bicPlaid)
drawHeatmap(dataMatrix, bicResult=bicPlaid, number=1)
drawHeatmap2(dataMatrix, bicResult=bicPlaid, number=1)


#######################################
# BICLUSTERING: Spectral (BCSpectral) #
#######################################

# Biclustering Spectral
bicSpectral <- biclust(dataMatrix, method=BCSpectral(), numberOfEigenvalues=3, 
                       minr=2, minc=2, withinVar=1)

# Visualización
bicluster(dataMatrix, bicSpectral)
drawHeatmap(dataMatrix, bicResult=bicSpectral, number=1)
drawHeatmap2(dataMatrix, bicResult=bicSpectral, number=1)


