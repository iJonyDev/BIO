dir <- getwd() # Directorio de trabajo
# Lee el archivo tomando la primera fila como encabezado.
dataset_p2 <- read.table("/Users/jonathanquishpe/Documents/BIO/GSE290268_RNAseq_FPKM_data.txt", 
                   sep = "\t", 
                   header = TRUE)  

# Operaciones basicas sobre dataframe

dim(dataset) # Numero de filas y columnas
colnames(dataset) # Nombres de las columnas
is.matrix(dataset) # Verificar si es una matriz
is.data.frame(dataset) # Verificar si es un data frame
str(dataset_p2) # Estructura del dataframe

# Estadisticas básicas

summary(dataset_p2[,2]) # Resumen de la columna 2
sd(dataset_p2[,2]) # Desviación estándar de la columna 2
var(dataset_p2[,2]) # Varianza de la columna 2
range(dataset_p2[,2]) # Rango de la columna 2
quantile(dataset_p2[,2]) # Cuantiles de la columna 2
IQR(dataset_p2[,2]) # Rango intercuartil de la columna 2


dataset[1:5,2] # Primeras 5 filas de la columna 2

# Histograma
hist(dataset_p2$UT1[1:500], main = "Histograma", xlab = "UT1", ylab = "Frequency", breaks = 20)

# Densidad
plot(density(dataset_p2$UT1[1:50]), main = "Densidad", xlab = "UT1", ylab = "Density")

# Q-Q plot
qqnorm(dataset_p2[1:1000,2])
