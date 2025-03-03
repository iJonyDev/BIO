# Lee el archivo omitiendo la primera fila, usando la segunda como encabezados
temp <- read.table("/Users/jonathanquishpe/JoniDev/BIO/GSE62301_non-normalized.txt", 
                   sep = "\t", 
                   header = FALSE,  
                   skip = 1)        # Omite la primera fila

# Asigna los nombres de las columnas
names(temp) <- as.character(unlist(temp[1, ]))

# Elimina la fila que contiene los nombres de columnas (ahora es la primera fila)
dataset <- temp[-1, ]

str(dataset) # Estructura del dataframe


# Identificar el número de columnas
num_cols <- ncol(dataset)

# Convertir todas las columnas excepto la primera y la última a numérico
for (i in 2:(num_cols-1)) {
  dataset[,i] <- as.numeric(as.character(dataset[,i]))
}

# Operaciones basicas sobre dataframe

dim(dataset) # Numero de filas y columnas
colnames(dataset) # Nombres de las columnas
is.matrix(dataset) # Verificar si es una matriz
is.data.frame(dataset) # Verificar si es un data frame
str(dataset) # Estructura del dataframe

# Estadisticas básicas

summary(dataset[,2]) # Resumen de la columna 2
sd(dataset[,2]) # Desviación estándar de la columna 2
var(dataset[,2]) # Varianza de la columna 2
range(dataset[,2]) # Rango de la columna 2
quantile(dataset[,2]) # Cuantiles de la columna 2
IQR(dataset[,2]) # Rango intercuartil de la columna 2


dataset[1:5,2] # Primeras 5 filas de la columna 2

# Histograma
hist(dataset[1:500,2], main = "Histograma", xlab = "Variable 2", ylab = "Frequency", breaks = 20)

# Densidad
plot(density(dataset[1:100,2]), main = "Densidad", xlab = "Variable 2", ylab = "Density")

# Q-Q plot
qqnorm(dataset[1:100,2])





