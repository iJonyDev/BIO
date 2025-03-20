# ----------Instalacion de paquetes y carga de datos-----------
# Instalacion de paquetes necesarios para el preprocesamiento
install.packages('caret')
install.packages("fastICA")
install.packages("mlbench")

# -------------------- Cargar el dataset --------------------
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
str(dataset) # Estructura del dataframe despues de la conversion




#-------------- Seleccionar variables de interés --------------

# Crear un data frame con las variables de interés
datos <- data.frame(
  var1 = dataset[,2],
  var2 = dataset[,3],
  var3 = dataset[,4]
)

# Histogramas de las variables de interés
# Configurar el panel para mostrar 3 gráficos en una fila
par(mfrow=c(1,3))

# Crear los histogramas
hist(datos$var1, 
     main="Histograma var1", 
     xlab="Valor", 
     ylab="Frecuencia",
     col="skyblue")

hist(datos$var2, 
     main="Histograma var2", 
     xlab="Valor", 
     ylab="Frecuencia",
     col="salmon")

hist(datos$var3, 
     main="Histograma var3", 
     xlab="Valor", 
     ylab="Frecuencia",
     col="lightgreen")

# Restaurar la configuración del panel
par(mfrow=c(1,1))


# Boxplot de las variables de interés
boxplot(datos, main="Boxplot",
               col=c("skyblue", "salmon", "lightgreen"),
               xlab="Variables", ylab="Valor",
               names=names(datos))

# Resumen con estadisticas de las variables de interés
summary(datos) 


# ----- Identificar valores atípicos (outliers) -----
# Obtener cuantitativamente el rango intercuantílico (IQR) 
# y los límites atípicos que determinan si un valor es un outlier o no
# Para ello, el IQR es la distancia entre Q3 y Q1.
var1_iqr <- IQR(datos$var1)
var2_iqr <- IQR(datos$var2)
var3_iqr <- IQR(datos$var3)
var1.Q1 <- quantile(datos$var1, 0.25)
var2.Q1 <- quantile(datos$var2, 0.25)
var3.Q1 <- quantile(datos$var3, 0.25)
var1.Q3 <- quantile(datos$var1, 0.75)
var2.Q3 <- quantile(datos$var2, 0.75)
var3.Q3 <- quantile(datos$var3, 0.75)
# Un valor se considerará outlier leve si el valor dista más de 1,5 veces el rango intercuantílico
# por debajo de Q1 o por encima de Q3.
leveDown1 <- var1.Q1 - 1.5*var1_iqr
leveUp1 <- var1.Q3 + 1.5*var1_iqr
leveDown2 <- var2.Q1 - 1.5*var2_iqr
leveUp2 <- var2.Q3 + 1.5*var2_iqr
leveDown3 <- var3.Q1 - 1.5*var3_iqr
leveUp3 <- var3.Q3 + 1.5*var3_iqr
# Un valor outlier extremo se considera cuando su valor dista más de 3 veces el IQR por debajo
# de Q1 o por encima de Q3.
extremeDown1 <- var1.Q1 - 3*var1_iqr
extremeUp1 <- var1.Q3 + 3*var1_iqr
extremeDown2 <- var2.Q1 - 3*var2_iqr
extremeUp2 <- var2.Q3 + 3*var2_iqr
extremeDown3 <- var3.Q1 - 3*var3_iqr
extremeUp3 <- var3.Q3 + 3*var3_iqr


# ------------- Recodificar los valores atípicos(outliers)  ------------
variable1 <- datos$var1
is.na(variable1) <- variable1 >= extremeUp1
is.na(variable1) <- variable1 >= leveUp1

variable2 <- datos$var2
is.na(variable2) <- variable2 >= extremeUp2
is.na(variable2) <- variable2 >= leveUp2

variable3 <- datos$var3
is.na(variable3) <- variable3 >= extremeUp3
is.na(variable3) <- variable3 >= leveUp3

# Crear una lista con las variables recodificadas
datos_recodificados <- list(
  variable1 = variable1,
  variable2 = variable2,
  variable3 = variable3
)

# Asignar los nombres originales de las columnas
names(datos_recodificados) <- colnames(datos)[1:3]

# Configurar el panel para mostrar 3 gráficos en una fila
par(mfrow=c(1,3))

# Crear los histogramas
hist(datos_recodificados[[1]], 
     main=paste(names(datos_recodificados)[1], " recodificada"), 
     xlab="Valor", 
     ylab="Frecuencia",
     col="skyblue")

hist(datos_recodificados[[2]], 
     main=paste(names(datos_recodificados)[2], " recodificada"), 
     xlab="Valor", 
     ylab="Frecuencia",
     col="salmon")

hist(datos_recodificados[[3]], 
     main=paste(names(datos_recodificados)[3], " recodificada"), 
     xlab="Valor", 
     ylab="Frecuencia",
     col="lightgreen")

# Restaurar la configuración del panel
par(mfrow=c(1,1))

# Boxplot de variables recodificadas
boxplot(datos_recodificados,
        main="Boxplot de las variables recodificadas",
        col=c("skyblue", "salmon", "lightgreen"),
        xlab="Variables", 
        ylab="Valor",
        names=names(datos_recodificados))
#Imprimimos el resumen de las variables recodificadas
for (i in 1:3) {
  cat("\n----- Resumen ",names(datos_recodificados)[i], " Recodificado -----\n")
  print(summary(datos_recodificados[[i]]))
}




# ------------- Eliminar outliers de las variables de interés ------------
for (i in 1:3) {
  # Eliminar valores NA de la variable i
  cat("\n----- Resumen ",names(datos_recodificados)[i], " sin outliers -----\n")
  datos_recodificados[[i]] <- datos_recodificados[[i]][!is.na(datos_recodificados[[i]])]
  print(summary(datos_recodificados[[i]]))
}
