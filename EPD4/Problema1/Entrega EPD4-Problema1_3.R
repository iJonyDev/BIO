# --------------- Resumen Antes y Despues de Transformar --------------------- 
# Extraer los datos transformados
preprocessed_data <- lapply(preprocessed_results, function(result) {
  return(result$transformed)
})

# Extraer los parámetros si se necesitan para futuras predicciones
preprocess_params <- lapply(preprocessed_results, function(result) {
  return(result$params)
})

# ---- Histogramas Antes y Despues de Transformar ----

# Determinar cuántas variables hay
n_vars <- length(preprocessed_data)
colors <- c("skyblue", "salmon", "lightgreen")

# Configurar el panel para mostrar los histogramas
par(mfrow=c(n_vars, 2) # 2 columnas por variable (original y preprocesado)
    , mar=c(4,4,2,1)) # Margenes de las gráficas (abajo, izquierda, arriba, derecha)
# Representacion visual de la configuración del panel para cada variable
# +---------------------+
#   |                     |
#   |    2 líneas (top)   |
#   |                     |
#   +---------------------+
#   |  |               |  |
#   |  |               |  |
#   |  |               |  |
#   |  |               |  |
#   |4 |    GRÁFICO    |1 |
#   |  |               |  |
#   |lín|               |lín|
#   |  |               |ea |
#   |  |               |  |
#   +---------------------+
#   |                     |
#   |   4 líneas (bottom) |
#   |                     |
#   +---------------------+

# Para cada variable mostrar histogramas
for (i in 1:n_vars) {
  # Obtener nombre de variable o asignar uno si no tiene
  var_name <- ifelse(!is.null(names(preprocessed_data)), 
                     names(preprocessed_data)[i], 
                     paste("var", i))
  
  # Color para esta variable
  col_i <- colors[i %% length(colors) + 1]
  
  # Histograma de datos no procesados
  hist(data_no_preprocessed[[i]], 
       col=col_i, 
       main=paste("No Preprocesado:", var_name), 
       xlab="Valor", 
       ylab="Frecuencia")
  
  # Histograma de datos preprocesados
  hist(preprocessed_data[[i]], 
       col=col_i, 
       main=paste("Preprocesado",tipo_preprocesamiento,":", var_name), 
       xlab="Valor", 
       ylab="Frecuencia")
}
# Restaurar configuración de panel por defecto
par(mfrow=c(1,1))

# ---- Boxplots Antes y Despues de Transformar ----
# Configurar el panel para mostrar los boxplots
par(mfrow=c(n_vars, 2), mar=c(2,4,2,1))

# Para cada variable, mostrar boxplots
for (i in 1:n_vars) {
  # Obtener nombre de variable o asignar uno si no tiene
  var_name <- ifelse(!is.null(names(preprocessed_data)), 
                     names(preprocessed_data)[i], 
                     paste("var", i))
  
  # Color para esta variable
  col_i <- colors[i %% length(colors) + 1]
  
  # Boxplot de datos no procesados
  boxplot(data_no_preprocessed[[i]], 
          col=col_i, 
          main=paste("No Preprocesado:", var_name), 
          xlab="Valor", 
          ylab="Valor")
  
  # Boxplot de datos preprocesados
  boxplot(preprocessed_data[[i]], 
          col=col_i, 
          main=paste("Preprocesado",tipo_preprocesamiento,":", var_name), 
          xlab="Valor", 
          ylab="Valor")
  
  # Imprimir resúmenes
  cat("\n--- ", var_name, " ---\n")
  cat("No Preprocesado: \n")
  print(summary(data_no_preprocessed[[i]]))
  cat("Preprocesado",tipo_preprocesamiento, ":\n")
  print(summary(preprocessed_data[[i]]))
}

# Restaurar configuración de panel por defecto
par(mfrow=c(1,1))


