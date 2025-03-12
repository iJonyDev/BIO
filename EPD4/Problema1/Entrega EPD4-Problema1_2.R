# ---- Preprocesamiento de los datos ----

# Crear una lista con las variables de interés sin preprocesar
data_no_preprocessed <- datos_recodificados


#####################################################
#### Procesado Tipo: Transformación Escala ##########
#####################################################
library(caret)
tipo_preprocesamiento <- "Transfomacion Escala"
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("scale")) # Transmacion por escala
  transformed <- predict(params, var_df) # Transformar los datos 
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

#####################################################
#### Procesado Tipo: Transformación Central #########
#####################################################
library(caret)
tipo_preprocesamiento <- "Transf.Central"
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("center")) # Transmacion central
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
#### Procesado Tipo: Transformación de Estandarizacion #########
################################################################
library(caret)
tipo_preprocesamiento <- "Transf.Estandarizacion "
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("center", "scale")) # Transmacion por estandarizacion
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
#### Procesado Tipo: Transformación de Normalizacion ###########
################################################################
library(caret)
tipo_preprocesamiento <- "Transf.Normalizacion "
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("range")) # Transmacion de normalizacion
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
######## Procesado Tipo: Transformación de Box-Cox #############
################################################################
library(mlbench)
library(caret)
tipo_preprocesamiento <- "Transf.Box-Cox "
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("BoxCox")) # Transmacion de BoxCox
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
###### Procesado Tipo: Transformación de Yeo-Johnson ###########
################################################################
library(mlbench)
library(caret)
tipo_preprocesamiento <- "Transf.Yeo-Johnson "
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("YeoJohnson")) # Transmacion de Yeo-Johnson
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
########## Procesado Tipo: Transformación PCA ##################
################################################################
library(mlbench)
tipo_preprocesamiento <- "Transf.PCA " # Análisis de Componentes Principales
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("center", "scale", "pca")) # Transmacion de PCA
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

################################################################
########## Procesado Tipo: Transformación ICA ##################
################################################################
library(mlbench)
library(caret)
library(fastICA)
tipo_preprocesamiento <- "Transf.ICA " # Analisis de Componentes Independientes
# Procesar todas las variables a la vez usando lapply
preprocessed_results <- lapply(data_no_preprocessed, function(variable) {
  var_df <- data.frame(x = variable) # Crear un data frame con la variable de interés
  # La transformación scale calcula la desviación estándar de un atributo
  # y divide cada valor por esa desviación estándar.
  params <- preProcess(var_df, method=c("center", "scale", "ica"), n.comp=5) # Transmacion de ICA
  transformed <- predict(params, var_df) # Transformar los datos
  # Devolver los parámetros y los datos transformados
  return(list(
    params = params,
    transformed = transformed$x
  ))
})

