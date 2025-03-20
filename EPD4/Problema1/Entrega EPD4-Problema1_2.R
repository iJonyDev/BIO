# ---------- Identificación de Outliers en RNA-Seq -----------
# Calcular el rango intercuartílico (IQR) para cada gen
iqr_values <- apply(dgeObj$counts, 1, IQR)

# Calcular los cuartiles (Q1 y Q3) para cada gen
q1_values <- apply(dgeObj$counts, 1, quantile, probs = 0.25)
q3_values <- apply(dgeObj$counts, 1, quantile, probs = 0.75)

# Definir los límites para outliers leves y extremos
leveDown <- q1_values - 1.5 * iqr_values
leveUp <- q3_values + 1.5 * iqr_values
extremeDown <- q1_values - 3 * iqr_values
extremeUp <- q3_values + 3 * iqr_values

# Identificar outliers en los recuentos de genes
outliers_leves <- dgeObj$counts < leveDown | dgeObj$counts > leveUp
outliers_extremos <- dgeObj$counts < extremeDown | dgeObj$counts > extremeUp

# Resumen de outliers
cat("Número de outliers leves por gen:\n")
print(rowSums(outliers_leves[1:10,])) # Mostrar solo los primeros 10 genes

cat("Número de outliers extremos por gen:\n")
print(rowSums(outliers_extremos[1:10,])) # Mostrar solo los primeros 10 genes

# ------------- Recodificar los valores atípicos (outliers)  ------------
# Crear una copia de los recuentos para recodificar los outliers
counts_recodificados <- dgeObj$counts

# Recodificar outliers leves como NA
counts_recodificados[outliers_leves] <- NA

# Recodificar outliers extremos como NA
counts_recodificados[outliers_extremos] <- NA

# Actualizar el objeto DGEList con los recuentos recodificados
dgeObjRecodificado <- dgeObj
dgeObjRecodificado$counts <- counts_recodificados

# Resumen de los recuentos recodificados
cat("\n----- Resumen de Recuentos Recodificados -----\n")
print(summary(counts_recodificados))

# ------------- Eliminar outliers de los recuentos de genes ------------
# Eliminar valores NA de los recuentos recodificados
counts_sin_outliers <- na.omit(dgeObjRecodificado$counts)

# Actualizar el objeto DGEList con los recuentos sin outliers
dgeObjSinOutliers <- dgeObj
dgeObjSinOutliers$counts <- counts_sin_outliers

# Resumen de los recuentos sin outliers
cat("\n----- Resumen de Recuentos sin Outliers -----\n")
print(summary(counts_sin_outliers))