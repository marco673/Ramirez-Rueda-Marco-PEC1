## Cargamos la base y las librerias
library(readr)
# Cargar el archivo de texto
library(readr)
df <- read_delim("metabolito.csv",
                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
##-------------------------------------------------------------------------------------
## Crear metadatos para las columnas
col_data <- data.frame(
METABOLOMICS.WORKBENCH.STUDY_ID = rep(df$M_T[1], ncol(df) - 1),
PROJECT.PROJECT_TITLE = rep("Fisiología del factor inducible por hipoxia-1a (HIF1a)",
                              ncol(df) - 1),)
SUBJECT.SUBJECT_TYPE = rep("Mouse", ncol(df) - 1)
  ##-------------------------------------------------------------------------------------
  
library(IRanges)
library(SummarizedExperiment)
library(S4Vectors)
library(MatrixGenerics)
library(matrixStats)
library(GenomicRanges)
library(stats4)
library(Biobase)
library(MatrixGenerics)
library(matrixStats)
  # Convertir los datos en un objeto de tipo SummarizedExperiment
contenedor <- SummarizedExperiment(
assays = list(counts = as.matrix(df[, -1])),
rowData = DataFrame(Gene_ID = rownames(df)),
colData = col_data)

library(dplyr)
 # Definir las variables de interés
variables_interes <- c("HIFFloxFPBS1-C18Neg",
                           "HIFFloxFPBS2-C18Neg",
                           "HIFFloxFPBS5-C18Neg",
                           "HIFmsdFPBS1-C18Neg",
                           "HIFmsdFPBS2-C18Neg")
# Crear un subconjunto del data frame solo con las variables seleccionadas
subset_df <- df %>%
select(M_T, all_of(variables_interes))
head(subset_df)
    ##-------------------------------------------------------------------------------------
library(corrplot)
numeric_subset <- subset_df %>% select(where(is.numeric))
cor_matrix <- cor(numeric_subset, use = "complete.obs")
corrplot(cor_matrix, method = "circle", type = "upper",
tl.col = "black", tl.srt = 45,addCoef.col = "white")
    ##-------------------------------------------------------------------------------------
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
for (var in variables_interes) {
hist(numeric_subset[[var]],
include.lowest = TRUE,
col = rainbow(10),
main = paste("Histograma de", var, sep = "\n"), xlab = var)}
    ##-------------------------------------------------------------------------------------
par(mfrow = c(2, 3))
for (var in variables_interes) {
boxplot(numeric_subset[[var]],
      main = paste("Boxplot de", var, sep = "\n"),
      ylab = var,
      col = rainbow(10))
    }

library(purrr)
library(cluster)
    # Seleccionar las primeras 100 observaciones
subset_numeric <- numeric_subset[1:100, ]
# Estandarizar variables
clusters <- scale(subset_numeric)
# Matriz de similaridades
simil <- dist(x = clusters, method = "euclidean", diag = TRUE)

enla.simple<-hclust(d = simil, method = "single")
#gráfico del dendrograma
plot(enla.simple,cex=0.6,hang=-1,main = "Dendrograma de agrupamiento jeráquico con enlace simple"
     ,xlab="Observaciones", ylab="Distancia Euclideana")

rect.hclust(enla.simple, k=10,
            border = c("blue","purple","yellow","green"))

enla.completo<-hclust(d = simil, method = "complete")
#gráfico del dendrograma
plot(enla.completo,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con enlace completo",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(enla.completo, k=5,
            border = c("blue","purple","yellow","green"))

enla.promedio<-hclust(d = simil, method = "average")
#gráfico del dendrograma
plot(enla.promedio,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con enlace promedio",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(enla.promedio, k=5,
            border = c("blue","purple","yellow","green"))

##-------------------------------------------------------------------------------------
wardD<-hclust(d = simil, method = "ward.D")
#gráfico del dendrograma
plot(wardD,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con el método de wardD",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(wardD, k=4,
            border = c("blue","purple","yellow","green"))
##-------------------------------------------------------------------------------------
wardD2<-hclust(d = simil, method = "ward.D2")
#gráfico del dendrograma
plot(wardD2,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con el método de wardD2",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(wardD2, k=5,
            border = c("blue","red","yellow","green"))

##-------------------------------------------------------------------------------------
centroide<-hclust(d = simil, method = "centroid")
#gráfico del dendrograma
plot(centroide,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con el método del centroide",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(centroide, k=5,
            border = c("blue","purple","yellow","green"))
##-------------------------------------------------------------------------------------
mediana<-hclust(d = simil, method = "median")
#gráfico del dendrograma
plot(mediana,
     cex=0.6,
     hang=-1,
     main = "Dendrograma de agrupamiento jeráquico con el método de la mediana",
     xlab="Observaciones", ylab="Distancia Euclideana")
rect.hclust(mediana, k=6,
            border = c("blue","purple","yellow","green"))

