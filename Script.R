# Instalar paquete
#install.packages("readr")

# Llamar librerias
library(readr)
library(knitr)

# Definir la URL del archivo
url = "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv"

# Leer el archivo CSV
data = read_csv(url)

# Verificar dataframe
kable(head(data))

# Instalar paquetes
#install.packages("BiocManager")
#BiocManager::install("SummarizedExperiment")

# Llamar libreria
library(SummarizedExperiment)

# Crear clase SummarizedExperiment

# Crear argumento assay (datos de expresion de metabolitos)
ensayo = t(as.matrix(data[,!(colnames(data) %in% c("Patient ID", "Muscle loss"))]))
colnames(ensayo) = data$`Patient ID` # Asignar nombre

# Crear argumento rowdata (metadatos de los metabolitos)
rowdata = data.frame(metabolito=tail(colnames(data),-2))

# Crear argumento coldata (metadatos de las muestras)
coldata = data.frame(grupo=data$`Muscle loss`, row.names=data$`Patient ID`)

# Añadir metadatos adicionales sobre el experimento
metadatos = list(experiment_date = "2010-08-01", protocol = "detección de metabolitos en orina a través del espectrómetro de RMN Varian INOVA de 600 MHz, sonda de triple eje con gradiente de 5-mm HCN; e identificados a través del software Chenomx NMRSuite 4.6", researchers = "Roman Eisner, Cynthia Stretch, Thomas Eastman, Jianguo Xia, David Hau, Sambasivarao Damaraju, Russell Greiner, David S. Wishart & Vickie E. Baracos", study_name = "Learning to predict cancer-associated skeletal muscle wasting from 1H-NMR profiles of urinary metabolites")

# Crear SummarizedExperiment
experimento = SummarizedExperiment(assay=list(ensayo=ensayo), rowData=rowdata, colData=coldata, metadata=metadatos)

# Visualizar
print(experimento)

# Instalar paquetes
#install.packages("kableExtra")

# Llamar librerias
library(dplyr)

# Definir tabla con la informacion del atributo metadata
infoDF = data.frame(matrix(rep(NA, 2 * length(metadatos)), ncol = 2))
for (i in 1:length(metadatos)) {
  infoDF[i, 1] = names(metadatos)[i]
  infoDF[i, 2] = metadatos[i]
}
colnames(infoDF) = c("Campo", "Descripción")

# Imprimir
infoDF %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# Imprimir primeras 6 filas de atributo assay
head(assay(experimento)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# Imprimir primeras 6 filas de atributo colData
head(colData(experimento)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# Imprimir ultimass 6 filas de atributo colData
tail(colData(experimento)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# Imprimir primeras 6 filas de atributo rowData
head(rowData(experimento)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# Obtener el numero de valores faltantes en la matriz
sum(is.nan(assay(experimento)))

# Realizar resumen estadistico descriptivo de los datos por metabolito
summary(t(assay(experimento)))

# Obtener medias y medianas para cada metabolito
medias = apply(t(assay(experimento)),2,mean)
medianas = apply(t(assay(experimento)),2,median)

# Calcular las diferencias
diferencias = medias-medianas

# Ordenar de menor a mayor
diferencias = diferencias[order(diferencias)]

# Imprimir
kable(diferencias)

# Calcular la desviacion estandar para cada metabolito
variabilidad = apply(t(assay(experimento)), 2, sd)

# Graficar diagrama de barras de la variabilidad de los primeros 31 metabolitos
barplot(head(variabilidad,31), main="Variabilidad de los primeros 31 metabolitos", col="blue", las=2, ylim=c(0,6000))
# Graficar diagrama de barras de la variabilidad de los ultimos 32 metabolitos
barplot(tail(variabilidad,-32), main="Variabilidad de los ultimos 32 metabolitos", col="blue", las=2, ylim=c(0,6000))

# Instalar paquete
# install.packages("pheatmap")

# Llamar libreria
library(pheatmap)

# Normalizar los datos
data_norm = scale(t(assay(experimento)))

# Obtener mapa de calor
pheatmap(data_norm, main="Mapa de calor de los metabolitos y muestras", show_rownames=FALSE, show_colnames=FALSE)

# Obtener la matriz de correlacion
cor_matrix = cor(t(assay(experimento)))

# Crear mapa de calor de las correlaciones
pheatmap(cor_matrix, main="Mapa de calor de correlación entre metabolitos", clustering_distance_rows="correlation", clustering_distance_cols="correlation", , show_rownames=FALSE, show_colnames=FALSE)

# Llamar libreria
library(ggplot2)

# Obtener los grupos de muestras
grupos = colData(experimento)

# Realizar analisis de componentes principales
pca = prcomp(t(assay(experimento)), scale.=TRUE)

# Obtener datos del PCA en data frame
pca_df = data.frame(pca$x)

# Graficar PCA (primeros dos componentes)
p1 = ggplot(pca_df, aes(PC1, PC2, color=grupos$grupo)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC1", y="PC2") +
  theme_minimal()
# Imprimir
print(p1)

# Graficar PCA (primero y tercer componente)
p2 = ggplot(pca_df, aes(PC1, PC3, color=grupos$grupo)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC1", y="PC3") +
  theme_minimal()
# Imprimir
print(p2)

# Graficar PCA (segundo y tercer componente)
p3 = ggplot(pca_df, aes(PC2, PC3, color=grupos$grupo)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC2", y="PC3") +
  theme_minimal()
# Imprimir
print(p3)

# Obtener los nombres de los metabolitos
metabolitos = rownames(rowData(experimento))

# Definir categorias de metabolitos
categorias = list(
  "Metabolismo de Carbohidratos" = c("Glucose", "Pyruvate", "Lactate", "Acetate", "Citrate", "Fumarate", "Sucrose", "Xylose", "cis-Aconitate", "trans-Aconitate"),
  
  "Metabolismo de Aminoacidos" = c("Alanine", "Glutamine", "Glycine", "Serine", "Leucine", "Isoleucine", "Lysine", "Valine", "Threonine", "Histidine", "2-Aminobutyrate", "3-Aminoisobutyrate", "Guanidoacetate", "Methylguanidine", "Asparagine", "Betaine", "Ethanolamine", "Glycolate", "N,N-Dimethylglycine", "O-Acetylcarnitine", "Pyroglutamate", "Taurine", "Trigonelline", "Tryptophan", "Tyrosine", "pi-Methylhistidine", "tau-Methylhistidine", "2-Hydroxyisobutyrate", "3-Hydroxyisovalerate"),
  
  "Metabolismo de Lipidos" = c("Carnitine", "Acetone", "Adipate", "Trimethylamine N-oxide"),
  
  "Metabolismo de Acidos Organicos y Ciclo de Krebs" = c("Succinate", "Citrate", "Fumarate", "2-Oxoglutarate", "Acetate", "Acetone", "Formate", "Fucose", "cis-Aconitate", "trans-Aconitate", "Tartrate"),
  
  "Metabolismo de Nucleotidos y Bases Nitrogenadas" = c("Hypoxanthine", "Uracil", "Quinolinate"),
  
  "Vias de Detoxificacion y Microbioma" = c("Indoxyl sulfate", "Hippurate", "3-Indoxylsulfate", "4-Hydroxyphenylacetate"),
  
  "Otras Vias Metabolicas" = c("1,6-Anhydro-beta-D-glucose", "Dimethylamine", "Methylamine", "Pantothenate", "myo-Inositol", "1-Methylnicotinamide"),
  
  "Biomarcadores de Estres y Energia" = c("Creatine", "Creatinine", "3-Hydroxybutyrate")
)

# Encontrar la categoría correspondiente para cada metabolito
categoria_df = data.frame(
  metabolito=metabolitos,
  categoria=sapply(metabolitos, function(metabolito) {
    category=NA
    for (cat in names(categorias)) {
      if (metabolito %in% categorias[[cat]]) {
        category=cat
        break
      }
    }
    return(category)
  })
)

# Realizar analisis de componentes principales
pca_metabolitos = prcomp(assay(experimento), scale.=TRUE)

# Obtener datos del PCA en data frame
pca_df_metabolitos = data.frame(pca_metabolitos$x)

# Graficar PCA (primeros dos componentes)
p1 = ggplot(pca_df_metabolitos, aes(PC1, PC2, color=categoria_df$categoria)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC1", y="PC2") +
  theme_minimal()
# Imprimir
print(p1)

# Graficar PCA (primer y tercer componente)
p2 = ggplot(pca_df_metabolitos, aes(PC1, PC3, color=categoria_df$categoria)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC1", y="PC3") +
  theme_minimal()
# Imprimir
print(p2)

# Graficar PCA (segundo y tercer componente)
p3 = ggplot(pca_df_metabolitos, aes(PC2, PC3, color=categoria_df$categoria)) +
  geom_point() +
  labs(title="PCA de metabolitos y muestras", x="PC2", y="PC3") +
  theme_minimal()
# Imprimir
print(p3)

