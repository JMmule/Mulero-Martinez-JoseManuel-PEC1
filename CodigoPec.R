#Descarga de los datos
url<-"https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv"
dataset<-"human_cachexia.csv"
download.file(url,dataset)
datos<-read.csv(dataset)
colnames(datos)