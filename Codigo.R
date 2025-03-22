#Descarga de los datos
url<-"https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv"
dataset<-"human_cachexia.csv"
download.file(url,dataset)
datos<-read.csv(dataset)
colnames(datos)

#Visualización del dataset
head(datos)

#Transformación de los datos
assayData<-as.matrix(datos[,-c(1:2)])
assayData<-t(assayData)
rownames(assayData)<-colnames(datos)[-c(1:2)]
colnames(assayData)<-datos$Patient.ID

#Verificación de los cambios
head(assayData)
dim(assayData)      
rownames(assayData) 
colnames(assayData) 

#Creación del DataFrame con los colData
columnData<-DataFrame(Muscle.loss=datos$Muscle.loss)
rownames(columnData)<-datos$Patient.ID

#Verificar que se ha creado correctamente
head(columnData)

#Creación del objeto de clase SummarizedExperiment
se<-SummarizedExperiment(assays=list(counts=assayData), colData=columnData)

#Creación de archivos
save(se,file="summarized_experiment_cachexia.Rda")
write.csv(assay(se),file="datos_cachexia.csv",row.names=TRUE)
write.csv(colData(se),file="metadatos_cachexia.csv",row.names=TRUE)

#Análisis exploratorio
#Estructura del objeto
head(se)
dim(se)
rownames(se)
colData(se)

#Valores faltantes
sum(is.na(assay(se)))

#Resumen estadístico
#Para que nos devuelva el resumen estadístico por metabolito tenemos que usar apply() ya que si usamos
#directamente summary() nos devolverá el resumen estadístico por columnas (que son los pacientes).
apply(assay(se),1,summary)  

#Heatmap con expresión de metabolitos entre pacientes y controles
matriz<-as.matrix(assay(se)) #Complexheatmap solo admite matrices asi que convertimos los datos en una matriz.
escalado<-t(scale(t(matriz))) #Escalamos los valores por filas (metabolitos) para que cada uno tenga su propia escala y ver claramente las diferencias.
colores<-colorRamp2(c(-2,0,2),c("blue","white","red"))
heatmap<-Heatmap(
  escalado,name="Expresión de metabolitos", 
  col=colores,
  show_row_names=TRUE,
  show_column_names=TRUE,
  row_names_gp=gpar(fontsize=6,rot=90),
  column_names_gp=gpar(fontsize=6,rot=45),                        
  cluster_rows=TRUE,  
  cluster_columns=TRUE,  
  top_annotation = HeatmapAnnotation(
    df=data.frame(Group=colData(se)$Muscle.loss),
    col=list(Group=c("cachexic"="yellow","control"="green")),
    annotation_legend_param=list(title="Grupo")),
  column_split=colData(se)$Muscle.loss,
  column_title="Muestras", 
  column_names_side="top",
  heatmap_legend_param=list(title="Z-score",title_position="topcenter") 
)
pdf("heatmap.pdf",width=12,height=10)
draw(heatmap)
dev.off()

#Análisis de Componentes Principales
pca<-prcomp(t(assay(se)),scale.=TRUE)
varianza<-summary(pca)$importance[2,1:2]*100
dfpca<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],Muscle.loss=colData(se)$Muscle.loss,Sample=colnames(assay(se)))
pcaplot<-ggplot(pca_df,aes(x=PC1,y=PC2,color=Muscle.loss,label=Sample))+
  geom_point(size=3,shape=1)+
  geom_text(aes(label=Sample),size=2,vjust=-1)+
  scale_color_manual(values=c("cachexic"="red","control"="blue"))+
  labs(title="PCA de las muestras",x=paste("PC1(",round(varianza[1],2),"%)",sep=""),y=paste("PC2 (",round(varianza[2],2),"%)",sep=""))+
  theme_minimal()+theme(legend.title = element_blank())
ggsave("pcaplot.png",plot=pcaplot,width=12,height=10,dpi=300)

#5 metabolitos con más contribución a PC1
pc1carga<- pca$rotation[, 1]
pc1ordenada<-sort(pc1carga,decreasing=TRUE)
top5<-names(pc1ordenada)[1:5]

#Cluster dendrogram
x<-assay(se)
clust.euclid.average<-hclust(dist(t(x)),method="average")
png("dendrograma_muestras.png",width=1000,height=800)
plot(clust.euclid.average, hang=-1,+
       main="Dendrograma de las muestras",+
       xlab="Muestras",+
       ylab="Distancia Euclidiana")
dev.off()

#Comparación de la expresión entre grupos (test de Wilcoxon)
control<-which(colData(se)$Muscle.loss=="control")
cachexic<-which(colData(se)$Muscle.loss=="cachexic")
wilcoxon<-function(x){
  tt=wilcox.test(x[control],x[cachexic])
  return(c(tt$statistic,tt$p.value,mean(x[control])-mean(x[cachexic])))}
ans<-apply(assay(se),1,ttest_wilcoxon)
tv<-ans[1,] 
pvals<-ans[2,]  
fc<-ans[3,]  
par(mar=c(8,4,4,2)) 
barplot(tv,names.arg=metabolitos,las=2,main="Estadístico de Wilcoxon por metabolito entre control y cachexia",cex.names=0.6)
