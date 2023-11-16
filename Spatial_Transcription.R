library(Seurat)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(clustree)


#read
SC <- Load10X_Spatial(data.dir = "./sc",slice = "sc")
dim(SC)

#SCT
SC <- SCTransform(SC,assay = "Spatial")  
SC <- RunPCA(SC)
ElbowPlot(SC,ndims = 50) 
pc.num=1:20
SC <- FindNeighbors(SC,dims = pc.num)
for(i in seq(0.4,2,by=0.1)){
  SC <- FindClusters(SC,resolution =i)
}
for (i in seq(0.4,2,by=0.1)) {
  SC <- SetIdent(SC,value = SC@meta.data[[paste0('SCT_snn_res.',i)]])
  pdf(file=paste0('Spatial_',i,'.pdf'))
  print(SpatialDimPlot(SC))
  dev.off()
}
SC <- SetIdent(SC,value = SC@meta.data[["SCT_snn_res.0.6"]])
SC <- RunTSNE(SC,dims = pc.num,label=T)
SC <- RunUMAP(SC,dims = pc.num,label=T)

## Export information for loupe
data.tsne <- SC@reductions[["tsne"]]@cell.embeddings
write.csv(data.tsne,"data.tsne.csv",row.names = T,quote = T)
write.csv(SC[['SCT_snn_res.0.6']],'SCT_snn_res.0.6.csv',row.names = T,quote = T)

#load name from loupe
SC_0.6 <- read.table('SCT_snn_res.0.6.csv',header = T,row.names = 1,sep = ",")
colnames(SC_0.6) <- 'celltype'
celltype <-SC_0.6

#to Metadata
SC <- AddMetaData(SC,celltype,col.name = "celltype")
SC@meta.data[["celltype"]] <- as.factor(SC@meta.data[["celltype"]])
SC <- SetIdent(SC,value=SC@meta.data$celltype)

#Find marker genes
#Normalization
SC@active.assay <- 'Spatial'
SC <- NormalizeData(SC, normalization.method = "LogNormalize", scale.factor = 10000)
SC <- FindVariableFeatures(SC, selection.method = "vst", nfeatures = 2000)
SC <- ScaleData(SC,vars.to.regress = "percent.mt") 
dif<-FindAllMarkers(SC,only.pos = T,assay = "RNA")
write.table(dif,"SC_dif.xls",row.names = T,col.names = T,quote = F,sep = "\t")
sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 5,wt = avg_log2FC)
write.table(sig.dif,"SC_sig.dif.xls",row.names = T,col.names = T,quote = F,sep = "\t")
#save
save(SC,file = "SC.Rda")








