#加载包
library(Seurat)
library(Rcpp)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(clustree)

#DoubletFinder
Doublet_sham <- read.table(file.choose(),header = T,row.names = 1,sep = "\t")
rownames(Doublet_sham) <- paste('sham',rownames(Doublet_sham),sep='_')
Doublet_SNI <- read.table(file.choose(),header = T,row.names = 1,sep = "\t")
rownames(Doublet_SNI) <- paste('SNI',rownames(Doublet_SNI),sep='_')
#merge
Doublet_sham_SNI <- rbind(Doublet_sham,Doublet_SNI)

#Data import
#sham
data_dir <- "./SC_sham/"  
list.files(data_dir)  
sham <- Read10X(data.dir = data_dir) 
colnames(sham) <- paste('sham',colnames(sham),sep='_')

#SNI
data_dir <- "./SC_SNI/"
list.files(data_dir)
SNI <- Read10X(data.dir = data_dir) 
colnames(SNI) <- paste('SNI', colnames(SNI),sep = '_')

#merge
sham <- CreateSeuratObject(counts = sham, project = "sham", min.cells = 3, min.features = 200)
SNI <- CreateSeuratObject(counts = SNI, project = "SNI", min.cells = 3, min.features = 200)
seurat_object <- merge(sham,SNI)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("nFeature_nCount_percent.mt.pdf", width = 28, height = 25, units = "cm")
Doublet_sham_SNI <- select(Doublet_sham_SNI,-DoubletFinder_score)
seurat_object <- AddMetaData(seurat_object,Doublet_sham_SNI,col.name = "Doublet_sham_SNI")

#filter
seurat_object <- subset(seurat_object, subset = nFeature_RNA < 7500 & percent.mt < 5 & Doublet_sham_SNI=="singlet" )
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("plot1+plot2.pdf", width = 56, height = 25, units = "cm")

#Normalization
seurat_object <- SCTransform(seurat_object,vars.to.regress="percent.mt")
seurat_object <- RunPCA(seurat_object)
P1 <- DimPlot(seurat_object, reduction = "pca")
ggsave("PCA.pdf", width = 56, height = 25, units = "cm")
ElbowPlot(seurat_object,ndims = 50)
pc.num=1:20
seurat_object <- RunTSNE(seurat_object, dims =pc.num)
seurat_object <- RunUMAP(seurat_object, dims =pc.num)
seurat_object <- FindNeighbors(seurat_object,dims=pc.num)%>%FindClusters(resolution=1.5)

#tsne plot
tsneplot1 <- DimPlot(seurat_object, reduction = "tsne", group.by = "orig.ident", pt.size = 1.5,label = TRUE)
ggsave("tsne1.5.pdf", width = 35, height = 20, units = "cm")
#save and Wait for the subgroup to be named
save(seurat_object,file = "waite_ids.Rda")

#clustree
for(i in seq(0.4,3,by=0.1)){
  seurat_object <- FindClusters(seurat_object,resolution =i)
  yourfilename=paste("FC",i,".pdf",sep = "")
  p1 <- DimPlot(seurat_object,reduction="tsne",label = T,label.size = 3)
  pdf(file=yourfilename)
  print(p1)
  dev.off()
}

clustree(seurat_object)
ggsave("clustree.pdf",width = 35, height = 60, units = "cm")

#Prepare the form for loupe naming
data.tsne <- seurat_object@reductions[["tsne"]]@cell.embeddings
data.cluster <- seurat_object[["SCT_snn_res.1.5"]]
data.tsne <- data.tsne %>% as.data.frame() %>% mutate(ID = rownames(data.tsne))
data.cluster <- data.cluster %>% as.data.frame() %>% mutate(ID = rownames(data.cluster))
data.tsne_sham <- data.tsne %>% filter(type=="sham")%>%select(tSNE_1,tSNE_2,type,barcode)
data.tsne_SNI <- data.tsne %>% filter(type=="SNI")%>%select(tSNE_1,tSNE_2,type,barcode)
data.cluster_sham <- data.cluster %>% filter(type=="sham")%>%select(SCT_snn_res.1.5,type,barcode)
data.cluster_SNI <- data.cluster %>% filter(type=="SNI")%>%select(SCT_snn_res.1.5,type,barcode)
row.names(data.tsne_sham) <- data.tsne_sham$barcode
row.names(data.tsne_SNI) <- data.tsne_SNI$barcode
row.names(data.cluster_sham) <- data.cluster_sham$barcode
row.names(data.cluster_SNI) <- data.cluster_SNI$barcode
data.tsne_sham <- data.tsne_sham %>%select(tSNE_1,tSNE_2)
data.tsne_SNI <- data.tsne_SNI%>%select(tSNE_1,tSNE_2)
data.cluster_sham <- data.cluster_sham %>%select(SCT_snn_res.1.5)
data.cluster_SNI <- data.cluster_SNI %>%select(SCT_snn_res.1.5)
#save Forms
write.csv(data.tsne_sham,"data.tsne_sham.csv",row.names = T,quote = T)
write.csv(data.tsne_SNI,"data.tsne_SNI.csv",row.names = T,quote = T)
write.csv(data.cluster_sham,"data.cluster1.5_sham.csv",row.names = T,quote = T)
write.csv(data.cluster_SNI,"data.cluster1.5_SNI.csv",row.names = T,quote = T)
#Naming
new.cluster.ids <- c('Oligo.','Inhibitory/Excitatory neuron','Oligo.','Oligo.','Oligo.','Excitatory neuron','Oligo.','Oligo.','Excitatory neuron','Inhibitory/Excitatory neuron','Microglia',
                     'Inhibitory neuron','Astrocyte','Excitatory neuron','Inhibitory/Excitatory neuron','Microglia','Oligo.','Excitatory neuron','OPC','Oligo.','Inhibitory neuron',
                     'Inhibitory neuron','Excitatory neuron','Excitatory neuron','Vascular','Inhibitory neuron','Excitatory neuron','Inhibitory neuron','Menin.','Excitatory neuron','Excitatory neuron',
                     'Excitatory neuron','Excitatory neuron','Inhibitory neuron','Inhibitory neuron','Microglia') 
new.cluster.ids_num <- c('0 Oligo.','1 Inhibitory/Excitatory neuron','2 Oligo.','3 Oligo.','4 Oligo.','5 Excitatory neuron','6 Oligo.','7 Oligo.','8 Excitatory neuron','9 Inhibitory/Excitatory neuron','10 Microglia',
                         '11 Inhibitory neuron','12 Astrocyte','13 Excitatory neuron','14 Inhibitory/Excitatory neuron','15 Microglia','16 Oligo.','17 Excitatory neuron','18 OPC','19 Oligo.','20 Inhibitory neuron',
                         '21 Inhibitory neuron','22 Excitatory neuron','23 Excitatory neuron','24 Vascular','25 Inhibitory neuron','26 Excitatory neuron','27 Inhibitory neuron','28 Menin.','29 Excitatory neuron','30 Excitatory neuron',
                         '31 Excitatory neuron','32 Excitatory neuron','33 Inhibitory neuron','34 Inhibitory neuron','35 Microglia') 

levels(seurat_object@meta.data[["SCT_snn_res.1.5"]]) <- new.cluster.ids_num
levels(seurat_object@meta.data[["seurat_clusters"]]) <- new.cluster.ids
seurat_object@meta.data[["back_num"]]<- seurat_object@meta.data[["SCT_snn_res.1.5"]]
seurat_object@meta.data[["back"]]<- seurat_object@meta.data[["seurat_clusters"]]

save(seurat_object,file = "SC_seurat.Rda")

#Subset of neuron
neuron_seu <- subset(seurat_object,subset = back_num=='1 Inhibitory/Excitatory neuron'|back_num=='5 Excitatory neuron'|back_num=='8 Excitatory neuron'|back_num=='9 Inhibitory/Excitatory neuron'|back_num==
                       '11 Inhibitory neuron'|back_num=='13 Excitatory neuron'|back_num=='14 Inhibitory/Excitatory neuron'|back_num=='17 Excitatory neuron'|back_num=='20 Inhibitory neuron'|back_num==
                       '21 Inhibitory neuron'|back_num=='22 Excitatory neuron'|back_num=='23 Excitatory neuron'|back_num=='25 Inhibitory neuron'|back_num=='26 Excitatory neuron'|back_num=='27 Inhibitory neuron'|back_num=='29 Excitatory neuron'|back_num=='30 Excitatory neuron'|back_num==
                       '31 Excitatory neuron'|back_num=='32 Excitatory neuron'|back_num=='33 Inhibitory neuron'|back_num=='34 Inhibitory neuron')
n_neuron_seu <- subset(seurat_object,subset = !(back_num=='1 Inhibitory/Excitatory neuron'|back_num=='5 Excitatory neuron'|back_num=='8 Excitatory neuron'|back_num=='9 Inhibitory/Excitatory neuron'|back_num==
                                                  '11 Inhibitory neuron'|back_num=='13 Excitatory neuron'|back_num=='14 Inhibitory/Excitatory neuron'|back_num=='17 Excitatory neuron'|back_num=='20 Inhibitory neuron'|back_num==
                                                  '21 Inhibitory neuron'|back_num=='22 Excitatory neuron'|back_num=='23 Excitatory neuron'|back_num=='25 Inhibitory neuron'|back_num=='26 Excitatory neuron'|back_num=='27 Inhibitory neuron'|back_num=='29 Excitatory neuron'|back_num=='30 Excitatory neuron'|back_num==
                                                  '31 Excitatory neuron'|back_num=='32 Excitatory neuron'|back_num=='33 Inhibitory neuron'|back_num=='34 Inhibitory neuron'))


seurat_object <- neuron_seu
#Normalization
seurat_object <- SCTransform(seurat_object,vars.to.regress="percent.mt")#SCTransform is designed to eliminate the influence of sequencing depth. The nCount_RNA is used to construct the sct mode, so you don't need to put it into vars.to.regress.
seurat_object <- RunPCA(seurat_object)
DimPlot(seurat_object, reduction = "pca")
ElbowPlot(seurat_object,ndims = 50)
pc.num=1:20
seurat_object <- RunTSNE(seurat_object, dims =pc.num)
seurat_object <- FindNeighbors(seurat_object,dims=pc.num)%>%FindClusters(resolution=0.8)
seurat_object <- SetIdent(seurat_object, value = seurat_object@meta.data$SCT_snn_res.0.8)
#Prepare the form for loupe naming
data.tsne <- seurat_object@reductions[["tsne"]]@cell.embeddings
data.cluster <- seurat_object[["SCT_snn_res.0.8"]]
data.tsne <- data.tsne %>% as.data.frame() %>% mutate(ID = rownames(data.tsne))
data.cluster <- data.cluster %>% as.data.frame() %>% mutate(ID = rownames(data.cluster))
data.tsne <- data.tsne %>% separate(ID,into = c("type","barcode"),sep = "_",remove = F)
data.cluster <- data.cluster %>% separate(ID,into = c("type","barcode"),sep = "_",remove = F)
data.tsne_sham <- data.tsne %>% filter(type=="sham")%>%select(tSNE_1,tSNE_2,type,barcode)
data.tsne_SNI <- data.tsne %>% filter(type=="SNI")%>%select(tSNE_1,tSNE_2,type,barcode)
data.cluster_sham <- data.cluster %>% filter(type=="sham")%>%select(SCT_snn_res.0.8,type,barcode)
data.cluster_SNI <- data.cluster %>% filter(type=="SNI")%>%select(SCT_snn_res.0.8,type,barcode)
row.names(data.tsne_sham) <- data.tsne_sham$barcode
row.names(data.tsne_SNI) <- data.tsne_SNI$barcode
row.names(data.cluster_sham) <- data.cluster_sham$barcode
row.names(data.cluster_SNI) <- data.cluster_SNI$barcode
data.tsne_sham <- data.tsne_sham %>%select(tSNE_1,tSNE_2)
data.tsne_SNI <- data.tsne_SNI%>%select(tSNE_1,tSNE_2)
data.cluster_sham <- data.cluster_sham %>%select(SCT_snn_res.0.8)
data.cluster_SNI <- data.cluster_SNI %>%select(SCT_snn_res.0.8)

write.csv(data.tsne_sham,"data.tsne_sham.csv",row.names = T,quote = T)
write.csv(data.tsne_SNI,"data.tsne_SNI.csv",row.names = T,quote = T)
write.csv(data.cluster_sham,"data.cluster0.8_sham.csv",row.names = T,quote = T)
write.csv(data.cluster_SNI,"data.cluster0.8_SNI.csv",row.names = T,quote = T)

#clustree
for(i in seq(0.4,1.5,by=0.1)){
  seurat_object <- FindClusters(seurat_object,resolution =i)
  yourfilename=paste("FC",i,".pdf",sep = "")
  p1 <- DimPlot(seurat_object,reduction="tsne",label = T,label.size = 3)
  pdf(file=yourfilename)
  print(p1)
  dev.off()
}
library(clustree)
clustree(seurat_object)
ggsave("clustree.pdf",width = 35, height = 60, units = "cm")

#Naming
new.cluster.ids <- c('Inhibitory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Excitatory_Neuron','Inhibitory/Excitatory neuron','Excitatory_Neuron','Inhibitory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Inhibitory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Excitatory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Inhibitory_Neuron','Excitatory_Neuron') 
new.cluster.ids_num <- c('0_Inhibitory_Neuron','1_Excitatory_Neuron','2_Excitatory_Neuron','3_Excitatory_Neuron','4_Excitatory_Neuron','5_Excitatory_Neuron','6_Inhibitory_Neuron','7_Inhibitory_Neuron','8_Excitatory_Neuron','9_Inhibitory/Excitatory neuron','10_Excitatory_Neuron','11_Inhibitory_Neuron','12_Excitatory_Neuron','13_Excitatory_Neuron','14_Excitatory_Neuron','15_Inhibitory_Neuron','16_Inhibitory_Neuron','17_Excitatory_Neuron','18_Excitatory_Neuron','19_Excitatory_Neuron','20_Inhibitory_Neuron','21_Excitatory_Neuron','22_Excitatory_Neuron','23_Excitatory_Neuron','24_Excitatory_Neuron','25_Inhibitory_Neuron','26_Inhibitory_Neuron','27_Inhibitory_Neuron','28_Inhibitory_Neuron','29_Inhibitory_Neuron','30_Excitatory_Neuron') 
levels(seurat_object@meta.data[["SCT_snn_res.0.8"]]) <- new.cluster.ids_num
levels(seurat_object@meta.data[["seurat_clusters"]]) <- new.cluster.ids
levels(seurat_object@meta.data[["SCT_snn_res.0.8"]]) <- gsub('Inhibitory/Excitatory neuron','Inhibitory Excitatory neuron',levels(seurat_object@meta.data[["SCT_snn_res.0.8"]]))
levels(seurat_object@meta.data[["SCT_snn_res.0.8"]]) <- gsub('_',' ',levels(seurat_object@meta.data[["SCT_snn_res.0.8"]]))

#Find marker genes
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object,vars.to.regress = "percent.mt") 
seurat_object <- SetIdent(seurat_object, value = seurat_object@meta.data$SCT_snn_res.0.8)
dif<-FindAllMarkers(seurat_object,only.pos = T,assay = "RNA")
write.table(dif,"neu_dif_num.xls",row.names = T,col.names = T,quote = F,sep = "\t")
sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 3,wt = avg_log2FC)
write.table(sig.dif,"neu_sig.dif_num.xls",row.names = T,col.names = T,quote = F,sep = "\t")
#save
save(seurat_object,file = "SC_neuron_seurat.Rda")


#load merged name (neuron_seu,n_neuron_seu) from loupe
load("./SC_seurat.Rda")
back_sham <- read.csv('SCT_snn_res.0.8_sham_annotion_nnum.csv',header = T,row.names = 1)
rownames(back_sham) <- paste('sham',rownames(back_sham),sep='_')
back_SNI <- read.csv('SCT_snn_res.0.8_SNI_annotion_nnum.csv',header = T,row.names = 1)
rownames(back_SNI) <- paste('SNI',rownames(back_SNI),sep='_')
back_num_sham <- read.csv('SCT_snn_res.0.8_sham_annotion.csv',header = T,row.names = 1)
rownames(back_num_sham) <- paste('sham',rownames(back_num_sham),sep='_')
back_num_SNI <- read.csv('SCT_snn_res.0.8_SNI_annotion.csv',header = T,row.names = 1)
rownames(back_num_SNI) <- paste('SNI',rownames(back_num_SNI),sep='_')
back_sham_SNI <- rbind(back_sham,back_SNI)
colnames(back_sham_SNI) <- 'back'
back_num_sham_SNI <- rbind(back_num_sham,back_num_SNI)
colnames(back_num_sham_SNI) <- 'back_num'
n_neuron_clusters <- n_neuron_seu@meta.data[["back"]]%>%as.data.frame()
row.names(n_neuron_clusters) <- row.names(n_neuron_seu@meta.data)
colnames(n_neuron_clusters) <- 'back'
n_neuron_clusters_num <- n_neuron_seu@meta.data[["back_num"]]%>%as.data.frame()
row.names(n_neuron_clusters_num) <- row.names(n_neuron_seu@meta.data)
colnames(n_neuron_clusters_num) <- 'back_num'
back <- rbind(n_neuron_clusters,back_sham_SNI)
back_num <- rbind(n_neuron_clusters_num,back_num_sham_SNI)
seurat_object <- AddMetaData(seurat_object,back,col.name = "back")
seurat_object <- AddMetaData(seurat_object,back_num,col.name = "back_num")
levels(seurat_object@meta.data[["back"]])
levels(seurat_object@meta.data[["back_num"]])
save(seurat_object,file = "SC_seurat.Rda")

#tsne plot
cols=c('Astrocyte'="#1F77B4",'Excitatory neuron'="#FF7F0F",'Inhibitory neuron'="#4DAF4A",'Inhibitory/Excitatory neuron'="#E41A1C",'Menin.'="#9467BD",'Microglia'="#F781BF",'Oligo.'="#7F7F7F",'OPC'="#A65628",'Vascular'='#BCBD22')
tsneplot2<-TSNEPlot(seurat_object,pt.size = 0.5,cols=cols,label = FALSE) 
pdf('tsne.pdf',width=7.7,height = 5)
tsneplot2
dev.off()

seu_sham <- subset(seurat_object,subset=orig.ident=='sham')
seu_SNI <- subset(seurat_object,subset=orig.ident=='SNI')

pdf('tsne_sham.pdf',width=7.7,height = 5)
TSNEPlot(seu_sham,pt.size = 0.5,cols=cols,label = FALSE)
dev.off()
pdf('tsne_SNI.pdf',width=7.7,height = 5)
TSNEPlot(seu_SNI,pt.size = 0.5,cols=cols,label = FALSE)
dev.off()

#Find marker genes
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object,vars.to.regress = "percent.mt") 
dif<-FindAllMarkers(seurat_object,only.pos = T,assay = "RNA")
write.table(dif,"dif_num.xls",row.names = T,col.names = T,quote = F,sep = "\t")
sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 3,wt = avg_log2FC)
#save
write.table(sig.dif,"sig.dif.xls",row.names = T,col.names = T,quote = F,sep = "\t")
save(seurat_object,file = "SC_seurat.Rda")


