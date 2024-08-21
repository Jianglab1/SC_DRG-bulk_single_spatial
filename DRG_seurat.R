library(Seurat)
library(Rcpp)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(tidyr)

#DoubletFinder
Doublet_sham <- read.table(file.choose(),header = T,row.names = 1,sep = "\t")
rownames(Doublet_sham) <- paste('sham',rownames(Doublet_sham),sep='_')
Doublet_SNI <- read.table(file.choose(),header = T,row.names = 1,sep = "\t")
rownames(Doublet_SNI) <- paste('SNI',rownames(Doublet_SNI),sep='_')
Doublet_sham_SNI <- rbind(Doublet_sham,Doublet_SNI)

#Data import
#DRG_sham
data_dir <- "./DRG_sham/" 
list.files(data_dir)  
sham <- Read10X(data.dir = data_dir) 
colnames(sham) <- paste('sham',colnames(sham),sep='_')

#DRG_SNI
data_dir <- "./DRG_SNI/"
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

#Filter
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
tsneplot <- DimPlot(seurat_object, reduction = "tsne", group.by = "orig.ident", pt.size = 1.5,label = TRUE)
ggsave("tsne_plot.pdf", width = 35, height = 20, units = "cm")

#save and Wait for the subgroup to be named
save(seurat_object,file = "waite_ids.Rda")

#Prepare the form for loupe naming
data.tsne <- seurat_object@reductions[["tsne"]]@cell.embeddings
data.cluster <- seurat_object[["SCT_snn_res.1.5"]]
data.tsne <- data.tsne %>% as.data.frame() %>% mutate(ID = rownames(data.tsne))
data.cluster <- data.cluster %>% as.data.frame() %>% mutate(ID = rownames(data.cluster))
data.tsne <- data.tsne %>% separate(ID,into = c("type","barcode"),sep = "_",remove = F)
data.cluster <- data.cluster %>% separate(ID,into = c("type","barcode"),sep = "_",remove = F)
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
new.cluster.ids_num <- c("0_Satglia glia", "1_Satglia glia", "2_Satglia glia", "3_Fibroblast_Dcn", "4_Schwann_N","5_Fibroblast_Dcn", "6_Vascular", "7_Satglia glia", "8_Vascular","9_Schwann_M","10_NP","11_Schwann_M","12_NF2",'13_Fibroblast_Dcn','14_Satglia glia','15_Satglia glia','16_NF1','17_NP','18_PEP1','19_PEP1',
                         '20_NF3','21_NF2','22_PNI_N1','23_PNI_N2','24_Macrophage','25_PEP1','26_cLTMR','27_Vascular','28_PEP2','29_Schwann_M','30_Satglia glia','31_cLTMR','32_Fibroblast_Dcn','33_NP','34_T cell','35_NF3','36_Vascular','37_SST','38_Vascular') 
new.cluster.ids <- c("Satglia glia", "Satglia glia", "Satglia glia", "Fibroblast_Dcn", "Schwann_N","Fibroblast_Dcn", "Vascular", "Satglia glia", "Vascular","Schwann_M","NP","Schwann_M","NF2",'Fibroblast_Dcn','Satglia glia','Satglia glia','NF1','NP','PEP1','PEP1',
                     'NF3','NF2','PNI_N1','PNI_N2','Macrophage','PEP1','cLTMR','Vascular','PEP2','Schwann_M','Satglia glia','cLTMR','Fibroblast_Dcn','NP','T cell','NF3','Vascular','SST','Vascular') 
levels(seurat_object@meta.data[["SCT_snn_res.1.5"]]) <- new.cluster.ids_num
levels(seurat_object@meta.data[["seurat_clusters"]]) <- new.cluster.ids

#tsne plot
seurat_object <- SetIdent(seurat_object, value = seurat_object@meta.data$SCT_snn_res.1.5)
tsneplot<-TSNEPlot(seurat_object,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot
ggsave("tsne_sham_SNI_num.pdf",tsneplot2,width = 25, height = 15, units = "cm")

seurat_object <- SetIdent(seurat_object, value = seurat_object@meta.data$seurat_clusters)
tsneplot<-TSNEPlot(seurat_object,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot
ggsave("tsne_sham_SNI.pdf",tsneplot,width = 25, height = 15, units = "cm")

#subset
sham <- subset(seurat_object,subset=orig.ident=='sham')
SNI <- subset(seurat_object,subset=orig.ident=='SNI')

#sham
sham <- SetIdent(sham, value = sham@meta.data$SCT_snn_res.1.5)
tsneplot2<-TSNEPlot(sham,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot2
ggsave("tsne_sham_num.pdf",tsneplot2,width = 25, height = 15, units = "cm")

sham <- SetIdent(sham, value = sham@meta.data$seurat_clusters)
tsneplot2<-TSNEPlot(sham,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot2
ggsave("tsne_sham.pdf",tsneplot2,width = 25, height = 15, units = "cm")

#SNI
SNI <- SetIdent(SNI, value = SNI@meta.data$SCT_snn_res.1.5)
tsneplot2<-TSNEPlot(SNI,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot2
ggsave("tsne_SNI_num.pdf",tsneplot2,width = 25, height = 15, units = "cm")

SNI <- SetIdent(SNI, value = SNI@meta.data$seurat_clusters)
tsneplot2<-TSNEPlot(SNI,label = TRUE, pt.size = 0.5,repel=T) 
tsneplot2
ggsave("tsne_SNI.pdf",tsneplot2,width = 25, height = 15, units = "cm")

#Find marker genes
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object,vars.to.regress = "percent.mt") 
dif<-FindAllMarkers(seurat_object,only.pos = T,assay = "RNA")
write.table(dif,"dif_num.xls",row.names = T,col.names = T,quote = F,sep = "\t")
sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 3,wt = avg_log2FC)
#save
write.table(sig.dif,"sig.dif.xls",row.names = T,col.names = T,quote = F,sep = "\t")
save(seurat_object,file = "DRG_seurat.Rda")


