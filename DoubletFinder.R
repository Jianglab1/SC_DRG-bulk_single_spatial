library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
library(readr)

#Data import
#SC_sham
data_dir <- "./SC_sham/"  
list.files(data_dir)  
SC_sham <- Read10X(data.dir = data_dir) 

#SC_SNI
data_dir <- "./SC_SNI/"
list.files(data_dir)   
SC_SNI <- Read10X(data.dir = data_dir) 

#DRG_sham
data_dir <- "./DRG_sham/"  
list.files(data_dir)  
DRG_sham <- Read10X(data.dir = data_dir) 

#DRG_SNI
data_dir <- "./DRG_SNI/"  
list.files(data_dir)  
DRG_SNI <- Read10X(data.dir = data_dir)  

#CreateSeuratObject
seurat_object <- CreateSeuratObject(counts = SC_sham, project = "SC_sham", min.cells = 3, min.features = 200)
seurat_object <- CreateSeuratObject(counts = SC_SNI, project = "SC_SNI", min.cells = 3, min.features = 200)
seurat_object <- CreateSeuratObject(counts = DRG_sham, project = "DRG_sham", min.cells = 3, min.features = 200)
seurat_object <- CreateSeuratObject(counts = DRG_SNI, project = "DRG_SNI", min.cells = 3, min.features = 200)

#Follow the steps below for each sample
#Pre-process Seurat object (sctransform) 
seurat_object <- SCTransform(seurat_object)
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object,ndims = 50)
pc.num=1:20
seurat_object <- RunTSNE(seurat_object, dims =pc.num)
seurat_object <- FindNeighbors(seurat_object,dims=pc.num)%>%FindClusters(resolution=1.5)

# find pK
sweep.res.list <- paramSweep_v3(seurat_object, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Exclude same-origin doublets that cannot be detected, and optimize the expected number of doublets.
DoubletRate = ncol(seurat_object)*7.6*1e-6
homotypic.prop <- modelHomotypic(seurat_object$seurat_clusters)     
nExp_poi <- round(DoubletRate*ncol(seurat_object)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seurat_object <- doubletFinder_v3(seurat_object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)

#Result display
DoubletFinder <-DimPlot(seurat_object,reduction="tsne",group.by = paste0("DF.classifications_0.25_",pK_bcmvn,"_",nExp_poi.adj))
ggsave("DoubletFinder.pdf",DoubletFinder,width = 20, height = 15, units = "cm")

#Multi-package data storage and export
doublets <- as.data.frame(
  cbind(
    colnames(seurat_object), 
    seurat_object@meta.data[,grepl("pANN_0.25_",colnames(seurat_object@meta.data))], 
    seurat_object@meta.data[,grepl("DF.classifications_0.25_",colnames(seurat_object@meta.data))]))
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)
out <- getwd()
write_delim(doublets, file = paste0(out,"/DoubletFinder_doublets_singlets.tsv"), delim = "\t")

# Calculate number of doublets and singlets
summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(out,"/DoubletFinder_doublet_summary.tsv"), "\t")






