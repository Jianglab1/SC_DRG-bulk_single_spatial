library(tidyverse)
library(ggtree)
library(ape)

load("./SC_neuron_seurat.Rda")
seurat_object@active.assay <- 'RNA'
# Calculate the mean expression per gene per cluster
mean_expression <- as.data.frame(AverageExpression(object = seurat_object, group.by = "SCT_snn_res.0.8",assays = 'RNA'))
colnames( mean_expression) <- gsub('RNA.','',colnames( mean_expression))
colnames( mean_expression) <- gsub('Inhibitory.Neuron','I',colnames( mean_expression))
colnames( mean_expression) <- gsub('Excitatory.Neuron','E',colnames( mean_expression))
colnames( mean_expression) <- gsub('Inhibitory.Excitatory.neuron','I/E',colnames( mean_expression))
colnames( mean_expression) <- gsub('.I',' I',colnames( mean_expression))
colnames( mean_expression) <- gsub('.E',' E',colnames( mean_expression))
colnames( mean_expression) <- gsub('I E','I/E',colnames( mean_expression))


mean_expression_filtered <- mean_expression[which(rowMeans(mean_expression) > 0), ]
distance_metrics <- dist(t(mean_expression_filtered), method = "euclidean")
hc <- hclust(distance_metrics, method = "complete")
tree =as.phylo(hc)

cols<- read.csv('cols.csv',header = T,row.names = 1)
rownames(cols) <- gsub('Inhibitory.Neuron','I',rownames(cols))
rownames(cols) <- gsub('Excitatory.Neuron','E',rownames(cols))
rownames(cols) <- gsub('Inhibitory Excitatory.neuron','I/E',rownames(cols))
tip_col=c('#E41A1C',
          '#8C4C6A',
          '#F58499',
          '#4384BE',
          '#1F5575',
          '#1C8D76',
          '#2AA478',
          '#76C582',
          '#B7DB7A',
          '#C3D257',
          '#606592',
          '#3C8A9B',
          '#47A363',
          '#4EAD4B',
          '#619462',
          '#747C78',
          '#9A509E',
          '#C8734D',
          '#DE8524',
          '#F38E3A',
          '#E589C2',
          '#D690C6',
          '#C797C9',
          '#BA98C2',
          '#B58698',
          '#AA6245',
          '#A05C35',
          '#8B7364',
          '#768994',
          '#4FA5DB',
          '#2A4386')
cols=c('black',
       'black',
       'black',
       '#8C4C6A',
       '#B58698',
       'black',
       'black',
       '#B7DB7A',
       'black',
       '#619462',
       '#F38E3A',
       'black',
       '#DE8524',
       'black',
       '#C797C9',
       'black',
       '#4EAD4B',
       '#C8734D',
       'black',
       'black',
       '#2AA478',
       'black',
       '#8B7364',
       '#768994',
       'black',
       'black',
       '#4384BE',
       '#D690C6',
       'black',
       'black',
       '#C3D257',
       '#4FA5DB',
       'black',
       'black',
       '#E589C2',
       'black',
       '#9A509E',
       '#AA6245',
       'black',
       'black',
       '#747C78',
       '#A05C35',
       'black',
       '#76C582',
       'black',
       'black',
       '#8C4C6A',
       'black',
       '#E41A1C',
       '#1F5575',
       'black',
       '#3C8A9B',
       'black',
       '#606592',
       '#47A363',
       'black',
       '#BA98C2',
       'black',
       '#1C8D76',
       '#2A4386'
)

plot(tree, type="unrooted",label.offset=0.01,
     rotate.tree=180, no.margin=FALSE, cex=1, edge.width=4,edge.color=cols)       

pdf('tree.pdf',width = 8,height = 8)
plot(tree, type="unrooted",label.offset=20,
     rotate.tree=180, no.margin=FALSE, cex=1, edge.width=4,edge.color=cols,tip.color=tip_col)
dev.off()