#install.packages('ComplexUpset')
library(ComplexUpset)
library(UpSetR)
library(paletteer)
library(CellChat)

#load the DEGs in bulk-seq and the marker in ST
dif <- read.csv('dif.csv',header = T)
up_in_DRG <- read.csv('bulk_DEG.csv',header = T)
up_in_DRG <- up_in_DRG$SYMBOL[up_in_DRG$sig=='up']
SC_R1 <- dif$gene[dif$cluster=='SC_R1']
SC_R2 <- dif$gene[dif$cluster=='SC_R2']
SC_R3 <- dif$gene[dif$cluster=='SC_R3']
SC_R4 <- dif$gene[dif$cluster=='SC_R4']
SC_R5 <- dif$gene[dif$cluster=='SC_R5']
SC_R6 <- dif$gene[dif$cluster=='SC_R6']
SC_R7 <- dif$gene[dif$cluster=='SC_R7']
SC_R8 <- dif$gene[dif$cluster=='SC_R8']

data_upset <-list(up_in_DRG = up_in_DRG,SC_R1 = SC_R1,SC_R2 = SC_R2,SC_R3 = SC_R3,SC_R4 = SC_R4,SC_R5 = SC_R5,SC_R6 = SC_R6,SC_R7 = SC_R7,SC_R8 = SC_R8)
pdf(file = 'upset1.pdf',width = 14.4,height =7.2)

upset(fromList(data_upset),sets=c('up_in_DRG','SC_R8','SC_R7','SC_R6','SC_R5','SC_R4','SC_R3','SC_R2','SC_R1'),
      nsets = 9,
      nintersects=NA,
      intersections=list(
        list('up_in_DRG','SC_R8'),
        list('up_in_DRG','SC_R7'),
        list('up_in_DRG','SC_R6'),
        list('up_in_DRG','SC_R5'),
        list('up_in_DRG','SC_R4'),
        list('up_in_DRG','SC_R3'),
        list('up_in_DRG','SC_R2'),
        list('up_in_DRG','SC_R1'),
        list('up_in_DRG','SC_R8','SC_R7'),
        list('up_in_DRG','SC_R3','SC_R7'),
        list('up_in_DRG','SC_R8','SC_R1'),
        list('up_in_DRG','SC_R1','SC_R6'),
        list('up_in_DRG','SC_R4','SC_R8'),
        list('up_in_DRG','SC_R5','SC_R8'),
        list('up_in_DRG','SC_R8','SC_R6'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R3'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R1'),
        list('up_in_DRG','SC_R7','SC_R2','SC_R3'),
        list('up_in_DRG','SC_R8','SC_R2','SC_R3'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R4'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R5'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R3','SC_R5'),
        list('up_in_DRG','SC_R8','SC_R7','SC_R4','SC_R1')),
      #show.numbers = "no",
      set_size.angles = 15,
      text.scale = c(1.5, 1.5, 1.5, 1.4, 1.5, 1.5),
      mb.ratio = c(0.55, 0.45),
      point.size = 3.5, line.size = 1,
      #order.by = "freq",
      empty.intersections = "on",
      matrix.color = "#000000",
      sets.bar.color = c("#F8766D","#F8766D","#000000","#000000", "#000000",  "#000000", "#000000","#000000", "#000000"),
      main.bar.color = "#000000",
      sets.x.label = "Gene Number",
      keep.order = TRUE,mainbar.y.label = "Intersection Gene Number",
      queries = list(list(query = intersects, params =list('up_in_DRG','SC_R8'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R1'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R4'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R5'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R7'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R6'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R7','SC_R1'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R7','SC_R3'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R7','SC_R4'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R2','SC_R3'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R8','SC_R7','SC_R5'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R1','SC_R4','SC_R7','SC_R8'),color = "#F8766D", active = T),
                     list(query = intersects, params =list('up_in_DRG','SC_R3','SC_R5','SC_R7','SC_R8'),color = "#F8766D", active = T)))
dev.off()

#
up_in_DRG <- read.csv('bulk_DEG.csv',header = T)
up_in_DRG <- up_in_DRG$SYMBOL[up_in_DRG$sig=='up']
inter_sect <- intersect(up_in_DRG,SC_R8)
up_in_DRG <- read.csv('bulk_DEG.csv',header = T)
intersect_R8 <- subset(up_in_DRG,subset = up_in_DRG$SYMBOL %in% inter_sect)
intersect_R8 <- intersect_R8[order(intersect_R8$logFC,decreasing = TRUE),]
write.csv(intersect_R8,"intersect_R8.csv")

