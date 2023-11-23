##################single cell data integration############
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
files<-list.files(path = "D:/pan_single/exp/",
           pattern = "*.csv",
           full.names = T)            
scRNAlist <- list()
cancer <- c("BRCA", "GBM", "LIHC", "LUAD", "OV", "PRAD", "RCC", "SARC", "SKCM")
for(i in 1:length(files)){
    counts <- read.table(files[i], stringsAsFactors = FALSE, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
    scRNAlist[[i]] <- CreateSeuratObject(counts, project = cancer[i], min.cells = 3, min.genes = 0)
    scRNAlist[[i]]@meta.data$orig.ident = cancer[i]
}
save(scRNAlist, file = "D:/pan_single/scRNAlist_new.Rdata")

######################### integration ########################
setwd("D:/Lei/pan_single/")
load("scRNAlist_new.Rdata")
scRNAlist
memory.limit(100000000)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst", nfeatures = 11000)#nfeatures = 2000
}
features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 10000)
length(features)
# [1] 7210
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features)
save(scRNA.anchors, file = "D:/Lei/pan_single/scRNA_anchors.Rdata")

load("scRNA_anchors.Rdata")
memory.limit(1000000000)
scRNA4 <- IntegrateData(anchorset = scRNA.anchors)
dim(scRNA4)
# [1]  7210 25664
save(scRNA4, file = "D:/pan_single/scRNA_9cancer_integrate.Rdata")
# table(scRNA4$orig.ident) 
# scRNA4@assays$RNA@counts[1:4,1:4]
# scRNA4@assays$integrated[1:4,1:4]
# View(as.data.frame(scRNA4@assays$integrated[1:4,1:4]))

load("scRNA_9cancer_integrate.Rdata")
exp_inte = as.data.frame(scRNA4@assays$integrated[1:7210, 1:25664])
# > dim(exp_inte)
# [1]  7210 25664
write.csv(exp_inte, "scRNA_9cancer_integrate.csv")

######################### dimension reduction and cluster #########################
length(VariableFeatures(scRNA4))
# [1] 7210
memory.limit(1000000000)
scRNA5 <- ScaleData(scRNA4, features = VariableFeatures(scRNA4))
scRNA5 <- RunPCA(scRNA5, features = VariableFeatures(scRNA5))
plot1 <- DimPlot(scRNA5, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA5, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("integrate/pca.png", plot = plotc, width = 8, height = 4)

pc.num=1:30
scRNA5 <- FindNeighbors(scRNA5, dims = pc.num) 
scRNA5 <- FindClusters(scRNA5, resolution = 0.1)
table(scRNA5@meta.data$seurat_clusters)
metadata <- scRNA5@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'integrate_res.1/cell_cluster.csv',row.names = F)
table(scRNA5$orig.ident, scRNA5@meta.data$seurat_clusters) ####res.1

## tSNE
scRNA5 = RunTSNE(scRNA5, dims = pc.num)
embed_tsne <- Embeddings(scRNA5, 'tsne')
write.csv(embed_tsne,'integrate_res.1/embed_tsne.csv')
## group_by_cluster
plot1 = DimPlot(scRNA5, reduction = "tsne", label=T) 
ggsave("integrate_res.1/tSNE.png", plot = plot1, width = 8, height = 7)
## group_by_sample
plot2 = DimPlot(scRNA5, reduction = "tsne", group.by='orig.ident') 
ggsave("integrate_res.1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
## combination
plotc <- plot1+plot2
ggsave("integrate_res.1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

## UMAP
scRNA5 <- RunUMAP(scRNA5, dims = pc.num)
embed_umap <- Embeddings(scRNA5, 'umap')
write.csv(embed_umap,'integrate_res.1/embed_umap.csv') 
## group_by_cluster
plot3 = DimPlot(scRNA5, reduction = "umap", label=T) 
ggsave("integrate_res.1/UMAP.png", plot = plot3, width = 8, height = 7)
## group_by_sample
plot4 = DimPlot(scRNA5, reduction = "umap", group.by='orig.ident')
ggsave("integrate_res.1/UMAP_sample.png", plot = plot4, width = 8, height = 7)
## combination
plotc <- plot3+plot4
ggsave("integrate_res.1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("integrate_res.1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
