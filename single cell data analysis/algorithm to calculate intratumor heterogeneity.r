## Developing an algorithm to calculate intratumor heterogeneity ##
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
setwd("D:/project/SC_ITH/")
##########################extract expression matrix of each patient#######################
## extract malignant cells
exp_tumor <- read.csv("SARC_tumor_scrna_exp.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE,row.names = 1)
View(head(exp_tumor))
dim(exp_tumor)

sample = read.csv("SARC_GSE131309_annotations.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)[,1:6]
dim(sample)
View(head(sample))

sample_tumor = subset(sample, sample[,6] == "Malignant")
dim(sample_tumor)

## gene expression info of each patient
patient = unique(sample_tumor[,5]) 
table(sample_tumor[,5])

## define normalization of gene exp matrix in each cluster
normfun <- function(data){
    xmax = max(data)
    xmin = min(data)
    y = (data-xmin)/(xmax-xmin)
    return(y)
}

ITH_all <- c()
cluster_num_all <- c()
for(p in c(1:4,6:length(patient))){
    geneSet <- subset(sample_tumor, sample_tumor[,5] == patient[p])
    exp_geneSet <- exp_tumor[ , colnames(exp_tumor) %in% geneSet[, 2]]
    write.csv(exp_geneSet, paste0("exp_", patient[p], ".csv"))

    #############cluster##############
    dir.create(file.path(getwd(), patient[p]))
    geneSet1 = read.table(paste0("exp_", patient[p], ".csv"), header=TRUE, row.names = 1, stringsAsFactors=FALSE, sep=",")
    geneSet1 <- CreateSeuratObject(counts = geneSet1, min.cells = 3, min.features = 200)
    all.genes <- rownames(geneSet1)

    ## Scaling the data
    geneSet1 <- ScaleData(geneSet1, features = all.genes) #z-score

    ## Perform linear dimensional reduction
    geneSet1 <- FindVariableFeatures(object = geneSet1) #default nfeatures=2000
    geneSet1 <- RunPCA(geneSet1, features = VariableFeatures(object = geneSet1))
    plot_pca = ElbowPlot(geneSet1)
    ggsave(paste0(patient[p],"/", patient[p],"_ElbowPlot.png"), plot = plot_pca, width = 8, height = 4)

    geneSet1 <- FindNeighbors(geneSet1, dims = 1:10)
    geneSet1 <- FindClusters(geneSet1, resolution = 0.5) 
    geneSet1 <- RunUMAP(geneSet1, dims = 1:10)
    plot_umap = DimPlot(geneSet1, reduction = "umap")
    ggsave(paste0(patient[p],"/", patient[p],"_umap.png"), plot = plot_umap, width = 8, height = 4)
    cluster_num = length(levels(geneSet1$RNA_snn_res.0.5))

    ## gene expression of each cluster
    for(num in 1:cluster_num){
    geneSet1_clu <- subset(x = geneSet1, RNA_snn_res.0.5 == levels(geneSet1$RNA_snn_res.0.5)[num])
    exp_subcluster = as.matrix(geneSet1_clu@assays$RNA@counts)
    write.csv(exp_subcluster, paste0(patient[p],"/exp_cluster",num,".csv"))
    }

    #########normalization and algorithm######
    geneset_mean = c()
    for(i in 1:cluster_num){
        geneSet_cluster <- read.table(paste0(patient[p],"/exp_cluster", i, ".csv"), stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, row.names = 1, sep = ",")
        geneSet_cluster <- as.matrix(geneSet_cluster)
        ## normalize gene exp in each cell
        geneSet_cluster.normalize <- apply(geneSet_cluster,2,function(x) normfun(x))
        ## calculate the mean value of each gene
        geneSet_cluster.mean <- apply(geneSet_cluster.normalize, 1, mean) %>% as.data.frame() 
        colnames(geneSet_cluster.mean) <- paste0("cluster", i, "_gene.mean")
        geneSet_cluster.mean <- as.matrix(geneSet_cluster.mean)
        geneset_mean = cbind(geneset_mean, geneSet_cluster.mean)
    }
    geneset_mean <- t(geneset_mean)

    ## calculate the ITH score
    euc <- dist(geneset_mean, method="euclidean")
    median.result <- median(euc)
    ITH_score <- round(median.result * cluster_num, 2)
    ITH_score <- as.data.frame(ITH_score)
    rownames(ITH_score) <- patient[p]
    ITH_all <- rbind(ITH_all, ITH_score)
    cluster_num <- as.data.frame(cluster_num)
    rownames(cluster_num) <- patient[p]
    cluster_num_all <- rbind(cluster_num_all, cluster_num)
}
colnames(ITH_all) <- "ITH_score"
colnames(cluster_num_all) <- "cluster_number"
result <- cbind(cluster_num_all, ITH_all)
write.csv(result, "SARC_ITH_Cluster.csv")
