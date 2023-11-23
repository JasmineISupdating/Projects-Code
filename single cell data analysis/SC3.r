## 1.3 Clustering analysis using R package "SC3" ##
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("SC3")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)

rm(list=ls())
rna = c(paste("scrna_seq", 1:19, sep=""))
k <- c()
for(i in 1:19){
    pbmc <- read.table(paste0("F:/LUAD/SC3/group_info/", rna[i], "_cluster/SCE_L_exp.csv"),stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, row.names=1, sep=",")
    ann <- read.table(paste0("F:/LUAD/SC3/ann/", rna[i], "_SCE_L.csv"), header=T, row.names=1, sep=",")
    sce <- SingleCellExperiment(
        assays = list(
            counts = as.matrix(pbmc),
            logcounts = log2(as.matrix(pbmc) + 1)
        ),
        colData=ann
    )
    ## define feature names in feature_symbol column
    rowData(sce)$feature_symbol <- rownames(sce)
    # remove features with duplicated names
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

    ## cluster
    sce <- sc3_estimate_k(sce)
    ks <- metadata(sce)$sc3$k_estimation
    k = rbind(k,ks)
}
rownames(k) = c(paste("scrna_seq",1:19,sep=""))
write.csv(k,"F:/LUAD/SC3/result/k_SCE_L.csv")