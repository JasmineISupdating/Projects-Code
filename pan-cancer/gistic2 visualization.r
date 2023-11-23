################visualization#################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)
library(tidyverse)
options(stringsAsFactors = F)

st = c("st1","st2","st3","st4")
for(i in 1:4){
    setwd(paste0("D:/project/pan_cancer/pancan/GISTIC2/",st[i],"_results/"))
    ## read GISTIC files
    cancer.gistic = readGistic(gisticAllLesionsFile = "all_lesions.conf_90.txt",
                            gisticAmpGenesFile = "amp_genes.conf_90.txt", 
                            gisticDelGenesFile = "del_genes.conf_90.txt", 
                            gisticScoresFile = "scores.gistic", 
                            isTCGA = T)
    cancer.gistic
    getSampleSummary(cancer.gistic)
    getGeneSummary(cancer.gistic)
    getCytobandSummary(cancer.gistic)
    write.GisticSummary(gistic=cancer.gistic, basename="cancer_gistic2")

    ### visualize the results ###
    ## genome plot
    a <- getCytobandSummary(cancer.gistic)$Cytoband
    pdf(file= paste0("D:/project/pan_cancer/pancan/GISTIC2/visualize/g_score_", st[i], ".pdf"),width=7,height=4)
    gisticChromPlot(gistic = cancer.gistic, ref.build = "hg38", fdrCutOff = 0.05, markBands = a[1:3], y_lims = c(-1,1))
    dev.off()
}
