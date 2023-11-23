####### Human and mouse gene symbol transformation ########
install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
genes = c("Tmx2", "Trp53", "Zfp286")
transMG <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    gs = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                values = x, mart = mouse,
                attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)
return(gs)
}
transMG(genes)