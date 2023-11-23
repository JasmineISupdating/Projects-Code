####################RUNinferCNV######################
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")
library(infercnv)

## Create an array matrix
count = read.table("expr.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, check.names=FALSE, quote="",fill=TRUE)
class(count)
mode(count)
dim(count)
annotation = read.table("groupFiles.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, check.names=FALSE)
gene = read.table("geneFile.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, check.names=FALSE)

a = apply(count, 2, as.numeric)
anno = t(annotation)
row_names = gene[,1]
col_names = anno[1,] 
dimnames(a) = list(row_names,col_names)
str(a)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = a,
                                    annotations_file = "/groupFiles.txt",
                                    delim = "\t",
                                    gene_order_file= "/geneFile.txt",
                                    ref_group_names = c("noncancer_cell"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             #out_dir = tempfile(), 
                             out_dir="/output_dir",
                             cluster_by_groups = TRUE,
                             #cluster_references = FALSE,
                             denoise = TRUE,
                             HMM = FALSE)

ref <- read.table("/output_dir/infercnv.references.txt",sep=" ", stringsAsFactors=FALSE, header=T, check.names=FALSE, row.names=1, quote="")
dim(ref)
obv <- read.table("/output_dir/infercnv.observations.txt",sep=" ", stringsAsFactors=FALSE, header=T, check.names=FALSE, row.names=1, quote="")
dim(obv)
