## CNV score calculation using tool "GISTIC2" ##
#########################make GISTIC input files##############################
### segment file data download and procession ###
setwd("F:/pan_cancer/pancan/GISTIC2/cancerfile/")
library(SummarizedExperiment)
library(TCGAbiolinks)
TCGAbiolinks:::getGDCprojects()$project_id
## download data
name = read.csv("cancer_name.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE )
Name = toupper(name[,1])

result = c()
for(i in 1:33){
    query <- GDCquery(project = paste0("TCGA-",Name[i]),
                    data.category = "Copy Number Variation",
                    data.type = "Masked Copy Number Segment")
    GDCdownload(query,method = "api")
    CNV_download <- GDCprepare(query = query, save = TRUE, save.filename = paste0(Name[i],"_CNV_download.rda"))

    #data procession
    load(paste0(Name[i],"_CNV_download.rda"))
    tumorCNV <- data
    tumorCNV = tumorCNV[,2:7]
    tumorCNV = tumorCNV[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
    #View(head(tumorCNV))

    ## tumor sample
    tum.sam = substr(tumorCNV$Sample,14,16) == "01A"
    #table(tum.sam)
    tumorCNV = tumorCNV[tum.sam, ]
    result = rbind(result, tumorCNV)
    write.table(tumorCNV,file = paste0(Name[i],'_CNV_segemnt_file.txt'),sep = '\t',quote = F,row.names = F)
    write.table(result,file = "Total_CNV_segemnt_file.txt",sep = '\t',quote = F,row.names = F)
}
dim(result)

### marker file data download and procession ###
setwd("F:/pan_cancer/pancan/GISTIC2/GISTIC/")
library(data.table)
Marker = data.table::fread("snp6.na35.remap.hg38.subset.txt", data.table=F)
str(Marker)
fre = Marker$freqcnv == "FALSE"
table(fre)

Marker = Marker[fre, 1:3]
colnames(Marker) = c("Marker_Name","Chromsome","Position")
head(Marker)
#   Marker_Name Chromsome  Position
# 2 SNP_A-1780271        15  33103578
# 4 SNP_A-1780274        20  35320106
# 5 SNP_A-1780277        12  75270266
# 6 SNP_A-1780278         1 218717316
# 7 SNP_A-1780283         4 126709121
# 8 SNP_A-1780285         6  90209746
dim(Marker)
write.table(Marker,file = 'Marker_file.txt',sep = '\t',quote = F,row.names = F)

#########################segment data grouping############################
setwd("D:/project/pan_cancer/pancan/GISTIC2/input_file/")
Segment = read.table("Total_CNV_segemnt_file.txt", sep = '\t', header=TRUE)
dim(Segment)
Segment$sample = substr(Segment[,1],1,15)
View(head(Segment))
dim(Segment)
length(unique(Segment[,1]))

group <- read.csv("D:/project/pan_cancer/pancan/pancan_6_cluster4_st.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
cordata = merge(group, Segment, by.x="sample",by.y="sample")
dim(cordata)
# [1] 2616692       8
View(head(cordata))
length(unique(cordata[,1]))
# [1] 9500
length(unique(cordata[,3]))
# [1] 9675
# names(cordata)
# [1] "sample"       "group"        "Sample"(original sample name)       "Chromosome"  
# [5] "Start"        "End"          "Num_Probes"   "Segment_Mean"

### grouping ###
st1 = subset(cordata[,3:8], cordata[,2]=="St1")
View(head(st1))
dim(st1)
# [1] 1097668      6
length(unique(st1[,1]))

st2 = subset(cordata[,3:8], cordata[,2]=="St2")
View(head(st2))
dim(st2)
# [1] 599860      6
length(unique(st2[,1]))
# [1] 1736

st3 = subset(cordata[,3:8], cordata[,2]=="St3")
View(head(st3))
dim(st3)
# [1] 744811      6
length(unique(st3[,1]))
# [1] 4100

st4 = subset(cordata[,3:8], cordata[,2]=="St4")
View(head(st4))
dim(st4)
# [1] 174353      6
length(unique(st4[,1]))
# [1] 644

write.table(st1,file = "St1_segment_file.txt",sep = '\t',row.names = F, col.names = F) #,quote = F
write.table(st2,file = "St2_segment_file.txt",sep = '\t',row.names = F, col.names = F)
write.table(st3,file = "St3_segment_file.txt",sep = '\t',row.names = F, col.names = F)
write.table(st4,file = "St4_segment_file.txt",sep = '\t',row.names = F, col.names = F)