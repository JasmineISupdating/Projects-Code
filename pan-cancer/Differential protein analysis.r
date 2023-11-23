## Differential protein analysis
setwd("D:/project/pan_cancer/pancan/protein/CPTAC/DE_protein/")
rna_seq = read.csv("D:/project/pan_cancer/pancan/protein/CPTAC/proteme_data/normalize_TCGA_CPTAC_15datasets.csv",check.names = FALSE,stringsAsFactor=FALSE, header=TRUE, row.names = 1)
rna_seq = na.omit(rna_seq)
dim(rna_seq)
# [1] 2389 1679
cluster = read.csv("D:/project/pan_cancer/pancan/protein/CPTAC/cluster/protein_4cluster_st.csv",check.names = FALSE,stringsAsFactor=FALSE, header=TRUE)
## find up-regulated protein in each subtype
St = c("StC1","StC2","StC3","StC4") 
for(j in 1:3){
    for(k in (j+1):4){
        c1=subset(cluster[,1], cluster[,2]==St[j]) #H
        c2=subset(cluster[,1], cluster[,2]==St[k]) #L
        pos_H=match(c1,colnames(rna_seq))
        pos_L=match(c2,colnames(rna_seq))
        rna_seq = rna_seq[rowSums(rna_seq)!=0,]
        p.value = apply(rna_seq, 1, function(x){
                t.test(x[pos_H], x[pos_L], var.equal=TRUE, na.rm=TRUE)$p.value
            })
        log2FC = apply(rna_seq, 1, function(x){
            mean(x[pos_H], na.rm=TRUE) - mean(x[pos_L], na.rm=TRUE)
            })
        rna_seq = rna_seq[rowSums(rna_seq)!=0,]
        resC = cbind(rownames(rna_seq), p.value, log2FC)
        resC = resC[order(as.numeric(resC[,2])),]
        FDR = rep(1,dim(resC)[1]); for(i in 1:dim(resC)[1]) FDR[i] = as.numeric(resC[i,2]) * dim(resC)[1]/i
        resC = cbind(resC,FDR) #gene names,pv,log2FC,fdr
        sum_low=resC[which(as.numeric(resC[,3])<(-0.322) & as.numeric(resC[,4])<0.05),]
        write.csv(sum_low,paste0("fdr0.01/", St[j],"vs",St[k],"_low_score.csv"),row.names=FALSE)
            
        sum_high=resC[which(as.numeric(resC[,3])>0.322 & as.numeric(resC[,4])<0.05),]
        write.csv(sum_high,paste0("fdr0.01/", St[j],"vs",St[k],"_high_score.csv"),row.names=FALSE)
        }
}
