#######################Logistic#########################
library("QuantPsyc")
setwd("D:/project/pan_cancer/pancan/logistic/")
## immuneScore ##
imm = read.csv("estimate_combine_new.csv", header=T, stringsAsFactors=FALSE, check.names=FALSE)[,c(1,4)]
# colnames(imm)=c("sample","ImmuneScore")
# imm=imm[-1,]
imm_1 = imm[which(as.numeric(imm[,2])> quantile(as.numeric(imm[,2]), probs = seq(0, 1, 0.5))[[2]]),] #up
dim(imm_1)
imm_0 = imm[which(as.numeric(imm[,2]) <= quantile(as.numeric(imm[,2]), probs = seq(0, 1, 0.5))[[2]]),] #down
dim(imm_0)
imm_1[,2] = 1
imm_0[,2] = 0
imm_all = rbind(imm_1,imm_0)
dim(imm_all)

b = read.csv("pancan_HRD_combine.csv", stringsAsFactors=FALSE, check.names=FALSE) #HRD
c = read.csv("pancan_TMB_combine.csv", stringsAsFactors=FALSE, check.names=FALSE) #TMB
HSC = read.csv("st_ratio_4col.csv", stringsAsFactors=FALSE, check.names=FALSE)[,c(1,2)]
ESC = read.csv("st_ratio_4col.csv", stringsAsFactors=FALSE, check.names=FALSE)[,c(1,3)]
# imm[,1] = substr(imm[,1],1,15)
b = b[match(imm_all[,1],b$sample),]
c = c[match(imm_all[,1],c$sample),]
HSC = HSC[match(imm_all[,1], HSC$sample),]
ESC = ESC[match(imm_all[,1], ESC$sample),]

dat = as.data.frame(cbind(imm_all, b[,3], c[,3], HSC[,2], ESC[,2]), stringsAsFactors = FALSE)
dat[,2] = as.numeric(dat[,2])
dat[,3] = as.numeric(dat[,3])
dat[,4] = as.numeric(dat[,4])
dat[,5] = as.numeric(dat[,5])
dat[,6] = as.numeric(dat[,6])

colnames(dat)=c("sample","ImmuneScore","HRD","TMB", "HSC_score", "ESC_score")
dat$ImmuneScore[dat$ImmuneScore=="1"] = 1
dat$ImmuneScore[dat$ImmuneScore=="0"] = 0
View(head(dat))

fit.full = glm(ImmuneScore ~ HRD + TMB + HSC_score + ESC_score, family = binomial(), data = dat) #binomial()
standardized_regression_coeffients=c('',lm.beta(fit.full))
summary(fit.full)
standardized_regression_coeffients
res = summary(fit.full)$coefficients
res = cbind( res, standardized_regression_coeffients)
write.csv(res, "immuneScore_TMB_HRD.csv")