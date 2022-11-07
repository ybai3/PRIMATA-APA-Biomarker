library(stringr)
library(dplyr)
library(hash)
library(survival)
library(survminer)
library(glmnet)
library(lmtest)
library(gridExtra)
library(gtsummary)

setwd("/ihome/hpark/yub20/SUR3UTR/pythons/TCGA_FPKM/JNCI_TumorNormal/result_JNCI")
mRNA_dir <- "/ihome/hpark/yub20/SUR3UTR/data/TCGA_harmonized/FPKM"
surv_dir <- "/ihome/hpark/yub20/SUR3UTR/new/survival"
#APA4_dir <- "/ihome/hpark/yub20/SUR3UTR/new/data"
miR_dir <- "/ihome/hpark/yub20/SUR3UTR/data"

miR_Family_Info <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/refs/miR_Family_Info.txt")
miR_Family_Info <- miR_Family_Info[str_detect(miR_Family_Info$MiRBase.ID, "hsa-"),]
miR_Family_Info$match4marker <- str_remove(miR_Family_Info$MiRBase.ID, "hsa-")
miR_hash <- hash(keys= miR_Family_Info$MiRBase.Accession, values=miR_Family_Info$miR.family)

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

numBS_dir <- list.files()
numBS_dir <- numBS_dir[str_detect(numBS_dir, "APA4")]
numBS_name <- str_extract(numBS_dir, "_[:upper:]*")
numBS_name <- str_remove(numBS_name, "_")
data_list <- vector(mode = "list", length = length(numBS_name))
names(data_list) <- numBS_name
for (i in c(1:11,13)){
  numBS_tmp <- read.delim(numBS_dir[i],header = F)
  numBS_tmp <- header.true(numBS_tmp)
  patName <- numBS_tmp$patName
  numBS_tmp <- apply(numBS_tmp[,-1],2,as.numeric)
  numBS_tmp <- as.data.frame(numBS_tmp)
  numBS_tmp <- cbind(patName, numBS_tmp)
  patientID <- str_sub(numBS_tmp$patName, 1, 12)
  numBS_tmp <- cbind(patientID, numBS_tmp)
  numBS_match <- numBS_tmp[numBS_tmp$patientID %in% patientID[duplicated(patientID)],]
  sampleType <- str_sub(numBS_match$patName, 14, -1)
  numBS_match <- cbind(sampleType, numBS_match)
  numBS_match <- numBS_match[order(numBS_match$patientID, numBS_match$sampleType),]
  deltaBS <- data.frame()
  for(i.tmp in 4:ncol(numBS_match)){
    numBS_match[,i.tmp] <- log2(numBS_match[,i.tmp]+1)
  }
  for (j in 1:length(unique(numBS_match$patientID))){
    pat.tmp <- unique(numBS_match$patientID)[j]
    pat.numBS <- numBS_match[numBS_match$patientID == pat.tmp,]
    deltaBS <- rbind(deltaBS, pat.numBS[2,-(1:3)]-pat.numBS[1,-(1:3)]) #log2(tumor/normal)
  }
  deltaBS <- cbind(unique(numBS_match$patientID), deltaBS)
  colnames(deltaBS)[1] <- "patID"
  
  miRNA_tmp <- read.delim(paste0("/zfs1/hpark/SUR3UTR_immune/data/",numBS_name[i],"/miRNA_HiSeq_gene"))
  for(i.tmp in 2:ncol(miRNA_tmp)){
    miRNA_tmp[,i.tmp] <- 2^miRNA_tmp[,i.tmp]-1
  }
  miRNA_tmp <- miRNA_tmp[miRNA_tmp$sample %in% miR_Family_Info$MiRBase.Accession,]
  miRNA_tmp$familyID <- values(miR_hash, keys = miRNA_tmp$sample)
  sum(duplicated(miRNA_tmp$familyID))
  miRNA_tmp[is.na(miRNA_tmp)] <- 0
  
  miRNA_tmp <- miRNA_tmp %>% 
    group_by(familyID) %>% 
    summarise_at(vars(-sample),sum)
  miRNA_tmp <- as.data.frame(miRNA_tmp)
  rownames(miRNA_tmp) <- miRNA_tmp$familyID
  ## transfer back to log2(RPM+1)
  for(i.tmp in 2:ncol(miRNA_tmp)){
    miRNA_tmp[,i.tmp] <- log2(miRNA_tmp[,i.tmp]+1)
  }
  
  miRNA_tmp <- t(miRNA_tmp[,-1])
  sampleType <- str_sub(rownames(miRNA_tmp), 14,-1)
  patientID <- str_sub(rownames(miRNA_tmp), 1, 12)
  miRNA_tmp <- as.data.frame(miRNA_tmp)
  miRNA_tmp <- cbind(sampleType, patientID, miRNA_tmp)
  miRNA_tmp <- miRNA_tmp[miRNA_tmp$patientID %in% miRNA_tmp$patientID[duplicated(miRNA_tmp$patientID)],]
  miRNA_tmp <- miRNA_tmp[order(miRNA_tmp$patientID, miRNA_tmp$sampleType),]
  deltamiR <- data.frame()
  for (j in 1:length(unique(miRNA_tmp$patientID))){
    pat.tmp <- unique(miRNA_tmp$patientID)[j]
    pat.miR <- miRNA_tmp[miRNA_tmp$patientID == pat.tmp,]
    deltamiR <- rbind(deltamiR, pat.miR[1,-(1:2)]-pat.miR[nrow(pat.miR),-(1:2)]) #log2(tumor/normal)
  }
  deltamiR <- cbind(unique(miRNA_tmp$patientID), deltamiR)
  colnames(deltamiR)[1] <- "patID"
  
  #commonPat <- intersect(intersect(rownames(mRNA.tmp),deltaBS$patID),deltamiR$patID)
  commonPat <- intersect(deltaBS$patID,deltamiR$patID)
  #mRNATumor <- mRNA.tmp[rownames(mRNA.tmp) %in% commonPat,]
  deltaBS <- deltaBS[deltaBS$patID %in% commonPat,]
  deltamiR <- deltamiR[deltamiR$patID %in% commonPat,]
  deltaBS <- deltaBS[,colnames(deltaBS) %in% colnames(deltamiR)]
  deltamiR <- deltamiR[,colnames(deltamiR) %in% colnames(deltaBS)]
  ## due to large amount of NAs in miR expression table, only 200 miRs remained.
  ## if proportion of NAs in miR expression is not considered, we will have 571 miRs. 
  data_list[[i]] <- list(deltaBS=deltaBS, deltamiR=deltamiR) #mRNATumor
}
data_list <- data_list[c(1:2,4:9, 11,13)]
sapply(data_list, function(x){print(nrow(x[[1]]))})

#BLCA BRCA HNSC KICH KIRC KIRP LIHC LUAD PRAD STAD 
#19   73   40   23   51   30   48   12   52   23 

## BRCA
clinic <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/BRCA/BRCA_clinicalMatrix")
colnames(clinic)[str_detect(colnames(clinic), "stage") | str_detect(colnames(clinic), "Stage")]
table(clinic$gender)
## 5 NA, 13 male. match patients first, then, consider whether do we need to keep gender in the model
clinic <- cbind(as.data.frame(do.call(rbind, str_split(clinic$sampleID, "\\-"))), clinic)
clinic$V4 <- as.integer(clinic$V4)
clinic <- clinic[clinic$V4 < 10,]
patientID <- paste0(clinic$V1, ".",clinic$V2, ".",clinic$V3)
clinic <- cbind(patientID, clinic)
clinic <- clinic[!duplicated(clinic$patientID),]
clinic <- clinic[clinic$patientID %in% data_list$BRCA$deltaBS$patID,]
table(clinic$gender)
## all female
clinic <- clinic[,colnames(clinic) %in% c("patientID","OS_Time_nature2012", "OS_event_nature2012", "pathologic_stage","age_at_initial_pathologic_diagnosis")]
clinic <- clinic[!is.na(clinic$OS_Time_nature2012),]
clinic <- clinic[clinic$pathologic_stage!="",]
clinic <- clinic[clinic$OS_Time > 0,]
sum(clinic$patientID %in% data_list$BRCA$deltaBS$patID)
deltaBS <- data_list$BRCA$deltaBS
deltamiR <- data_list$BRCA$deltamiR
#mRNA <- data_list$BRCA$mRNATumor

clinic$pathologic_stage <- str_remove(clinic$pathologic_stage, "A")
clinic$pathologic_stage <- str_remove(clinic$pathologic_stage, "B")
clinic <- clinic[clinic$pathologic_stage!="Stage IV",]


deltaBS <- deltaBS[deltaBS$patID %in% clinic$patientID,]
deltamiR <- deltamiR[deltamiR$patID %in% clinic$patientID,]
deltamiR <- deltamiR[match(deltaBS$patID, deltamiR$patID),]
clinic <- clinic[match(deltaBS$patID, clinic$patientID),]


## clinical var alone
cox_mat <- cbind(OS_Time=clinic$OS_Time_nature2012, OS=clinic$OS_event_nature2012, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage)
rownames(cox_mat) <- clinic$patientID
#cox_mat <- cox_mat[,!(colnames(cox_mat) %in% c("SNAR-A2", "OR10G9"))]
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
clinic_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    clinic_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    clinic_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, clinic_strata)
colnames(clinic)[2:3] <- c("OS_Time", "OS")
BRCA_plt_clinic <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ clinic_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ clinic_strata, data = clinic)
# pval 0.09
# Select miRNA
identical(deltaBS$patID, rownames(cox_mat))
identical(deltaBS$patID, rownames(cox_mat))
lasso_tab <- cbind(cox_mat, deltaBS[,-1])
lasso_tab <- cbind(model.matrix(~-1+stage, data=cox_mat), lasso_tab[,-4])
select_feature <- character()
set.seed("02092022")
for(i in 1:100){
  cv.fit.tmp <- cv.glmnet(as.matrix(lasso_tab)[,-c(4:5)], Surv(time = lasso_tab$OS_Time, event = lasso_tab$OS), nfolds = 10,
                          family = "cox", type.measure="deviance", alpha=1,penalty.factor=c(rep(0,4), rep(1,(ncol(deltaBS)-1))))
  select_feature <- c(select_feature, rownames(coef(cv.fit.tmp, s = 'lambda.min'))[coef(cv.fit.tmp, s = 'lambda.min')[,1]!= 0])
  
}
setdiff(select_feature, colnames(lasso_tab)[1:6])
## 100 runs, miR-3659, 17
table(select_feature)

## deltaBS
identical(deltaBS$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, miR_3659=deltaBS$`miR-3659`)

rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_3659 <- as.numeric(cox_mat$miR_3659)
cox_mat$stage <- as.factor(cox_mat$stage)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -c(5:7)])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.0208

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltaBS_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltaBS_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltaBS_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltaBS_strata)
#colnames(clinic)[2:3] <- c("OS_Time", "OS")
BRCA_plt_deltaBS<- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic)
# log-rank P: 0.008

## selected miRNA numTS alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ miR_3659 , data = cox_mat)
summary(coxmodel)

## miR expression
identical(deltamiR$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, miR_3659=deltamiR$`miR-3659`)
rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_3659 <- as.numeric(cox_mat$miR_3659)
cox_mat$stage <- as.factor(cox_mat$stage)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -c(5:7)])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.6998
top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltamiR_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltamiR_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltamiR_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltamiR_strata)
#colnames(clinic)[2:3] <- c("OS_Time", "OS")
BRCA_plt_deltamiR <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic)
# 0.1
## selected miRNA expression alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ miR_3659 , data = cox_mat)
summary(coxmodel)

BRCA_plt_clinic
BRCA_plt_deltaBS
BRCA_plt_deltamiR

grid.arrange(grobs=list(BRCA_plt_clinic$plot,BRCA_plt_deltaBS$plot,BRCA_plt_deltamiR$plot),
             nrow = 1, ncol = 3)



## HNSC
clinic <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/HNSC/HNSC_clinicalMatrix")
colnames(clinic)[str_detect(colnames(clinic), "stage") | str_detect(colnames(clinic), "Stage")]
#View(clinic[,str_detect(colnames(clinic), "stage") | str_detect(colnames(clinic), "Stage")])
## pathologic_stage has 73 missing. clinical_stage has 14 missing. There is discrepancy between these two columns
table(clinic$gender)
clinic <- cbind(as.data.frame(do.call(rbind, str_split(clinic$sampleID, "\\-"))), clinic)
clinic$V4 <- as.integer(clinic$V4)
clinic <- clinic[clinic$V4 < 10,]
patientID <- paste0(clinic$V1, ".",clinic$V2, ".",clinic$V3)
clinic <- cbind(patientID, clinic)
clinic <- clinic[!duplicated(clinic$patientID),]
clinic <- clinic[clinic$patientID %in% data_list$HNSC$deltaBS$patID,]
table(clinic$gender)
## 12 female, 28 male
clinic <- clinic[,colnames(clinic) %in% c("patientID","OS.time", "OS", "pathologic_stage","age_at_initial_pathologic_diagnosis", "gender")]
clinic <- clinic[!is.na(clinic$OS.time),]
table(clinic$pathologic_stage)
clinic <- clinic[clinic$pathologic_stage!="",]
clinic <- clinic[clinic$OS.time > 0,]
sum(clinic$patientID %in% data_list$HNSC$deltaBS$patID)
deltaBS <- data_list$HNSC$deltaBS
deltamiR <- data_list$HNSC$deltamiR
#mRNA <- data_list$HNSC$mRNATumor

deltaBS <- deltaBS[deltaBS$patID %in% clinic$patientID,]
deltamiR <- deltamiR[deltamiR$patID %in% clinic$patientID,]
deltamiR <- deltamiR[match(deltaBS$patID, deltamiR$patID),]
clinic <- clinic[match(deltaBS$patID, clinic$patientID),]


## clinical var alone
cox_mat <- cbind(OS_Time=clinic$OS.time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender)
rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
clinic_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    clinic_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    clinic_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, clinic_strata)
colnames(clinic)[2:3] <- c("OS_Time", "OS")
HNSC_plt_clinic <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ clinic_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ clinic_strata, data = clinic)
# pval 0.1
# Select miRNA
identical(deltaBS$patID, rownames(cox_mat))
identical(deltaBS$patID, rownames(cox_mat))
lasso_tab <- cbind(cox_mat, deltaBS[,-1])
lasso_tab <- cbind(model.matrix(~-1+stage, data=cox_mat), lasso_tab[,-4])
lasso_tab$gender <- ifelse(lasso_tab$gender == "FEMALE", 0, 1)
select_feature <- character()
set.seed("02092022")
for(i in 1:100){
  cv.fit.tmp <- cv.glmnet(as.matrix(lasso_tab)[,-c(5:6)], Surv(time = lasso_tab$OS_Time, event = lasso_tab$OS), nfolds = 10,
                          family = "cox", type.measure="deviance", alpha=1,penalty.factor=c(rep(0,6), rep(1,(ncol(deltaBS)-1))))
  select_feature <- c(select_feature, rownames(coef(cv.fit.tmp, s = 'lambda.min'))[coef(cv.fit.tmp, s = 'lambda.min')[,1]!= 0])
  
}
setdiff(select_feature,colnames(lasso_tab)[1:8])
## 100 runs, "miR-663/663a/1908" 5  "miR-651"  3
table(select_feature)

## deltaBS
identical(deltaBS$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender,miR_663=deltaBS$`miR-663/663a/1908`, miR_651=deltaBS$`miR-651`)

rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_663 <- as.numeric(cox_mat$miR_663)
cox_mat$miR_651 <- as.numeric(cox_mat$miR_651)
cox_mat$stage <- as.factor(cox_mat$stage)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -c(6:7)])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.001041

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltaBS_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltaBS_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltaBS_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltaBS_strata)
HNSC_plt_deltaBS<- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic)
# log-rank P: 0.004

## selected miRNA numTS alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ `miR_663` , data = cox_mat)
summary(coxmodel)
coxmodel <- coxph(Surv((OS_Time), OS) ~ `miR_651` , data = cox_mat)
summary(coxmodel)


## miR expression
identical(deltamiR$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender,miR_663=deltamiR$`miR-663/663a/1908`, miR_651=deltamiR$`miR-651`)
rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_663 <- as.numeric(cox_mat$miR_663)
cox_mat$miR_651 <- as.numeric(cox_mat$miR_651)
cox_mat$stage <- as.factor(cox_mat$stage)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -c(6:7)])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.3951
top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltamiR_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltamiR_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltamiR_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltamiR_strata)
#colnames(clinic)[2:3] <- c("OS_Time", "OS")
HNSC_plt_deltamiR <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic)
# 0.04

## selected miRNA expression alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ miR_663 , data = cox_mat)
summary(coxmodel)
coxmodel <- coxph(Surv((OS_Time), OS) ~ miR_651 , data = cox_mat)
summary(coxmodel)


HNSC_plt_clinic
HNSC_plt_deltaBS
HNSC_plt_deltamiR

grid.arrange(grobs=list(BRCA_plt_clinic$plot,BRCA_plt_deltaBS$plot,BRCA_plt_deltamiR$plot,
                        HNSC_plt_clinic$plot,HNSC_plt_deltaBS$plot,HNSC_plt_deltamiR$plot),
             nrow = 2, ncol = 3)

## KIRC
clinic <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/KIRC/KIRC_clinicalMatrix")
colnames(clinic)[str_detect(colnames(clinic), "stage") | str_detect(colnames(clinic), "Stage")]
table(clinic$gender)
clinic <- cbind(as.data.frame(do.call(rbind, str_split(clinic$sampleID, "\\-"))), clinic)
clinic$V4 <- as.integer(clinic$V4)
clinic <- clinic[clinic$V4 < 10,]
patientID <- paste0(clinic$V1, ".",clinic$V2, ".",clinic$V3)
clinic <- cbind(patientID, clinic)
clinic <- clinic[!duplicated(clinic$patientID),]
clinic <- clinic[clinic$patientID %in% data_list$KIRC$deltaBS$patID,]
table(clinic$gender)
## 16 female, 35 male
clinic <- clinic[,colnames(clinic) %in% c("patientID","OS.time", "OS", "pathologic_stage","age_at_initial_pathologic_diagnosis", "gender")]
clinic <- clinic[!is.na(clinic$OS.time),]
table(clinic$pathologic_stage)
clinic <- clinic[clinic$pathologic_stage!="",]
clinic <- clinic[clinic$OS.time > 0,]
sum(clinic$patientID %in% data_list$KIRC$deltaBS$patID)
deltaBS <- data_list$KIRC$deltaBS
deltamiR <- data_list$KIRC$deltamiR
#mRNA <- data_list$KIRC$mRNATumor
deltaBS <- deltaBS[deltaBS$patID %in% clinic$patientID,]
deltamiR <- deltamiR[deltamiR$patID %in% clinic$patientID,]
deltamiR <- deltamiR[match(deltaBS$patID, deltamiR$patID),]
#mRNA <- mRNA[match(deltaBS$patID, rownames(mRNA)),]
clinic <- clinic[match(deltaBS$patID, clinic$patientID),]


## clinical var alone
cox_mat <- cbind(OS_Time=clinic$OS.time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender)
rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
clinic_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    clinic_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    clinic_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, clinic_strata)
colnames(clinic)[2:3] <- c("OS_Time", "OS")
KIRC_plt_clinic <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ clinic_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ clinic_strata, data = clinic)
# pval 4e-05
# Select miRNA
identical(deltaBS$patID, rownames(cox_mat))
identical(deltaBS$patID, rownames(cox_mat))
lasso_tab <- cbind(cox_mat, deltaBS[,-1])
lasso_tab <- cbind(model.matrix(~-1+stage, data=cox_mat), lasso_tab[,-4])
lasso_tab$gender <- ifelse(lasso_tab$gender == "FEMALE", 0, 1)
select_feature <- character()
set.seed("02092022")
for(i in 1:100){
  cv.fit.tmp <- cv.glmnet(as.matrix(lasso_tab)[,-c(5:6)], Surv(time = lasso_tab$OS_Time, event = lasso_tab$OS), nfolds = 10,
                          family = "cox", type.measure="deviance", alpha=1,penalty.factor=c(rep(0,6), rep(1,(ncol(deltaBS)-1))))
  select_feature <- c(select_feature, rownames(coef(cv.fit.tmp, s = 'lambda.min'))[coef(cv.fit.tmp, s = 'lambda.min')[,1]!= 0])
  
}
setdiff(select_feature, colnames(lasso_tab)[1:8])
## 100 runs, "miR-3917"  76
table(select_feature)

## deltaBS
identical(deltaBS$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender,miR_3917=deltaBS$`miR-3917`)

rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_3917 <- as.numeric(cox_mat$miR_3917)

cox_mat$stage <- as.factor(cox_mat$stage)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -6])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.003926

top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltaBS_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltaBS_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltaBS_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltaBS_strata)
#colnames(clinic)[2:3] <- c("OS_Time", "OS")
KIRC_plt_deltaBS<- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltaBS_strata, data = clinic)
# log-rank P: 3e-07

## selected miRNA numTS alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ `miR_3917` , data = cox_mat)
summary(coxmodel)


## miR expression
identical(deltamiR$patID, clinic$patientID)
cox_mat <- cbind(OS_Time=clinic$OS_Time, OS=clinic$OS, age = clinic$age_at_initial_pathologic_diagnosis, stage=clinic$pathologic_stage, gender=clinic$gender,miR_3917=deltamiR$`miR-3917`)
rownames(cox_mat) <- clinic$patientID
cox_mat <- as.data.frame(cox_mat)
cox_mat$OS_Time <- as.integer(cox_mat$OS_Time)
cox_mat$OS <- as.integer(cox_mat$OS)
cox_mat$age <- as.integer(cox_mat$age)
cox_mat$miR_3917 <- as.numeric(cox_mat$miR_3917)
cox_mat$stage <- as.factor(cox_mat$stage)
cox_mat$gender <- as.factor(cox_mat$gender)
coxmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat)
predict.hazard <- predict(coxmodel, cox_mat[,-c(1:2)])
clinicmodel <- coxph(Surv((OS_Time), OS) ~ . , data = cox_mat[, -6])
lrtest(coxmodel, clinicmodel)
# LRT pval = 0.2696
top.pat <- names(predict.hazard)[predict.hazard >= median(predict.hazard)]
bot.pat <- names(predict.hazard)[predict.hazard < median(predict.hazard)]
deltamiR_strata <- character(length = nrow(clinic))
for (i in 1:nrow(clinic)){
  if (clinic$patientID[i] %in% top.pat){
    deltamiR_strata[i] <- "high"
  }else if(clinic$patientID[i] %in% bot.pat){
    deltamiR_strata[i] <- "low"
  }
}
clinic <- cbind(clinic, deltamiR_strata)
#colnames(clinic)[2:3] <- c("OS_Time", "OS")
KIRC_plt_deltamiR <- ggsurvplot(
  fit = survfit(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic), palette=c("red","blue"),
  xlab = "Days", 
  ylab = "Overall survival probability")
survdiff(Surv((OS_Time), OS) ~ deltamiR_strata, data = clinic)
# 3e-06
## selected miRNA expression alone
coxmodel <- coxph(Surv((OS_Time), OS) ~ `miR_3917` , data = cox_mat)
summary(coxmodel)

KIRC_plt_clinic
KIRC_plt_deltaBS
KIRC_plt_deltamiR


grid.arrange(grobs=list(BRCA_plt_clinic$plot,BRCA_plt_deltaBS$plot,BRCA_plt_deltamiR$plot,
                        HNSC_plt_clinic$plot,HNSC_plt_deltaBS$plot,HNSC_plt_deltamiR$plot,
                        KIRC_plt_clinic$plot,KIRC_plt_deltaBS$plot,KIRC_plt_deltamiR$plot),
             nrow = 3, ncol = 3)

clinicmodel %>% 
  gtsummary::tbl_regression(exp = TRUE) 

LRT_Ptab <- data.frame(Pval = c(0.0208,0.6998,0.001041,0.3951,0.003926,0.2696), Model=rep(c("clinical + deltaTS", "clinical + deltaExpr"), 3), Cancer = rep(c("BRCA","HNSC","KIRC"), each=2))
LRT_Ptab$Model <- factor(LRT_Ptab$Model, levels = c("clinical + deltaTS", "clinical + deltaExpr"))
ggplot(LRT_Ptab, aes(x=Cancer, y=-log10(Pval), fill = Model))+geom_bar(stat="identity", position=position_dodge()) +ylab("Additional Predictive Power \n -log10 (LRT P-value)") +theme_classic()+scale_fill_manual(values=c("red","blue"))#+
scale_y_continuous(trans='log10')#+scale_y_reverse()
