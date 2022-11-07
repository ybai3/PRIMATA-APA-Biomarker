library(survival)
library(survminer)
library(stringr)
library(glmnet)
library(ROCR)
library(survivalROC)
library(survcomp)

readin_clinical <- function(cancertype){
  mRNA <- read.delim(paste0(mRNA_dir,"/TCGA-",cancertype, "_FPKM.txt"), stringsAsFactors=FALSE)
  numBS <- read.delim(paste0(APA4_dir, "/APA4immuneProlif_", cancertype, ".txt"), stringsAsFactors=FALSE)
  miRNA <- read.delim(paste0(Ctrl4_dir, "/ctrl4immuneProlif_", cancertype, ".txt"), stringsAsFactors=FALSE)
  PDUI <- read.delim(paste0(PDUI_dir, "/TCGA_", cancertype,"_Combined_PDUIs.txt"))
  colnames(mRNA) <- str_replace_all(colnames(mRNA), pattern = "\\.", replacement = "-")
  mRNA_dup <- str_sub(colnames(mRNA)[-1], start = 6, end = -5)[which(duplicated(str_sub(colnames(mRNA)[-1], start = 6, end = -5)))]
  mRNA_remove <- character()
  if(length(mRNA_dup)>0){
    for(i in 1:length(mRNA_dup)){
      candidate <- colnames(mRNA)[str_detect(colnames(mRNA), mRNA_dup[i])]
      if (sum(str_sub(candidate, 14,15)=="01")==1){
        mRNA_remove <- c(mRNA_remove, candidate[which(str_sub(candidate, 14,15)!="01")])
      }else if(sum(str_sub(candidate, 14,15)=="01")==0){
        mRNA_remove <- c(mRNA_remove, candidate[-1])
      }else{
        mRNA_remove <- c(mRNA_remove, candidate[c(which(str_sub(candidate, 14,15)!="01"),which(str_sub(candidate, 14,15)=="01")[-1])])
      }
    }
    mRNA <- mRNA[,!(colnames(mRNA) %in% mRNA_remove)]
  }
  colnames(mRNA)[-1] <- str_sub(colnames(mRNA)[-1], start = 6, end = -5)
  mRNA_gene <- mRNA$geneName
  mRNA_t <- as.data.frame(t(mRNA[,-1]))
  colnames(mRNA_t) <- mRNA_gene
  numBS <- numBS[numBS$patName %in% rownames(mRNA_t),]
  numBS <- numBS[,-c((ncol(numBS)-1):ncol(numBS))]
  miRNA <- miRNA[miRNA$patName %in% rownames(mRNA_t),]
  miRNA <- miRNA[,-c((ncol(miRNA)-1):ncol(miRNA))]
  mRNA_t <- mRNA_t[rownames(mRNA_t) %in% numBS$patName,]
  mRNA_t$patID <- rownames(mRNA_t)
  mRNA_t <- mRNA_t[,c(ncol(mRNA_t), 1:(ncol(mRNA_t)-1))]
  mRNA_t <- mRNA_t[match(numBS$patName, mRNA_t$patID),]
  rownames(miRNA) <- miRNA$patName
  rownames(numBS) <- numBS$patName
  
  for(i_tmp in 2:ncol(mRNA_t)){
    mRNA_t[,i_tmp] <- log2(mRNA_t[,i_tmp] + 1)
  }
  for(i_tmp in 2:ncol(numBS)){
    numBS[,i_tmp] <- log2(numBS[,i_tmp] + 1)
  }
  
  PDUI <- PDUI[,-c(2:3)]
  PDUI_tmp <- str_sub(colnames(PDUI), start = 6, end = 12)
  PDUI_dup <- PDUI_tmp[duplicated(PDUI_tmp)]
  PDUI_remove <- character()
  if(length(PDUI_dup)>0){
    for(i in 1:length(PDUI_dup)){
      candidate <- colnames(PDUI)[str_detect(colnames(PDUI), PDUI_dup[i])]
      if (sum(str_sub(candidate, 14)=="01")==1){
        PDUI_remove <- c(PDUI_remove, candidate[which(str_sub(candidate, 14)!="01")])
      }else if(sum(str_sub(candidate, 14)=="01")==0){
        PDUI_remove <- c(PDUI_remove, candidate[-1])
      }else{
        PDUI_remove <- c(PDUI_remove, candidate[c(which(str_sub(candidate, 14)!="01"),which(str_sub(candidate, 14)=="01")[-1])])
      }
    }
    PDUI <- PDUI[,!(colnames(PDUI) %in% PDUI_remove)]
  }
  colname.tmp <- str_sub(colnames(PDUI), start = 6, end = 12)
  colnames(PDUI) <- str_replace(colname.tmp, "\\.", "-")
  rownames(PDUI) <- PDUI[,1]
  PDUI <- PDUI[ ,colnames(PDUI) %in% numBS$patName]
  PDUI <- t(PDUI)
  PDUI <- as.data.frame(PDUI)
  PDUI$patName <- rownames(PDUI)
  PDUI <- PDUI[, c(ncol(PDUI), 1:(ncol(PDUI)-1))]
  PDUI <- PDUI[match(numBS$patName, PDUI$patName),]
  PDUI <- PDUI[,c(TRUE, colSums(is.na(PDUI[,-1]))==0)]
  
  return(list(mRNA=mRNA_t, numBS=numBS, miRNA=miRNA,PDUI = PDUI))
}

