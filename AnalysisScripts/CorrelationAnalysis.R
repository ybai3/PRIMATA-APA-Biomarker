## Figure 1.B and Supp. Figure 2

library(stringr)
library(dplyr)
library(hash)
library(ggplot2)

setwd("/ihome/hpark/yub20/SUR3UTR/pythons/TCGA_FPKM/JNCI_TumorNormal/result_JNCI")
mRNA_dir <- "/ihome/hpark/yub20/SUR3UTR/data/TCGA_harmonized/FPKM"
surv_dir <- "/ihome/hpark/yub20/SUR3UTR/new/survival"
#APA4_dir <- "/ihome/hpark/yub20/SUR3UTR/new/data"
miR_dir <- "/ihome/hpark/yub20/SUR3UTR/data"
miRNA_conversion <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/miRNA_conversion.txt")
miRNA_conversion <- miRNA_conversion[str_detect(miRNA_conversion$ID, "hsa-"),]

miR_Family_Info <- read.delim("/zfs1/hpark/SUR3UTR_immune/data/refs/miR_Family_Info.txt")
miR_Family_Info <- miR_Family_Info[str_detect(miR_Family_Info$MiRBase.ID, "hsa-"),]
miR_Family_Info$match4marker <- str_remove(miR_Family_Info$MiRBase.ID, "hsa-")
miR_name <- str_replace_all(str_remove(c(miRNA_conversion$Mature1_ID, miRNA_conversion$Mature2_ID), "^[:alpha:]*-"), "-", "\\.")
#hash_tab <- data.frame(keys= c(miRNA_conversion$Mature1_Acc, miRNA_conversion$Mature2_Acc), values=miR_name)
#hash_tab <- hash_tab[hash_tab$keys!="",]
miR_hash <- hash(keys= miR_Family_Info$MiRBase.Accession, values=miR_Family_Info$miR.family)

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

numBS_name <- c("HNSC","KIRC","LUAD","OV","STAD","UCEC",
                "LGG","BRCA","LUSC","SKCM")
numBS_dir <- paste0("APA4immuneProlif_", numBS_name,".txt")
data_list <- vector(mode = "list", length = length(numBS_name))
names(data_list) <- numBS_name
for (i in c(1:3,5,8:9)){
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
  for(i_tmp in 4:ncol(numBS_match)){
    numBS_match[,i_tmp] <- log2(numBS_match[,i_tmp] +1)
  }
  
  miRNA_tmp <- read.delim(paste0("/zfs1/hpark/SUR3UTR_immune/data/",numBS_name[i],"/miRNA_HiSeq_gene"))
  miRNA_tmp <- miRNA_tmp[miRNA_tmp$sample %in% miR_Family_Info$MiRBase.Accession,]
  miRNA_tmp$familyID <- values(miR_hash, keys = miRNA_tmp$sample)
  sum(duplicated(miRNA_tmp$familyID))
  miRNA_tmp <- aggregate(. ~ familyID, data=miRNA_tmp[,-1], FUN=sum) ## by default, na.rm =T
  #miRNA_tmp <- miRNA_tmp[!duplicated(miRNA_tmp$familyID),]
  rownames(miRNA_tmp) <- miRNA_tmp$familyID
  miRNA_tmp <- miRNA_tmp[,-ncol(miRNA_tmp)]
  miRNA_tmp <- t(miRNA_tmp[,-1])
  
  sampleType <- str_sub(rownames(miRNA_tmp), 14,-1)
  patientID <- str_sub(rownames(miRNA_tmp), 1, 12)
  miRNA_tmp <- as.data.frame(miRNA_tmp)
  miRNA_tmp <- cbind(sampleType, patientID, miRNA_tmp)
  miRNA_tmp <- miRNA_tmp[miRNA_tmp$patientID %in% miRNA_tmp$patientID[duplicated(miRNA_tmp$patientID)],]
  miRNA_tmp <- miRNA_tmp[order(miRNA_tmp$patientID, miRNA_tmp$sampleType),]
  
  numBS_match=numBS_match[,-((ncol(numBS_match)-1):ncol(numBS_match))]
  numBS_match <- numBS_match[,colnames(numBS_match) %in% colnames(miRNA_tmp)]
  miRNA_tmp <- miRNA_tmp[,colnames(miRNA_tmp) %in% colnames(numBS_match)]
  numBS_match <- numBS_match[numBS_match$patientID %in% miRNA_tmp$patientID,]
  miRNA_tmp <- miRNA_tmp[miRNA_tmp$patientID %in% numBS_match$patientID,]
  
  numBS_tumor <- numBS_match[numBS_match$sampleType=="Tumor",]
  numBS_normal <- numBS_match[numBS_match$sampleType=="Normal",]
  
  miRNA_tumor <- miRNA_tmp[as.integer(miRNA_tmp$sampleType)<10,]
  miRNA_normal <- miRNA_tmp[as.integer(miRNA_tmp$sampleType)>=10,]
  
  #data_list[[i]] <- list(numBS_tumor=numBS_tumor, numBS_normal=numBS_normal, miRNA_tumor=miRNA_tumor, miRNA_normal=miRNA_normal)
  
  numBS_tumor <- numBS_tumor[match(miRNA_tumor$patientID, numBS_tumor$patientID),]
  numBS_normal <- numBS_normal[match(miRNA_normal$patientID, numBS_normal$patientID),]
  numBS_tumor <- numBS_tumor[,match(colnames(miRNA_tumor), colnames(numBS_tumor))]
  numBS_normal <- numBS_normal[,match(colnames(miRNA_normal), colnames(numBS_normal))]
  
  numBS_tumor_sd <- apply(numBS_tumor[,-c(1:2)], 2, sd)/apply(numBS_tumor[,-c(1:2)], 2, mean)
  numBS_normal_sd <- apply(numBS_normal[,-c(1:2)], 2, sd)/apply(numBS_normal[,-c(1:2)], 2, mean)
  miRNA_tumor_sd <- apply(miRNA_tumor[,-c(1:2)], 2, sd)/apply(miRNA_tumor[,-c(1:2)], 2, mean)
  miRNA_normal_sd <- apply(miRNA_normal[,-c(1:2)], 2, sd)/apply(miRNA_normal[,-c(1:2)], 2, mean)
  
  sd_value = c(numBS_tumor_sd, numBS_normal_sd, miRNA_tumor_sd, miRNA_normal_sd)
  sd_tab <- data.frame(sd=sd_value, cancertype=rep(numBS_name[i], length(sd_value)),
                       measure = c(rep("numBS", length(c(numBS_tumor_sd, numBS_normal_sd))), rep("Expr", length(c(miRNA_tumor_sd, miRNA_normal_sd)))),
                       sampletype= c(rep("tumor",length(numBS_tumor_sd)), rep("normal", length(numBS_normal_sd)), rep("tumor", length(miRNA_tumor_sd)), rep("normal", length(miRNA_normal_sd))))
  #ave_value = c(numBS_tumor_ave, numBS_normal_ave, miRNA_tumor_ave, miRNA_normal_ave)
  #ave_tab <- data.frame(ave=ave_value, cancertype=rep(numBS_name[i], length(ave_value)),
  #                      measure = c(rep("numBS", length(c(numBS_tumor_ave, numBS_normal_ave))), rep("Expr", length(c(miRNA_tumor_ave, miRNA_normal_ave)))),
  #                      sampletype= c(rep("tumor",length(numBS_tumor_ave)), rep("normal", length(numBS_normal_ave)), rep("tumor", length(miRNA_tumor_ave)), rep("normal", length(miRNA_normal_ave))))
  cor.tumor <- data.frame()  
  cor.normal <- data.frame()
  for (j in 3:ncol(numBS_tumor)){
    cor.tmp <- cor.test(numBS_tumor[,j], miRNA_tumor[,j], method = "spearman")
    cor.tumor <- rbind(cor.tumor, c(miRNA = colnames(numBS_tumor)[j], rho = cor.tmp$estimate, pval = cor.tmp$p.value))
    cor.tmp <- cor.test(numBS_normal[,j], miRNA_normal[,j], method = "spearman")
    cor.normal <- rbind(cor.normal, c(miRNA = colnames(numBS_normal)[j], rho = cor.tmp$estimate, pval = cor.tmp$p.value))
  }
  colnames(cor.tumor) <- c("miRNA", "rho", "pval")
  colnames(cor.normal) <- c("miRNA", "rho", "pval")
  cor.tumor$FDR <- p.adjust(cor.tumor$pval, "BH")
  cor.normal$FDR <- p.adjust(cor.normal$pval, "BH")
  
  data_list[[i]] <- list(cor.tumor=cor.tumor, cor.normal=cor.normal, 
                         sd_tab = sd_tab)#, ave_tab = ave_tab)
}

data_list <- data_list[c(1:3,5,8:9)]
sapply(data_list, function(x){print(nrow(x[[1]]))})
#saveRDS(data_list,"delta_data.RDS")
#data_list <- data_list[c(1:2,5)]

#BRCA HNSC KIRC LUAD LUSC STAD 
#73   40   51   12    7   30 
plot_tab_sd <- data.frame()
#plot_tab_ave <- data.frame()
for(i in 1:length(data_list)){
  plot_tab_sd <- rbind(plot_tab_sd, data_list[[i]]$sd_tab)
  #plot_tab_ave <- rbind(plot_tab_ave, data_list[[i]]$ave_tab)
}
plot_tab_sd$measure <- factor(plot_tab_sd$measure, levels = c("numBS", "Expr"))
ggplot(plot_tab_sd, aes(x=cancertype,y=sd, fill=measure)) + geom_boxplot()+facet_wrap(vars(sampletype))+theme_bw()+ylab("Standard Coefficient")+scale_fill_manual(values=c("red", "blue"))
ggplot(plot_tab_ave, aes(x=cancertype,y=ave, fill=measure)) + geom_boxplot()+facet_wrap(vars(sampletype))

plot_tab <- data.frame()
for(i in 1:length(data_list)){
  plot_tab <- rbind(plot_tab, cbind(data_list[[i]]$cor.tumor, cancer=rep(names(data_list)[i], nrow(data_list[[i]]$cor.tumor)), type=rep("tumor", nrow(data_list[[i]]$cor.tumor))))
  plot_tab <- rbind(plot_tab, cbind(data_list[[i]]$cor.normal, cancer=rep(names(data_list)[i], nrow(data_list[[i]]$cor.normal)), type=rep("normal", nrow(data_list[[i]]$cor.normal))))
}
plot_tab$rho <- as.numeric(plot_tab$rho)
plot_tab$pval <- as.numeric(plot_tab$pval)
ggplot(plot_tab,aes(x=cancer, y=-log10(pval), color=type))+geom_boxplot()+theme_bw()+geom_hline(yintercept=-log10(0.05), color="red",linetype = "dashed")
positive_plt <- ggplot(plot_tab[plot_tab$rho>0,],aes(x=cancer, y=-log10(pval), fill=type))+geom_boxplot()+theme_bw()
negative_plt <- ggplot(plot_tab[plot_tab$rho<0,],aes(x=cancer, y=-log10(pval), fill=type))+geom_boxplot()+theme_bw()
library(gridExtra)
grid.arrange(positive_plt , negative_plt, ncol=1)
ggplot(plot_tab,aes(x=cancer, y=rho, fill=type))+geom_boxplot()
ggplot(plot_tab,aes(x=cancer, y=rho, fill=type))+geom_violin(alpha=0.5)

## plot R2 as box plot
plot_tab$R2 <- plot_tab$rho ^2
ggplot(plot_tab,aes(x=cancer, y=R2, fill=type))+geom_boxplot()+theme_bw()
ggplot(plot_tab,aes(x=cancer, y=R2, fill=type))+geom_violin()+theme_bw()
library(ggpubr)
plot_tab$cancer <- factor(plot_tab$cancer, levels = c("BRCA", "KIRC", "HNSC","STAD", "LUAD", "LUSC"))
ggplot(plot_tab,aes(x=cancer, y=R2, fill=type)) + ylab(expression("R"^2))+
  geom_boxplot(alpha=0.8, width=0.8) + scale_fill_manual(values=c("blue", "red"))+# geom_point(position = position_jitterdodge(seed=1,dodge.width = 1), size =0.5) +
  stat_compare_means(aes(group=type), label="p.format", method = "wilcox.test") + #ylim(c(0,0.85))+
  theme(text=element_text(size=10))+theme_classic()
ggplot(plot_tab[plot_tab$cancer %in% c("BRCA", "KIRC", "HNSC", "STAD"),],aes(x=cancer, y=R2, fill=type)) + ylab(expression("R"^2))+
  geom_boxplot(alpha=0.8, width=0.8) + scale_fill_manual(values=c("blue", "red"))+# geom_point(position = position_jitterdodge(seed=1,dodge.width = 1), size =0.5) +
  stat_compare_means(aes(group=type), label="p.format", method = "wilcox.test") + #ylim(c(0,0.85))+
  theme(text=element_text(size=10))+theme_classic()




## Figure 1.C, Supp Fig.1

setwd("C:/Users/byl12/OneDrive - University of Pittsburgh/Graduate/Pitt/ParkLab/IndependentProject/SUR3UTR/data/GDC_TCGA_FPKM/result_TC3A")
source("C:/Users/byl12/OneDrive - University of Pittsburgh/Graduate/Pitt/ParkLab/IndependentProject/SUR3UTR/GDC_FPKM/Independence/functions_logremoved.R")

par(mfrow=c(1,1))
group<-NULL
all_p<-NULL
all_q<-NULL
all_v<-NULL
cancerType1=c("BRCA","SKCM","HNSC","KIRC","LGG","LUAD","LUSC","OV","STAD","UCEC")

#Data_option<-"immune"
Data_option<-"proliferation"
Num<-rep(0,10)
Num1<-rep(0,10)
sampleSize<-rep(0,10)
names(sampleSize) <- cancerType1
for(ind_cancer in 1:length(cancerType1)){
  cancertype=cancerType1[ind_cancer]
  data_all<-process_data(cancertype = cancertype,option="Data_option")
  df1=data_all$df1
  df1<-df1[,1:(ncol(df1)-1)]
  Num[ind_cancer]<-ncol(df1)
  sampleSize[ind_cancer]<-nrow(df1)
  df.ctrl1=data_all$df.ctrl1
  df.ctrl1=df.ctrl1[,1:(ncol(df.ctrl1)-1)]
  index1<-which(apply(df1,2,sum)==0)
  index2<-which(apply(df.ctrl1,2,sum)==0)
  index<-c(index1,index2)
  df1<-df1[,-index]
  Num1[ind_cancer]<-ncol(df1)
  df.ctrl1<-df.ctrl1[,-index]
  
  group<-c(group,rep(cancertype,dim(df1)[2]))
  set.seed(12315)
  p.cor<-rep(0,dim(df1)[2])
  v.cor<-rep(0,dim(df1)[2])
  for(i in 1:dim(df1)[2]){
    temp<-cor.test(df1[,i],df.ctrl1[,i], method = "spearman")
    p.cor[i]<-temp$p.value
    v.cor[i]<-temp$estimate
  }
  all_p<-c(all_p,p.cor)
  q.temp<-rep(NA,length(p.cor))
  q.temp[which(!is.na(p.cor))]<-p.adjust(p.cor[!is.na(p.cor)],method = "BH")
  all_q<-c(all_q,q.temp)
  all_v<-c(all_v,v.cor)
}
library(hash)
samplesize <- hash(keys=names(sampleSize), values= sampleSize)



data<-data.frame(group,all_p,all_q,all_v)
Num_DE<-rep(0,10)
Num_nonDE<-rep(0,10)
percentage_DE<-rep(0,10)
percentage_NonDE<-rep(0,10)
for(ind_cancer in 1:length(cancerType1)){
  cancertype=cancerType1[ind_cancer]
  temp<-data[data$group==cancertype,]
  Num_DE[ind_cancer]<- sum(temp$all_q<0.05)
  Num_nonDE[ind_cancer]<-sum(temp$all_q>=0.05)
  percentage_DE[ind_cancer]<-round(sum(temp$all_q<0.05)/nrow(temp),2)
  percentage_NonDE[ind_cancer]<-1-percentage_DE[ind_cancer]
}
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
#addline_format(cancerType2)
cancerType2<-rep("0",10)
Num2<-rep(0,10)
for(i in 1:10){
  Num2[i]<-paste("(",sampleSize[i],")",sep="")
}
#Num_DE<-round(percentage_DE*Num1)
#Num_nonDE<-Num1-Num_DE
Num<-c(Num_DE,Num_nonDE)
# Num<-rep(0,20)
# for(i in 1:10){
# Num[2*i-1]<-Num_DE[i]
# Num[2*i]<-Num_nonDE[i]    
# 
# }
data_Percentage<-data.frame(c(cancerType1,cancerType1),c(percentage_DE,percentage_NonDE),
                            c(rep("Correlated",10),rep("Not Correlated",10)))
colnames(data_Percentage)<-c("group","ratio","DE")
class(data_Percentage$group)
class(data_Percentage$DE)
class(data_Percentage$ratio)
data_Percentage<-cbind(data_Percentage,Num)
data_Percentage$sampleSize <- values(samplesize, keys=data_Percentage$group)
data_Percentage <- data_Percentage[order(data_Percentage$sampleSize, decreasing = T),]

for(i in 1:10){
  cancerType2[i]<-paste(names(sampleSize)[i],paste("(",sampleSize[i],")",sep=""),sep=" ")
}
data_Percentage$group <- factor(data_Percentage$group, levels = unique(data_Percentage$group))
library(ggplot2)
ggplot(data=data_Percentage, aes(x=group, y=ratio,fill=DE)) +
  geom_bar(stat="identity",position="dodge")+scale_y_continuous(name="Percentage")+ylab("Proportion of miRNAs")+
  scale_x_discrete(name="Cancer type (sample size)",breaks=cancerType1, 
                   labels=addline_format(cancerType2))+#ggtitle("Correlation between numTS and expression")+
  geom_text(aes(label=Num), position=position_dodge(width=0.9), vjust=-0.25)+theme_classic()+
  scale_fill_manual(values=c("black","darkgrey"))+ scale_y_continuous(breaks = seq(0, 1, by = 0.2))

ggplot(data=data_Percentage[data_Percentage$DE == "Not Correlated",], aes(x=group, y=ratio,fill=DE)) +
  geom_bar(stat="identity",position="dodge")+ylab("Percentage of miRNAs")+
  scale_x_discrete(name="Cancer type (sample size)",breaks=cancerType1, 
                   labels=addline_format(cancerType2))+#ggtitle("Correlation between numTS and expression")+
  geom_text(aes(label=Num), position=position_dodge(width=0.9), vjust=-0.25)+theme_classic()+
  scale_fill_manual(values=c("darkgrey"))+ scale_y_continuous(breaks = seq(0, 1, by = 0.2))#+ scale_y_sqrt()

ggplot(data=data_Percentage[data_Percentage$DE == "Correlated",], aes(x=group, y=ratio,fill=DE)) +
  geom_bar(stat="identity",position="dodge")+ylab("Percentage of miRNAs")+
  scale_x_discrete(name="Cancer type (sample size)",breaks=cancerType1, 
                   labels=addline_format(cancerType2))+#ggtitle("Correlation between numTS and expression")+
  geom_text(aes(label=Num), position=position_dodge(width=0.9), vjust=-0.25)+theme_classic()+
  scale_fill_manual(values=c("black"))+ scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))#+ scale_y_sqrt()

data
colnames(data_Percentage)[c(1,3)] <- c("Cancer type", "Correlation between numTS and Expr.")
ggplot(data_Percentage, aes(x=sampleSize, y=Num, color=`Cancer type`, shape=`Correlation between numTS and Expr.`))+geom_point()+theme_bw()+xlab("Sample Size")+ylab("Number of miRNAs")
ggplot(data_Percentage[data_Percentage$`Correlation between numTS and Expr.` =="Correlated",], aes(x=sampleSize, y=Num, color=`Cancer type`))+geom_point(size = 3)+theme_bw()+xlab("Sample Size")+ylab("Number of miRNAs with correlation")


library(dplyr)
library(tidyr)
cnt_num <- data %>% filter(all_q<0.05) %>% group_by(group) %>% dplyr::count(all_v > 0) %>% 
  spread(`all_v > 0`, n) 
cnt_num$number <- paste0("Pos:",cnt_num$`TRUE`, " Neg:", cnt_num$`FALSE`)
library(ggrepel)
ggplot(cnt_num, aes(x=`TRUE`, y = `FALSE`, col=group, label = number))+geom_point(size=3)+theme_bw()+xlab("#Positively correlated")+ylab("#Negatively correlated")+
  geom_abline(intercept = 0, slope = 1)+geom_text_repel()#+scale_x_continuous(limits = c(0, 90))+scale_y_continuous(limits = c(0, 60))

cnt_plt <- data %>% filter(all_q<0.05) %>% group_by(group) %>% dplyr::count(all_v > 0)
colnames(cnt_plt)[2] <- "Correlation"
cnt_plt$Correlation <- ifelse(cnt_plt$Correlation, "Positive", "Negative")
cnt_plt <- rbind(cnt_plt, data.frame(group="KIRC", Correlation="Negative", n=0))
cnt_plt <- cnt_plt[order(cnt_plt$group),]
cnt_plt$sampleSize <- values(samplesize, keys=cnt_plt$group)



cnt_plt <- cnt_plt[order(cnt_plt$sampleSize, decreasing = T),]
cnt_plt$group <- factor(cnt_plt$group, levels = unique(cnt_plt$group))
cnt_plt$Correlation <- factor(cnt_plt$Correlation, levels = c("Positive", "Negative"))

cnt_num <- cnt_num[match(unique(cnt_plt$group),cnt_num$group),]
Num <- integer()
for(i in 1:10){
  Num <- c(Num, c(cnt_num$`FALSE`[i], cnt_num$`TRUE`[i]))
}
Num[is.na(Num)] <- 0

Num2<-rep(0,10)
for(i in 1:10){
  Num2[i]<-paste("(",unique(cnt_plt$sampleSize)[i],")",sep="")
}

cancerType2<-rep("0",10)
for(i in 1:10){
  cancerType2[i]<-paste(unique(cnt_plt$group)[i],Num2[i],sep=" ")
}


ggplot(data=cnt_plt, aes(x=group, y=n,fill=Correlation)) +
  geom_bar(stat="identity",position="dodge")+scale_y_continuous(name="Number of miRNAs with significant correlation")+
  scale_x_discrete(name="Cancer type (sample size)",breaks=unique(cnt_plt$group), 
                   labels=addline_format(cancerType2))+
  geom_text(aes(label=cnt_plt$n), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values=c("black","darkgrey"))+theme_classic()#+scale_fill_manual(values=c("skyblue","red"))


### Figure 1.D

setwd("C:/Users/byl12/OneDrive - University of Pittsburgh/Graduate/Pitt/ParkLab/IndependentProject/SUR3UTR/data/GDC_TCGA_FPKM/result_TC3A")
source("C:/Users/byl12/OneDrive - University of Pittsburgh/Graduate/Pitt/ParkLab/IndependentProject/SUR3UTR/GDC_FPKM/Independence/functions_logremoved.R")

par(mfrow=c(1,1))
group<-NULL
all_p<-NULL
all_q<-NULL
all_v<-NULL
all_miR <- NULL
cancerType1=c("BRCA","SKCM","HNSC","KIRC","LGG","LUAD","LUSC","OV","STAD","UCEC")

## data_option doesn't matter here, since Y we don't use Y.
Data_option<-"immune"
#Data_option<-"proliferation"
Num<-rep(0,10)
Num1<-rep(0,10)
sampleSize<-rep(0,10)
names(sampleSize) <- cancerType1
for(ind_cancer in 1:length(cancerType1)){
  cancertype=cancerType1[ind_cancer]
  data_all<-process_data(cancertype = cancertype,option=Data_option)
  df1=data_all$df1
  df1<-df1[,1:(ncol(df1)-1)]
  Num[ind_cancer]<-ncol(df1)
  df.ctrl1=data_all$df.ctrl1
  df.ctrl1=df.ctrl1[,1:(ncol(df.ctrl1)-1)]
  index1<-which(apply(df1,2,sum)==0)
  index2<-which(apply(df.ctrl1,2,sum)==0)
  index<-c(index1,index2)
  df1<-df1[,-index]
  Num1[ind_cancer]<-ncol(df1)
  df.ctrl1<-df.ctrl1[,-index]
  sampleSize[ind_cancer]<-nrow(df1)
  
  group<-c(group,rep(cancertype,dim(df1)[2]))
  set.seed(12315)
  p.cor<-rep(0,dim(df1)[2])
  v.cor<-rep(0,dim(df1)[2])
  for(i in 1:dim(df1)[2]){
    temp<-cor.test(df1[,i],df.ctrl1[,i])
    p.cor[i]<-temp$p.value
    v.cor[i]<-temp$estimate
  }
  all_miR <- c(all_miR, colnames(df1))
  all_p<-c(all_p,p.cor)
  q.temp<-rep(NA,length(p.cor))
  q.temp[which(!is.na(p.cor))]<-p.adjust(p.cor[!is.na(p.cor)],method = "BH")
  all_q<-c(all_q,q.temp)
  all_v<-c(all_v,v.cor)
}
library(hash)
samplesize <- hash(keys=names(sampleSize), values= sampleSize)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

data<-data.frame(all_miR, group,all_p,all_q,all_v)
for(i in 1:nrow(data)){
  if(data$all_q[i] > 0.05){
    data$all_v[i] <- 0
  }
}
library(tidyr)
data_wide <- spread(data[,c(1,2,5)], group, all_v)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
Heatmap(na.omit(data_wide[,-1]), col = colorRamp2(c(0.25, 0, -0.25), brewer.pal(n=3, name="RdBu")), 
        show_row_names = T, show_column_names = F)

order_miR <- data_wide$all_miR[order(rowSums(data_wide[,-1]>0, na.rm = T), rowSums(data_wide[,-1]!=0, na.rm = T),decreasing = F)]
setdiff(order_miR, data$all_miR)
setdiff(data$all_miR, order_miR)
data$all_miR <- factor(data$all_miR, levels = order_miR)
library(ggplot2)
ggplot(data, aes(y= all_miR, x=group, fill=all_v))+geom_tile()+scale_fill_gradient2(low="navy", mid="white", high="red", 
                                                                                    midpoint=0)+ theme(axis.title.y=element_blank(),axis.text.y=element_blank())+ggtitle("Pearson's Correlation Coefficient")
library(tidyr)
#data<-data.frame(all_miR, group,all_p,all_q,all_v)
data<-data.frame(all_miR, group,all_p,all_q,all_v)
data_wide <- spread(data[,c(1,2,4)], group, all_q)
order_miR <- data_wide$all_miR[order(rowSums(data_wide[,-1], na.rm = T), decreasing = T)]
setdiff(order_miR, data$all_miR)
setdiff(data$all_miR, order_miR)
data$all_miR <- factor(data$all_miR, levels = order_miR)
ggplot(data, aes(y= all_miR, x=group, fill=all_q))+geom_tile()+scale_fill_gradient2(low="navy", mid="white", high="white", 
                                                                                    midpoint=0.05)+ theme(axis.title.y=element_blank(),axis.text.y=element_blank())+ggtitle("Pearson's Correlation P value")

plt_bar <- data.frame(miR=data_wide$all_miR, CorrelatedCount=rowSums(data_wide[,-1]<0.05,na.rm=T))
plt_bar <- plt_bar %>% group_by(CorrelatedCount) %>% count()
plt_bar <- rbind(plt_bar, data.frame(CorrelatedCount=c(10,9,8), n=c(0,0,0)))
plt_bar$CorrelatedCount <- factor(plt_bar$CorrelatedCount, levels = 0:10)
ggplot(plt_bar, aes(x=CorrelatedCount, y=n))+geom_bar(stat="identity")+ylab("Number of miRs")+xlab("Number of cancer types where Expr. significantly correlated with numTS")+
  geom_text(aes(label=plt_bar$n), position=position_dodge(width=0.9), vjust=-0.25)+theme_classic()

