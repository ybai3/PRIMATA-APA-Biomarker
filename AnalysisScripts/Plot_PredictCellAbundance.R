library(ggplot2)
library(stringr)

group_function<-function(x,num){
  res<-NULL
  for(i in 1:length(x)){
    res<-c(res,rep(x[i],num))
  }
  return(res)
}


########################################RMSE

setwd("/ihome/hpark/yub20/SUR3UTR/GDC_TCGA_FPKM_analysis/timer_predict/eachTypeAlone/rep300")
cancerType1=c("BRCA","SKCM","LUAD","HNSC","KIRC","LGG","LUSC","OV","STAD","UCEC")
num_simulation<-300
dir_name <- list.files()
dir_name <- dir_name[str_detect(dir_name, ".RData")]
cell_type <- str_remove(dir_name, "_TIMER_result_list.RData")

plt_tab <- data.frame(RMSE=NA, cancer=NA, type=NA, color=NA, celltype=NA)
#colnames(plt_tab) <- c("RMSE", "cancer", "type","color","celltype")
for (i in 1:length(dir_name)){
  load(dir_name[i])
  result_Rmse_case<-Data$result_Rmse_case
  result_Rmse_ctrl<-Data$result_Rmse_ctrl
  cancer<-group_function(cancerType1,num_simulation)
  type<-c(rep("numTS",num_simulation*length(cancerType1)),rep("expr",num_simulation*length(cancerType1)))
  color<-c(rep("numTS",num_simulation*length(cancerType1)),rep("expr",num_simulation*length(cancerType1)))
  data_Rcase<-cbind(as.vector(result_Rmse_case),cancer)
  data_Rctrl<-cbind(as.vector(result_Rmse_ctrl),cancer)
  Data_R<-rbind(data_Rcase,data_Rctrl)
  Data_R<-cbind(Data_R,type)
  Data_R<-cbind(Data_R,color)
  colnames(Data_R)<-c("RMSE","cancer","type","color")
  Data_R<-as.data.frame(Data_R)
  class(Data_R$value)
  class(Data_R$cancer)
  class(Data_R$type)
  class(Data_R$color)
  Data_R$RMSE<-as.numeric(as.character(Data_R$RMSE))
  
  Data_R$color<-as.character(Data_R$color)
  Data_R$celltype <- rep(cell_type[i], nrow(Data_R))
  plt_tab <- rbind(plt_tab, Data_R)
}
plt_tab <- plt_tab[-1,]
plt_tab$cancer <- factor(plt_tab$cancer, levels = c("BRCA","LGG", "OV", "LUAD", "UCEC", "HNSC", "SKCM", "KIRC", "STAD", "LUSC"))
plt_tab$type <- factor(plt_tab$type, levels = c("numTS", "expr"))
ggplot(data = plt_tab, aes(x = cancer, y = RMSE,fill=color)) + geom_boxplot()+ facet_wrap(~celltype)+
  labs(x=NULL,y="RMSE")+theme_bw()+theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.background = element_blank(),
                                         panel.border = element_rect(colour = "black")) + theme(axis.text.x = element_text(angle = 45,  hjust = 1))+theme_classic()
plt_tab.tmp <- plt_tab[plt_tab$celltype=="T.cell.CD8.",] 
ggplot(data = plt_tab.tmp, aes(x = cancer, y = RMSE, fill=type)) + geom_boxplot(alpha=0.8)+
  labs(x=NULL,y="RMSE")+theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1))+theme_classic()+scale_fill_manual(values=c("red","blue"))

library(dplyr)
library(ggpubr)
forest_tab <- plt_tab %>% 
  group_by(cancer, type, celltype) %>% 
  summarize(avg=mean(RMSE), sd=sd(RMSE), median=median(RMSE),low=mean(RMSE)-sd(RMSE), high=mean(RMSE)+sd(RMSE))

ggplot(forest_tab,aes(x = avg, y = cancer, xmin = low, xmax = high, group=type,color=type)) +
  #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
  geom_pointrange(position = position_dodge(width = 1),fatten = .5, size=1) +
  stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
  #geom_vline(xintercept = mean(median), linetype = 3)+ 
  coord_flip() + facet_wrap(~celltype) +theme_classic()

forest_list <- vector(mode = "list", length = 6)
for(i in 1:length(unique(forest_tab$celltype))){
  celltype.tmp <- unique(forest_tab$celltype)[i]
  plt_tab.tmp <- plt_tab[plt_tab$celltype==celltype.tmp,] 
  forest_tab.tmp <- forest_tab[forest_tab$celltype==celltype.tmp,]
  yintercept <- forest_tab.tmp %>% group_by(type) %>% summarise(vintercept=mean(median))
  #yintercept <- plt_tab.tmp %>% group_by(type) %>% summarise(vintercept=median(RMSE))
  test.tmp <- wilcox.test(plt_tab.tmp$RMSE[plt_tab.tmp$type=="numTS"],plt_tab.tmp$RMSE[plt_tab.tmp$type=="expr"])
  forest_list[[i]] <- ggplot(forest_tab.tmp,aes(x = avg, y = cancer, xmin = low, xmax = high, group=type,color=type)) +
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
    geom_pointrange(position = position_dodge(width = 0.5),fatten = 3, size=1) +
    scale_color_manual(values = c("blue", "red")) +
    #stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
    geom_vline(xintercept = yintercept$vintercept[1], linetype = 2, color="blue")+ 
    geom_vline(xintercept = yintercept$vintercept[2], linetype = 2, color="red")+ 
    coord_flip()+theme_classic()+xlab("RMSE")+ggtitle(celltype.tmp)
  if(F){
    forest_tab.tmp$id <- paste0(forest_tab.tmp$cancer,forest_tab.tmp$type,forest_tab.tmp$celltype)
    plt_tab.tmp$id <- paste0(plt_tab.tmp$cancer,plt_tab.tmp$type,plt_tab.tmp$celltype)
    join_tab.tmp <- plt_tab.tmp %>% left_join(forest_tab.tmp, by="id")
    a <- ggplot(forest_tab.tmp,aes(x = avg, y = cancer, xmin = low, xmax = high, group=type,color=type)) +
      #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
      geom_pointrange(position = position_dodge(width = 0.5),fatten = 3, size=1) +
      scale_color_manual(values = c("blue", "red")) +
      #stat_compare_means(aes(x= RMSE, group=cancer),comparisons = list(c("numTS","expr")),method = "wilcox.test")+
      geom_vline(xintercept = yintercept$vintercept[1], linetype = 2, color="blue")+ 
      geom_vline(xintercept = yintercept$vintercept[2], linetype = 2, color="red")+ 
      coord_flip()+theme_classic()+xlab("RMSE")+ggtitle(celltype.tmp)
    stat.test <- compare_means(RMSE ~ type, data = plt_tab.tmp, group.by="cancer", method="wilcox.test")
    stat.test <- stat.test %>% mutate(y.position = rep(max(celltypeCount$percent)*1.2),41)
    p + theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.7))
    a+stat_pvalue_manual(stat.test,label="p.format",size=2)
    
    
    ## merge
    ggplot(join_tab.tmp, aes(y=cancer.x)) +
      #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
      geom_pointrange(aes(x = avg, xmin = low, xmax = high, group=type,color=type),position = position_dodge(width = 0.5),fatten = 3, size=1) +
      scale_color_manual(values = c("blue", "red")) +
      stat_compare_means(aes(x=RMSE, group=type),comparisons = list(c("numTS","expr")),method = "wilcox.test")+
      geom_vline(xintercept = yintercept$vintercept[1], linetype = 2, color="blue")+ 
      geom_vline(xintercept = yintercept$vintercept[2], linetype = 2, color="red")+ 
      coord_flip()+theme_classic()+xlab("RMSE")+ggtitle(celltype.tmp)
  }
  
}


forest_tab <- plt_tab %>% 
  group_by(cancer, type, celltype) %>% 
  summarize(avg=mean(RMSE), runs= n(),sd=sd(RMSE), se=sd(RMSE)/sqrt(n()), median=median(RMSE),low=mean(RMSE)-2*(sd(RMSE)/sqrt(n())), high=mean(RMSE)+2*(sd(RMSE)/sqrt(n())))

forest_list <- vector(mode = "list", length = 10)
for(i in 1:length(unique(forest_tab$cancer))){
  cancertype.tmp <- unique(forest_tab$cancer)[i]
  plt_tab.tmp <- plt_tab[plt_tab$cancer==cancertype.tmp,] 
  forest_tab.tmp <- forest_tab[forest_tab$cancer==cancertype.tmp,]
  forest_tab.tmp$celltype <- factor(forest_tab.tmp$celltype, levels = c("T.cell.CD4.", "T.cell.CD8.", "B.cell", "Neutrophil","Macrophage","Myeloid.dendritic.cell"))
  axis_stick <- integer(length = nrow(forest_tab.tmp))
  axis_label <- character(length = nrow(forest_tab.tmp))
  for(j in 1:nrow(forest_tab.tmp)){
    if(forest_tab.tmp$celltype[j] == "T.cell.CD4."){
      axis_stick[j] <- 2
      axis_label[j] <- "CD4 T cell"
    }else if(forest_tab.tmp$celltype[j] == "T.cell.CD8."){
      axis_stick[j] <- 4
      axis_label[j] <- "CD8 T cell"
    }else if(forest_tab.tmp$celltype[j] == "B.cell"){
      axis_stick[j] <- 6
      axis_label[j] <- "B cell"
    }else if(forest_tab.tmp$celltype[j] == "Neutrophil"){
      axis_stick[j] <- 8
      axis_label[j] <- "Neutrophil"
    }else if(forest_tab.tmp$celltype[j] == "Macrophage"){
      axis_stick[j] <- 10
      axis_label[j] <- "Macrophage"
    }else if(forest_tab.tmp$celltype[j] == "Myeloid.dendritic.cell"){
      axis_stick[j] <- 12
      axis_label[j] <- "Myeloid DC"
    }
  }
  forest_tab.tmp$axis_stick <- axis_stick
  forest_tab.tmp$axis_label <- axis_label
  forest_tab.tmp$type <- as.character(forest_tab.tmp$type)
  forest_tab.tmp$type <- ifelse(forest_tab.tmp$type == "numTS", "miRNA numTS", "miRNA expression")
  #yintercept <- forest_tab.tmp %>% group_by(type) %>% summarise(vintercept=mean(median))
  yintercept <- plt_tab.tmp %>% group_by(type) %>% summarise(vintercept=median(RMSE))
  #test.tmp <- wilcox.test(plt_tab.tmp$RMSE[plt_tab.tmp$type=="numTS"],plt_tab.tmp$RMSE[plt_tab.tmp$type=="expr"])
  # forest plot
  if(F){
    forest_list[[i]] <- ggplot(forest_tab.tmp) + 
      geom_rect(data=data.frame(ymin = c(3,7,11), ymax = c(5,9,13), xmin = -Inf, xmax = Inf),
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5) +
      #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
      geom_pointrange(aes(x = avg, y = axis_stick, xmin = low, xmax = high, group=type,color=type),position = position_dodge(width = 0.8),fatten = 1.2, size=1) + 
      scale_y_continuous(breaks = axis_stick, labels=forest_tab.tmp$axis_label) +
      scale_color_manual(values = c("blue", "red")) +
      #stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
      geom_vline(xintercept = yintercept$vintercept[1], linetype = 2, color="blue")+ 
      geom_vline(xintercept = yintercept$vintercept[2], linetype = 2, color="red")+ 
      coord_flip()+theme_classic()+xlab("RMSE")+ggtitle(cancertype.tmp)+ylab("Cell type")+theme(axis.text=element_text(size=14, color="black"),
                                                                                                axis.title=element_text(size=16,color="black"),
                                                                                                plot.title = element_text(size=16,color="black"))
    
  }
  forest_tab.tmp$type <- factor(forest_tab.tmp$type, levels = c("miRNA numTS", "miRNA expression"))
  forest_list[[i]] <- ggplot(forest_tab.tmp) + 
    geom_rect(data=data.frame(xmin = c(3,7,11), xmax = c(5,9,13), ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5) +
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
    #geom_pointrange(aes(x = avg, y = axis_stick, xmin = low, xmax = high, group=type,color=type),position = position_dodge(width = 0.8),fatten = 1.2, size=1) + 
    geom_point(aes(x = axis_stick, y=avg, group=type,color=type),position=position_dodge(0.8),size=0.8) + 
    geom_errorbar(aes(x = axis_stick, y=avg, ymin=low, ymax=high, group=type,color=type), width=.2,position=position_dodge(0.8))+
    scale_x_continuous(breaks = axis_stick, labels=forest_tab.tmp$axis_label) +
    scale_color_manual(values = c("red", "blue")) +
    #stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
    geom_hline(yintercept = yintercept$vintercept[1], linetype = 2, color="red")+ 
    geom_hline(yintercept = yintercept$vintercept[2], linetype = 2, color="blue")+ labs(color='Model type') +
    theme_classic()+ylab("RMSE")+ggtitle(cancertype.tmp)+xlab("Cell type")+theme(axis.text=element_text(size=14, color="black"),
                                                                                 axis.title=element_text(size=16,color="black"),
                                                                                 plot.title = element_text(size=16,color="black"),
                                                                                 legend.text=element_text(size=12,color="black"),
                                                                                 legend.title=element_text(size=14,color="black"))
}

library(gridExtra)
library(grid)
grid.arrange(grobs=forest_list, nrow=5)

forest_list <- vector(mode = "list", length = 6)
for(i in 1:length(unique(forest_tab$celltype))){
  celltype.tmp <- c("T.cell.CD4.", "T.cell.CD8.", "B.cell", "Neutrophil","Macrophage","Myeloid.dendritic.cell")[i]
  plt_tab.tmp <- plt_tab[plt_tab$celltype==celltype.tmp,] 
  forest_tab.tmp <- forest_tab[forest_tab$celltype==celltype.tmp,]
  #forest_tab.tmp$celltype <- factor(forest_tab.tmp$celltype, levels = c("T.cell.CD4.", "T.cell.CD8.", "B.cell", "Neutrophil","Macrophage","Myeloid.dendritic.cell"))
  axis_stick <- integer(length = nrow(forest_tab.tmp))
  axis_label <- character(length = nrow(forest_tab.tmp))
  celltype.tmp <- c("CD4 T cell.", "CD8 T cell.", "B cell", "Neutrophil","Macrophage","Myeloid DC")[i]
  for(j in 1:nrow(forest_tab.tmp)){
    if(forest_tab.tmp$cancer[j] == "BRCA"){
      axis_stick[j] <- 2
      axis_label[j] <- "BRCA"
    }else if(forest_tab.tmp$cancer[j] == "LGG"){
      axis_stick[j] <- 4
      axis_label[j] <- "LGG"
    }else if(forest_tab.tmp$cancer[j] == "OV"){
      axis_stick[j] <- 6
      axis_label[j] <- "OV"
    }else if(forest_tab.tmp$cancer[j] == "LUAD"){
      axis_stick[j] <- 8
      axis_label[j] <- "LUAD"
    }else if(forest_tab.tmp$cancer[j] == "UCEC"){
      axis_stick[j] <- 10
      axis_label[j] <- "UCEC"
    }else if(forest_tab.tmp$cancer[j] == "HNSC"){
      axis_stick[j] <- 12
      axis_label[j] <- "HNSC"
    }else if(forest_tab.tmp$cancer[j] == "SKCM"){
      axis_stick[j] <- 14
      axis_label[j] <- "SKCM"
    }else if(forest_tab.tmp$cancer[j] == "KIRC"){
      axis_stick[j] <- 16
      axis_label[j] <- "KIRC"
    }else if(forest_tab.tmp$cancer[j] == "STAD"){
      axis_stick[j] <- 18
      axis_label[j] <- "STAD"
    }else if(forest_tab.tmp$cancer[j] == "LUSC"){
      axis_stick[j] <- 20
      axis_label[j] <- "LUSC"
    }
  }
  forest_tab.tmp$axis_stick <- axis_stick
  forest_tab.tmp$axis_label <- axis_label
  forest_tab.tmp$type <- as.character(forest_tab.tmp$type)
  forest_tab.tmp$type <- ifelse(forest_tab.tmp$type == "numTS", "miRNA numTS", "miRNA expression")
  #yintercept <- forest_tab.tmp %>% group_by(type) %>% summarise(vintercept=mean(median))
  yintercept <- plt_tab.tmp %>% group_by(type) %>% summarise(vintercept=median(RMSE))
  #test.tmp <- wilcox.test(plt_tab.tmp$RMSE[plt_tab.tmp$type=="numTS"],plt_tab.tmp$RMSE[plt_tab.tmp$type=="expr"])
  # forest plot
  if(F){
    forest_list[[i]] <- ggplot(forest_tab.tmp) + 
      geom_rect(data=data.frame(ymin = c(3,7,11), ymax = c(5,9,13), xmin = -Inf, xmax = Inf),
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5) +
      #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
      geom_pointrange(aes(x = avg, y = axis_stick, xmin = low, xmax = high, group=type,color=type),position = position_dodge(width = 0.8),fatten = 1.2, size=1) + 
      scale_y_continuous(breaks = axis_stick, labels=forest_tab.tmp$axis_label) +
      scale_color_manual(values = c("blue", "red")) +
      #stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
      geom_vline(xintercept = yintercept$vintercept[1], linetype = 2, color="blue")+ 
      geom_vline(xintercept = yintercept$vintercept[2], linetype = 2, color="red")+ 
      coord_flip()+theme_classic()+xlab("RMSE")+ggtitle(celltype.tmp)+ylab("Cell type")+theme(axis.text=element_text(size=14, color="black"),
                                                                                              axis.title=element_text(size=16,color="black"),
                                                                                              plot.title = element_text(size=16,color="black"))
    
  }
  forest_tab.tmp$type <- factor(forest_tab.tmp$type, levels = c("miRNA numTS", "miRNA expression"))
  forest_list[[i]] <- ggplot(forest_tab.tmp) + 
    geom_rect(data=data.frame(xmin = c(3,7,11,15,19), xmax = c(5,9,13,17,21), ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5) +
    #geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
    #geom_pointrange(aes(x = avg, y = axis_stick, xmin = low, xmax = high, group=type,color=type),position = position_dodge(width = 0.8),fatten = 1.2, size=1) + 
    geom_point(aes(x = axis_stick, y=avg, group=type,color=type),position=position_dodge(0.8),size=0.8) + 
    geom_errorbar(aes(x = axis_stick, y=avg, ymin=low, ymax=high, group=type,color=type), width=.2,position=position_dodge(0.8))+
    scale_x_continuous(breaks = axis_stick, labels=forest_tab.tmp$axis_label) +
    scale_color_manual(values = c("red", "blue")) +
    #stat_compare_means(comparisons = list(c("numTS","expr")),method = "wilcox.test")+
    geom_hline(yintercept = yintercept$vintercept[1], linetype = 2, color="red")+ 
    geom_hline(yintercept = yintercept$vintercept[2], linetype = 2, color="blue")+ labs(color='Model type') +
    theme_classic()+ylab("RMSE")+ggtitle(celltype.tmp)+xlab("Cell type")+theme(axis.text=element_text(size=14, color="black"),
                                                                               axis.title=element_text(size=16,color="black"),
                                                                               plot.title = element_text(size=16,color="black"),
                                                                               legend.text=element_text(size=12,color="black"),
                                                                               legend.title=element_text(size=14,color="black"))
}

grid.arrange(grobs=forest_list, nrow=3)

