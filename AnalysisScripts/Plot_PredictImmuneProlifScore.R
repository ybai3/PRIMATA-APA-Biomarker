setwd("/ihome/hpark/yub20/SUR3UTR/GDC_TCGA_FPKM_analysis/logTransNumTS")
#source("D:/OneDrive - University of Pittsburgh/Graduate/Pitt/ParkLab/IndependentProject/SUR3UTR/YJPart//prediction comparision between FPM and binding sites/functions_YLupdate_gridSearch_removeLogtransform.R")
library(ggplot2)
num_simulation<-300
cancerType1=c("HNSC","KIRC","LUAD","OV","STAD","UCEC",
              "LGG","BRCA","LUSC","SKCM")

group_function<-function(x,num){
  res<-NULL
  for(i in 1:length(x)){
    res<-c(res,rep(x[i],num))
  }
  return(res)
}


########################################RMSE
#load("immune_rep300_result_list.Rdata")
load("proliferation_rep300_result_list.Rdata")
#Data_option<-"immune"
Data_option<-"proliferation"
#cancer<-group_function(cancerType1,num_simulation)
result_Rmse_case<-Data$result_Rmse_case
result_Rmse_ctrl<-Data$result_Rmse_ctrl
cancer<-group_function(cancerType1,num_simulation)
type<-c(rep("#target site",num_simulation*length(cancerType1)),rep("expression",num_simulation*length(cancerType1)))
color<-c(rep("#target site",num_simulation*length(cancerType1)),rep("expression",num_simulation*length(cancerType1)))
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
Data_R$color <- ifelse(Data_R$color == "#target site", "numTS", "expression")
Data_R$type <- ifelse(Data_R$type == "#target site", "numTS", "expression")
Data_R$type <- factor(Data_R$type, levels=c("numTS", "expression"))
library(tidyr)
library(dplyr)
meanTab <- Data_R %>% group_by(cancer, type) %>% summarise(
  value = mean(RMSE)
)
colnames(meanTab)[3] <- "mean"
library(tidyr)
meanTab_wide <- spread(meanTab, type, mean)
mean((meanTab_wide$expression-meanTab_wide$numTS)/meanTab_wide$numTS)
Data_R$cancer <- factor(Data_R$cancer, levels = c("BRCA", "LGG", "OV", "LUAD", "UCEC", "HNSC", "SKCM", "KIRC", "STAD", "LUSC"))
ggplot(data = Data_R, aes(x = type, y = RMSE)) + geom_boxplot(aes(fill=type),width=0.6)+
  #geom_text(data=meanTab, aes(label = round(mean,2), x=type, y=mean))+ 
  facet_wrap(~cancer, nrow=2)+
  labs(x=NULL,y="RMSE")+theme_bw()+theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.background = element_blank(),
                                         panel.border = element_rect(colour = "black"))+
  scale_fill_manual(values = c("red","blue"))
first5<- ggplot(data = Data_R[Data_R$cancer %in% c("BRCA", "LGG", "OV", "LUAD", "UCEC"),], aes(x = type, y = RMSE)) + geom_boxplot(aes(fill=type),width=0.6)+
  #geom_text(data=meanTab, aes(label = round(mean,2), x=type, y=mean))+ 
  facet_wrap(~cancer, nrow=1)+
  labs(x=NULL,y="RMSE")+theme_bw()+theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.background = element_blank(),
                                         panel.border = element_rect(colour = "black"))+
  scale_fill_manual(values = c("red","blue"))
second5<- ggplot(data = Data_R[Data_R$cancer %in% c("HNSC", "SKCM", "KIRC", "STAD", "LUSC"),], aes(x = type, y = RMSE)) + geom_boxplot(aes(fill=type),width=0.6)+
  #geom_text(data=meanTab, aes(label = round(mean,2), x=type, y=mean))+ 
  facet_wrap(~cancer, nrow=1)+
  labs(x=NULL,y="RMSE")+theme_bw()+theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.background = element_blank(),
                                         panel.border = element_rect(colour = "black"))+
  scale_fill_manual(values = c("red","blue"))
library(gridExtra)
library(grid)
#ggarrange(list(first5, second5), nrow = 2)
grid.arrange(first5, second5, nrow = 2)

plt_list <- vector("list", length = length(cancerType1))
names(plt_list) <- cancerType1
for(i in 1:length(cancerType1)){
  plt_list[[i]] <- ggplot(data = Data_R[Data_R$cancer==cancerType1[i],], aes(x = type, y = RMSE,fill=type)) + geom_boxplot(width=0.6)+ 
    labs(x=NULL,y="RMSE")+ theme_bw()+theme(legend.position = "none",
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            strip.background = element_blank(),
                                            panel.border = element_rect(colour = "black"))+scale_fill_manual(values=c("red","blue"))
}
library(ggpubr)
ggarrange(plotlist = plt_list, nrow=2, ncol=5)