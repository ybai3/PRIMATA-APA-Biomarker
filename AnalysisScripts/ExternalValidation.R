library(glmnet)
library(ROCR)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

KLUAD_mRNAexpr_dir = "/zfs1/hpark/SUR3UTR_immune/new/K_LUAD"
TCGA_bindingsite_dir = "/zfs1/hpark/SUR3UTR_immune/pythons/TCGA_FPKM/TC3A/result_TC3A"
STable_dir = "/zfs1/hpark/SUR3UTR_immune/new/K_LUAD/STable7_Davoli.txt"
KLUAD_bindingsite_dir = "/zfs1/hpark/SUR3UTR_immune/new/K_LUAD"
gene4immune <- c("CD247","CD2","CD3E","GZMH","NKG7","PRF1","GZMK")
gene4Prolif <- c("CENPE","CCNA2","CCNB2","MCM6","CCNF","BUB1","CDC20","CDC6","CDK1","PLK1")
cancerType1 = c("LUAD")
num_simulation<-200

boots.data <- function(data, train.portion) {
  data <- data[sample(rownames(data),size = nrow(data), replace = T),]
  data <- apply(data, 2, as.double)
  n <- nrow(data)
  ntr <- floor(n*train.portion)
  trainid <- sample(seq(n), size = ntr)
  ntrainid <- !seq(n)%in%trainid
  x.train <- data[trainid, -ncol(data)]
  y.train <- data[trainid, ncol(data)]
  x.test <- data[ntrainid, -ncol(data)]
  y.test <- data[ntrainid, ncol(data)]
  p <- ncol(data) - 1
  list(x.train = x.train, y.train = y.train, x.test = x.test, y.test = y.test, p = p)
}

featureSelection <- function(Data, num_simulation, trainportion=0.75){
  p=ncol(Data)-1
  penalty<-rep(1,ncol(Data)-1)
  nonzero.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
  coef.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
  R_squared.case1 <- rep(NA, num_simulation)
  Rmse.case1 <- rep(NA, num_simulation)
  #source("boots.data.r")
  for(i in 1:num_simulation) {
    div.data <- boots.data(Data, trainportion)
    x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
    
    
    foldid <- sample(rep(1:10, length.out = length(y.train)))
    cv.fit.CCLE <- list()
    length(cv.fit.CCLE) <- 11
    names(cv.fit.CCLE) <- paste0("alpha_",seq(0,1,0.1))
    searchResult <- matrix(nrow = length(cv.fit.CCLE), ncol = 3)
    colnames(searchResult) <- c("MSE","lambda.min","alpha")
    for (j in 1:length(cv.fit.CCLE)){
      cv.fit.tmp <- cv.glmnet(x.train, y.train,type.measure="mse", alpha=seq(0,1,0.1)[j],penalty.factor=penalty, foldid = foldid)
      searchResult[j,1] <- cv.fit.tmp$cvm[cv.fit.tmp$lambda == cv.fit.tmp$lambda.1se]
      searchResult[j,2] <- cv.fit.tmp$lambda.1se
      searchResult[j,3] <- seq(0,1,0.1)[j]
      cv.fit.CCLE[[j]] <- cv.fit.tmp
    }
    selectedAlpha <- searchResult[searchResult[,1] == min(searchResult[,1]), ]
    cv.fit <- cv.fit.CCLE[[which(searchResult[,1] == min(searchResult[,1]))]]
    
    pred <-predict(cv.fit, newx=x.test,s="lambda.1se")
    var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
    nonzero.case1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
    coef.case1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
    R_squared.case1[i]=cor(y.test, pred)**2 
    Rmse.case1[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
  }
  MiRNA_case1_percentage<-apply(nonzero.case1,2,function(x){return(sum(x)/num_simulation)})
  MiRNA_case1_percentage<-as.data.frame(MiRNA_case1_percentage)
  rownames(MiRNA_case1_percentage)<-colnames(Data)[1:(ncol(Data)-1)]
  
  positive_case1<-apply(coef.case1,2,function(x){return(sum(x>0)/num_simulation)})
  negative_case1<-apply(coef.case1,2,function(x){return(sum(x<0)/num_simulation)})
  indicator_case1<-rep(-1,ncol(coef.case1))
  for(i in 1:length((indicator_case1))){
    if(positive_case1[i]>negative_case1[i]){
      indicator_case1[i]<-1
    }else if(positive_case1[i]<negative_case1[i]){
      indicator_case1[i]<--1
    }else{
      indicator_case1[i]<-0
    }
  }
  MiRNA_case1<-cbind(MiRNA_case1_percentage,indicator_case1)
  result_case1<-list(MIRNA=MiRNA_case1,R_square=R_squared.case1,RMSE=Rmse.case1)
  return(result_case1)
}



Data_option<-"immune"
#Data_option<-"proliferation"
result_Rsquare_case<-matrix(NA,num_simulation,length(cancerType1))
result_Rsquare_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_Rmse_case<-matrix(NA,num_simulation,length(cancerType1))
result_Rmse_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_corR_case<-matrix(NA,num_simulation,length(cancerType1))
result_corR_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_corP_case<-matrix(NA,num_simulation,length(cancerType1))
result_corP_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_corRsp_case<-matrix(NA,num_simulation,length(cancerType1))
result_corRsp_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_corPsp_case<-matrix(NA,num_simulation,length(cancerType1))
result_corPsp_ctrl<-matrix(NA,num_simulation,length(cancerType1))
result_miRNA_case<-list()
length(result_miRNA_case) <- length(cancerType1)
result_miRNA_ctrl<-list()
length(result_miRNA_ctrl) <- length(cancerType1)
for(i1 in 1:length(cancerType1)){
  print(i1)
  cancertype=cancerType1[i1]
  data_all<-process_data(cancertype = cancertype,option=Data_option)
  df1=data_all$df1
  df.ctrl1=data_all$df.ctrl1
  #############
  set.seed(12315)
  result_case1_NoPrior=case1(Data=df1,option="NoPrior",Num_simulation = num_simulation)
  result_Rsquare_case[,i1]<-result_case1_NoPrior$R_square
  result_Rmse_case[,i1]<-result_case1_NoPrior$RMSE
  result_corR_case[,i1]<-result_case1_NoPrior$cor_R
  result_corP_case[,i1]<-result_case1_NoPrior$cor_P
  result_corRsp_case[,i1]<-result_case1_NoPrior$cor_Rsp
  result_corPsp_case[,i1]<-result_case1_NoPrior$cor_Psp
  result_miRNA_case[[i1]]<-result_case1_NoPrior$MIRNA
  result_ctrl1_NoPrior=ctrl1(Data=df.ctrl1,option="NoPrior",Num_simulation = num_simulation)
  result_Rsquare_ctrl[,i1]<-result_ctrl1_NoPrior$R_square
  result_Rmse_ctrl[,i1]<-result_ctrl1_NoPrior$RMSE
  result_corR_ctrl[,i1]<-result_ctrl1_NoPrior$cor_R
  result_corP_ctrl[,i1]<-result_ctrl1_NoPrior$cor_P
  result_corRsp_ctrl[,i1]<-result_ctrl1_NoPrior$cor_Rsp
  result_corPsp_ctrl[,i1]<-result_ctrl1_NoPrior$cor_Psp
  result_miRNA_ctrl[[i1]]<-result_ctrl1_NoPrior$MIRNA
}

colnames(result_Rsquare_case)<-cancerType1
colnames(result_Rmse_case)<-cancerType1
colnames(result_corR_case)<-cancerType1
colnames(result_corP_case)<-cancerType1
colnames(result_corRsp_case)<-cancerType1
colnames(result_corPsp_case)<-cancerType1
names(result_miRNA_case)<-cancerType1
colnames(result_Rsquare_ctrl)<-cancerType1
colnames(result_Rmse_ctrl)<-cancerType1
colnames(result_corR_ctrl)<-cancerType1
colnames(result_corP_ctrl)<-cancerType1
colnames(result_corRsp_ctrl)<-cancerType1
colnames(result_corPsp_ctrl)<-cancerType1
names(result_miRNA_ctrl)<-cancerType1

Data<-list(result_miRNA_case=result_miRNA_case,result_Rsquare_case=result_Rsquare_case,result_Rmse_case=result_Rmse_case,result_corR_case=result_corR_case,result_corP_case=result_corP_case,result_corRsp_case=result_corRsp_case,result_corPsp_case=result_corPsp_case,
           result_miRNA_ctrl=result_miRNA_ctrl,result_Rsquare_ctrl=result_Rsquare_ctrl,result_Rmse_ctrl=result_Rmse_ctrl,result_corR_ctrl=result_corR_ctrl,result_corP_ctrl=result_corP_ctrl,result_corRsp_ctrl=result_corRsp_ctrl,result_corPsp_ctrl=result_corPsp_ctrl)

save(Data,file=paste0(Data_option,"_BT_LUAD_immune.Rdata"))
