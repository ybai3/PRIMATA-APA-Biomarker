## B cell abundance as an example. Same code was ran on all other cell types. 

setwd("/ihome/hpark/yub20/SUR3UTR/pythons/TCGA_FPKM/TC3A/result_TC3A")

library(glmnet)
library(ROCR)
library(stringr)
divide.data <- function(data, train.portion) {
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




auc <- function(predic, y) {
  pred <- prediction(predic, y)
  auc <- as.numeric(performance(pred, measure = "auc", x.measure = "cutoff")@y.values)
  auc
}


process_data<-function(cancertype, option){
  MATin=paste("APA4immuneProlif_",cancertype,".txt",sep="")
  df=read.delim(MATin)
  #df=read.delim(MATin)
  df<-df[,-((ncol(df)-1):ncol(df))]
  for (tmp_i in 2:ncol(df)){
    df[,tmp_i] <- log2(df[,tmp_i]+1)
  }
  
  #MATin="APA4survival_BRCA.txt"
  MATin=paste("ctrl4immuneProlif_",cancertype,".txt",sep="")
  df.ctrl=read.delim(MATin)
  
  
  df.ctrl<-df.ctrl[,-((ncol(df.ctrl)-1):ncol(df.ctrl))]
  
  infiltration_estimation_for_tcga <- read.csv("infiltration_estimation_for_tcga.csv")
  timer <- infiltration_estimation_for_tcga[,c(1:7)]
  timer <- timer[as.integer(str_sub(timer$cell_type, 14,15))<10,]
  timer$cell_type <- str_sub(timer$cell_type, 6,12)
  timer <- timer[timer$cell_type %in% df$patName,]
  timer <- timer[!duplicated(timer$cell_type),]  
  com_pat <- intersect(timer$cell_type, df$patName)
  timer <- timer[timer$cell_type %in% com_pat,]    
  df <- df[df$patName %in% com_pat,]
  df.ctrl <- df.ctrl[df.ctrl$patName %in% com_pat,]
  timer <- timer[match(df$patName, timer$cell_type),]
  df <- cbind(df, timer[,option])
  df <- df[,-1]
  df.ctrl <- cbind(df.ctrl, timer[,option])
  df.ctrl <- df.ctrl[,-1]
  return(list(df1=df,df.ctrl1=df.ctrl))
}



case1<-function(Data,option,Num_simulation,Prior_list=NULL){
  p=ncol(Data)-1
  if(option=="NoPrior"){
    penalty<-rep(1,ncol(Data)-1)
    
  }else if(option=="IncludePrior"){
    MiRNA_vector<-colnames(Data)
    MiRNA_vector[which(MiRNA_vector%in%immune_BRCA)]
    penalty<-rep(1,ncol(Data))
    penalty[which(MiRNA_vector%in%immune_BRCA)]<-0
    penalty<-penalty[-length(penalty)]
  }else{
    MiRNA_vector<-colnames(Data)
    MiRNA_vector[which(MiRNA_vector%in%immune_BRCA)]
    Data<-Data[,c(which(MiRNA_vector%in%immune_BRCA),ncol(Data))]
    p=ncol(Data)-1
    #nonzero.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
    #coef.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
    R_squared.case1 <- rep(NA, num_simulation)
    Rmse.case1 <- rep(NA, num_simulation)
    #source("divide.data.r")
    for(i in 1:num_simulation) {
      div.data <- divide.data(Data, c(0.75))
      x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
      glm_data<-data.frame(x.train,y.train)
      res.glm<-glm(y.train~.,family="gaussian",data=glm_data)
      pred <-predict(res.glm, newdata=as.data.frame(x.test))
      #var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
      #nonzero.case1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
      #coef.case1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
      
      R_squared.case1[i]=cor(y.test, pred)**2 
      Rmse.case1[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
    }
    
    result_case1<-list(R_square=R_squared.case1,RMSE=Rmse.case1)
    return(result_case1)
  }
  
  
  nonzero.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
  coef.case1 <- matrix(NA, nrow=num_simulation, ncol=p)
  R_squared.case1 <- rep(NA, num_simulation)
  cor.case1 <- rep(NA, num_simulation)
  cor.p.case1 <- rep(NA, num_simulation)
  Rmse.case1 <- rep(NA, num_simulation)
  #source("divide.data.r")
  for(i in 1:num_simulation) {
    div.data <- divide.data(Data, c(0.75))
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
    cor.tmp <- cor.test(y.test, pred, method = "pearson")
    cor.case1[i] <- cor.tmp$estimate
    cor.p.case1[i] <- cor.tmp$p.value
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
  result_case1<-list(MIRNA=MiRNA_case1,R_square=R_squared.case1,RMSE=Rmse.case1,cor_R=cor.case1,cor_P=cor.p.case1)
  return(result_case1)
}

ctrl1<-function(Data,option,Num_simulation,Prior_list=NULL){
  p=ncol(Data)-1
  if(option=="NoPrior"){
    penalty<-rep(1,ncol(Data)-1)
    
  }else if(option=="IncludePrior"){
    MiRNA_vector<-colnames(Data)
    MiRNA_vector[which(MiRNA_vector%in%immune_BRCA)]
    penalty<-rep(1,ncol(Data))
    penalty[which(MiRNA_vector%in%immune_BRCA)]<-0
    penalty<-penalty[-length(penalty)]
  }else{
    MiRNA_vector<-colnames(Data)
    MiRNA_vector[which(MiRNA_vector%in%immune_BRCA)]
    Data<-Data[,c(which(MiRNA_vector%in%immune_BRCA),ncol(Data))]
    p=ncol(Data)-1
    #nonzero.ctrl1 <- matrix(NA, nrow=num_simulation, ncol=p)
    #coef.ctrl1 <- matrix(NA, nrow=num_simulation, ncol=p)
    R_squared.ctrl1 <- rep(NA, num_simulation)
    Rmse.ctrl1 <- rep(NA, num_simulation)
    #source("divide.data.r")
    for(i in 1:num_simulation) {
      div.data <- divide.data(Data, c(0.75))
      x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
      glm_data<-data.frame(x.train,y.train)
      res.glm<-glm(y.train~.,family="gaussian",data=glm_data)
      pred <-predict(res.glm, newdata=as.data.frame(x.test))

      
      R_squared.ctrl1[i]=cor(y.test, pred)**2 
      Rmse.ctrl1[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
    }
    
    result_ctrl1<-list(R_square=R_squared.ctrl1,RMSE=Rmse.ctrl1)
    return(result_ctrl1)
  }
  
  
  nonzero.ctrl1 <- matrix(NA, nrow=num_simulation, ncol=p)
  coef.ctrl1 <- matrix(NA, nrow=num_simulation, ncol=p)
  R_squared.ctrl1 <- rep(NA, num_simulation)
  Rmse.ctrl1 <- rep(NA, num_simulation)
  cor.ctrl1 <- rep(NA, num_simulation)
  cor.p.ctrl1 <- rep(NA, Num_simulation)
  #source("divide.data.r")
  for(i in 1:num_simulation) {
    div.data <- divide.data(Data, c(0.75))
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
    nonzero.ctrl1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
    coef.ctrl1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
    cor.tmp<-cor.test(y.test, pred, method = "pearson")
    cor.ctrl1[i] <- cor.tmp$estimate
    cor.p.ctrl1[i] <- cor.tmp$p.value
    R_squared.ctrl1[i]=cor(y.test, pred)**2 
    Rmse.ctrl1[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
  }
  MiRNA_ctrl1_percentage<-apply(nonzero.ctrl1,2,function(x){return(sum(x)/num_simulation)})
  MiRNA_ctrl1_percentage<-as.data.frame(MiRNA_ctrl1_percentage)
  rownames(MiRNA_ctrl1_percentage)<-colnames(Data)[1:(ncol(Data)-1)]
  
  positive_ctrl1<-apply(coef.ctrl1,2,function(x){return(sum(x>0)/num_simulation)})
  negative_ctrl1<-apply(coef.ctrl1,2,function(x){return(sum(x<0)/num_simulation)})
  indicator_ctrl1<-rep(-1,ncol(coef.ctrl1))
  for(i in 1:length((indicator_ctrl1))){
    if(positive_ctrl1[i]>negative_ctrl1[i]){
      indicator_ctrl1[i]<-1
    }else if(positive_ctrl1[i]<negative_ctrl1[i]){
      indicator_ctrl1[i]<--1
    }else{
      indicator_ctrl1[i]<-0
    }
  }
  MiRNA_ctrl1<-cbind(MiRNA_ctrl1_percentage,indicator_ctrl1)
  result_ctrl1<-list(MIRNA=MiRNA_ctrl1,R_square=R_squared.ctrl1,RMSE=Rmse.ctrl1,cor_R=cor.ctrl1,cor_P=cor.p.ctrl1)
  return(result_ctrl1)
}



num_simulation<-300
cancerType1=c("BRCA","SKCM","LUAD","HNSC","KIRC","LGG","LUSC","OV","STAD","UCEC")
Data_option<-"B.cell_TIMER"
#Data_option<-"proliferation"
result_Rsquare_case<-matrix(NA,num_simulation,10)
result_Rsquare_ctrl<-matrix(NA,num_simulation,10)
result_Rmse_case<-matrix(NA,num_simulation,10)
result_Rmse_ctrl<-matrix(NA,num_simulation,10)
result_corR_case<-matrix(NA,num_simulation,10)
result_corR_ctrl<-matrix(NA,num_simulation,10)
result_corP_case<-matrix(NA,num_simulation,10)
result_corP_ctrl<-matrix(NA,num_simulation,10)
for(i1 in 1:10){
  print(i1)
  cancertype=cancerType1[i1]
  data_all<-process_data(cancertype = cancertype,option=Data_option)
  df1=data_all$df1
  df.ctrl1=data_all$df.ctrl1
  #############
  set.seed(220916)
  result_case1_NoPrior=case1(Data=df1,option="NoPrior",Num_simulation = num_simulation)
  result_Rsquare_case[,i1]<-result_case1_NoPrior$R_square
  result_Rmse_case[,i1]<-result_case1_NoPrior$RMSE
  result_corR_case[,i1]<-result_case1_NoPrior$cor_R
  result_corP_case[,i1]<-result_case1_NoPrior$cor_P
  result_ctrl1_NoPrior=ctrl1(Data=df.ctrl1,option="NoPrior",Num_simulation = num_simulation)
  result_Rsquare_ctrl[,i1]<-result_ctrl1_NoPrior$R_square
  result_Rmse_ctrl[,i1]<-result_ctrl1_NoPrior$RMSE
  result_corR_ctrl[,i1]<-result_ctrl1_NoPrior$cor_R
  result_corP_ctrl[,i1]<-result_ctrl1_NoPrior$cor_P
}

colnames(result_Rsquare_case)<-cancerType1
colnames(result_Rmse_case)<-cancerType1
colnames(result_corR_case)<-cancerType1
colnames(result_corP_case)<-cancerType1
colnames(result_Rsquare_ctrl)<-cancerType1
colnames(result_Rmse_ctrl)<-cancerType1
colnames(result_corR_ctrl)<-cancerType1
colnames(result_corP_ctrl)<-cancerType1

Data<-list(result_Rsquare_case=result_Rsquare_case,result_Rmse_case=result_Rmse_case,result_corR_case=result_corR_case,result_corP_case=result_corP_case,
           result_Rsquare_ctrl=result_Rsquare_ctrl,result_Rmse_ctrl=result_Rmse_ctrl,result_corR_ctrl=result_corR_ctrl,result_corP_ctrl=result_corP_ctrl)

setwd("/ihome/hpark/yub20/SUR3UTR/GDC_TCGA_FPKM_analysis/timer_predict/eachTypeAlone/rep300")
save(Data,file=paste0(Data_option,"_result_list.RData"))
