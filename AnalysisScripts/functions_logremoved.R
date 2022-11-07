
library(glmnet)
library(ROCR)

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


process_data<-function(cancertype,option){
  if(option=="immune"){
    MATin=paste("APA4immuneProlif_",cancertype,".txt",sep="")
    df=read.delim(MATin)[,-1]
    #df=read.delim(MATin)
    df<-df[,1:(ncol(df)-1)]
    
    #MATin="APA4survival_BRCA.txt"
    MATin=paste("ctrl4immuneProlif_",cancertype,".txt",sep="")
    df.ctrl=read.delim(MATin)[,-1]
    
    
    df.ctrl<-df.ctrl[,1:(ncol(df.ctrl)-1)]
    df.ctrl1<-df.ctrl
    
    ##################case I, take log to (sum of MiRs across or the genes)
    df1<-df
    for (i in 1:ncol(df1)){
      df1[,i] <- log2(df1[,i] + 1)
    }
    ###################
    #######################case III, log(FPM) and log(MiRs) separate coefficient
    df3<-cbind(df1[,1:(ncol(df1)-1)],df.ctrl1[,1:(ncol(df.ctrl1)-1)],df1[,ncol(df1)])
  }else{
    MATin=paste("APA4immuneProlif_",cancertype,".txt",sep="")
    df=read.delim(MATin)[,-1]
    #df=read.delim(MATin)
    df<-df[,-(ncol(df)-1)]
    
    #MATin="APA4survival_BRCA.txt"
    MATin=paste("ctrl4immuneProlif_",cancertype,".txt",sep="")
    df.ctrl=read.delim(MATin)[,-1]
    
    
    df.ctrl<-df.ctrl[,-(ncol(df.ctrl)-1)]
    df.ctrl1<-df.ctrl
    
    df1<-df
    #######################case III, log(FPM) and log(MiRs) separate coefficient
    df3<-cbind(df1[,1:(ncol(df1)-1)],df.ctrl1[,1:(ncol(df.ctrl1)-1)],df1[,ncol(df1)])
  }
  return(list(df1=df1,df3=df3,df.ctrl1=df.ctrl1))
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
      #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
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
  Rmse.case1 <- rep(NA, num_simulation)
  #source("divide.data.r")
  for(i in 1:num_simulation) {
    div.data <- divide.data(Data, c(0.75))
    x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
    
    cv.fit <- cv.glmnet(x.train, y.train,type.measure="mse", alpha=0.5,penalty.factor=penalty)
    id<-which(cv.fit$lambda == cv.fit$lambda.min)
    #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
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


case3=function(Data,option,Num_simulation,Prior_list=NULL){
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
    R_squared.case3 <- rep(NA, num_simulation)
    Rmse.case3 <- rep(NA, num_simulation)
    source("divide.data.r")
    for(i in 1:num_simulation) {
      div.data <- divide.data(Data, c(0.75))
      x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
      glm_data<-data.frame(x.train,y.train)
      res.glm<-glm(y.train~.,family="gaussian",data=glm_data)
      #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
      pred <-predict(res.glm, newdata=as.data.frame(x.test))
      #var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
      #nonzero.case1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
      #coef.case1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
      
      R_squared.case3[i]=cor(y.test, pred)**2 
      Rmse.case3[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
    }
    
    result_case3<-list(R_square=R_squared.case3,RMSE=Rmse.case3)
    return(result_case3)
  }
  
  p1<-ncol(Data)-1
  nonzero.case3 <- matrix(NA, nrow=num_simulation, ncol=p1)
  coef.case3 <- matrix(NA, nrow=num_simulation, ncol=p1)
  R_squared.case3 <- rep(NA, num_simulation)
  Rmse.case3 <- rep(NA, num_simulation)
  source("divide.data.r")
  for(i in 1:num_simulation) {
    div.data <- divide.data(Data, c(0.75))
    x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
    
    cv.fit <- cv.glmnet(x.train, y.train,type.measure="mse", alpha=0.5,penalty.factor=penalty)
    id<-which(cv.fit$lambda == cv.fit$lambda.min)
    #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
    pred <-predict(cv.fit, newx=x.test,s="lambda.1se")
    var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
    nonzero.case3[i,]<-as.numeric(1:p1 %in% as.numeric(var_sel[[1]]))
    coef.case3[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
    R_squared.case3[i]=cor(y.test, pred)**2 
    Rmse.case3[i] <- sqrt(sum((y.test-pred)^2)/length(y.test))
  }
  MiRNA_case3_percentage_bind<-apply(nonzero.case3[,1:(ncol(nonzero.case3)/2)],2,function(x){return(sum(x)/num_simulation)})
  MiRNA_case3_percentage_bind<-as.data.frame(MiRNA_case3_percentage_bind)
  rownames(MiRNA_case3_percentage_bind)<-colnames(Data)[1:(ncol(Data)-1)]
  MiRNA_case3_percentage_FPM<-apply(nonzero.case3[,(ncol(nonzero.case3)/2+1):ncol(nonzero.case3)],2,function(x){return(sum(x)/num_simulation)})
  MiRNA_case3_percentage_FPM<-as.data.frame(MiRNA_case3_percentage_FPM)
  rownames(MiRNA_case3_percentage_FPM)<-colnames(Data)[1:(ncol(Data)-1)]
  MIRNA_case3<-data.frame(bind=MiRNA_case3_percentage_bind,FPM=MiRNA_case3_percentage_FPM)
  result_case3<-list(MIRNA=MIRNA_case3,R_squared=R_squared.case3,RMSE=Rmse.case3)
  return(result_case3)
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
      #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
      pred <-predict(res.glm, newdata=as.data.frame(x.test))
      #var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
      #nonzero.ctrl1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
      #coef.ctrl1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
      
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
  #source("divide.data.r")
  for(i in 1:num_simulation) {
    div.data <- divide.data(Data, c(0.75))
    x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
    
    cv.fit <- cv.glmnet(x.train, y.train,type.measure="mse", alpha=0.5,penalty.factor=penalty)
    id<-which(cv.fit$lambda == cv.fit$lambda.min)
    #pred <-predict(cv.fit, newx=x.test,s="lambda.min")
    pred <-predict(cv.fit, newx=x.test,s="lambda.1se")
    var_sel <- predict(cv.fit,newx=x.test,  s="lambda.1se", type="nonzero")
    nonzero.ctrl1[i,]<-as.numeric(1:p %in% as.numeric(var_sel[[1]]))
    coef.ctrl1[i,]<- as.numeric(coef(cv.fit, s="lambda.1se"))[-1]
    
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
  result_ctrl1<-list(MIRNA=MiRNA_ctrl1,R_square=R_squared.ctrl1,RMSE=Rmse.ctrl1)
  return(result_ctrl1)
}






