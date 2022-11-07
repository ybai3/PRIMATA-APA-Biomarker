#' Generate a bootstrap sample and divide data to train and test sets
#'
#' @description  \code{boots.data} is called in \code{featureSelection} function to generate a bootstrap sample and divide data into train and test set.
#'
#' @param data \code{featureSelection} would give the input. Basically, a data frame consisting numTS matrix and outcome score (in the last column).
#' @param train.portion Portion of train set.
#' @return A list consisting 1) train set predictor matrix, 2) train set outcome matrix, 3) test set predictor matrix, 4) test set outcome matrix, 5) number of features.
#'
#' @details This is a function called in \code{featureSelection} function. In each iterative run of model training, it generates a bootstrap sample and splits whole dataset into train and test set, so that model trained in train set can be evaluated in the completely unseen test set. It allows us to do feature selection as well as evaluate model performance.
#'
#'
#' @examples
#' div.data  <- divide.data(df.numTS, train.portion = 0.75)
#'
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

#' Feature selection using Elastic net regression conjugated with bootstrapping
#'
#' @description \code{featureSelection} generates bootstrap sample and splits data into train set and test set. It trains model using elastic net regression and evaluate the performance by calculating R square and RMSE on test set.
#'
#' @param Data A data frame consisting predictor numTS matrix and outcome score (in the last column). Please check our example data (Example_TCGA_numTS.txt and Example_Seo_numTS.txt) to make sure there is no formatting issue.
#' @param num_simulation Desired number of bootstrap samples.
#' @param trainportion Desired portion of train set. Default to 0.75.
#' @return A list containing 1) table of miRNA with the selected frequency and direction of association. 2) R-square estimated in the test set. 3) RMSE estimated in the test set.
#'
#' @details \code{trainModel} train predictive models using elastic net regression conjugated with bootstrapping. It conducts elastic net regression conjugated with bootstrapping to select stable sets of predictive miRNAs. It first creates \code{num_simulation} bootstrap samples from the original dataset by resampling with replacement. Then, with each bootstrap sample, it selects miRNAs using elastic net regression. This step would generate \code{num_simulation} sets of selected miRNAs from \code{num_simulation} bootstrap samples. Then, we calculated the frequency of being selected for each miRNA to indicate its predictive power.
#'
#' @examples
#' Example_TCGA_numTS <- read.delim("./data/Example_TCGA_numTS.txt")
#' resultExample <- featureSelection(Data = Example_TCGA_numTS, num_simulation = 200)
#'
#'
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
