#' Read in output of PRIMATA-APA, process the data for model training
#'
#' @description Read in the output of PRIMATA-APA, ensure the miRNA expression and numTS features are matched, take log transformation for numTS features, merge outcome scores with predictor matrices.
#'
#' @param numTSdir Path to the numTS table output from PRIMATA-APA. Please check our example data (Example_numTS.txt) to make sure there is no formatting issue.
#' @param Exprdir Path to the miR expression table output from PRIMATA-APA. Please check our example data (Example_EXPR.txt) to make sure there is no formatting issue.
#' @param hallmark User provided outcome features to be predicted. Please check our example data (Example_hallmark.txt) to make sure there is no formatting issue.
#' @return The return will be a list of two data frames consisting numTS predictors with outcome score and expression predictors with outcome score.
#'
#' @details \code{process_data} will read in the output of PRIMATA-APA programs (provided on the github) and process the data for predictive model training.
#'
#' @examples
#' dataList  <- process_data(numTSdir="Path/to/numTSfile", EXPRdir="Path/to/Exprfile", hallmark="Path/to/outcomeScorefile")
#'
#'
process_data <- function(numTSdir, EXPRdir, hallmark){
  df=read.delim(numTSdir)
  df.ctrl=read.delim(EXPRdir)
  df.hallmark=read.delim(hallmark)

  ## These lines are just to ensure we only keep patients with all info available
  ## Theoratically, PRIMATA-APA will match miRNA expression and miRNA numTS, so the input should be matched already.
  commonPat <- Reduce(intersect, list(df$patName, df.ctrl$patName, df.hallmark$patName))
  df <- df[df$patName %in% commonPat,]
  df <- df[match(df$patName, commonPat),]
  df.ctrl <- df.ctrl[df.ctrl$patName %in% commonPat,]
  df.ctrl <- df.ctrl[match(df.ctrl$patName, commonPat),]
  df.hallmark <- df.hallmark[df.hallmark$patName %in% commonPat,]
  df.hallmark <- df.hallmark[match(df.hallmark$patName, commonPat),]
  ## remove the patient name, no longer need them as records were matched by patient ID.
  df <- df[, -1]
  df.ctrl <- df.ctrl[,-1]
  df.hallmark <- df.hallmark[,-1]
  ## miRNA expression is in log(RPM+1), we do the same to miRNA numTS
  df1<-df
  for(i in 1:ncol(df1)){
    df1[,i] <- log2(df1[,i] + 1)
  }
  ## We combine predictor matrices with outcome matrix. We put the outcome score as the last column of the predictor matrices.
  df1 <- cbind(df1, df.hallmark)
  df.ctrl <- cbind(df.ctrl, df.hallmark)
  return(list(df=df1, df.ctrl=df.ctrl))
}
