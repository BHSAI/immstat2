#load("C:\\plu\\project\\Kenya Study\\Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")

# convert .RData -> .rdb/.rdx
e = local({load("/home/plu/data/kenya/Pf1000_USAMRU-K_ACTs_Dataxfer.rdata"); environment()})
tools:::makeLazyLoadDB(e, "New")
lazyLoad("New")

library(caret)
library(doMC)
# library(gtools)
library(stringr)
library(data.table)

#identify and remove features with near zero variance and zero variance
q <- data.frame(t(wide.df))
h <- nearZeroVar(q, saveMetrics=TRUE)
q <- q[,(colnames(q) %in% rownames(h[h[,"nzv"]==FALSE | h[,"zeroVar"]==TRUE,]))]
q <- q[,(colnames(q) %in% rownames(h[h[,"zeroVar"]==FALSE,]))]


#identify and remove highly correlated features
# rho = cor(q, method="pearson", use="pairwise.complete.obs")
# iiRemove <- findCorrelation(rho, cutoff=0.75) # indices to remove to reduce pair-wise corr
# X <- q[, -iiRemove]
# dim(X)

#corrplot(rho, method = "color")

# if (anyNA(q)){
#   #check for features (columns) and samples (rows) where more than 5% of the data is missing
#   pMiss <- function(x){sum(is.na(x))/length(x)*100}
#   #check samples (rows)
#   rmiss <- apply(t(wide.df),1,pMiss)
#   hist(rmiss)
#   #check features (columns)
#   cmiss <- apply(t(wide.df),2,pMiss)
#   hist(cmiss)
#   
#   #immputation/missing value predictation
#   registerDoMC(cores = 6)
#   preProc_Imput <- preProcess(method="bagImpute", q)
#   complete_q <- predict(preProc_Imput, q)
# }
# 
# write.csv(complete_q , "/home/plu/imp_kenya_20190315.csv")

myData <- read.csv("/home/plu/data/kenya/imp_kenya_20190315.csv",row.names = 1, header = TRUE)
#myData <- read.csv("C:\\plu\\project\\Kenya Study\\data\\imp_kenya_20190315.csv",row.names = 1, header = TRUE)

if(anyNA(myData)){
  myData <- myData[ , colSums(is.na(myData)) == 0]
}


RF <- function(X, Y, repeats){
  
  # auc_obj <- c()
  # Sens <- c()
  # Spec <- c(0)
  # nfold = 5
  # ntime = 10
  # 
  # set.seed(112358)
  # idx <- createMultiFolds(Y,k = nfold, times = ntime)
  # #sapply(idx, length)
  #aucroc <- c()
  
  # set.seed(112358)
  # results <- foreach(i = 1:(nfold*ntime), .combine = cbind) %dopar% {
  #repeats <- 100
  set.seed(112358)
  
  control <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = repeats,
    verboseIter = FALSE,
    returnData = FALSE,
    allowParallel = TRUE,
    sampling = "up",
    selectionFunction = "oneSE",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    savePredictions = T
  )
  
  set.seed(112358)
  rfUnbalance <- train(
    x = data.frame(X),
    y = Y,
    method = "rf",
    ntree = 1500,
    tuneLength = 10,
    metric = "ROC",
    preProcess = c("center", "scale"),
    trControl = control)
  
  #print(rfUnbalance)
  
  #print(rfUnbalance$finalModel)
  
  #imp <- varImp(rfUnbalance)
  
  #confusionMatrix(rfUnbalance$finalModel$predicted,Y, positive = 'Protected')
  #return(rfUnbalance$results$ROC[as.numeric(row.names(rfUnbalance$bestTune))])
  return(rfUnbalance)
}


permutationTEST <- function(X, Y, r){
  
  ls <- RF(X, Y, 100)
  df <- ls$resample
  df$Rep <- stringr::str_split_fixed(df$Resample, '\\.', 2)[,2]
  varImp <- varImp(ls)$importance
  varImp$Ann <- ann.df[row.names(varImp(ls)$importance),"Description"]
  
  write(data.table::setDT(df)[, mean(ROC), by=Rep]$V1, "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_aucroc", ncolumns = 1)
  write.csv(ls$finalModel$confusion, "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_finalModel.csv", row.names = TRUE, col.names = TRUE)
  write.csv(varImp[order(varImp$Overall,decreasing = TRUE),], "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_varImp.csv", row.names = TRUE, col.names = TRUE)
  prediction_data = ls$pred[ls$pred$mtry == ls$finalModel$mtry, ]
  write.csv(confusionMatrix(prediction_data$pred, prediction_data$obs)$table, "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_accuracyCV.csv",row.names = TRUE, col.names = TRUE)
  write(confusionMatrix(prediction_data$pred, prediction_data$obs)$overall[2], "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_Kappa")
  
  
  e <- mean(df$ROC)
  print(paste("e is ",e))
  c <- 0
  permutes <- c()
  holder <- c()
  ep <- c()
  
  
  for(i in 1:r){
    permutes <- as.character(sample(Y, replace = FALSE))
    holder <- rbind(holder, permutes)
  }
  
  
  for (i in 1:r) {
    print(paste(i,"permutation"))
    ep[i] <- mean(RF(X, holder[i,], 10)$resample$ROC)
    print(ep[i])
    if(ep[i] >= e){
      c <- c + 1
    }
  }
  
  #results <- foreach(i = 1:r, .combine = cbind) %dopar% {
  #   print(i)
  #   ep <- greedyRLS(X, holder[i,], k)
  #   print(ep)
  #   return(ep)
  # }
  # for(j in 1:r){
  #   if(results[j] >= e){
  #     c <- c + 1
  #   }
  #}
  p <- c/r
  return(list("P"= p, "aucroc" = e, "aucroc_permutation" = ep, "permutes" = holder))
}

# results_combine_trans <- read.csv("C:\\plu\\project\\ML_malaria_vaccine\\MAL027_R\\GSE89292_D59.csv",row.names = T)[-1]
# outcome <- read.csv("C:\\plu\\project\\ML_malaria_vaccine\\MAL027_R\\GSE89292_outcome.csv")[-1]

#results_combine_trans <- read.csv("/home/plu/data/MAL027/GSE89292_D59.csv")[-1]
#outcome <- read.csv("/home/plu/data/MAL027/GSE89292_outcome.csv")[-1]


registerDoMC(cores = 11)

r = 100

# X <- myData[which(!(is.na(pheno.df$ACPR42))),]
# names(X) <- make.names(names(X))
# Y <- pheno.df[which(!(is.na(pheno.df$ACPR42))),"ACPR42"]
# Y <- as.factor(Y)
# levels(Y) <- c("NP", "P")
# output <- permutationTEST(X, Y, r)

X <- myData[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28))),]
names(X) <- make.names(names(X))
Y <- pheno.df[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28))),"ACPR28"]
Y <- as.factor(Y)
levels(Y) <- c("NP", "P")
output <- permutationTEST(X, Y, r)

lapply(output, print)
write.csv(output$permutes, "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806.csv")
write(output$aucroc_permutation, "/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806_aucrocPerm", ncolumns = 1)

save.image("/home/plu/rf_up_statusDay28_kenya_complete_ASMQ_cv_20200806.RData")

#prediction_data = ls$pred[ls$pred$mtry == ls$finalModel$mtry, ]
#print(confusionMatrix(prediction_data$pred, prediction_data$obs)$table)

#prediction_data = rfUnbalance$pred[rfUnbalance$pred$mtry == rfUnbalance$finalModel$mtry, ]
#print(confusionMatrix(prediction_data$pred, prediction_data$obs))

#prediction_data = rfUnbalance2$pred[rfUnbalance2$pred$mtry == rfUnbalance2$finalModel$mtry, ]
#print(confusionMatrix(prediction_data$pred, prediction_data$obs))