#load("/home/plu/data/kenya/Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")

load("C:\\plu\\project\\Kenya Study\\Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")

q <- data.frame(t(wide.df))
# min(q, na.rm = TRUE)
# 
# #check negative vaules
# has.neg <- apply(q_raW, 2, function(x) any(x < 0))
# #which(has.neg)
# length(which(has.neg))
# 
# #check negative vaules
# has.neg <- apply(q, 2, function(x) any(x < 0))
# which(has.neg)
# length(which(has.neg))
# 
# sum_nneg <- 0
# sum_npos <- 0
# sum_nzero <- 0
# for(i in 1:ncol(q_raW)){
#   nneg <- length(q_raW[,i][q_raW[,i] < 0 & !(is.na(q_raW[,i]))])
#   npos <- length(q_raW[,i][q_raW[,i] > 0 & !(is.na(q_raW[,i]))])
#   nzero <- length(q_raW[,i][q_raW[,i] == 0 & !(is.na(q_raW[,i]))])
#   sum_nneg = sum_nneg + nneg
#   sum_npos = sum_npos + npos
#   sum_nzero = sum_nzero + nzero
# }
# 
# percent_nneg <- sum_nneg/(nrow(q_raW)*ncol(q_raW))
# percent_npos <- sum_npos/(nrow(q_raW)*ncol(q_raW))
# 
# #replace negative value with zero
# q[q<0 & !(is.na(q))] <- 0

# convert .RData -> .rdb/.rdx
# e = local({load("/home/plu/data/kenya/Pf1000_USAMRU-K_ACTs_Dataxfer.rdata"); environment()})
# tools:::makeLazyLoadDB(e, "New")
# lazyLoad("New")

library(caret)
library(doMC)
library(gtools)
library(stringr)
library(data.table)
library(gtools)

# #remove features with near zero variance and zero variance
# h <- nearZeroVar(q, saveMetrics=TRUE)
# q <- q[,(colnames(q) %in% rownames(h[h[,"nzv"]==FALSE | h[,"zeroVar"]==TRUE,]))]
# q <- q[,(colnames(q) %in% rownames(h[h[,"zeroVar"]==FALSE,]))]


comput.p <- function(X,Y,a){
  if(diff(range(X, na.rm = TRUE)) == 0 | diff(range(Y, na.rm = TRUE)) == 0 ){
    return(wilcox.test(X, Y, 
                       paired=FALSE, 
                       conf.level=0.95)$p.value)
  }
  else if(shapiro.test(X)$p.value > a & shapiro.test(Y)$p.value > a){
    return(t.test(X, Y, 
                  paired=FALSE, 
                  conf.level=0.95)$p.value)
  }else{
    return(wilcox.test(X, Y, 
                       paired=FALSE, 
                       conf.level=0.95)$p.value)
  }
}

# P <- q[which(pheno.df$Study.Arm == "AL" & pheno.df$ACPR28 == 1),]
# NP <- q[which(pheno.df$Study.Arm == "AL" &pheno.df$ACPR28 == 0),]

DEP <- function(P, NP, FDR){
  results <- data.frame(matrix(ncol = 5, nrow = 0))
  names(results) <- c("name","pvalue","log2FoldChange","fc","ann")
  
  for(i in names(P)){
    #print(i)
    p <- comput.p(P[,i], NP[,i],0.05)
    fc <- gtools::foldchange((mean(P[,i],na.rm = TRUE) + abs(min(q, na.rm = TRUE))),
                             (mean(NP[,i],na.rm = TRUE) + abs(min(q, na.rm = TRUE))))
    results[nrow(results) + 1,] = list(i,p,foldchange2logratio(fc, base=2),fc)
  }
  
  results$q = p.adjust(results$pvalue, method = "BH")
  results$ann <- ann.df[results$name,"Description"]
  
  #return(results$name[which(results$BH <= FDR)])
  #return(results[which(results$q <= FDR),])
  return(results)
}

DEP28_AL <- DEP(q[which(pheno.df$Study.Arm == "AL" & pheno.df$ACPR28 == 1),],
                q[which(pheno.df$Study.Arm == "AL" &pheno.df$ACPR28 == 0),], 0.1)

DEP28_ASMQ <- DEP(q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$ACPR28 == 1),],
                  q[which(pheno.df$Study.Arm == "ASMQ" &pheno.df$ACPR28 == 0),], 0.1)

DEP42_AL <- DEP(q[which(pheno.df$Study.Arm == "AL" & pheno.df$ACPR42 == 1),],
                q[which(pheno.df$Study.Arm == "AL" &pheno.df$ACPR42 == 0),], 0.1)

DEP42_ASMQ <- DEP(q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$ACPR42 == 1),],
                  q[which(pheno.df$Study.Arm == "ASMQ" &pheno.df$ACPR42 == 0),], 0.1)

write.csv(DEP28_ASMQ,"C:\\plu\\paper\\Kenya\\Figures&Tables\\DEP28_ASMQ.csv")
write.csv(DEP42_ASMQ,"C:\\plu\\paper\\Kenya\\Figures&Tables\\DEP42_ASMQ.csv")
  
volcanoPlot <- function(res){
  tiff("kenya.tiff",units = "in", width = 5, height = 5, res = 300)
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=NA,
                 xlim=c(-1,1), ylim=c(0,6), cex.lab = 1.3, cex.axis = 1.5
  ))
  # Add colored points: red if padj<0.1, orange of log2FC>1, green if both)
  with(subset(res, q<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  #with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
  #with(subset(res, BH<.1 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  dev.off()
}

volcanoPlot(DEP28_AL)
volcanoPlot(DEP28_ASMQ)
volcanoPlot(DEP42_AL)
volcanoPlot(DEP42_ASMQ) 
