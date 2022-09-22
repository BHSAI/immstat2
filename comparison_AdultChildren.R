load("C:\\plu\\project\\Kenya Study\\Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")

library("caret")
#remove features with near zero variance and zero variance

q <- data.frame(t(wide.df))
h <- nearZeroVar(q, saveMetrics=TRUE)
q <- q[,(colnames(q) %in% rownames(h[h[,"nzv"]==FALSE | h[,"zeroVar"]==TRUE,]))]
q <- q[,(colnames(q) %in% rownames(h[h[,"zeroVar"]==FALSE,]))]

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

DEP <- function(P, NP, FDR){
  results <- data.frame(matrix(ncol = 2, nrow = 0))
  names(results) <- c("name","pvalue")
  
  for(i in names(P)){
    #print(i)
    p <- comput.p(P[,i], NP[,i],0.05)
    results[nrow(results) + 1,] = list(i,p)
  }
  
  results$BH = p.adjust(results$pvalue, method = "BH")
  
  
  results$FC <- gtools::foldchange(colMeans(P),colMeans(NP))
  #return(results$name[which(abs(results$FC) > 1 & results$BH <= FDR)])
  
  return(results$name[which(results$BH <= FDR)])
}

DEP28_pooled <- DEP(q[which(pheno.df$ACPR28 == 1),],q[which(pheno.df$ACPR28 == 0),], 0.1)
#anyNA(q[which(!(is.na(pheno.df$ACPR28))),DEP28_pooled])

DEP28_AL <- DEP(q[which(pheno.df$Study.Arm == "AL" & pheno.df$ACPR28 == 1),],
                q[which(pheno.df$Study.Arm == "AL" &pheno.df$ACPR28 == 0),], 0.1)

DEP28_ASMQ <- DEP(q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$ACPR28 == 1),],
                  q[which(pheno.df$Study.Arm == "ASMQ" &pheno.df$ACPR28 == 0),], 0.1)

#write(DEP28_ASMQ, "C:\\plu\\project\\Kenya Study\\result\\FeatureRankCorrelation\\DEP28_ASMQ_277")


DEP42_pooled <- DEP(q[which(pheno.df$ACPR42 == 1),],q[which(pheno.df$ACPR42 == 0),], 0.1)
#anyNA(q[which(!(is.na(pheno.df$ACPR42))),DEP42_pooled])

DEP42_AL <- DEP(q[which(pheno.df$Study.Arm == "AL" & pheno.df$ACPR42 == 1),],
                q[which(pheno.df$Study.Arm == "AL" &pheno.df$ACPR42 == 0),], 0.1)

DEP42_ASMQ <- DEP(q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$ACPR42 == 1),],
                  q[which(pheno.df$Study.Arm == "ASMQ" &pheno.df$ACPR42 == 0),], 0.1)

if (anyNA(q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]) |
    anyNA(q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ])){
  #q <- read.csv("/home/plu/data/kenya/imp_kenya_20190315.csv",row.names = 1, header = TRUE)
  #q <- read.csv("C:\\Users\\P\\Downloads\\KS\\data\\imp_kenya_20190315.csv",row.names = 1, header = TRUE)
  q <- read.csv("C:\\plu\\project\\Kenya Study\\data\\imp_kenya_20190315.csv",row.names = 1, header = TRUE)
}

##Age.Cat1

#Day28
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "1" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_17 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "2" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "3" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]

df <- data.frame(matrix(ncol = 4, nrow = 78))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_17),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 17", "adult"), each = 26)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens") 

#Day42
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "1" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_17 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "2" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "3" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]


df <- data.frame(matrix(ncol = 4, nrow = 30))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_17),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 17", "adult"), each = 10)


library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens") 

##Age.Cat1
##All features except selected features by univariate analysis

#Day28

Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "1" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]
Age_17 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "2" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "3" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]

df <- data.frame(matrix(ncol = 4, nrow = 3183))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_17),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 17", "adult"), each = 1061)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens") 

#Day42
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "1" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]
Age_17 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "2" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat1 == "3" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]


df <- data.frame(matrix(ncol = 4, nrow = 3231))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_17),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 17", "adult"), each = 1077)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens")

##Age.Cat2

#Day28
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "1" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_15 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "2" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "3" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]

df <- data.frame(matrix(ncol = 4, nrow = 78))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_15),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 15", "elder"), each = 26)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens") 

#Day42
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "1" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_15 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "2" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat2 == "3" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]


df <- data.frame(matrix(ncol = 4, nrow = 30))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_15),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 15", "elder"), each = 10)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens")

##Age.Cat3

pheno.df$Age.Cat3 = pheno.df$Age.Yrs

for(i in 1:length(pheno.df$Age.Cat3)){
  if (pheno.df$Age.Cat3[i] < 5){
    pheno.df$Age.Cat3[i] = 1
  } else if(pheno.df$Age.Cat3[i] > 12){
    pheno.df$Age.Cat3[i] = 3
  } else {pheno.df$Age.Cat3[i] = 2}
}


#Day28
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "1" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_12 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "2" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "3" & !(is.na(pheno.df$ACPR28))),DEP28_ASMQ]

df <- data.frame(matrix(ncol = 4, nrow = 831))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_12),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 12", "12 or older"), each = 277)

df$Age <- factor(df$Age, levels = c("5 or younger", "5 to 12", "12 or older"))

#install.packages("car")
library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
tiff("kenya.tiff",units = "in", width = 5.5, height = 4, res = 300)
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.9), 
                 grid = FALSE,regLine = TRUE, smooth = FALSE,
                 ylim = c(-0.7, 2.6), 
                 xlim = c(-0.7, 2.6),
                 xlab="Mean Immune Response Intensity in Subjects Aged 12 or Older", ylab="Mean Immune Response Intensity", 
                 main="Day 28, 277 Immune Responses")
dev.off()

lm(df$Mean[which(df$Age == "5 or younger")] ~ df$Mean_adult[which(df$Age == "5 or younger")], data=df)
lm(df$Mean[which(df$Age == "5 to 12")] ~ df$Mean_adult[which(df$Age == "5 to 12")], data=df)
lm(df$Mean[which(df$Age == "12 or older")] ~ df$Mean_adult[which(df$Age == "12 or older")], data=df)

#Day42
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "1" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_12 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "2" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "3" & !(is.na(pheno.df$ACPR42))),DEP42_ASMQ]


df <- data.frame(matrix(ncol = 4, nrow = 30))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_12),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 12", "12 or older"), each = 10)

df$Age <- factor(df$Age, levels = c("5 or younger", "5 to 12", "12 or older"))

tiff("kenya.tiff",units = "in", width = 5.5, height = 4, res = 300)
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.9), 
                 grid = FALSE,regLine = TRUE, smooth = FALSE,
                 ylim = c(-0.7, 2.6), 
                 xlim = c(-0.7, 2.6),
                 xlab="Mean Immune Response Intensity in Subjects Aged 12 or Older", ylab="Mean Immune Response Intensity", 
                 main="Day 42, 10 Immune Responses")
dev.off()

lm(df$Mean[which(df$Age == "5 or younger")] ~ df$Mean_adult[which(df$Age == "5 or younger")], data=df)
lm(df$Mean[which(df$Age == "5 to 12")] ~ df$Mean_adult[which(df$Age == "5 to 12")], data=df)
lm(df$Mean[which(df$Age == "12 or older")] ~ df$Mean_adult[which(df$Age == "12 or older")], data=df)

##Age.Cat3
##All features except selected features by univariate analysis

pheno.df$Age.Cat3 = pheno.df$Age.Yrs

for(i in 1:length(pheno.df$Age.Cat3)){
  if (pheno.df$Age.Cat3[i] < 5){
    pheno.df$Age.Cat3[i] = 1
  } else if(pheno.df$Age.Cat3[i] > 12){
    pheno.df$Age.Cat3[i] = 3
  } else {pheno.df$Age.Cat3[i] = 2}
}

#Day28

Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "1" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]
Age_12 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "2" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "3" & !(is.na(pheno.df$ACPR28))),colnames(q)[which(!(colnames(q) %in% DEP28_ASMQ))]]

df <- data.frame(matrix(ncol = 4, nrow = 3183))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_12),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 12", "elder"), each = 1061)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens") 

#Day42
Age_5 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "1" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]
Age_12 <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "2" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]
Age_adult <- q[which(pheno.df$Study.Arm == "ASMQ" & pheno.df$Age.Cat3 == "3" & !(is.na(pheno.df$ACPR42))),colnames(q)[which(!(colnames(q) %in% DEP42_ASMQ))]]


df <- data.frame(matrix(ncol = 4, nrow = 3231))
x <- c("Mean_adult", "Mean", "Feature", "Age")
colnames(df) <- x

df$Mean_adult <- rep(colMeans(Age_adult), times=3)
df$Mean <- c(colMeans(Age_5),colMeans(Age_12),colMeans(Age_adult))
df$Feature <- rep(colnames(Age_adult), times=3)
df$Age <- rep(c("5 or younger", "5 to 12", "elder"), each = 1077)

library(car)
# Enhanced Scatterplot of MPG vs. Weight 
# by Number of Car Cylinders 
car::scatterplot(df$Mean ~ df$Mean_adult | Age, data=df,
                 legend=list(coords="topleft", cex=0.8), 
                 xlab="Mean Signal Intensity in Adults", ylab="Mean Signal Intensity", 
                 main="Mean Ab Signal Intensitities to Pf Antigens")