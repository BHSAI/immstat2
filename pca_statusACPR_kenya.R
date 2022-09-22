load("C:\\plu\\project\\Kenya Study\\Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")
qq <- data.frame(t(wide.df))

q <- read.csv("C:\\plu\\project\\Kenya Study\\data\\imp_kenya_20190315.csv",row.names = 1, header = TRUE)


ACPR28_AL <- q[which(pheno.df$Study.Arm == "AL" & !(is.na(pheno.df$ACPR28))),]
group_ACPR28_AL <- pheno.df$ACPR28[which(pheno.df$Study.Arm == "AL" & !(is.na(pheno.df$ACPR28)))]

ACPR28_ASMQ <- q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28))),]
group_ACPR28_ASMQ <- pheno.df$ACPR28[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28)))]

ACPR42_AL <- q[which(pheno.df$Study.Arm == "AL" & !(is.na(pheno.df$ACPR42))),]
group_ACPR42_AL <- pheno.df$ACPR42[which(pheno.df$Study.Arm == "AL" & !(is.na(pheno.df$ACPR42)))]

ACPR42_ASMQ <- q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR42))),]
group_ACPR42_ASMQ <- pheno.df$ACPR42[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR42)))]

#Variable importance
varImp_Day28_ASMQ <- read.csv("C:\\plu\\project\\Kenya Study\\result\\rf_statusDay28_kenya_complete_ASMQ_cv_20190327_varImp.csv")
ACPR28_ASMQ_varImp <- q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR28))),
                       as.character(varImp_Day28_ASMQ$X[which(varImp_Day28_ASMQ$Overall >= 50)])]

varImp_Day42_ASMQ <- read.csv("C:\\plu\\project\\Kenya Study\\result\\rf_status_kenya_complete_ASMQ_cv_20190325_varImp.csv")
ACPR42_ASMQ_varImp <- q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR42))),
                        as.character(varImp_Day42_ASMQ$X[which(varImp_Day42_ASMQ$Overall >= 50)])]

# ACPR42_ASMQ_varImp <- q[which(pheno.df$Study.Arm == "ASMQ" & !(is.na(pheno.df$ACPR42))),
#                         varImp_Day42_ASMQ$X[1:2]]


plotDataPCA <- ACPR42_ASMQ_varImp 
# X <- ACPR42_ASMQ
# preProc <- preProcess(X)
# plotDataPCA <- data.frame(predict(preProc, X))

group <- group_ACPR42_ASMQ

#plot the top 2 variables
dataPlot <- plotDataPCA
dataPlot$ACPR <- as.factor(group)

ggplot(dataPlot, aes(x = PF3D7_0317300_s3_454, y = PF3D7_1411600_383, colour = ACPR)) +
  geom_point() +
  stat_ellipse() +
  xlab("rifin (RIF)") +
  ylab("GTP binding protein, putative")

#Box plot
library(ggplot2)
# Basic box plot
p <- ggplot(dataPlot, aes(x=ACPR, y=PF3D7_0900200_e2s1_19))
p
  
#names(plotDataPCA) <- ann.df[names(plotDataPCA),"Description"]

#PCA
ir.pca <- prcomp(plotDataPCA,
                 center = TRUE,
                 scale. = TRUE) 
names(ir.pca)
# print method
print(ir.pca)

# plot method
par(mfrow=c(1,1))
plot(ir.pca, type = "l")

# summary method
summary(ir.pca)

#install.packages("devtools")
library(devtools)

#install_github("vqv/ggbiplot")
library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = as.factor(group), var.axes = FALSE,
              ellipse = TRUE, circle = TRUE)

g <- g + scale_color_discrete(name = 'ACPR42')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
#g <- g + scale_fill_discrete(name = "ACPR", labels = c("0", "1"))
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(g)

#ACPR28
#output a high resolution figure
tiff("test.tiff", units="in", width=5, height=3.3, res=300)
ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(group), var.axes = FALSE,
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = 'ACPR28') +
  theme(legend.direction = 'horizontal',
        legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))
dev.off()

#ACPR42
#output a high resolution figure
tiff("test.tiff", units="in", width=5, height=3.3, res=300)
ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(group), var.axes = FALSE,
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = 'ACPR42') +
  theme(legend.direction = 'horizontal',
        legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))
dev.off()

#compute standard deviation of each principal component
std_dev <- ir.pca$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var[1:10]

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

#scree plot
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
