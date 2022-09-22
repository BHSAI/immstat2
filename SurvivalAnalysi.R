load("C:\\plu\\project\\Kenya Study\\Pf1000_USAMRU-K_ACTs_Dataxfer.rdata")

pheno.df$subjectid <- paste(substr(pheno.df$Specimen.ID,2,4),substr(pheno.df$Specimen.ID,6,8),sep = "")

id <- pheno.df[,c("Sample.ID","subjectid")]

prediction <- read.csv("C:\\plu\\project\\Kenya Study\\result\\Prediciton_kenya_20190621.csv")

merged <- merge(id, prediction, by = "Sample.ID")

id_prediction <- merged[,c("Sample.ID","subjectid","ASMQ_D28_all","ASMQ_D42_all","ASMQ_LTF_2fs")]

# # installing/loading the latest installr package:
# install.packages("installr"); library(installr) # install+load installr
# 
# updateR() # updating R.


# Load required packages
library(survival)
#install.packages("survminer")
library(survminer)
library(dplyr)

#import data
df <- read.csv("C:\\plu\\project\\Kenya Study\\data\\outcomesetvivo1n2.csv")

#replace "MISSED" and "W/D" with NA
df[ df == "MISSED" | df == "W/D" ] <- NA

df$Survtime <- NA
df$Infection <- NA
nTime = 6

for (i in 1:nrow(df)){
  if (all(is.na(df[i,3:8]))){
    next
  }
  if (all(df[i,3:8] == "NEG") & !(anyNA(df[i,3:8]))){
    df$Infection[i] <- 0
    df$Survtime[i] <- nTime * 7 #every 7 days
    next
  }
  for(j in 1:nTime){
    if(is.na(df[i,(j+2)])){
      if(all(is.na(df[i,(j+2):(nTime+2)]))){
        df$Infection[i] <- 0
        df$Survtime[i] <- (j-1) * 7
      }
      break
    }else if(df[i,(j+2)] == "POS"){
      df$Infection[i] <- 1
      df$Survtime[i] <- j * 7
      break
    }
  }
}

df1 <- df[which(df$subjectid %in% id_prediction$subjectid),]

df2 <- merge(df1, id_prediction, by = "subjectid")

#remove rows with undetermined Survtime
df3 <- df2[!(is.na(df2$Survtime)),]

# #subset ASMQ arm
# df2 <- df1[df1$arm == "ASMQ" & df1$day28 != "N/A" & df1$day42 != "N/A",]

###############
#3-class model#
###############
df3$LTF <- NA

for(i in 1:nrow(df3)){
  if(df3$ASMQ_D28_all[i] == "P" & df3$ASMQ_D42_all[i] == "P"){
    df3$LTF[i] <- 2
  } else if(df3$ASMQ_D28_all[i] == "P" & df3$day42[i] != "P") {
    df3$LTF[i] <- 1
  } else {df3$LTF[i] <- 0}
}

df3$status <- NA
for(i in 1:nrow(df3)){
  if(df3$Infection[i] == "1"){
    df3$status[i] <- 2
  } else {df3$status[i] <- 1}
}

df3 <- df3[which(df3$arm == "ASMQ"),]

#Univariate Cox regression
res.cox <- coxph(Surv(Survtime, status) ~ LTF, data = df3)
res.cox
summary(res.cox)

km.fit <- survfit(Surv(Survtime, status) ~ LTF, data = df3)
km.surv <- ggsurvplot(km.fit, risk.table = TRUE)
km.surv

LTF_df <- with(df3,
               data.frame(LTF = c (1, 2,3)))


# Survival curves
fit <- survfit(res.cox, newdata = LTF_df)
png("surv.tiff", units="in", width=6, height=4, res=300, bg = "transparent")
ggsurvplot(fit, data = df3, conf.int = FALSE, legend.labs=c("Vulnerable","Susceptible", "Immune"),
           ggtheme = theme_minimal(), risk.table = FALSE,
           xlab = "Time(days)", ylab = "Disease free probability",
           font.x = c(14, "black"), font.y = c(14, "black"), font.tickslab = c(12,"black"))
dev.off()

ggsurvplot(fit, data = df3, conf.int = FALSE,
           ggtheme = theme_minimal(), risk.table = FALSE,
           xlab = "Time(days)", ylab = "Disease free probability",
           font.x = c(14, "black"), font.y = c(14, "black"), font.tickslab = c(12,"black"))

length(which(df3$LTF == 0 & df3$Survtime == 35))
length(which(df3$LTF == 1 & df3$Survtime == 35))
length(which(df3$LTF == 2 & df3$Survtime == 35))

###############
#2-class model#
###############
df4 <- df3
df4$LTF <- NA

for(i in 1:nrow(df4)){
  if(df4$ASMQ_LTF_2fs[i] == "P"){
    df4$LTF[i] <- 2
  } else if(df4$ASMQ_LTF_2fs[i] == "NP") {
    df4$LTF[i] <- 1
  }
}

df4$status <- NA
for(i in 1:nrow(df4)){
  if(df4$Infection[i] == "1"){
    df4$status[i] <- 2
  } else {df4$status[i] <- 1}
}

df4 <- df4[which(df4$arm == "ASMQ" & df4$day28 == "ACPR"),]

#Univariate Cox regression
res.cox <- coxph(Surv(Survtime, status) ~ LTF, data = df4)
res.cox
summary(res.cox)

LTF_df <- with(df4,
               data.frame(LTF = c (1,2)))


# Survival curves
fit <- survfit(res.cox, newdata = LTF_df)
ggsurvplot(fit, data = df4, conf.int = FALSE, legend.labs=c("Susceptible", "Immune"),
           ggtheme = theme_minimal(),
           xlab = "Time(Day)", ylab = "Disease free probability")



# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = df1$Survtime, event = df1$Infection)
surv_object 

# Fit the Kaplan-Meier curves



fit1 <- survfit(Surv(Survtime, Infection) ~ studyarm, data = df1)
summary(fit1)

#examine the corresponding survival curve
#plot the p-value of a log rank test, pval = TRUE argument
ggsurv <- ggsurvplot(fit1, data = df1, pval = TRUE,
           xlab = "Time(Day)", ylab = "Survival/Infection-free probability")
ggsurv$plot + geom_hline(yintercept = 0.5)

print(ggsurv)


#Pairwise comparisons using Log-Rank test 
res <- pairwise_survdiff(Surv(Survtime, Infection) ~ studyarm, data = df1)
res



# Import the ovarian cancer dataset and have a look at it
data(ovarian)
glimpse(ovarian)

# Dichotomize age and change data labels
ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))

# Data seems to be bimodal
hist(ovarian$age) 
#bi-modal distribution suggests a cutoff of 50 years
ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
ovarian$age_group <- factor(ovarian$age_group)

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 

# Fit the Kaplan-Meier curves
fit1 <- survfit(surv_object ~ rx, data = ovarian)
summary(fit1)

#examine the corresponding survival curve
#plot the p-value of a log rank test, pval = TRUE argument
ggsurvplot(fit1, data = ovarian, pval = TRUE)

