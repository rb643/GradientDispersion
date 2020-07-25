library(corrplot)
library(viridis)
library(ggplot2)
library(scico)
library(MASS)
library(sfsmisc)

setwd("~/GitHub/GradientDispersion/")

# load some data for linear models (includes dispersion and within connectivity for each Yeo network
# + clustering and pathlength) 
data <- read.csv('./Data/camcan_lms.csv', header =T)
library(MBESS)
library(corpcor)
library(yhat)
library(miscTools)
library(boot)

lm.out <- lm(cattell ~ age + DMN + DMN_w + motion + sex,data=data)
regrOut<-calc.yhat(lm.out)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm.out, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)

write.csv(result$combCIaps,'./Scratch/Visual.csv')


data$Age_z <- (data$age - mean(data$age))/sd(data$age)
data$Motion_z <- (data$motion - mean(data$motion))/sd(data$motion)
data$Cattell_z <- (data$cattell - mean(data$cattell))/sd(data$cattell)

data$DMN_z <- (data$DMN - mean(data$DMN))/sd(data$DMN)
data$Visual_z <- (data$Visual - mean(data$Visual))/sd(data$Visual)
data$FP_z <- (data$Fronto.Parietal - mean(data$Fronto.Parietal))/sd(data$Fronto.Parietal)
data$DAN_z <- (data$Dorsal.Attention_ - mean(data$Dorsal.Attention_))/sd(data$Dorsal.Attention_)
data$VAN_z <- (data$Ventral.Attention - mean(data$Ventral.Attention))/sd(data$Ventral.Attention)
data$Limbic_z <- (data$Limbic - mean(data$Limbic))/sd(data$Limbic)
data$SM_z <- (data$Sensory.Motor - mean(data$Sensory.Motor))/sd(data$Sensory.Motor)

data$DMN_Wz <- (data$DMN_w - mean(data$DMN_w))/sd(data$DMN_w)
data$Visual_Wz <- (data$Visual_w - mean(data$Visual_w))/sd(data$Visual_w)
data$FP_Wz <- (data$FrontoParietal_w - mean(data$FrontoParietal_w))/sd(data$FrontoParietal_w)
data$DAN_Wz <- (data$DAN_w - mean(data$DAN_w))/sd(data$DAN_w)
data$VAN_Wz <- (data$VAN_w - mean(data$VAN_w))/sd(data$VAN_w)
data$Limbic_Wz <- (data$Limbic_w - mean(data$Limbic_w))/sd(data$Limbic_w)
data$SM_Wz <- (data$Somatomotor_w - mean(data$Somatomotor_w))/sd(data$Somatomotor_w)


data$DMN_age <- data$Age_z * data$DMN_z
data$FP_age <- data$Age_z * data$FP_z
data$Limbic_age <- data$Age_z * data$Limbic_z
data$VAN_age <- data$Age_z * data$VAN_z
data$DAN_age <- data$Age_z * data$DAN_z
data$SM_age <- data$Age_z * data$SM_z
data$Visual_age <- data$Age_z * data$Visual_z

lm1 <- lm(Cattell_z ~ Age_z + DMN_z + DMN_Wz + Motion_z + sex,data=data)
lm2 <- lm(Cattell_z ~ Age_z + FP_z + FP_Wz + Motion_z + sex,data=data)
lm3 <- lm(Cattell_z ~ Age_z + Limbic_z +  Limbic_Wz + motion + sex,data=data)
lm4 <- lm(Cattell_z ~ Age_z + VAN_z + VAN_Wz + Motion_z + sex,data=data)
lm5 <- lm(Cattell_z ~ Age_z + DAN_z + DAN_Wz + Motion_z + sex,data=data)
lm6 <- lm(Cattell_z ~ Age_z + SM_z + SM_Wz + Motion_z + sex,data=data)
lm7 <- lm(Cattell_z ~ Age_z + Visual_z + Visual_Wz + Motion_z + sex,data=data)

regrOut<-calc.yhat(lm1)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm1, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/DMN_common_boot.csv')

regrOut<-calc.yhat(lm2)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm2, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/FP_common_boot.csv')

regrOut<-calc.yhat(lm3)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm3, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/Limbic_common_boot.csv')

regrOut<-calc.yhat(lm4)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm4, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/VAN_common_boot.csv')

regrOut<-calc.yhat(lm5)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm5, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/DAN_common_boot.csv')

regrOut<-calc.yhat(lm6)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm6, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/SM_common_boot.csv')

regrOut<-calc.yhat(lm7)
boot.out<- boot(data,boot.yhat,1000,lmOut=lm7, regrout0=regrOut)
result<-booteval.yhat(regrOut,bty= "perc",boot.out)
write.csv(result$combCIaps,'./Scratch/Visual_common_boot.csv')

lm1 <- lm(Cattell_z ~ Age_z + DMN_z + Age_z*DMN_z + Age_z*sex + sex*DMN_z + DMN_Wz + Motion_z + sex,data=data)
lm2 <- lm(Cattell_z ~ Age_z + FP_z + Age_z*FP_z + Age_z*sex + sex*FP_z + FP_Wz + Motion_z + sex,data=data)
lm3 <- lm(Cattell_z ~ Age_z + Limbic_z + Age_z*Limbic_z + Age_z*sex + sex*Limbic_z + Limbic_Wz + motion + sex,data=data)
lm4 <- lm(Cattell_z ~ Age_z + VAN_z + Age_z*VAN_z + Age_z*sex + sex*VAN_z + VAN_Wz + Motion_z + sex,data=data)
lm5 <- lm(Cattell_z ~ Age_z + DAN_z + Age_z*DAN_z + Age_z*sex + sex*DAN_z + DAN_Wz + Motion_z + sex,data=data)
lm6 <- lm(Cattell_z ~ Age_z + SM_z + Age_z*SM_z + Age_z*sex + sex*SM_z + SM_Wz + Motion_z + sex,data=data)
lm7 <- lm(Cattell_z ~ Age_z + Visual_z + Age_z*Visual_z + Age_z*sex + sex*Visual_z + Visual_Wz + Motion_z + sex,data=data)


