## moderated mediation analysis
library(tidyverse)
library(knitr)
library(lavaan)
library(mediation)
library(rockchalk) 
library(multilevel) 
library(bda) 
library(gvlma) 
library(stargazer)
library(psych)

setwd("~/GitHub/yesweCanCAM")

# load some data
df <- read.csv('./Data/camcan_lms.csv', header = TRUE)
to_keep <- apply(!is.na(df), 1, all)
df <-  df[to_keep,]

# inspect the data  
df %>% 
  dplyr::select(age, sex, motion, cattell, Fronto.Parietal, DMN) %>% 
  pairs.panels(scale = FALSE, pch = ".")

# center the variables
# centre or log-scale all variable that show some heteroscedascity
df$DMN_log <- log(df$DMN)
df$FP_log <- log(df$Fronto.Parietal)
df$VAN_log <- log(df$Ventral.Attention)
df$DAN_log <- log(df$Dorsal.Attention_)
df$cattell_log <- log(df$cattell) #Outcome variable; cattell
df$age_c   <- c(scale(df$age, center=TRUE, scale=FALSE)) #Centering IV; age

### Using the mediation package
library(mediation)
#M = FP dispersion, X = age , Y = cattell
fitM <- lm(FP_log ~ age_c + motion + sex, data=df) #path a
fitY <- lm(cattell_log ~ age_c + FP_log + motion + sex, data=df) #total effect
summary(fitM)
summary(fitY)
gvlma(fitM)
gvlma(fitY)

fitMedBoot <- mediation::mediate(fitM, fitY, boot = TRUE, sims = 999, treat="age_c", mediator="FP_log")
summary(fitMedBoot)
plot(fitMedBoot)

#plot
pdf(file = './Figures/subplots/Medfit_FP_Cattel.pdf', height = 3, width = 4)
  plot(fitMedBoot)
dev.off()

# print to text
sink("./Data/MediationResults/Medfit_FP_Cattel.txt")
print(summary(fitMedBoot))
sink()  

# plot
df$age_high_low <- (df$age > 50)
df$FP_log_highlow <- (df$FP_log > median(df$FP_log))
pdf(file = './Figures/subplots/Medfit_Age_Cattel.pdf', height = 5, width = 5)
ggplot(df,aes(x=age_c,y=cattell_log, fill = FP_log_highlow, colour = FP_log_highlow)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "loess") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
  )
dev.off()

#M = DMN dispersion, X = age , Y = cattell
fitM <- lm(DMN_log ~ age_c + motion + sex, data=df) #path a
fitY <- lm(cattell_log ~ age_c + DMN_log + motion + sex, data=df) #total effect
summary(fitM)
summary(fitY)
gvlma(fitM)
gvlma(fitY)

fitMedBoot <- mediation::mediate(fitM, fitY, boot = TRUE, sims = 999, treat="age_c", mediator="DMN_log")
summary(fitMedBoot)
plot(fitMedBoot)

#plot
pdf(file = './Figures/subplots/Medfit_DMN_Cattel.pdf', height = 3, width = 4)
plot(fitMedBoot)
dev.off()

# print to text
sink("./Data/MediationResults/Medfit_DMN_Cattel.txt")
print(summary(fitMedBoot))
sink()  

# plot
df$DMN_log_highlow <- (df$DMN_log > median(df$DMN_log))
pdf(file = './Figures/subplots/Medfit_Age_Cattel_DMN.pdf', height = 5, width = 5)
ggplot(df,aes(x=age_c,y=cattell_log, fill = DMN_log_highlow, colour = DMN_log_highlow)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "loess") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
  )
dev.off()

#M = VAN dispersion, X = age , Y = cattell
fitM <- lm(VAN_log ~ age_c + motion + sex, data=df) #path a
fitY <- lm(cattell_log ~ age_c + VAN_log + motion + sex, data=df) #total effect
summary(fitM)
summary(fitY)
gvlma(fitM)
gvlma(fitY)

fitMedBoot <- mediation::mediate(fitM, fitY, boot = TRUE, sims = 999, treat="age_c", mediator="VAN_log")
summary(fitMedBoot)
plot(fitMedBoot)

#plot
pdf(file = './Figures/subplots/Medfit_VAN_Cattel.pdf', height = 3, width = 4)
plot(fitMedBoot)
dev.off()

# print to text
sink("./Data/MediationResults/Medfit_VAN_Cattel.txt")
print(summary(fitMedBoot))
sink()  

# plot
df$VAN_log_highlow <- (df$VAN_log > median(df$VAN_log))
pdf(file = './Figures/subplots/Medfit_Age_Cattel_VAN.pdf', height = 5, width = 5)
ggplot(df,aes(x=age_c,y=cattell_log, fill = VAN_log_highlow, colour = VAN_log_highlow)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "loess") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
  )
dev.off()

#M = VAN dispersion, X = age , Y = cattell
fitM <- lm(DAN_log ~ age_c + motion + sex, data=df) #path a
fitY <- lm(cattell_log ~ age_c + DAN_log + motion + sex, data=df) #total effect
summary(fitM)
summary(fitY)
gvlma(fitM)
gvlma(fitY)

fitMedBoot <- mediation::mediate(fitM, fitY, boot = TRUE, sims = 999, treat="age_c", mediator="DAN_log")
summary(fitMedBoot)
plot(fitMedBoot)

#plot
pdf(file = './Figures/subplots/Medfit_DAN_Cattel.pdf', height = 3, width = 4)
plot(fitMedBoot)
dev.off()

# print to text
sink("./Data/MediationResults/Medfit_DAN_Cattel.txt")
print(summary(fitMedBoot))
sink()  

# plot
df$DAN_log_highlow <- (df$DAN_log > median(df$DAN_log))
pdf(file = './Figures/subplots/Medfit_Age_Cattel_DAN.pdf', height = 5, width = 5)
ggplot(df,aes(x=age_c,y=cattell_log, fill = DAN_log_highlow, colour = DAN_log_highlow)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "loess") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
  )
dev.off()
