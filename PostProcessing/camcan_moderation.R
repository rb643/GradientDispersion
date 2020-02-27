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


# center the variables
#Centering Data
Y     <- df$cattell_log #Outcome variable
Xc    <- c(scale(df$age, center=TRUE, scale=FALSE)) #Centering IV; age
Zc    <- c(scale(df$Dorsal.Attention_,  center=TRUE, scale=FALSE)) #Centering moderator;
motion <- df$motion
sex <- df$sex

fitMod <- lm(Y ~ Xc + Zc + Xc*Zc + motion + sex) #Model interacts IV & moderator + motion confounrd
summary(fitMod)

# inspect
coef(summary(fitMod))

# summarize
stargazer(fitMod,type="text", title = "FP")


# plotting
df$agebin <- cut(df$age, breaks = 3, labels = c("young","middle","old"))

# FP
FP <- ggplot(df,aes(x=Fronto.Parietal,y=cattell, fill = agebin)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "lm") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  labs(x = "FP Dispersion", y = "Cattell score") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
  )

# FP
DMN <- ggplot(df,aes(x=DMN,y=cattell, fill = agebin)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "lm") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  labs(x = "DMN Dispersion", y = "Cattell score") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
  )

# VAN
VAN <- ggplot(df,aes(x=Ventral.Attention,y=cattell, fill = agebin)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "lm") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  labs(x = "VAN Dispersion", y = "Cattell score") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
  )

# DAN
DAN <- ggplot(df,aes(x=Dorsal.Attention_,y=cattell, fill = agebin)) +
  geom_point(shape=21, alpha = 0.7, size = 1.5) + 
  geom_smooth(alpha = 0.6,method = "lm") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  labs(x = "DAN Dispersion", y = "Cattell score") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
  )

library(ggpubr)

ggarrange(DAN, VAN, FP, DMN, ncol = 4, common.legend = T, 
          labels = c("DAN", "VAN", "FP", "DMN") )
