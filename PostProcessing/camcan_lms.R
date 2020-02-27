library(tidyverse)
library(knitr)
library(corrplot)
library(viridis)
library(ggplot2)
library(scico)
library(MASS)
library(sfsmisc)

setwd("~/GitHub/yesweCanCAM")

# load some data for linear models (includes dispersion and within connectivity for each Yeo network
# + clustering and pathlength) 
df.corr <- read.csv('./Data/camcan_lms.csv', header =T)
  
m1 <- rlm(age ~ DMN + DMN_w + motion + sex, data = df.corr)
m2 <- rlm(age ~ Fronto.Parietal + FrontoParietal_w + motion + sex, data = df.corr)
m3 <- rlm(age ~ Sensory.Motor + Somatomotor_w + motion + sex, data = df.corr)
m4 <- rlm(age ~ Dorsal.Attention_ + DAN_w + motion + sex, data = df.corr)
m5 <- rlm(age ~ Ventral.Attention + VAN_w + motion + sex, data = df.corr)
m6 <- rlm(age ~ Limbic + Limbic_w + motion + sex, data = df.corr)
m7 <- rlm(age ~ Visual + Visual_w + motion + sex, data = df.corr)

# merge models
models <- list(m1, m2, m3, m4, m5, m6, m7)

# extract the relevant stats
F <- p <- t <- numeric()
for (m in 1:7){
  F[[m]] <- f.robftest(models[[m]], var = 2)$statistic[[1]]
  p[[m]] <- f.robftest(models[[m]], var = 2)$p.value
  t[[m]] <- as.data.frame(summary(models[[m]])$coefficients)$`t value`[2]
}

output <- as.data.frame(cbind(F,p,t))
output$model <- as.factor(c("DMN","FrontoParietal","Sensory.Motor","DAN","VAN","Limbic","Visual"))
output$pFDR <- p.adjust(output$p, method = "fdr")

## correcting for CT
df.corr <- df.corr[-71,] #one missing CT subject..
ct <- as.data.frame(t(read.csv('./Data/yeo_ct.csv',header = FALSE)))
colnames(ct) <- c("Visual_CT","Sensory.Motor_CT","DAN_CT","VAN_CT","Limbic_CT","FrontoParietal_CT","DMN_CT")
df.corr <- cbind(df.corr,ct)

c1 <- rlm(age ~ DMN + DMN_w + DMN_CT + motion + sex, data = df.corr)
c2 <- rlm(age ~ Fronto.Parietal + FrontoParietal_w + FrontoParietal_CT + motion + sex, data = df.corr)
c3 <- rlm(age ~ Sensory.Motor + Somatomotor_w + Sensory.Motor_CT + motion + sex, data = df.corr)
c4 <- rlm(age ~ Dorsal.Attention_ + DAN_w + DAN_CT + motion + sex, data = df.corr)
c5 <- rlm(age ~ Ventral.Attention + VAN_w + VAN_CT + motion + sex, data = df.corr)
c6 <- rlm(age ~ Limbic + Limbic_w + Limbic_CT + motion + sex, data = df.corr)
c7 <- rlm(age ~ Visual + Visual_w + Visual_CT + motion + sex, data = df.corr)

ctmodels <- list(c1, c2, c3, c4, c5, c6, c7)

F2 <- p2 <- t2 <- numeric()
for (m in 1:7){
  F2[[m]] <- f.robftest(ctmodels[[m]], var = 2)$statistic[[1]]
  p2[[m]] <- f.robftest(ctmodels[[m]], var = 2)$p.value
  t2[[m]] <- as.data.frame(summary(ctmodels[[m]])$coefficients)$`t value`[2]
}

output2 <- as.data.frame(cbind(F2,p2,t2))
output2$model <- as.factor(c("DMN","FrontoParietal","Sensory.Motor","DAN","VAN","Limbic","Visual"))
output2$pFDR <- p.adjust(output2$p, method = "fdr")


## plot some of the residuals
cmap_yeo = c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E")

Rage <- lm(age ~ motion + sex + Limbic_w, data = df.corr)
Rdmn <- lm(Limbic ~ motion + sex +  Limbic_w, data = df.corr) 

df.res <- as.data.frame(cbind(Rage$residuals, Rdmn$residuals))
colnames(df.res) <- c("age","Limbic")

pdf('./Figures/subplots/Limbic_CT_res.pdf', height = 6.6 , width = 4.2)
ggplot(df.res,aes(x=age,y=Limbic)) +
  geom_point(shape=21, alpha = 0.7,fill = cmap_yeo[5],size = 1.5) + 
  geom_smooth(colour = 'green', fill = 'green',alpha = 0.6,method = "loess") +
  geom_smooth(colour = 'blue', fill = 'blue',alpha = 0.6,method = "lm") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
        )
dev.off()
