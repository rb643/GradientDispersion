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

## get all the age residuals
m7 <- lm(age ~ DMN_w + motion + sex, data = df.corr)
m6 <- lm(age ~ FrontoParietal_w + motion + sex, data = df.corr)
m2 <- lm(age ~ Somatomotor_w + motion + sex, data = df.corr)
m3 <- lm(age ~ DAN_w + motion + sex, data = df.corr)
m4 <- lm(age ~ VAN_w + motion + sex, data = df.corr)
m5 <- lm(age ~ Limbic_w + motion + sex, data = df.corr)
m1 <- lm(age ~ Visual_w + motion + sex, data = df.corr)

## get all the model residuals
rm7 <- lm(DMN ~ DMN_w + motion + sex, data = df.corr)
rm6 <- lm(Fronto.Parietal ~ FrontoParietal_w + motion + sex, data = df.corr)
rm2 <- lm(Sensory.Motor ~ Somatomotor_w + motion + sex, data = df.corr)
rm3 <- lm(Dorsal.Attention_ ~ DAN_w + motion + sex, data = df.corr)
rm4 <- lm(Ventral.Attention ~ VAN_w + motion + sex, data = df.corr)
rm5 <- lm(Limbic ~ Limbic_w + motion + sex, data = df.corr)
rm1 <- lm(Visual ~ Visual_w + motion + sex, data = df.corr)

## put all models in list and create a list of network names
models <- list(m1, m2, m3, m4, m5, m6, m7)
models2 <- list(rm1, rm2, rm3, rm4, rm5, rm6, rm7)
names <- c("Visual","SM","DAN","VAN","Limbic","FP","DMN")

## collect all residuals and put them in dataframe
df <- data.frame()
for (m in 1:7){
  df.res <- as.data.frame(cbind(models[[m]]$residuals, models2[[m]]$residuals))
  colnames(df.res) <- c("age","value")
  df.res$network <- names[m]
  
  df <- rbind(df,df.res)
}

## sort the network to make sure ggplot doesn't automatically put them in alphabetical order
df$network <- factor(df$network, levels=unique(df$network))

## get the Yeo colormap (in the order of the names)
cmap_yeo = c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E")

## plot all
pdf('./Figures/subplots/Figure_2B.pdf', height = 5 , width = 11.5)
ggplot(df,aes(x=age,y=value)) +
  geom_point(aes(y=value,x=age, fill=network),size=1.5,alpha=0.7,shape=21) +
  geom_smooth(colour = 'black', fill = 'grey',alpha = 0.5,method = "loess") +
  #geom_smooth(colour = 'black', fill = 'black',alpha = 0.6,method = "lm") +
  scale_fill_manual(values = cmap_yeo) + 
  scale_colour_manual(values = cmap_yeo) + 
  #geom_smooth(colour = 'blue', fill = 'blue',alpha = 0.6,method = "lm") +
  #scale_color_scico(palette = "lajolla", direction = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1)
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()
  ) + 
  facet_wrap(~network, ncol = 7)
dev.off()
  
