library(tidyverse)
library(knitr)
library(corrplot)
library(viridis)
library(ggplot2)
library(scico)
library(MASS)
library(ggridges)
library(reshape2)

setwd("~/GitHub/yesweCanCAM")
cmap_yeo = c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E")

df <- read.csv('./Data/yeo_SSD_age.csv', header = FALSE)
real <- df[1,]
real <- melt(real)
colnames(real) <- c("variable","real")
null <- df[2:nrow(df),]
null <- melt(null)
colnames(null) <- c("variable","null")
full <- merge(null, real, by = "variable")
full$type <- "dispersion"
full$type <- as.factor(full$type)

### within connecticity
null2 <- read.csv('./Data/yeo_perm_t_within_meanregress.csv', header = FALSE)
null2 <- melt(null2)
colnames(null2) <- c("variable","null")

real2 <- read.csv('./Data/yeo_within_t_meanregress.csv',header = FALSE)
real2 <- melt(real2)
colnames(real2) <- c("variable","real")

full2 <- merge(null2, real2, by = "variable")
full2$type <- "within network connectivity"
full2$type <- as.factor(full2$type)

full3 <- rbind(full, full2)
full3$type <- as.factor(full3$type)

pdf('./Figures/subplots/Figure_2A.pdf', height = 5 , width = 11.5)
ggplot(data=full3, aes(y=type,x=null,fill=variable)) +
  geom_density_ridges(alpha = 0.5) +
  geom_point(aes(y=type,x=real, fill=variable),size=5,shape=21) +
  scale_fill_manual(values = cmap_yeo) + 
  scale_colour_manual(values = cmap_yeo) + 
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        title = element_blank()) +
  theme(
    #axis.text.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    panel.spacing = unit(2, "lines")
  ) +
  coord_flip() +
  facet_wrap(~variable, ncol = 7)
dev.off()

