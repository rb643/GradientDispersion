## load some libraries
library(tidyverse)
library(knitr)
library(ggplot2)
library(viridis)
library(ggpmisc)
library(reshape2)
library(ggridges)
source("https://bit.ly/2q4XQ66")

setwd("~/GitHub/yesweCanCAM")

df <- read.csv('./Data/variance.csv')
df <- melt(df,id.vars = c("Sub","Age","Motion"), measure.vars = c("Var_G1","Var_G2","Var_G3"))


ggplot(data = df, aes(y = value, x = variable, colour = Age)) +
  geom_point(position = position_jitter(width = .15), 
            size = 1.5, alpha = 0.5) +
  scale_colour_viridis(option = "C", direction = 1) +
  geom_boxplot(outlier.shape = NA,alpha = 0.5) +
  #guides(fill = FALSE) +
  #guides(color = FALSE) +
  xlab("Gradient") +
  ylab("Variance explained") +
  theme(legend.position = "bottom") +
  theme_cowplot()
