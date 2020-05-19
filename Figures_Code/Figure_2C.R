library(ggplot2)
library(tibble)
library(reshape)
library(ggtextures)

setwd("~/GitHub/yesweCanCAM")
cmap_yeo = c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E")

df <- read.csv('./Data/yeo_SSD_age.csv', header = FALSE)
colnames(df) <- c("Visual","SM","DAN","VAN","Limbic","FP","DMN")
real <- df[1,]
real <- melt(real)
colnames(real) <- c("variable","real")
real$type <- "dispersion"

df2 <- read.csv('./Data/seg.csv',header = FALSE)
colnames(df2) <- paste("Seg",c("Visual","SM","DAN","VAN","Limbic","FP","DMN"),sep = "_")
main <- read.csv('./Data/camcan_lms.csv')
df2 <- cbind(df2,main)

m1 <- lm(age ~ Seg_Visual + motion + sex, data = df2)
m2 <- lm(age ~ Seg_SM + motion + sex, data = df2)
m3 <- lm(age ~ Seg_DAN + motion + sex, data = df2)
m4 <- lm(age ~ Seg_VAN + motion + sex, data = df2)
m5 <- lm(age ~ Seg_Limbic + motion + sex, data = df2)
m6 <- lm(age ~ Seg_FP + motion + sex, data = df2)
m7 <- lm(age ~ Seg_DMN + motion + sex, data = df2)

models <- list(m1, m2, m3, m4, m5, m6, m7)

t <- numeric()
for (m in 1:7){
  t[[m]] <- as.data.frame(summary(models[[m]])$coefficients)$`t value`[2]
}

output <- as.data.frame(t)
colnames(output) <- "real"
output$variable <- as.factor(c("Visual","SM","DAN","VAN","Limbic","FP","DMN"))
output$type <- "segmentation"

df <- rbind(output,real)
df$type <- factor(df$type, levels=unique(df$type))
df$variable <- as.character(df$variable)
df$variable <- factor(df$variable, levels=unique(df$variable))
df$image <- ifelse(df$type == "dispersion","./Utilities/hideout.svg","./Utilities/temple.svg")

seg <- read.csv('./Data/segregation_tstat.csv')
seg <- subset(seg, meas == 4)
seg$variable <- as.factor(c("Visual","SM","DAN","VAN","Limbic","FP","DMN"))
seg$yeo <- NULL
seg$meas <- NULL
colnames(seg) <- c("real","variable")
seg$type <- "clustering"
seg$image <- "./Utilities/formal-invitation.svg"

df <- rbind(df,seg)
df$type <- factor(df$type, levels = c("dispersion","segmentation","clustering"))

pdf('./Figures/subplots/Figure_2C.pdf', height = 5 , width = 11.5)
ggplot(df, aes(type, real, image = image, fill = variable)) +
  geom_textured_col(img_width = unit(0.5, "null")) +
  geom_bar(stat="identity",alpha = 0.6) +
  scale_fill_manual(values = cmap_yeo) + 
  scale_colour_manual(values = cmap_yeo) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom") +
  labs(y="T") +
  facet_wrap(~variable, ncol = 7)
dev.off()

