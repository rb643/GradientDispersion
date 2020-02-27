library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(scico)
library(dplyr)

setwd("~/GitHub/yesweCanCAM")

cmap_yeo = as.data.frame(c("#781286","#4682B4","#00760E","#C43AFA","#DCF8A4","#E69422","#CD3E4E"))
colnames(cmap_yeo) <- "color"
cmap_yeo$network <- as.factor(c("Visual","Sensory-Motor","Dorsal Attention",
                      "Ventral Attention","Limbic","Fronto-Parietal",
                      "DMN"))
cmap_yeo$color <- as.character(cmap_yeo$color)

Dist_age <- t(read.csv('./Data/yeo_dist_age.csv',header = FALSE))
Dist_age_p <- t(read.csv('./Data/yeo_dist_age_p.csv',header = FALSE))

p= matrix(NA, 7, 7)
p[lower.tri(p, diag=FALSE)] <- Dist_age_p
p <- 1-p
p[p<0.95] <- NA
p[p>0.95] <- 1

pt= matrix(NA, 7, 7)
pt[lower.tri(pt, diag=FALSE)] <- Dist_age[,1]
pt <- pt*p

colnames(p) <- rownames(p) <- colnames(pt) <- rownames(pt) <- c("Visual","Sensory-Motor","Dorsal Attention",
                                "Ventral Attention","Limbic","Fronto-Parietal",
                                "DMN")

library(circlize)
df <- melt(pt)
df <- df[complete.cases(df),]
## map values to colours
colors <- scico(length(df$value),palette = "oslo", direction = 1)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.05), start.degree = 90, gap.degree = 4)

pdf('./Figures/subplots/Figure_2D1.pdf', height = 5 , width = 5)
chordDiagram(df, annotationTrack = c("grid","name"), transparency = 0.3,
             col = colors,
             #col = ifelse(df$value > 0, "red", colors),
             order=cmap_yeo$network, grid.col = setNames(cmap_yeo$color, cmap_yeo$network))
dev.off()

## separate plot for the legend since circlize doesn't have that option
plotdefault2 <- data.frame(freq = seq(min(df$value),max(df$value), length.out=length(df$value)),
                           y = as.factor(1))  

pdf('./Figures/subplots/Figure_2D1_cbar.pdf',height = 8, width = 1)
ggplot(plotdefault2,aes(freq,y)) +
  geom_tile(aes(fill=freq)) + 
  scale_fill_scico(palette = "oslo", direction = 1) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  coord_flip()
dev.off()


######## for within
Dist_age <- t(read.csv('./Data/yeo_between_t_meanregress.csv',header = FALSE))
Dist_age_p <- t(read.csv('./Data/yeo_between_p_meanregress.csv',header = FALSE))

p= matrix(NA, 7, 7)
p[lower.tri(p, diag=FALSE)] <- Dist_age_p
p <- 1-p
p[p<0.95] <- NA
p[p>0.95] <- 1

pt= matrix(NA, 7, 7)
pt[lower.tri(pt, diag=FALSE)] <- Dist_age[,1]
pt <- pt*p

colnames(p) <- rownames(p) <- colnames(pt) <- rownames(pt) <- c("Visual","Sensory-Motor","Dorsal Attention",
                                                                "Ventral Attention","Limbic","Fronto-Parietal",
                                                                "DMN")

library(circlize)
df <- melt(pt)
df <- df[complete.cases(df),]

## map values to colours
colors <- scico(length(df$value),palette = "oslo", direction = 1)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.05), start.degree = 90, gap.degree = 4)

pdf('./Figures/subplots/Figure_2D2.pdf', height = 5 , width = 5)
chordDiagram(df, annotationTrack = c("grid","name"), transparency = 0.3,
             col = colors,
             #col = ifelse(df$value > 0, "red", "green"),
             order=cmap_yeo$network, grid.col = setNames(cmap_yeo$color, cmap_yeo$network))
dev.off()

## separate plot for the legend since circlize doesn't have that option
plotdefault2 <- data.frame(freq = seq(min(df$value),max(df$value), length.out=length(df$value)),
                           y = as.factor(1))  

pdf('./Figures/subplots/Figure_2D2_cbar.pdf',height = 8, width = 1)
ggplot(plotdefault2,aes(freq,y)) +
  geom_tile(aes(fill=freq)) + 
  scale_fill_scico(palette = "oslo", direction = 1) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  coord_flip()
dev.off()







#### if we also want to plot non-significant ones in a different color
df_non <- melt(pt)
df_non <- df_non[!complete.cases(df_non),]
df_non$value <- 0.1

colors_non <- rep("black",nrow(df_non))

df_full <- rbind(df,df_non)
color_full <- unlist(list(colors,colors_non))

par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.05), start.degree = 90, gap.degree = 4)
chordDiagram(df_full, annotationTrack = c("grid","name"), transparency = 0.3,
             col = 'red',
             #col = ifelse(df$value > 0, "red", colors),
             order=cmap_yeo$network, grid.col = setNames(cmap_yeo$color, cmap_yeo$network))
