#### Script to visualize the residuals from sex and motion regression on various
#### gradient descriptors
setwd("~/GitHub/yesweCanCAM")

library(ggplot2)
library(ggpubr)
library(lme4)
library(viridis)
library(ggpmisc)
library(scico)

G1 <- read.csv('./Data/bimodality_values_gradient1.csv')
G2 <- read.csv('./Data/bimodality_values_gradient2.csv')
G3 <- read.csv('./Data/bimodality_values_gradient3.csv')

## regress out confounds
G1_GLM <- data.frame(G1$age)
G1_GLM$r1 <- lm(r1 ~  sex + motion, data = G1)$residuals
G1_GLM$dip <- lm(dip ~ sex + motion, data = G1)$residuals

G2_GLM <- data.frame(G2$age)
G2_GLM$r1 <- lm(r1 ~  sex + motion, data = G2)$residuals
G2_GLM$dip <- lm(dip ~ sex + motion, data = G2)$residuals

G3_GLM <- data.frame(G3$age)
G3_GLM$r1 <- lm(r1 ~  sex + motion, data = G3)$residuals
G3_GLM$dip <- lm(dip ~ sex + motion, data = G3)$residuals

## reshape and plot the correlations between age and residuals
library(reshape2)
G1_GLM <- melt(G1_GLM, id.vars = "G1.age")
# standardize fill variable
G1_GLM <- group_by(G1_GLM,variable)
G1_GLM <- mutate(G1_GLM, z1 = (value-mean(value))/sd(value))

P1 <- ggplot(G1_GLM, aes(x=G1.age, y=value, fill = z1)) +
  geom_point(stat = "identity",shape = 21, size = 1.5, show.legend = FALSE)  +
  geom_smooth(method=lm) +
  scale_fill_scico(palette = "lajolla") +
  xlab("Age") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G1_GLM$variable)))

D1 <- ggplot(G1_GLM, aes(value, fill = variable)) +
  geom_density(alpha = 0.6)  +
  #scale_fill_viridis(option = "C", discrete = TRUE) +
  scale_fill_manual(values = scico(2, begin = 0.4, palette = "lajolla")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G1_GLM$variable)))

G2_GLM <- melt(G2_GLM, id.vars = "G2.age")
G2_GLM <- group_by(G2_GLM,variable)
G2_GLM <- mutate(G2_GLM, z1 = (value-mean(value))/sd(value))
P2 <- ggplot(G2_GLM, aes(x=G2.age, y=value, fill = z1)) +
  geom_point(stat = "identity",shape = 21, size = 1.5, show.legend = FALSE)  +
  geom_smooth(method=lm) +
  scale_fill_scico(palette = "lajolla") +
  xlab("Age") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G2_GLM$variable)))

D2 <- ggplot(G2_GLM, aes(value, fill = variable)) +
  geom_density(alpha = 0.6)  +
  scale_fill_manual(values = scico(2, begin = 0.4, palette = "lajolla")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G2_GLM$variable)))

G3_GLM <- melt(G3_GLM, id.vars = "G3.age")
G3_GLM <- group_by(G3_GLM,variable)
G3_GLM <- mutate(G3_GLM, z1 = (value-mean(value))/sd(value))
P3 <- ggplot(G3_GLM, aes(x=G3.age, y=value, fill = z1)) +
  geom_point(stat = "identity",shape = 21, size = 1.5, show.legend = FALSE)  +
  geom_smooth(method=lm) +
  scale_fill_scico(palette = "lajolla") +
  xlab("Age") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G3_GLM$variable)))

D3 <- ggplot(G3_GLM, aes(value, fill = variable)) +
  geom_density(alpha = 0.6)  +
  scale_fill_manual(values = scico(2, begin = 0.4, palette = "lajolla")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free", nrow = length(unique(G3_GLM$variable)))

## combined the plots of each gradient into one
finalFigure <- ggarrange(P1,D1,P2,D2,P3,D3, nrow = 1, ncol = 6,
                         widths = c(3,1.5,3,1.5,3,1.5),
                         labels = c("G1","G1_Dens","G2","G2_Dens","G3","G3_Dens"))

pdf('./Figures/supplemental/SFigure_1.pdf',height = 8, width = 12)
plot(finalFigure)
dev.off()
