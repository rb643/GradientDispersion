setwd("~/GitHub/yesweCanCAM")

library(ggplot2)
library(ggpubr)
library(lme4)
library(viridis)
library(ggpmisc)
library(reshape2)
library(ggridges)
library(scico)

G1 <- read.csv("./Data/deciles_gradient1.txt", header = FALSE)
colnames(G1) <- c("18-27", "28-37", "38-47", "48-57", "58-67", "68-77", "78-87")
G <- melt(G1)
colnames(G)[2] <- "G1"

G2 <- read.csv("./Data/deciles_gradient2.txt", header = FALSE)
colnames(G2) <- c("18-27", "28-37", "38-47", "48-57", "58-67", "68-77", "78-87")
G$G2 <- melt(G2)$value

G3 <- read.csv("./Data/deciles_gradient3.txt", header = FALSE)
colnames(G3) <- c("18-27", "28-37", "38-47", "48-57", "58-67", "68-77", "78-87")
G$G3 <- melt(G3)$value

G4 <- melt(G)
colnames(G4) <- c("variable", "gradient", "value")

ggplot(G4, aes(x = value, y = variable, fill = ..x..)) +
  geom_density_ridges_gradient(ascale = 2, gradient_lwd = 1.)  +
  scale_fill_viridis(option = "C", direction = 1) +
  #scale_fill_scico(palette = "imola", direction = 1) +
  scale_y_discrete(expand = c(0.2, 0)) +   # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0))  +
  facet_wrap(~gradient, scales = "free_x", shrink = FALSE) +
  theme_ridges() +
  labs(y = "Age bin") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(3, "lines"))
