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
library(R.matlab)

setwd("~/GitHub/GradientDispersion")

# load some data
df <- read.csv('./Data/camcan_lms.csv', header = TRUE)
subs <- read.csv('./Data/sublist514.csv', header = FALSE)
df$CCID <- subs$V1
demos <- read.csv('./Data/camcan_demos.csv')

df <- merge(df,demos,by = 'CCID')
df$agebin <- cut(df$age, breaks = c(18,28,38,48,58,68,78,90))

hand <- read.table('~/Downloads/participants.tsv', header = T)

df <- merge(df,hand, by = 'participant_id', all.x = TRUE)

library(tableone)
library(boot)
library(htmlTable)


table1::table1(~ factor(sex) + age.x + motion + hand +factor(qual)  + factor(SC_r) + factor (SEG_r)| agebin, data=df)

holdout <- read.table('./Data/holdout_subjects.csv')
colnames(holdout) <- 'participant_id'
demos$participant_id <- paste0("sub-",demos$CCID)

dh <- merge(holdout,demos, by = 'participant_id')
dh <- merge(dh,hand, by = 'participant_id', all.x = T)
dh$agebin <- cut(dh$age, breaks = c(18,28,38,48,58,68,78,90))


motion_hold <- matrix(NA,106,1)
for (i in 1:length(unique(dh$participant_id))){
  sub <- as.character(dh$participant_id)[i]
  motion_hold[i] <- mean(read.table(paste0('~/Downloads/vol2surf_motion_refrms/',sub,'_task-Rest_bold_metric_REFRMS.1D'))$V1)
  
}

dh$motion <- motion_hold
table1::table1(~ factor(gender_code) + age  + motion +hand +factor(qual)  + factor(SC_r) + factor (SEG_r)| agebin, data=dh)



df.select <- df[,c("age.x", "sex", "motion", "cattell", "Fronto.Parietal", "DMN", "Visual", "Sensory.Motor", "Dorsal.Attention_", "Ventral.Attention", "Limbic")]
df.select <- melt(df.select, measure.vars = c("Fronto.Parietal", "DMN", "Visual", "Sensory.Motor", "Dorsal.Attention_", "Ventral.Attention", "Limbic"))
to_keep <- apply(!is.na(df.select), 1, all)
df.select <-  df.select[to_keep,]


library(ggplot2)
library(viridis)

ggplot(df.select, aes(y = value, x = variable, fill = variable)) +
  geom_jitter(aes(fill = variable),shape = 21, alpha = .5) +
  geom_boxplot(alpha = .7) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  facet_wrap(~agebin, ncol =7, nrow = 1, scales = "free_x") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(fill=FALSE) +
  theme_minimal()


library(scico)
fake_scico <- scico(3, palette = "vik")
ggstatsplot::ggcorrmat(data=df.select,type="robust",matrix.type = "upper", colors = c(fake_scico[1],fake_scico[2], fake_scico[3]), p.adjust.method = "BH", sig.level = 0.01)
