# Repeated measures ANOVA

library(tidyverse)
library(lmerTest)
library(here)

dat <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv")

mod_front <- lmer(LPP_front ~ block + (1|pid), data = dat)
av_front <- anova(mod_front)
df_num_front <- av_front$NumDF
df_den_front <- sprintf("%.2f", av_front$DenDF)
f_front <- sprintf("%.2f", av_front[["F value"]])

mod_lpp <- lmer(LPP ~ block + (1|pid), data = dat)
av_lpp <- anova(mod_lpp)
df_num_lpp <- av_lpp$NumDF
df_den_lpp <- sprintf("%.2f", av_lpp$DenDF)
f_lpp <- sprintf("%.2f", av_lpp[["F value"]])

mod_epn <- lmer(EPN ~ block + (1|pid), data = dat)
av_epn <- anova(mod_epn)
df_num_epn <- av_epn$NumDF
df_den_epn <- round(av_epn$DenDF, digits = 0)
f_epn <- sprintf("%.2f", av_epn[["F value"]])

mod_n170 <- lmer(N170 ~ block + (1|pid), data = dat)
av_n170 <- anova(mod_n170)
df_num_n170 <- av_n170$NumDF
df_den_n170 <- round(av_n170$DenDF, digits = 0)
f_n170 <- sprintf("%.2f", av_n170[["F value"]])

mod_ar <- lmer(arousal ~ block + (1|pid), data = dat)
av_ar <- anova(mod_ar)
df_num_ar <- av_ar$NumDF
df_den_ar <- round(av_ar$DenDF, digits = 0)
f_ar <- sprintf("%.2f", av_ar[["F value"]])

mod_val <- lmer(valence ~ block + (1|pid), data = dat)
av_val <- anova(mod_val)
df_num_val <- av_val$NumDF
df_den_val <- round(av_val$DenDF, digits = 0)
f_val <- sprintf("%.2f", av_val[["F value"]])

mod_diff <- lmer(difficulty ~ block + (1|pid), data = dat)
av_diff <- anova(mod_diff)
df_num_diff <- av_diff$NumDF
df_den_diff <- round(av_diff$DenDF, digits = 0)
f_diff <- sprintf("%.2f", av_diff[["F value"]])

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_repeated_measures_ANOVA_analysis-data", ".RData"))
