# violin plots for paper 3
## read in packages
library(tidyverse)
library(gridExtra)
library(patchwork)
library(eegUtils)
library(here)
library(grid)
library(lmerTest)
library(emmeans)
library(here)
library(performance)
## load in data
# read in data
dat <- read_csv(here("data", "paper_three", "dat_for_analyses_2021-04-02.csv"))

# RC2 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC2_block_mod <- lmer(RC2 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC2_block_mod)
## check results
summary(RC2_block_mod)
## all contrasts
emmeans(RC2_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "Tukey")

# RC3 analyses
## relevel block factor so that neutral watch is the reference
RC3_block_mod <- lmer(RC3 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC3_block_mod)
## check results
summary(RC3_block_mod)
## all contrasts
emmeans(RC3_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "Tukey")


RC3_block_aov <- anova(RC3_block_mod) # significant anova
df_num_rc3 <- RC3_block_aov$NumDF
df_den_rc3 <- sprintf("%.2f", RC3_block_aov$DenDF)
f_rc3 <- sprintf("%.2f", RC2_block_aov[["F value"]])

# RC5 analyses
## relevel block factor so that neutral watch is the reference
RC5_block_mod <- lmer(RC5 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC5_block_mod)
## check results
summary(RC5_block_mod)
## all contrasts
emmeans(RC5_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "Tukey")

# RC8 analyses
## relevel block factor so that neutral watch is the reference
RC8_block_mod <- lmer(nRC8 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC8_block_mod)
## check results
summary(RC8_block_mod)
## all contrasts
emmeans(RC8_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "Tukey")

# plot it
## plotting dataset
# read in data
dat <- read_csv(here("data", "paper_three", "dat_for_analyses_2021-04-02.csv")) %>%
  select(block, group, RC5, RC2, RC3, nRC8) %>%
  pivot_longer(cols = c(RC5, RC2, RC3, nRC8),
               names_to = "comp",
               values_to = "mv") %>%
  mutate(comp = factor(comp,
                       levels = c("RC5",
                                  "nRC8",
                                  "RC2",
                                  "RC3"),
                       labels = c("134 ms Comp",
                                  "257 ms Comp",
                                  "371 ms Comp",
                                  "736 ms Comp")),
         block = factor(block,
                        levels = c("Neg_Watch",
                                   "Neu_Watch",
                                   "Pos_Watch"),
                        labels = c("Negative",
                                   "Neutral",
                                   "Positive")))

# plotting code
ggplot(dat, aes(x = comp, y = mv)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Negative` = "purple",
                               `Neutral` = "grey",
                               `Positive` ="blue")) +
  geom_boxplot(aes(group = interaction(block, comp)), outlier.size = 2,
               width = 0.25, fill = "white", position = position_dodge(width = .75)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 16)
  ) +
  labs(x = "Component",
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  annotate(geom = "segment", # 134 ms peak annotations
           x = c(0.75),
           xend = c(1),
           y = 6.75,
           yend = 6.75,
           color = "black") +
  annotate(geom = "segment",
           x = c(0.75, 1),
           xend = c(0.75, 1),
           y = 6.75,
           yend = c(6.3, 5.1),
           color = "black") +
  annotate(geom = "text",
           x = c(0.875),
           y = c(6.9),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 257 ms negative peak annotations
           x = c(1.75, 2.025, 1.75),
           xend = c(1.975, 2.25, 2.25),
           y = c(2.6, 2.6, -4),
           yend = c(2.6, 2.6, -4),
           color = "black") +
  annotate(geom = "segment",
           x = c(1.75, 1.975, 2.025, 2.25, 1.75, 2.25),
           xend = c(1.75, 1.975, 2.025, 2.25, 1.75, 2.25),
           y = c(2.6, 2.6, 2.6, 2.6, -4, -4),
           yend = c(1.9, 1.8, 1.8, 2.3, -3.8, -3.4),
           color = "black") +
  annotate(geom = "text",
           x = c(1.875, 2, 2.125),
           y = c(2.7, -4.5, 2.7),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 371 ms positive peak annotations
           x = c(2.75, 3.025, 2.75),
           xend = c(2.975, 3.25, 3.25),
           y = c(4.85, 4.85, -1.5),
           yend = c(4.85, 4.85, -1.5),
           color = "black")+
  annotate(geom = "segment",
           x = c(2.75, 2.975, 3.025, 3.25, 2.75, 3.25),
           xend = c(2.75, 2.975, 3.025, 3.25, 2.75, 3.25),
           y = c(4.85, 4.85, 4.85, 4.85, -1.5, -1.5),
           yend = c(4.65, 3.5, 3.5, 4.5, -0.2, -0.5),
           color = "black") +
  annotate(geom = "text",
           x = c(2.8625, 3.1375, 3),
           y = c(5.05, 5.05, -1.9),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 736 ms positive peak annotations
           x = c(3.75, 4.025, 3.75),
           xend = c(3.975, 4.25, 4.25),
           y = c(3.65, 3.65, -1.3),
           yend = c(3.65, 3.65, -1.3),
           color = "black") +
  annotate(geom = "segment",
           x = c(3.75, 3.975, 4.025, 4.25, 3.75, 4.25),
           xend = c(3.75, 3.975, 4.025, 4.25, 3.75, 4.25),
           y = c(3.65, 3.65, 3.65, 3.65, -1.3, -1.3),
           yend = c(3.45, 2, 2, 2.3, -0.5, -0.6),
           color = "black") +
  annotate(geom = "text",
           x = c(3.8625, 4.1375, 4),
           y = c(3.75, 3.75, -1.7),
           label = "*",
           size = 5)
