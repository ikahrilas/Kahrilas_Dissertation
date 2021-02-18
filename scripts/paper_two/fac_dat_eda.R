# load packages
library(tidyverse)
library(here)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(patchwork)

# read in factor score data
temp_score_dat <- read_csv(here("data", "paper_two", "temp_fac_score_dat.csv"))
temp_spat_score_dat <- read_csv(here("data", "paper_two", "temp_spat_fac_score_dat.csv"))

# read in quesionnaire data
per_sr_dat <- read_csv(here("data", "paper_two", "archive", "data_from_PER_R_project", "created_data", "per_measures.csv"))

# merge the data sets and omit variables
dat <- full_join(temp_score_dat, per_sr_dat, by = c("pid", "block")) %>%
  select(pid, age, sex, Race, block:valence, anticipating:masq_aa, tmms_repair:depression)

# define block as a factor variable
dat$block <- factor(dat$block, levels = c("Neg_Inc", "Neg_Dec", "Neg_Watch",
                                          "Neu_Watch",
                                          "Pos_Watch", "Pos_Dec", "Pos_Inc"))

# observe histograms for each of the factors
comps <- c("RC2", "RC3", "RC5", "RC11", "RC12", "RC12")

# define electrodes
rc2_elec <- c("A29", "B26")
rc3_elec <- c("A29", "B26", "A26", "B23",
              "B28", "A30", "B27", "A25", "B22")
rc5_elec <- c("A29", "B26")
rc11_elec <- c("A29", "B26")
rc12_elec <- c("B21", "B28")
pos_rc12_elec <- c("A29", "B26")

## list for function
elec_list <- list(rc2_elec,
                  rc3_elec,
                  rc5_elec,
                  rc11_elec,
                  rc12_elec,
                  pos_rc12_elec)

map2(comps, elec_list, ~ {
  dat %>%
    filter(elec %in% .y) %>%
  ggplot(., aes(.data[[.x]])) +
    geom_histogram(color = "black", fill = "white") +
    facet_wrap(~ block) +
    theme_classic()
})

# observe means and standard deviations
map2(comps, elec_list, ~ {
dat %>%
  filter(elec %in% .y) %>%
  group_by(block) %>%
  summarise(mean(.data[[.x]], na.rm = TRUE),
            sd(.data[[.x]], na.rm = TRUE))
})

# boxplots
map2(comps, elec_list, ~ {
dat %>%
  filter(elec %in% .y) %>%
ggplot(aes(block, .data[[.x]])) +
  geom_boxplot() +
  theme_classic()
})

# all looks fine, so create data set that has factor scores calculated based on appropriate electrodes
# and save the file

fac_score_dat <- dat %>%
  filter(elec %in% rc2_elec) %>%
  group_by(pid, block) %>%
  summarise(RC2 = mean(RC2, na.rm = TRUE)) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc3_elec) %>%
      group_by(pid, block) %>%
      summarise(RC3 = mean(RC3, na.rm = TRUE))
            ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc5_elec) %>%
      group_by(pid, block) %>%
      summarise(RC5 = mean(RC5, na.rm = TRUE))
            ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc11_elec) %>%
      group_by(pid, block) %>%
      summarise(RC11 = mean(RC11, na.rm = TRUE))
      ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc12_elec) %>%
      group_by(pid, block) %>%
      summarise(RC12 = mean(RC12, na.rm = TRUE))
  ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% pos_rc12_elec) %>%
      group_by(pid, block) %>%
      summarise(pos_RC12 = mean(RC12, na.rm = TRUE))
  ) %>%
  select(pid...1, block...2, RC2, RC3, RC5, RC11, RC12, pos_RC12) %>%
  rename("pid" = pid...1,
         "block" = block...2) %>%
  full_join(per_sr_dat, by = c("pid", "block")) %>%
  select(pid:valence, ethnicity, Race, anticipating:masq_aa, tmms_repair:depression)


write_csv(x = fac_score_dat,
          file = paste0("data/paper_two/created_data/temp_fac_score_dat_analyses", Sys.Date(), ".csv"))


glimpse(fac_score_dat)

fac_score_dat_long <- fac_score_dat %>%
  pivot_longer(cols = RC2:pos_RC12,
               names_to = "component",
               values_to = "fac_score") %>%
  relocate(pid, block, component, fac_score, everything())

fac_score_dat_long$component <- factor(fac_score_dat_long$component,
                                       levels = c("RC5", "RC11", "RC12", "pos_RC12", "RC3", "RC2"))

levels(fac_score_dat_long$component) <- c("RC5", "RC11", "neg_RC12", "RC12", "RC3", "RC2")

levels(fac_score_dat_long$component) <- c("125 ms Peak",
                                          "170 ms Peak",
                                          "250 ms Negative Peak",
                                          "250 ms Positive Peak",
                                          "375 ms Peak",
                                          "800 ms Peak")

# envisioning 3 grouped box plots, one with watch conditions, one with positive conditions
# and one with negative conditions, with simple topo plots with electrode regions highlighted.
# use this stackoverflow post: https://stackoverflow.com/questions/29263046/how-to-draw-the-boxplot-with-significant-level
# for significance bars
fac_score_dat_long$block <- as.factor(fac_score_dat_long$block)
levels(fac_score_dat_long$block) <- c("Negative Decrease",
                                      "Negative Increase",
                                      "Negative Watch",
                                      "Neutral Watch",
                                      "Positive Decrease",
                                      "Positive Increase",
                                      "Positive Watch")

watch_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Watch")) %>%
  mutate(facet = 1,
         block = factor(block, levels = c("Negative Watch", "Neutral Watch", "Positive Watch")))
pos_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Positive")) %>%
  mutate(facet = 2,
         block = factor(block, levels = c("Positive Decrease", "Positive Watch", "Positive Increase")))
neg_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Negative")) %>%
  mutate(facet = 3,
         block = factor(block, levels = c("Negative Decrease", "Negative Watch", "Negative Increase")))



# first plot the neutral conditions
p_1 <-
  ggplot(watch_cases, aes(x = component, y = fac_score)) +
    geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
    scale_fill_manual(values = c(`Negative Watch` = "red",
                                 `Neutral Watch` = "grey",
                                 `Positive Watch` ="blue")) +
    geom_boxplot(aes(group = interaction(block, component)),
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
    facet_wrap(~ facet, ncol = 1) +
    theme_classic() +
    theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
    labs(x = NULL,
         y = expression(paste("Amplitude (",mu,"V)")),
         fill = "Block") +
    annotate(geom = "segment", # 170 ms peak annotations
             x = c(1.75, 2.025),
             xend = c(1.975, 2.25),
             y = 5.5,
             yend = 5.5,
             color = "black") +
    annotate(geom = "segment",
             x = c(1.75, 1.975, 2.025, 2.25),
             xend = c(1.75, 1.975, 2.025, 2.25),
             y = 5.5,
             yend = c(5, 4, 4, 5),
             color = "black") +
    annotate(geom = "text",
             x = c(1.8625, 2.1375),
             y = 5.6,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 250 ms negative peak annotations
             x = c(2.75),
             xend = c(3),
             y = 2.5,
             yend = 2.5,
             color = "black") +
    annotate(geom = "segment",
             x = c(2.75, 3),
             xend = c(2.75, 3),
             y = 2.5,
             yend = c(2.25, 1.5),
             color = "black") +
    annotate(geom = "text",
             x = c(2.875),
             y = 2.6,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 250 ms positive peak annotations
             x = c(3.75),
             xend = c(4),
             y = 3.9,
             yend = 3.9,
             color = "black") +
    annotate(geom = "segment",
             x = c(3.75, 4),
             xend = c(3.75, 4),
             y = 3.9,
             yend = c(3.6, 3.6),
             color = "black") +
    annotate(geom = "text",
             x = c(3.875),
             y = 4,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 375 ms peak annotations
             x = c(4.75, 5.025),
             xend = c(4.975, 5.25),
             y = 3.35,
             yend = 3.35,
             color = "black") +
    annotate(geom = "segment",
             x = c(4.75, 4.975, 5.025, 5.25),
             xend = c(4.75, 4.975, 5.025, 5.25),
             y = c(3.35, 3.35, 3.35, 3.35),
             yend = c(3.15, 1.8, 1.8, 2)) +
    annotate(geom = "text",
             x = c(4.8625, 5.1375),
             y = 3.45,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 800 ms peak annotations
             x = c(5.75, 6.025),
             xend = c(5.975, 6.25),
             y = 5,
             yend = 5,
             color = "black") +
    annotate(geom = "segment",
             x = c(5.75, 5.975, 6.025, 6.25),
             xend = c(5.75, 5.975, 6.025, 6.25),
             y = c(5, 5, 5, 5),
             yend = c(4.6, 2.85, 2.85, 3.6)) +
    annotate(geom = "text",
             x = c(5.8625, 6.1375),
             y = 5.1,
             label = "*",
             size = 5)

p_2 <-
  ggplot(pos_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Positive Decrease` = "cadetblue1",
                               `Positive Watch` = "blue",
                               `Positive Increase` = "purple")) +
  geom_boxplot(aes(group = interaction(block, component)),
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
  ) +
  labs(x = NULL,
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = c(3),
           xend = c(3.25),
           y = 1.75,
           yend = 1.75,
           color = "black") +
  annotate(geom = "segment",
           x = c(3, 3.25),
           xend = c(3, 3.25),
           y = 1.75,
           yend = c(1.5, 1.5),
           color = "black") +
  annotate(geom = "text",
           x = c(3.125),
           y = 1.85,
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 375 ms peak annotations
           x = c(5),
           xend = c(5.25),
           y = 2.65,
           yend = 2.65,
           color = "black") +
  annotate(geom = "segment",
           x = c(5, 5.25),
           xend = c(5, 5.25),
           y = c(2.65, 2.65),
           yend = c(2.05, 2.3)) +
  annotate(geom = "text",
           x = c(5.1275),
           y = 2.75,
           label = "*",
           size = 5) +

p_3 <-
  ggplot(neg_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Negative Decrease` = "coral",
                               `Negative Watch` = "red",
                               `Negative Increase` = "magenta")) +
  geom_boxplot(aes(group = interaction(block, component)),
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_text(size = 16)
      ) +
  labs(x = "Component",
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  annotate(geom = "segment", # 125 ms peak annotations
           x = c(1),
           xend = c(1.25),
           y = 5.85,
           yend = 5.85,
           color = "black") +
  annotate(geom = "segment",
           x = c(1, 1.25),
           xend = c(1, 1.25),
           y = 5.85,
           yend = c(4.4, 5.5),
           color = "black") +
  annotate(geom = "text",
           x = c(1.125),
           y = 5.95,
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = c(3),
           xend = c(3.25),
           y = 2.65,
           yend = 2.65,
           color = "black") +
  annotate(geom = "segment",
           x = c(3, 3.25),
           xend = c(3, 3.25),
           y = 2.65,
           yend = c(2.3, 1.0),
           color = "black") +
  annotate(geom = "text",
           x = c(3.125),
           y = 2.75,
           label = "*",
           size = 5)

grid.arrange(p_1, p_2, p_3)

p_1 / p_2 / p_3









ggplot(fac_score_dat_long, aes(x = component, y = fac_score, fill = block)) +
  geom_violin(width = 1.4,
              position = position_dodge(width = 1)) +
  geom_boxplot(aes(group = interaction(block, component)),
               width = 0.1,
               color="grey",
               alpha = 0.2,
               position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum()


ggplot(fac_score_dat_long, aes(x = component, y = fac_score, fill = block))+
  geom_violin(position = position_dodge(width = 1)) +
  geom_boxplot(aes(col = block), fill = "white",
               position = position_dodge(width = 1), width = 0.3, outlier.shape = NA)
  geom_boxplot(position = position_dodge(width = 1), alpha = 0, width = 0.3)
