# Principal components analysis

# load packages
library(tidyverse)
library(here)
library(factoextra)
library(sjstats)
library(psych)

# read in average referenced data set
## using average reference as per Dein (2012)
dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
## mastoid referenced data
dat_mast <- read_csv(here("data", "paper_two", "created_data", "erp_mast.csv"))

# rearrange data for temporal PCA with timepoints as columns
## average reference
dat_temp <- dat %>%
  pivot_longer(A1:EXG2, names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)
## mastoid reference
dat_mast_temp <- dat_mast %>%
  pivot_longer(A1:EXG2, names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)
## just n170/epn electrodes and cut off at 1500 ms
dat_1500_epn_n170_elec <- dat %>%
  filter(ms < 1500) %>%
  select(prop_trials:pid, A29, B26) %>%
  pivot_longer(c(A29, B26), names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)
### now with 1000 ms cutoff
dat_1000_epn_n170_elec <- dat %>%
  filter(ms < 1000) %>%
  select(prop_trials:pid, A29, B26) %>%
  pivot_longer(c(A29, B26), names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)
### now with 500 ms cutoff
dat_500_epn_n170_elec <- dat %>%
  filter(ms < 500) %>%
  select(prop_trials:pid, A29, B26) %>%
  pivot_longer(c(A29, B26), names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)

# pca
## whole thing avr
dat_pca <- prcomp(na.omit(dat_temp[,-c(1:5)]), center = TRUE, scale. = FALSE)
scree_all <- fviz_eig(dat_pca)
ggsave(filename = here("images", "paper_2", "pca_images", "scree_all.png"), plot = scree_all, device = "png")
## 1500 ms cut off with limited n170/epn electrodes
dat_pca_1500_n170_epn_elec <- prcomp(na.omit(dat_1500_epn_n170_elec[,-c(1:5)]), center = TRUE, scale. = FALSE)
scree_1500_n170_epn_elec <- fviz_eig(dat_pca_1500_n170_epn_elec)
ggsave(filename = here("images", "paper_2", "pca_images", "scree_1500_n170_epn_elec.png"), plot = scree_1500_n170_epn_elec, device = "png")
### 1000 ms cut off
dat_pca_1000_n170_epn_elec <- prcomp(na.omit(dat_1000_epn_n170_elec[,-c(1:5)]), center = TRUE, scale. = FALSE)
scree_1000_n170_epn_elec <- fviz_eig(dat_pca_1000_n170_epn_elec)
ggsave(filename = here("images", "paper_2", "pca_images", "scree_1000_n170_epn_elec.png"), plot = scree_1000_n170_epn_elec, device = "png")
### 500 ms cut off
dat_pca_500_n170_epn_elec <- prcomp(na.omit(dat_500_epn_n170_elec[,-c(1:5)]), center = TRUE, scale. = FALSE)
scree_500_n170_epn_elec <- fviz_eig(dat_pca_500_n170_epn_elec)
ggsave(filename = here("images", "paper_2", "pca_images", "scree_500_n170_epn_elec.png"), plot = scree_500_n170_epn_elec, device = "png")

# pca with promax rotation
## on average referenced data with promax rotation
dat_pca_promax <- principal(dat_temp[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
## without rotation
dat_pca_nr <- principal(dat_temp[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", missing = TRUE)
# pca with promax rotation on mastoid referenced data
dat_mast_pca_promax <- principal(dat_mast_temp[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
## without rotation
dat_mast_pca_nr <- principal(dat_mast_temp[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", missing = TRUE)
## with 1500 ms cutoff and just n170 and epn electrodes and promax
dat_mast_pca_promax_1500_epn_n170_elec <- principal(dat_1500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
### with 1000 ms cutoff and just n170 and epn electrodes and promax
dat_mast_pca_promax_1000_epn_n170_elec <- principal(dat_1000_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
### with 500 ms cutoff
dat_pca_promax_500_epn_n170_elec <- principal(dat_500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)

# make visualization
## avr
cov_loadings_mat <- matrix(dat_pca_promax$loadings, nrow = 1638)
ms_vec <- dat %>%
  filter(pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_df <- data.frame(cov_loadings_mat)
names(cov_loadings_df) <- paste0("RC", c(1:10))
cov_loadings_df <- cov_loadings_df %>%
  mutate(ms = ms_vec)
## avr no rotation
cov_loadings_mat_nr <- matrix(dat_pca_nr$loadings, nrow = 1638)
ms_vec <- dat %>%
  filter(pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mat_nr <- data.frame(cov_loadings_mat_nr)
names(cov_loadings_mat_nr) <- paste0("RC", c(1:10))
cov_loadings_mat_nr <- cov_loadings_mat_nr %>%
  mutate(ms = ms_vec)
## mastoid
cov_loadings_mat_mast <- matrix(dat_mast_pca_promax$loadings, nrow = 1638)
ms_vec <- dat_mast %>%
  filter(pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mast_df <- data.frame(cov_loadings_mat_mast)
names(cov_loadings_mast_df) <- paste0("RC", c(1:10))
cov_loadings_mast_df <- cov_loadings_mast_df %>%
  mutate(ms = ms_vec)
## 1500 ms cutoff with n170 and epn electrodes
cov_loadings_mat_trimmed <- matrix(dat_mast_pca_promax_1500_epn_n170_elec$loadings, nrow = 870)
ms_vec <- dat %>%
  filter(ms < 1500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_trimmed_df <- data.frame(cov_loadings_mat_trimmed)
names(cov_loadings_trimmed_df) <- paste0("RC", c(1:10))
cov_loadings_trimmed_df <- cov_loadings_trimmed_df %>%
  mutate(ms = ms_vec)
### 1000 ms cut off
cov_loadings_mat_trimmed_1000 <- matrix(dat_mast_pca_promax_1000_epn_n170_elec$loadings, nrow = 614)
ms_vec <- dat %>%
  filter(ms < 1000,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_trimmed_1000_df <- data.frame(cov_loadings_mat_trimmed_1000)
names(cov_loadings_trimmed_1000_df) <- paste0("RC", c(1:10))
cov_loadings_trimmed_1000_df <- cov_loadings_trimmed_1000_df %>%
  mutate(ms = ms_vec)
### 500 ms cutoff
cov_loadings_mat_trimmed_500 <- matrix(dat_pca_promax_500_epn_n170_elec$loadings, nrow = 359)
ms_vec <- dat %>%
  filter(ms < 500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_trimmed_500_df <- data.frame(cov_loadings_mat_trimmed_500)
names(cov_loadings_trimmed_500_df) <- paste0("RC", c(1:10))
cov_loadings_trimmed_500_df <- cov_loadings_trimmed_500_df %>%
  mutate(ms = ms_vec)

# long form so components are factor variables and make visualizations
## avr
cov_loadings_long <- pivot_longer(cov_loadings_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long$component <- factor(cov_loadings_long$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
avr_loadings_all <- ggplot(cov_loadings_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "avr_loadings_all.png"), plot = avr_loadings_all, device = "png")
## avr no rotation
cov_loadings_long_nr <- pivot_longer(cov_loadings_mat_nr, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long_nr$component <- factor(cov_loadings_long_nr$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
avr_loadings_all_nr <- ggplot(cov_loadings_long_nr, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "avr_loadings_all_nr.png"), plot = avr_loadings_all_nr, device = "png")
## mast
cov_loadings_long_mast <- pivot_longer(cov_loadings_mast_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long_mast$component <- factor(cov_loadings_long_mast$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
mast_loadings_all <- ggplot(cov_loadings_long_mast, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "mast_loadings_all.png"), plot = mast_loadings_all, device = "png")
## 1500 ms cutoff and epn/n170 electrodes only
cov_loadings_trimmed_df_long <- pivot_longer(cov_loadings_trimmed_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_trimmed_df_long$component <- factor(cov_loadings_trimmed_df_long$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
trimmed_loadings <- ggplot(cov_loadings_trimmed_df_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "trimmed_loadings.png"), plot = trimmed_loadings, device = "png", width = 10)
### 1000 ms cutoff
cov_loadings_trimmed_df_1000_long <- pivot_longer(cov_loadings_trimmed_1000_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_trimmed_df_1000_long$component <- factor(cov_loadings_trimmed_df_1000_long$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
trimmed_loadings_1000 <- ggplot(cov_loadings_trimmed_df_1000_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "trimmed_loadings_1000.png"), plot = trimmed_loadings_1000, device = "png", width = 14)
### 500 ms cutoff
cov_loadings_trimmed_df_500_long <- pivot_longer(cov_loadings_trimmed_500_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_trimmed_df_500_long$component <- factor(cov_loadings_trimmed_df_500_long$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
trimmed_loadings_500 <- ggplot(cov_loadings_trimmed_df_500_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "trimmed_loadings_500.png"), plot = trimmed_loadings_500, device = "png", width = 14)
