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
## mastoid reference with 1500 ms cutoff and Pz electrode
dat_mast_temp <- dat_mast %>%
  filter(ms < 1500) %>%
  select(prop_trials:pid, B28) %>%
  pivot_longer(B28, names_to = "elec") %>%
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


# pca analyses with principal function
## on average referenced data with promax rotation
dat_pca_promax <- principal(dat_temp[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
## without rotation
dat_pca_nr <- principal(dat_temp[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", missing = TRUE)

## pca with promax rotation on mastoid referenced data
dat_mast_pca_promax <- principal(dat_mast_temp[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
## without rotation
dat_mast_pca_nr <- principal(dat_mast_temp[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", missing = TRUE)
## varimax rotation
dat_mast_pca_var <- principal(dat_mast_temp[,-c(1:5)], nfactors = 10, rotate = "varimax", cor = "cov", missing = TRUE)

## with 1500 ms cutoff and just n170 and epn electrodes and promax
dat_mast_pca_promax_1500_epn_n170_elec <- principal(dat_1500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
### with 1000 ms cutoff and just n170 and epn electrodes and promax
dat_mast_pca_promax_1000_epn_n170_elec <- principal(dat_1000_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
### with 500 ms cutoff
dat_pca_promax_500_epn_n170_elec <- principal(dat_500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
### 500 ms cutoff with no rotation
dat_pca_nr_500_epn_n170_elec <- principal(dat_500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", covar = TRUE, missing = TRUE)
### varimax rotation
dat_pca_varimax_500_epn_n170_elec <- principal(dat_500_epn_n170_elec[,-c(1:5)], nfactors = 10, rotate = "varimax", cor = "cov", missing = TRUE)

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
## mastoid 1500 ms promax
cov_loadings_mat_mast <- matrix(dat_mast_pca_promax$loadings, nrow = 870)
ms_vec <- dat_mast %>%
  filter(ms < 1500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mast_df <- data.frame(cov_loadings_mat_mast)
names(cov_loadings_mast_df) <- paste0("RC", c(1:10))
cov_loadings_mast_df <- cov_loadings_mast_df %>%
  mutate(ms = ms_vec)

## mastoid 1500 ms nr
cov_loadings_mat_mast_nr <- matrix(dat_mast_pca_nr$loadings, nrow = 870)
ms_vec <- dat_mast %>%
  filter(ms < 1500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mast_df_nr <- data.frame(cov_loadings_mat_mast_nr)
names(cov_loadings_mast_df_nr) <- paste0("RC", c(1:10))
cov_loadings_mast_df_nr <- cov_loadings_mast_df_nr %>%
  mutate(ms = ms_vec)

## mastoid 1500 ms var
cov_loadings_mat_mast_var <- matrix(dat_mast_pca_var$loadings, nrow = 870)
ms_vec <- dat_mast %>%
  filter(ms < 1500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mast_df_var <- data.frame(cov_loadings_mat_mast_var)
names(cov_loadings_mast_df_var) <- paste0("RC", c(1:10))
cov_loadings_mast_df_var <- cov_loadings_mast_df_var %>%
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
### 500 ms nr
cov_loadings_mat_trimmed_500_nr <- matrix(dat_pca_nr_500_epn_n170_elec$loadings, nrow = 359)
ms_vec <- dat %>%
  filter(ms < 500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_trimmed_500_df_nr <- data.frame(cov_loadings_mat_trimmed_500_nr)
names(cov_loadings_trimmed_500_df_nr) <- paste0("RC", c(1:10))
cov_loadings_trimmed_500_df_nr <- cov_loadings_trimmed_500_df_nr %>%
  mutate(ms = ms_vec)
### 500 ms varimax
cov_loadings_mat_trimmed_500_var <- matrix(dat_pca_varimax_500_epn_n170_elec$loadings, nrow = 359)
ms_vec <- dat %>%
  filter(ms < 500,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_df_trimmed_500_var <- data.frame(cov_loadings_mat_trimmed_500_var)
names(cov_loadings_df_trimmed_500_var) <- paste0("RC", c(1:10))
cov_loadings_df_trimmed_500_var <- cov_loadings_df_trimmed_500_var %>%
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

## mast promax
cov_loadings_long_mast <- pivot_longer(cov_loadings_mast_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long_mast$component <- factor(cov_loadings_long_mast$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
mast_loadings_all <- ggplot(cov_loadings_long_mast, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "mast_loadings_promax.png"), plot = mast_loadings_all, device = "png", width = 14)

## mast varimax
cov_loadings_long_mast_var <- pivot_longer(cov_loadings_mast_df_var, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long_mast_var$component <- factor(cov_loadings_long_mast_var$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
mast_loadings_all_var <- ggplot(cov_loadings_long_mast_var, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "mast_loadings_varimax.png"), plot = mast_loadings_all_var, device = "png", width = 14)


## mast nr
cov_loadings_long_mast_nr <- pivot_longer(cov_loadings_mast_df_nr, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_long_mast_nr$component <- factor(cov_loadings_long_mast$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
mast_loadings_all_nr <- ggplot(cov_loadings_long_mast_nr, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "mast_loadings_nr.png"), plot = mast_loadings_all_nr, device = "png", width = 14)


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
### 500 ms cutoff nr
cov_loadings_trimmed_df_500_long_nr <- pivot_longer(cov_loadings_trimmed_500_df_nr, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_trimmed_df_500_long_nr$component <- factor(cov_loadings_trimmed_df_500_long_nr$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
trimmed_loadings_500_nr <- ggplot(cov_loadings_trimmed_df_500_long_nr, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "trimmed_loadings_500_nr.png"), plot = trimmed_loadings_500_nr, device = "png", width = 14)
## varimax
cov_loadings_df_trimmed_500_var_long <- pivot_longer(cov_loadings_df_trimmed_500_var, cols = RC1:RC10, names_to = "component", values_to = "mv")
cov_loadings_df_trimmed_500_var_long$component <- factor(cov_loadings_df_trimmed_500_var_long$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
trimmed_loadings_500_var <- ggplot(cov_loadings_df_trimmed_500_var_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 2)
ggsave(filename = here("images", "paper_2", "pca_images", "trimmed_loadings_500_var.png"), plot = trimmed_loadings_500_var, device = "png", width = 14)


###############################
####Derive Factor Scores######
###############################
# EPN/N170
dat_long_elec_of_int <- dat %>%
  filter(ms < 500) %>%
  mutate(epn_n170_elec = (A29 + B26) / 2) %>%
  select(prop_trials:pid, epn_n170_elec)

epn_n170_loadings <- cov_loadings_trimmed_500_df %>%
  select(ms, RC3, RC4) %>%
  rename(epn_loading = RC4,
         n170_loading = RC3)

dat_w_loadings <- full_join(dat_long_elec_of_int, epn_n170_loadings, by = "ms") %>%
  mutate(epn_factor_score = epn_n170_elec * epn_loading,
         n170_factor_score = epn_n170_elec * n170_loading)
#LPP
lpp_loadings <- cov_loadings_mast_df %>%
  select(ms, RC2, RC7)

lpp_dat_w_loadings <- dat_mast %>%
  filter(ms < 1500) %>%
  select(prop_trials:pid, B28) %>%
  full_join(lpp_loadings, by = "ms") %>%
  mutate(lpp_rc2_score = B28 * RC2,
         lpp_rc7_score = B28 * RC7)

# now plot it to see what it looks like

########################
#######LPP PLOTS#############
########################

###ALL BLOCKS###
lpp_dat_w_loadings %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc2_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

lpp_dat_w_loadings %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc7_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

###WATCH BLOCKS###
lpp_dat_w_loadings %>%
  filter(block %in% c("Pos_Watch", "Neu_Watch", "Neg_Watch")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc2_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

lpp_dat_w_loadings %>%
  filter(block %in% c("Pos_Watch", "Neu_Watch", "Neg_Watch")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc7_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

###POS BLOCKS###
lpp_dat_w_loadings %>%
  filter(block %in% c("Pos_Watch", "Pos_Dec", "Pos_Inc")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc2_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

lpp_dat_w_loadings %>%
  filter(block %in% c("Pos_Watch", "Pos_Dec", "Pos_Inc")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc7_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

###NEG BLOCKS###
lpp_dat_w_loadings %>%
  filter(block %in% c("Neg_Watch", "Neg_Dec", "Neg_Inc")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc2_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

lpp_dat_w_loadings %>%
  filter(block %in% c("Neg_Watch", "Neg_Dec", "Neg_Inc")) %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score),
            lpp_rc7_score = mean(lpp_rc7_score)) %>%
  ggplot(., aes(ms, lpp_rc7_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####LPP plot all blocks with varimax rotation####
lpp_loadings_var <- cov_loadings_mast_df_var %>%
  select(ms, RC2)

lpp_dat_w_loadings_var <- dat_mast %>%
  filter(ms < 1500) %>%
  select(prop_trials:pid, B28) %>%
  full_join(lpp_loadings_var, by = "ms") %>%
  mutate(lpp_rc2_score = B28 * RC2)

lpp_dat_w_loadings_var %>%
  group_by(ms, block) %>%
  summarize(lpp_rc2_score = mean(lpp_rc2_score)) %>%
  ggplot(., aes(ms, lpp_rc2_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

########################
#####EPN/N170##########
########################

####PLOT ALL BLOCKS EPN####
dat_w_loadings %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
ggplot(., aes(ms, epn_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT WATCH BLOCKS EPN###
dat_w_loadings %>%
  filter(block %in% c("Neu_Watch", "Pos_Watch", "Neg_Watch")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, epn_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT NEGATIVE BLOCKS EPN###
dat_w_loadings %>%
  filter(block %in% c("Neg_Dec", "Neg_Inc", "Neg_Watch")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, epn_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT POSITIVE BLOCKS EPN###
dat_w_loadings %>%
  filter(block %in% c("Pos_Dec", "Pos_Inc", "Pos_Watch")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, epn_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT ALL BLOCKS N170####
dat_w_loadings %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, n170_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT WATCH BLOCKS N170####
dat_w_loadings %>%
  filter(block %in% c("Neu_Watch", "Pos_Watch", "Neg_Watch")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, n170_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT NEGATIVE BLOCKS N170####
dat_w_loadings %>%
  filter(block %in% c("Neg_Watch", "Neg_Dec", "Neg_Inc")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, n170_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()

####PLOT POSITIVE BLOCKS N170####
dat_w_loadings %>%
  filter(block %in% c("Pos_Watch", "Pos_Dec", "Pos_Inc")) %>%
  group_by(block, ms) %>%
  summarize(epn_factor_score = mean(epn_factor_score),
            n170_factor_score = mean(n170_factor_score)) %>%
  ggplot(., aes(ms, n170_factor_score)) +
  geom_line(aes(color = block)) +
  theme_classic()


# unused code for conducting PCAs by regulation block
# # PCA for ERPs grouped by regulation block
#
# ## split data frame with 500 ms cut off by block
# dat$block <- as.factor(dat$block)
#
# dat_500 <- dat %>%
#   filter(ms < 500) %>%
#   select(prop_trials:pid, A29, B26) %>%
#   pivot_longer(c(A29, B26), names_to = "elec") %>%
#   pivot_wider(names_from = ms, values_from = value)
#
# dat_list <- split(dat_500, dat_500$block)
# block_names <- list(names(dat_list))
#
# # iterative and get scree plots
# map2(dat_list, names(dat_list), ~ {
# pca <- prcomp(na.omit(.x[,-c(1:5)]), center = TRUE, scale. = FALSE)
# scree_pca <- fviz_eig(pca)
# ggsave(filename = here("images", "paper_2", "pca_images", paste0(.y, "_scree_500_n170_epn_elec.png")), plot = scree_pca, device = "png")
# }
# )
#
# # conduct rotated PCAs make plots
#
# ## make ms vector
# ms_vec <- dat %>%
#   filter(ms < 500,
#          pid == 206201832,
#          block == "Pos_Inc") %>%
#   select(ms) %>%
#   pull()
#
# ## map function for pcas
# promax_list <- map2(dat_list, names(dat_list), ~ {
# ## promax rotation
# pca_promax <- principal(.x[,-c(1:5)], nfactors = 10, rotate = "promax", cor = "cov", missing = TRUE)
#
# ## promax df
# cov_loadings_mat_promax <- matrix(pca_promax$loadings, nrow = 359)
# cov_loadings_mat_promax_df <- data.frame(cov_loadings_mat_promax)
# names(cov_loadings_mat_promax_df) <- paste0("RC", c(1:10))
# cov_loadings_mat_promax_df <- cov_loadings_mat_promax_df %>%
#   mutate(ms = ms_vec)
#
# ## long form
# cov_loadings_mat_promax_df <- pivot_longer(cov_loadings_mat_promax_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
# cov_loadings_mat_promax_df$component <- factor(cov_loadings_mat_promax_df$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
# cov_loadings_mat_promax_df <- cov_loadings_mat_promax_df %>%
#   mutate(block = .y)
# }
# )
# # bind list together
# promax_df <- bind_rows(promax_list)
#
# promax_all_blocks <- ggplot(promax_df, aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# # save image
# ggsave(filename = here("images", "paper_2", "pca_images", "all_blocks_promax_500.png"), plot = promax_all_blocks, device = "png", width = 14)
#
#
# # unrotated all blocks
# nr_list <- map2(dat_list, names(dat_list), ~ {
# ## unrotated
# pca_nr <- principal(.x[,-c(1:5)], nfactors = 10, rotate = "none", cor = "cov", missing = TRUE)
#
# ## nr df
# cov_loadings_mat_nr <- matrix(pca_nr$loadings, nrow = 359)
# cov_loadings_mat_nr_df <- data.frame(cov_loadings_mat_nr)
# names(cov_loadings_mat_nr_df) <- paste0("RC", c(1:10))
# cov_loadings_mat_nr_df <- cov_loadings_mat_nr_df %>%
#   mutate(ms = ms_vec)
#
# ## long form
# cov_loadings_mat_nr_df <- pivot_longer(cov_loadings_mat_nr_df, cols = RC1:RC10, names_to = "component", values_to = "mv")
# cov_loadings_mat_nr_df$component <- factor(cov_loadings_mat_nr_df$component, levels = c("RC1", "RC2", "RC3", "RC4", "RC5", "RC6", "RC7", "RC8", "RC9", "RC10"))
# cov_loadings_mat_nr_df <- cov_loadings_mat_nr_df %>%
#   mutate(block = .y)
# })
#
# # bind list together for nr
# nr_df <- bind_rows(nr_list)
#
# nr_all_blocks <- ggplot(nr_df, aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# # save image
# ggsave(filename = here("images", "paper_2", "pca_images", "all_blocks_nr_500.png"), plot = nr_all_blocks, device = "png", width = 14)
#
# # now, max cov loading plots for each individual block
# nr_df$block <- as.factor(nr_df$block)
# promax_df$block <- as.factor(promax_df$block)
# block_vector <- levels(promax_df$block)
#
# map(block_vector, ~ {
# nr_plot <- nr_df %>%
#   filter(block == .x) %>%
# ggplot(., aes(ms, mv)) +
#   geom_line() +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0(.x, "_nr_500.png")), plot = nr_plot, device = "png", width = 14)
#
# promax_plot <- promax_df %>%
#   filter(block == .x) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line() +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0(.x, "_promax_500.png")), plot = promax_plot, device = "png", width = 14)
# })
#
# # now make factor loading plots for those contditions that you will contrast
# # within the same statistical model
#
# # unrotated
# watch_nr_plots <- nr_df %>%
#   filter(block %in% c("Neu_Watch", "Neg_Watch", "Pos_Watch")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("watch_nr_500.png")), plot = watch_nr_plots, device = "png", width = 14)
#
# pos_nr_plots <- nr_df %>%
#   filter(block %in% c("Pos_Watch", "Pos_Inc", "Pos_Dec")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("pos_nr_500.png")), plot = pos_nr_plots, device = "png", width = 14)
#
# neg_nr_plots <- nr_df %>%
#   filter(block %in% c("Neg_Watch", "Neg_Inc", "Neg_Dec")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("neg_nr_500.png")), plot = neg_nr_plots, device = "png", width = 14)
#
# # promax rotated
#
# watch_promax_plots <- promax_df %>%
#   filter(block %in% c("Neu_Watch", "Neg_Watch", "Pos_Watch")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("watch_promax_500.png")), plot = watch_promax_plots, device = "png", width = 14)
#
# pos_promax_plots <- promax_df %>%
#   filter(block %in% c("Pos_Watch", "Pos_Inc", "Pos_Dec")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("pos_promax_500.png")), plot = pos_promax_plots, device = "png", width = 14)
#
# neg_promax_plots <- promax_df %>%
#   filter(block %in% c("Neg_Watch", "Neg_Inc", "Neg_Dec")) %>%
#   ggplot(., aes(ms, mv)) +
#   geom_line(aes(color = block)) +
#   facet_wrap(~ component, nrow = 2)
#
# ggsave(filename = here("images", "paper_2", "pca_images", paste0("neg_promax_500.png")), plot = neg_promax_plots, device = "png", width = 14)
#
