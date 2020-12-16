# 2-step PCA

# load packages
library(tidyverse)
library(here)
library(factoextra)
library(sjstats)
library(psych)
library(devtools)
library(eegUtils) # remotes::install_github("craddm/eegUtils")
library(patchwork)

# read in average referenced data set
## average referenced data
dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
## mastoid referenced data
dat_mast <- read_csv(here("data", "paper_two", "created_data", "erp_mast.csv"))

# make dataset with timepoints as variables and 2000 ms cutoff for early and later components
dat_2000 <- dat %>%
  filter(ms < 2000) %>%
  pivot_longer(c(A1:EXG2), names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)

# conduct parallel analysis
parallel <- fa.parallel(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
                        fm = "ml",
                        fa = "both",
                        n.iter = 50,
                        quant = .95,
                        SMC = TRUE)
## results of parallel analysis suggest 43 factors and 22 components
# cleaner ggplot style scree plot
## create dataframe
obs <- data.frame(parallel$fa.values)
obs$type <- c('Observed Data')
obs$num <- c(row.names(obs))
obs$num <- as.numeric(obs$num)
colnames(obs) <- c('eigenvalue', 'type', 'num')

## create quantiles for simulated eigenvalues
percentile <- apply(parallel$values,2,function(x) quantile(x,.95))
min <- as.numeric(nrow(obs))
min <- (4*min) - (min-1)
max <- as.numeric(nrow(obs))
max <- 4*max
percentile1 <- percentile[min:max]

## create data frame for simulated data
sim <-  data.frame(percentile)
sim$type <-  c('Simulated Data (95th %ile)')
sim$num <-  c(row.names(obs))
sim$num <-  as.numeric(sim$num)
colnames(sim) <-  c('eigenvalue', 'type', 'num')

## merge the dataframes
eigendat <- rbind(obs, sim)

## apa theme for scree plot
apatheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='Arial'),
        legend.title=element_blank(),
        legend.position=c(.7,.8),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'))

## now plot it
p <-  ggplot(filter(eigendat, num < 43), aes(x=num, y=eigenvalue, shape=type)) +
  #Add lines connecting data points
  geom_line()+
  #Add the data points.
  geom_point(size=4)+
  #Label the y-axis 'Eigenvalue'
  scale_y_continuous(name='Eigenvalue')+
  #Label the x-axis 'Factor Number', and ensure that it ranges from 1-max # of factors, increasing by one with each 'tick' mark.
  scale_x_continuous(name='Factor Number', breaks=min(eigendat$num):max(eigendat$num))+
  #Manually specify the different shapes to use for actual and simulated data, in this case, white and black circles.
  scale_shape_manual(values=c(16,1)) +
  #Add vertical line indicating parallel analysis suggested max # of factors to retain
  geom_vline(xintercept = parallel$nfact, linetype = 'dashed')+
  #Apply our apa-formatting theme
  apatheme
p

# perform temporal PCA with covariance matrix and no rotation
# with 43 factors, which is informed by the parallel analysis
dat_pca_nr <- principal(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
                        nfactors = 43,
                        rotate = "none",
                        cor = "cov",
                        missing = TRUE)
# perform promax rotation on kaiser-normalized factor loadings
dat_pca_promax <- kaiser(dat_pca_nr, rotate = "Promax")

# create data frame with covariance loadings
ms_vec <- dat %>%
  filter(ms < 2000,
         pid == 206201832,
         block == "Pos_Inc") %>%
  select(ms) %>%
  pull()
cov_loadings_mat <- matrix(dat_pca_promax$loadings, nrow = length(ms_vec))
cov_loadings_df <- data.frame(cov_loadings_mat)
names(cov_loadings_df) <- paste0("RC", c(1:43))
cov_loadings_df <- cov_loadings_df %>%
  mutate(ms = ms_vec)

# long form for plotting
cov_loadings_long <- pivot_longer(cov_loadings_df, cols = RC1:RC43, names_to = "component", values_to = "mv")
cov_loadings_long$component <- factor(cov_loadings_long$component, levels = c(paste0("RC", 1:43)))

# plot it
temp_loadings_plot <- ggplot(cov_loadings_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 8)
temp_loadings_plot
## based on these loading plots, should retain 1-3, 5, 6, 9, 10, 11, 12-32, 34, 36-38, 41

# this code extracts the covariance loadings for each rotated factor along with its corresponding timepoint
# in a data frame, such that it can be merged with the raw data to derive factor scores.
## name of electrodes to be used in tidyeval
elec_vec <- c(paste0("A", 1:32), paste0("B", 1:32), "EXG1", "EXG2")
## covariance loading by component extraction
loadings_lst <- map(paste0("RC", 1:43), ~ {
  select(cov_loadings_df, all_of(.x), ms) %>%
    rename("component" = .x)
})
## 43 data frames with weighted ERPs for each component
weighted_dfs_lst <- map(loadings_lst, ~ {
  full_join(filter(dat, ms < 2000), all_of(.x), by = "ms") %>%
    mutate(across(.cols = all_of(elec_vec), .fns = ~ all_of(.x) * component))
})

# read in EEG coordinate data
elec_loc <- read_csv(here("data", "paper_two", "Equidistant Layout.csv"))
elec_loc <- elec_loc %>%
  rename("channel" = `channel name`) %>%
  filter(channel != "CMS", channel != "DRL") %>%
  select(-number)

elec_loc$radian_phi <- pi/180 * elec_loc$phi

elec_loc <- elec_loc %>%
  mutate(x = theta * cos(radian_phi),
         y = theta * sin(radian_phi))

# start writing code for topo plotting
long_loading <- pivot_longer(weighted_dfs_lst[[3]],
                             cols = all_of(elec_vec),
                             names_to = "elec",
                             values_to = "mv")

long_loading_loc <- left_join(long_loading, elec_loc, by = c("elec" = "channel"))

# try plotting
long_loading_loc %>%
  filter(!is.na(mv), block %in% c("Neu_Watch", "Pos_Watch", "Neg_Watch")) %>%
  group_by(block, elec, x, y) %>%
  dplyr::summarize(mv = mean(mv, na.rm = TRUE)) %>%
  mutate(block = if_else(block == "Neu_Watch", "Neutral",
                         if_else(block == "Pos_Watch", "Positive", "Negative")),
         block = factor(block, levels = c("Neutral", "Positive", "Negative"))) %>%
  ggplot(., aes(x = x, y = y, fill = mv)) +
  stat_scalpmap() +
  geom_mask(scale_fac = 1.7) +
  geom_head() +
  geom_channels(size = 0.125) +
  scale_fill_viridis_c(limits = c(-1, 1), oob = scales::squish) +
  scale_color_manual(breaks = c("black", "white"),
                     values = c("black", "white"),
                     guide = FALSE) +
  labs(fill = expression(paste("Weighted ", mu,"V"))) +
  coord_equal() +
  theme_void() +
  facet_wrap(~ block, ncol = 3) +
  theme(strip.text.x = element_text(size = 12, vjust = 1),
        strip.text.y = element_text(size = 12, angle = 90))


plot_butterfly(rename(long_loading_loc, "time" = "ms", "electrode" = "elec", "amplitude" = "mv"))


erp_scalp(rename(long_loading_loc, "time" = "ms", "electrode" = "elec", "amplitude" = "mv"),
          montage = "biosemi64alpha",
          size = 0.5,
          color = as.factor(.data$block))


save.image(file = paste0("scripts/paper_two/", Sys.Date(), "_two_step_pca", ".RData"))
