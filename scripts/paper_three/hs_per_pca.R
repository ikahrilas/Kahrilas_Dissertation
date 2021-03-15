# PCA for merged headspace/per data

# load packages
library(paran)
library(tidyverse)
library(here)
library(factoextra)
library(sjstats)
library(psych)
library(devtools)
library(eegUtils) # remotes::install_github("craddm/eegUtils")
library(patchwork)
library(GPArotation)
library(Hmisc)

# read in average referenced data set
hs_dat <- read_csv(here("data", "paper_three", "erp_dat.csv")) %>%
  filter(ms < 2000,
         block %in% c("Pos_Watch", "Neu_Watch", "Neg_Watch")) %>%
  pivot_longer(cols = -c(pid:ms, group),
               names_to = "elec",
               values_to = "mv") %>%
  pivot_wider(names_from = "ms",
              values_from = "mv") %>%
  filter(pid != 22585577)

per_dat <- read_csv(here("data", "paper_two", "pre_pca_dat.csv")) %>%
  filter(block %in% c("Pos_Watch", "Neu_Watch", "Neg_Watch")) %>% #only passive watch conditions
  select(-prop_trials) %>%
  mutate(group = "per")

# merge them together
dat <- bind_rows(hs_dat, per_dat)

# The parallel code below is commented out since it is computationally intensive. Uncomment and run
# to see results of parallel analysis, though results are below in the the comments.
# conduct parallel analyses on temporal data to determine number of components to retain for temporal PCA
# dat_pa <- dat_2000 %>% select(-c(pid, block, elec, n_trials, prop_trials)) # filter out variables
## run the parallel analysis - this takes a long time.
paran(dat %>% select(-c(pid, block, n_trials, group, elec)),
      centile = 95,
      iterations = 100, # the default is 30p (p = # of columns) which takes really long
      status = TRUE,    # i've seen little variation in using 1 - 100, so settled on 100
      graph = TRUE,
      cfa = FALSE
      )
## results suggest that 22 components should be retained

# perform temporal PCA with covariance matrix and promax rotation
# with 14 components, which is informed by the parallel analysis
## promax rotation with kappa = 3, tends to give best results for ERPs and is the default for SAS
## covariance matrix (mean corrected)
## derive factor scores using "Harman" method, which finds weights based upon so-called "idealized" variables
dat_pca_promax <- principal(select(dat, -c(pid, block, elec, n_trials, group)),
                            nfactors = 22,
                            rotate = "promax", # SPSS seems to do a Kaiser normalization before doing
                            m = 3,             ## Promax, this is done here by the call to "promax"
                            cor = "cov",       ## which does the normalization before calling Promax in GPArotation.
                            missing = TRUE,
                            scores = TRUE,
                            method = "Harman")

# create data frame with covariance loadings
cov_loadings_mat <- as.matrix(unclass(dat_pca_promax$loadings))
cov_loadings_df <- as_tibble(cov_loadings_mat) %>%
  mutate(ms = as.numeric(dimnames(cov_loadings_mat)[[1]])) %>%
  relocate(ms, everything()) # components in order of % of explained variance

# long form for plotting
cov_loadings_long <- pivot_longer(cov_loadings_df, cols = -ms, names_to = "component", values_to = "mv")
cov_loadings_long$component <- factor(cov_loadings_long$component)

# plot it to analyze time course
temp_loadings_plot <- ggplot(cov_loadings_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 5)
temp_loadings_plot

# components to explore based on time course: 2-5, 7, 9 ,11, & 12
comp_to_retain <- paste0("RC", c(2, 3, 5, 7, 12, 8, 17))

# replot with just components of interest and save image
cov_loadings_long %>%
  filter(component %in% comp_to_retain) %>%
  ggplot(., aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 4)
ggsave(here("images", "paper_3", "cov_loadings_comp_of_interest.png"), device = "png", width = 10)

cov_loadings_df <- cov_loadings_df %>% select(ms, all_of(comp_to_retain))

# find peaks for components of interest
map_chr(comp_to_retain, ~{
  max_ms <- cov_loadings_df %>%
    filter(cov_loadings_df[[.x]] == max(cov_loadings_df[[.x]])) %>%
    select(ms) %>%
    pull()
  return(paste("The maximum timepoint for", .x, "is", max_ms))
})

# The topography of each factor is encoded by the mean amplitude of its factor scores at each site.
# One can use this information to reproduce the portion of an observation's waveform represented by
# a given factor by multiplying the time point factor loadings by the observation's factor score and
# then multiplying each time point by its standard deviation (Dien, 1998a).

# extract factor scores from PCA
factor_scores_df <- data.frame(dat_pca_promax$scores)

# merge with original data that has block type and electrode variables
dat_fac_scores <- bind_cols(dat %>% select(pid:elec),
                                    factor_scores_df %>%
                                   select(comp_to_retain)) # save the data set only with components you want to
# run the analyses on
# save data set to work space for statistical analyses
write_csv(dat_fac_scores,
          here("data", "paper_three", "hs_per_fac_score_dat.csv"))

# read in EEG coordinate data
elec_loc <- read_csv(here("data", "paper_two", "Equidistant Layout.csv"))
elec_loc <- elec_loc %>%
  rename("channel" = `channel name`) %>%
  filter(channel != "CMS", channel != "DRL") %>%
  select(-number)

elec_loc$radian_phi <- pi/180 * elec_loc$phi

## data frame with electrode coordinates
elec_loc <- elec_loc %>%
  mutate(x = theta * cos(radian_phi),
         y = theta * sin(radian_phi))

# create data frame with valence and regulation variables and merge with electrode
# coordinate data and merge with factor score data
topo_dat <- dat_fac_scores %>%
  group_by(block, elec) %>%
  summarise(across(all_of(comp_to_retain), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(
    valence = case_when(
      str_detect(block, "Pos") ~ "Positive",
      str_detect(block, "Neg") ~ "Negative",
      str_detect(block, "Neu") ~ "Neutral"
    ),
    regulation = case_when (
      str_detect(block, "Watch") ~ "Watch",
      str_detect(block, "Inc") ~ "Increase",
      str_detect(block, "Dec") ~ "Decrease"
    )) %>%
  left_join(elec_loc, by = c("elec" = "channel"))

# create faceted topoplots
topo_facet <- function(component) {
  p <- ggplot(topo_dat,aes(x = x, y = y, fill = topo_dat[[component]], label = elec)) +
    geom_topo(grid_res = 300,
              interp_limit = "head",
              chan_markers = "text",
              chan_size = 2) +
    scale_fill_distiller(palette = "RdBu") +
    theme_void() +
    coord_equal() +
    labs(fill = expression(paste("Factor Score"))) +
    facet_grid(regulation ~ valence, switch = "both") +
    theme(strip.text.x = element_text(size = 12, vjust = 1),
          strip.text.y = element_text(size = 12, angle = 90)) +
    ggtitle(component) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16))
}
# iterate the function over each component
map(comp_to_retain, ~ {
  topo_facet(.x)
# save the plot
ggsave(here("images", "paper_3", "component_topos", paste0(.x, ".png")),
       plot = last_plot(),
       device = "png",
       width = 10,
       dpi = "retina")
})


# plot ERPs for components of interest
## iterate over each component and multiply covariance loadings by factor scores for each
## observation.
## create "comp" variable as well.
temp_raw_df <- map_df(1:length(names(dat_pca_promax$R2)), ~ {
  scores_matrix <- as.matrix(dat_pca_promax$scores)
  loadings_matrix <- t(as.matrix(unclass(dat_pca_promax$loadings)))

  tmp <- matrix(loadings_matrix[.x,],
                nrow = nrow(scores_matrix),
                ncol = ncol(loadings_matrix),
                byrow = TRUE,
                dimnames = list(NULL, dimnames(loadings_matrix)[[2]]))

  tmp_scores_mat <- matrix(rep(scores_matrix[,.x],
                               times = ncol(loadings_matrix)),
                           ncol = ncol(loadings_matrix))

  raw_mat <- tmp * tmp_scores_mat

  raw_dat <- as_tibble(raw_mat) %>%
    mutate(comp = names(dat_pca_promax$R2)[.x]) %>%
    bind_cols(dat %>%
                select(pid:elec)) %>%
    relocate(pid:elec, comp)
})

# write data set with "raw" factor ERPs to work space
comp_for_erp <- c("RC2", "RC3", "RC5", "RC7", "RC8", "RC17")

temp_raw_df %>%
  filter(comp %in% comp_for_erp) %>%
  pivot_longer(-c(pid:comp),
               names_to = "ms",
               values_to = "mv") %>%
  mutate(ms = as.numeric(ms)) %>%
  write_csv(here("data", "paper_three", "hs_per_temp_fac_score_erp.csv"))

# derive ERP plots in mV units for each temporal component
## define list of electrodes of interest for each component based on visual
## inspection of topo plots

###################################################
# temporal components that change across conditions
## RC2 @ c(A29, B26)
## RC3 @ c(A29, B26, A26, B23, B28)
## RC5 @ c(A29, B26)
## RC11 @ c(A29, B26)
## RC12 @ c(B21, B28) for negativity
## RC12 @ c(A29, A31, B30, B26) for positivity

elec_selections <- list(c("A29", "B26"),
                        c("A29", "B26", "A26", "B23", "B28"),
                        c("A29", "B26"),
                        c("A29", "B26"),
                        c("B21", "B28")) # site of negative activity for RC12

# list of component sites for map function
component_list <- list("RC2",
                       "RC3",
                       "RC5",
                       "RC11",
                       "RC12")

# iterate over each component and electrode selection to create ERP plots
map2(component_list, elec_selections, ~ {
  watch <- temp_raw_df %>%
    filter(elec %in% c(.y),
           comp == .x,
           block %in% c("Pos_Watch", "Neu_Watch", "Neg_Watch")) %>%
    pivot_longer(cols = -c(pid:comp),
                 names_to = "ms",
                 values_to = "mv") %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude (",mu,"V)")),
         title = paste("Average", .x, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))

  positive <- temp_raw_df %>%
    filter(elec %in% c(.y),
           comp == .x,
           block %in% c("Pos_Watch", "Pos_Inc", "Pos_Dec")) %>%
    pivot_longer(cols = -c(pid:comp),
                 names_to = "ms",
                 values_to = "mv") %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude (",mu,"V)")),
         title = paste("Average", .x, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))

  negative <- temp_raw_df %>%
    filter(elec %in% c(.y),
           comp == .x,
           block %in% c("Neg_Watch", "Neg_Inc", "Neg_Dec")) %>%
    pivot_longer(cols = -c(pid:comp),
                 names_to = "ms",
                 values_to = "mv") %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude (",mu,"V)")),
         title = paste("Average", .x, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
  # save each plot
  ggsave(plot = watch, here("images", "paper_2", "temporal_component_ERPs", paste0(.x, "_watch.png")),
         device = "png",
         width = 12,
         dpi = "retina")
  ggsave(plot = positive, here("images", "paper_2", "temporal_component_ERPs", paste0(.x, "_positive.png")),
         device = "png",
         width = 12,
         dpi = "retina")
  ggsave(plot = negative, here("images", "paper_2", "temporal_component_ERPs", paste0(.x, "_negative.png")),
         device = "png",
         width = 12,
         dpi = "retina")
})


# prepare factor score data from temporal PCA for 2nd step spatial PCA with electrodes as columns
## only bother with components of interest from temporal PCA
comp_to_retain <- c("RC2", "RC3", "RC5", "RC7", "RC8", "RC17")

pre_spatial_pca_dat <- bind_cols(select(dat_2000, pid:elec),
                                 select(factor_scores_df, all_of(comp_to_retain))) %>%
  pivot_longer(cols = all_of(comp_to_retain),
               names_to = "comp",
               values_to = "mv") %>%
  pivot_wider(names_from = elec,
              values_from = mv)

# conduct parallel analysis for each temporal factor
## define function that returns number of factors to retain
parallel_fun <- function(component) {
  test <- paran(pre_spatial_pca_dat %>%
                  filter(comp == component) %>%
                  select(-c(pid:comp)),
                centile = 95,
                status = FALSE,
                graph = FALSE)
  return(test$Retained)
}

component_vector <- map_dbl(comp_to_retain, ~ parallel_fun(.x))
names(component_vector) <- comp_to_retain
# the average number of components to retain for each component is 6.25, so
# each spatial PCA will retain 6 components

# define function that performs spatial PCA for each of the 8 components
# from the temporal PCA
spatial_pca_fun <- function(comp_num){
  principal_info(pre_spatial_pca_dat %>%
                   filter(comp == names(component_vector)[comp_num]) %>%
                   select(-c(pid:comp)),
                 nfactors = 6,
                 rotate = "infomaxQ",
                 maxit = 100000,
                 cor = "cov",
                 method = "Harman")
}

# run the function and return the results into a list of length 8
spatial_pca_lst <- map(1:length(component_vector), ~ spatial_pca_fun(.x))

# extract factor scores into a dataframe and write csv file

## temporospatial components to retain
temp_spat_comp_retain <- c("RC2-TC1", "RC2-TC6", "RC3-TC1", "RC3-TC2", "RC5-TC6", "RC5-TC5",
                           "RC5-TC2", "RC5-TC1", "RC5-TC4", "RC11-TC4", "RC11-TC1", "RC11-TC3",
                           "RC12-TC3", "RC12-TC2")

data.frame(spatial_pca_lst[[1]]$scores) %>%
  rename_with(.cols = everything(), .fn = ~ paste(names(component_vector)[1], .x, sep = "-")) %>%
  bind_cols(., pre_spatial_pca_dat %>%
              filter(comp == names(component_vector)[1]) %>%
              select(pid:prop_trials)) %>%
  relocate(pid:prop_trials) %>%
  bind_cols(data.frame(spatial_pca_lst[[2]]$scores) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste(names(component_vector)[2], .x, sep = "-"))) %>%
  bind_cols(data.frame(spatial_pca_lst[[3]]$scores) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste(names(component_vector)[3], .x, sep = "-"))) %>%
  bind_cols(data.frame(spatial_pca_lst[[4]]$scores) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste(names(component_vector)[4], .x, sep = "-"))) %>%
  bind_cols(data.frame(spatial_pca_lst[[5]]$scores) %>%
              rename_with(.cols = everything(),
                          .fn = ~ paste(names(component_vector)[5], .x, sep = "-"))) %>%
  select(pid:prop_trials, temp_spat_comp_retain) %>%
  write_csv(file = here("data", "paper_two", "temp_spat_fac_score_dat.csv"))

# extract factor scores from PCA
factor_scores_df <- data.frame(dat_pca_promax$scores)
names(factor_scores_df) <- comp_vector

# merge with original data that has block type and electrode variables
dat_2000_fac_scores <- bind_cols(dat_2000 %>%
                                   select(pid:elec),
                                 factor_scores_df)

# save data set to work space for statistical analyses
write_csv(dat_2000_fac_scores,
          here("data", "paper_two", "temp_fac_score_dat.csv"))


# run function to derive factor scores for each electrode for each participant in each block (for topo plots)
temp_component_names <- str_replace(names(component_vector), "R", "T")

temp_spat_pca_df <- map_df(1:length(temp_component_names), ~ {

  scores_matrix <- as.matrix(spatial_pca_lst[[.x]]$scores)
  loadings_matrix <- t(as.matrix(unclass(spatial_pca_lst[[.x]]$loadings)))
  spat_comp_names <- str_replace(dimnames(scores_matrix)[[2]], "T", "S")

  tmp <- map_df(1:length(spat_comp_names), ~ {
    tmp_loadings_mat <- matrix(loadings_matrix[.x,],
                               nrow = nrow(scores_matrix),
                               ncol = ncol(loadings_matrix),
                               byrow = TRUE,
                               dimnames = list(NULL, dimnames(loadings_matrix)[[2]]))

    tmp_scores_mat <- matrix(rep(scores_matrix[,.x],
                                 times = ncol(loadings_matrix)),
                             ncol = ncol(loadings_matrix))

    raw_mat <- tmp_loadings_mat * tmp_scores_mat

    raw_dat <- as_tibble(raw_mat) %>%
      mutate(comp = spat_comp_names[.x]) %>%
      bind_cols(pre_spatial_pca_dat %>%
                  filter(comp == "RC2") %>%
                  select(pid:prop_trials)) %>%
      relocate(pid:prop_trials, comp)
  })

  tmp %>%
    mutate(comp = paste(temp_component_names[.x], comp, sep = "-")) %>%
    pivot_longer(cols = c(A1:EXG2),
                 names_to = "elec",
                 values_to = "fac_score")
})

# prepare data frame with temporospatial factors for topo plotting
temp_spat_pca_topo_df <- temp_spat_pca_df %>%
  pivot_wider(names_from = "comp",
              values_from = "fac_score") %>%
  group_by(block, elec) %>%
  summarise(across(`TC2-SC1`:`TC12-SC4`, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(
    valence = case_when(
      str_detect(block, "Pos") ~ "Positive",
      str_detect(block, "Neg") ~ "Negative",
      str_detect(block, "Neu") ~ "Neutral"
    ),
    regulation = case_when (
      str_detect(block, "Watch") ~ "Watch",
      str_detect(block, "Inc") ~ "Increase",
      str_detect(block, "Dec") ~ "Decrease"
    )) %>%
  left_join(elec_loc, by = c("elec" = "channel"))

# define function for temporospatial topo plots - nearly identical to the earlier function except that
# the images are saved to a different location
topo_facet_spat <- function(component) {
  p <- ggplot(temp_spat_pca_topo_df,aes(x = x, y = y, fill = temp_spat_pca_topo_df[[component]], label = elec)) +
    geom_topo(grid_res = 300,
              interp_limit = "head",
              chan_markers = "text",
              chan_size = 2) +
    scale_fill_distiller(palette = "RdBu") +
    theme_void() +
    coord_equal() +
    labs(fill = expression(paste("Factor Score"))) +
    facet_grid(regulation ~ valence, switch = "both") +
    theme(strip.text.x = element_text(size = 12, vjust = 1),
          strip.text.y = element_text(size = 12, angle = 90)) +
    ggtitle(component) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16))
  # save the plot
  ggsave(here("images", "paper_2", "spat_component_topos", paste0(component, ".png")),
         plot = p,
         device = "png",
         width = 10,
         dpi = "retina")
}
# define string of temporospatial components
temp_spat_comp_names <- str_subset(names(temp_spat_pca_topo_df), pattern = "TC")

# iterate the function over each component
map(temp_spat_comp_names, ~ topo_facet_spat(.x))

##############################
# TS factors that change across conditions:
## 2-1 @ c(A31, B32)
## 2-6 @ B26
## 3-1 @ B32
## 3-2 @ c(B7, B8, B9, B14, B15)
## 5-6 @ c(A30, A31)
## 5-5 @ B21
## 5-2 @ B26
## 5-1 @ c(B30, B31, B32)
## 5-4 @ A27
## 11-4 @ A27
## 11-1 @ c(A31, B30, B32)
## 11-3 @ B26
## 12-3 @ c(A27, A29)
## 12-2 @ B29

# define vectors of temporal and spatial components for map function
temp_comp_to_retain <- str_replace_all(comp_to_retain, "R", "T")
spat_comp_vec <- paste0("SC", 1:6)

# multiply temporospatial factor scores and factor loadings from temporal PCA
# to create data frame for ERP waveforms

temp_spat_erp_df <- map_df(temp_comp_to_retain, ~ {

  loadings_retained_comp_mat <- as_tibble(unclass(dat_pca_promax$loadings)) %>%
    select(.data[[str_replace(.x, "T", "R")]]) %>%
    as.matrix()

  temp_spat_scores_df <- temp_spat_pca_df %>%
    pivot_wider(names_from = comp,
                values_from = fac_score) %>%
    select(contains(paste0(.x, "-")))

  loadings_retained_comp_mat <- matrix(t(loadings_retained_comp_mat),
                                       nrow = 24024,
                                       ncol = 1126,
                                       byrow = TRUE)
  final_tbl <- map_df(spat_comp_vec, ~ {
    tmp <- temp_spat_scores_df %>%
      select(contains(.x)) %>%
      as.matrix()

    scores_mat <- matrix(rep(tmp, times = 1126),
                         ncol = 1126)

    prod_tbl <- as_tibble(loadings_retained_comp_mat * scores_mat)

    names(prod_tbl) <- dimnames(dat_pca_promax$loadings)[[1]]

    prod_tbl %>%
      mutate(comp = .x) %>%
      bind_cols(dat_2000 %>%
                  select(pid:elec)) %>%
      relocate(pid:elec, comp)
  })

  final_tbl %>%
    mutate(comp = paste(.x, comp, sep = "-"))
})

# derive ERP plots in mV units for each temporospatial component
## define list of electrodes of interest for each component based on visual
## inspection of topo plots

elec_selections_temp_spat <- list(c("A31", "B32"),
                                  "B26",
                                  "B32",
                                  c("B7", "B8", "B9", "B14", "B15"),
                                  c("A30", "A31"),
                                  "B21",
                                  "B26",
                                  c("B30", "B31", "B32"),
                                  "A27",
                                  "A27",
                                  c("A31", "B30", "B32"),
                                  "B26",
                                  c("A27", "A29"),
                                  "B29")

# list of component sites for map function
component_list_temp_spat <- list("TC2-SC1",
                                 "TC2-SC6",
                                 "TC3-SC1",
                                 "TC3-SC2",
                                 "TC5-SC6",
                                 "TC5-SC5",
                                 "TC5-SC2",
                                 "TC5-SC1",
                                 "TC5-SC4",
                                 "TC11-SC4",
                                 "TC11-SC1",
                                 "TC11-SC3",
                                 "TC12-SC3",
                                 "TC12-SC2")

# iterate over each component and electrode selection to create ERP plots
map2(component_list_temp_spat, elec_selections_temp_spat, ~ {
  p <- temp_spat_erp_df %>%
    filter(elec %in% c(.y),
           comp == .x) %>%
    pivot_longer(cols = -c(pid:comp),
                 names_to = "ms",
                 values_to = "mv") %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude (",mu,"V)")),
         title = paste("Average", .x, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
  # save each plot
  ggsave(plot = p, filename = here("images", "paper_2", "temp_spat_ERPs", paste0(.x, ".png")),
         device = "png",
         width = 12)
})
