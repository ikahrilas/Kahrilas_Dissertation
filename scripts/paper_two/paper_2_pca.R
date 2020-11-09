# Principal components analysis

# load packages
library(tidyverse)
library(here)
library(factoextra)
library(sjstats)
library(psych)

# read in average referenced data set
# using average reference as per Dein (2012)
dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))

# rearrange data for temporal PCA with timepoints as columns
dat_temp <- dat %>%
  pivot_longer(A1:EXG2, names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)

# pca
dat_pca <- prcomp(na.omit(dat_temp[,-c(1:5)]), center = TRUE, scale. = FALSE)
promax(dat_pca)
fviz_eig(dat_pca_promax$values)

# pca with promax rotation
dat_pca_promax <- principal(dat_temp[,-c(1:5)], nfactors = 10, rotate = "promax", covar = TRUE, missing = TRUE)
dat_pca_promax
biplot.psych(dat_pca_promax)
scree_dat <- scree(dat_temp[,-c(1:5)])
