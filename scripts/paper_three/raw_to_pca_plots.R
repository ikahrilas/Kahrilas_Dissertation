# compose final erp plots with raw and component data
# load packages
library(tidyverse)
library(patchwork)
library(here)

# load plot data
topos <- readRDS("data/paper_three/topo_plots.rds")
erps <- readRDS("data/paper_three/erp_plots.rds")
raw_topos <- readRDS("data/paper_three/raw_topo_plots.rds")
raw_erps <- readRDS("data/paper_three/raw_erp_plots.rds")

# arrow
plot_arrow <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0.5,
           xend = 0.5,
           y = 1,
           yend = 0,
           size = 1,
           arrow = arrow(type = "closed",
                         length = unit (.1, "inches"))) +
  annotate(geom = "text",
           x = 0.45,
           y = 0.51,
           label = "PCA",
           angle = 90,
           size = 7)

# composition and saving
layout <- '
AAABBBCCC
DDDDDDDDD
EEEEEEEEE
FFFGGGHHH
IIIIIIIII
'

neg_rc8_component_plots <-
  raw_topos[[1]][[2]] +
  raw_topos[[1]][[1]] +
  raw_topos[[1]][[3]] +
  raw_erps[[1]] +
  plot_arrow +
  topos[[3]][[2]] +
  topos[[3]][[1]] +
  topos[[3]][[3]] +
  erps[[3]] +
  plot_layout(design = layout,
              heights = c(1.1, 1.4, 0.8, 1.1, 1.4),
              widths = 2,
              guides = "auto")
  # plot_annotation(title = "257 ms Component (EPN)",
  #                 theme = theme(plot.title = element_text(hjust = 0.5,
  #                                                         size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "neg_rc8_plots.png"),
       plot = neg_rc8_component_plots,
       height = 8,
       width = 8)

rc2_component_plots <-
  raw_topos[[2]][[2]] +
  raw_topos[[2]][[1]] +
  raw_topos[[2]][[3]] +
  raw_erps[[2]] +
  plot_arrow +
  topos[[5]][[2]] +
  topos[[5]][[1]] +
  topos[[5]][[3]] +
  erps[[5]] +
  plot_layout(design = layout,
              heights = c(1.1, 1.4, 0.8, 1.1, 1.4),
              widths = 2,
              guides = "auto")
  # plot_annotation(title = "371 ms Component (LPP)",
  #                 theme = theme(plot.title = element_text(hjust = 0.5,
  #                                                         size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc2_plots.png"),
       plot = rc2_component_plots,
       height = 8,
       width = 8)

rc3_component_plots <-
  raw_topos[[3]][[2]] +
  raw_topos[[3]][[1]] +
  raw_topos[[3]][[3]] +
  raw_erps[[3]] +
  plot_arrow +
  topos[[6]][[2]] +
  topos[[6]][[1]] +
  topos[[6]][[3]] +
  erps[[6]] +
  plot_layout(design = layout,
              heights = c(1.1, 1.4, 0.8, 1.1, 1.4),
              widths = 2,
              guides = "auto")
  # plot_annotation(title = "736 ms Component (LPP)",
  #                 theme = theme(plot.title = element_text(hjust = 0.5,
  #                                                         size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc3_plots.png"),
       plot = rc3_component_plots,
       height = 8,
       width = 8)
