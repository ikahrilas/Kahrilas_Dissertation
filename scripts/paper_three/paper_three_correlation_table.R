# Correlation table

# load packages
library(tidyverse)
library(here)
library(magick)
library(kableExtra)
library(GGally)

## peaks from hs_per_pca script:
# [1] "The maximum timepoint for RC2 is 370.800031605892"
# [2] "The maximum timepoint for RC3 is 736.346627189118"
# [3] "The maximum timepoint for RC5 is 134.269881522629"
# [4] "The maximum timepoint for RC7 is 194.868515015035"
# [5] "The maximum timepoint for RC12 is 102.993167462032"
# [6] "The maximum timepoint for RC8 is 257.421943136229"
# [7] "The maximum timepoint for RC17 is 1176.17541866626"

# load data
cor_dat <- read_csv("data/paper_three/dat_for_analyses_2021-03-19.csv") %>%
  filter(!is.na(RC2)) %>%
  rename("371 ms Comp" = "RC2", # rename components
         "736 ms Comp" = "RC3",
         "134 ms Comp" = "RC5",
         "195 ms Comp" = "RC7",
         "257 ms Comp" = "RC8",
         "1176 ms Comp" = "RC17")

glimpse(cor_dat)

# create correlation matrix
cor_list_fun <- function(block_type) {
  tab <- cor_dat %>%
    filter(block == block_type) %>%
    select(`371 ms Comp`:`1176 ms Comp`,
           phq_total,
           masq_pa:masq_ad,
           pswq_total,
           panas_pos_total:panas_neg_total,
           pss_total) %>%
    relocate(`134 ms Comp`,
             `195 ms Comp`,
             `257 ms Comp`,
             `371 ms Comp`,
             `736 ms Comp`,
             `1176 ms Comp`,
             everything())

  #Compute correlation matrix
  x <- as.matrix(tab)
  correlation_matrix <- Hmisc::rcorr(x)
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value

  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .05, "*", " ")

  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

  ## build a new matrix that includes the correlations with their appropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")

  ## remove upper triangle of correlation matrix
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)

  ## remove last column and return the correlation matrix
  tab <- cbind(Rnew[1:length(Rnew)-1])

  ## clean it up and prepare for latex rendering
  tab <- tab %>%
    as.data.frame() %>%
    map_df(., ~ {
      as.character(.x) %>%
        str_trim() %>%
        str_remove("^0+")
    }) %>%
    map_df(., ~ {
      if_else(str_detect(.x, "-"), paste0("-", str_remove(.x, "-0")), .x)
    })
  colnames(tab) <- paste0(c(1:15), ".")

  tab %>%
    mutate(" " = c("1. 134 ms Comp",
                   "2. 195 ms Comp",
                   "3. 257 ms Comp",
                   "4. 371 ms Comp",
                   "5. 736 ms Comp",
                   "6. 1176 ms Comp",
                   "7. Depression",
                   "8. MASQ PA",
                   "9. MASQ NA",
                   "10. MASQ AA",
                   "11. MASQ AD",
                   "12. Worry",
                   "13. PANAS PA",
                   "14. PANAS NA",
                   "15. Stress")) %>%
    select(" ", everything())
}

corr_list <- map(unique(cor_dat$block), ~ cor_list_fun(.x))

#save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_correlation_table-data", ".RData"))

## create tables
blocks <- c("Negative Watch",
            "Neutral Watch",
            "Positive Watch")

map2(corr_list, blocks, ~ {
  .x %>%
    kable(linesep = "", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 14)), caption = .y) %>%
    kable_styling(latex_options = "scale_down", font_size = 12) %>%
    column_spec(1:15, width = "3em") %>%
    footnote(general_title = "Note.",
             escape = FALSE,
             general = "*$p < .05$. PA = positive affectivity, NA = negative affectivity,
          AA = Anxious Arousal, AD = Anhedonic Depression.",
             threeparttable = TRUE,
             footnote_as_chunk = TRUE)
})

# create scatterplot matrix for each block

pos <-
  cor_dat %>%
  filter(block == "Pos_Watch") %>%
  select(`371 ms Comp`:`1176 ms Comp`,
         phq_total,
         masq_pa:masq_ad,
         pswq_total,
         panas_pos_total:panas_neg_total,
         pss_total) %>%
  relocate(`134 ms Comp`,
           `195 ms Comp`,
           `257 ms Comp`,
           `371 ms Comp`,
           `736 ms Comp`,
           `1176 ms Comp`,
           everything()) %>%
    ggpairs(title = "Positive Watch")

neg <-
  cor_dat %>%
  filter(block == "Neg_Watch") %>%
  select(`371 ms Comp`:`1176 ms Comp`,
         phq_total,
         masq_pa:masq_ad,
         pswq_total,
         panas_pos_total:panas_neg_total,
         pss_total) %>%
  relocate(`134 ms Comp`,
           `195 ms Comp`,
           `257 ms Comp`,
           `371 ms Comp`,
           `736 ms Comp`,
           `1176 ms Comp`,
           everything()) %>%
  ggpairs(title = "Negative Watch")

neu <-
  cor_dat %>%
  filter(block == "Neu_Watch") %>%
  select(`371 ms Comp`:`1176 ms Comp`,
         phq_total,
         masq_pa:masq_ad,
         pswq_total,
         panas_pos_total:panas_neg_total,
         pss_total) %>%
  relocate(`134 ms Comp`,
           `195 ms Comp`,
           `257 ms Comp`,
           `371 ms Comp`,
           `736 ms Comp`,
           `1176 ms Comp`,
           everything()) %>%
  ggpairs(title = "Neutral Watch")
