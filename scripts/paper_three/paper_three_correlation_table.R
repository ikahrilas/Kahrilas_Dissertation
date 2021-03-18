# Correlation table

# load packages
library(tidyverse)
library(here)
library(magick)
library(kableExtra)

## peaks from hs_per_pca script:
# [1] "The maximum timepoint for RC2 is 370.800031605892"
# [2] "The maximum timepoint for RC3 is 736.346627189118"
# [3] "The maximum timepoint for RC5 is 134.269881522629"
# [4] "The maximum timepoint for RC7 is 194.868515015035"
# [5] "The maximum timepoint for RC12 is 102.993167462032"
# [6] "The maximum timepoint for RC8 is 257.421943136229"
# [7] "The maximum timepoint for RC17 is 1176.17541866626"

# load data
cor_dat <- read_csv("data/paper_three/dat_for_analyses_2021-03-18.csv")%>%
  rename("381 ms Comp" = "RC2", # rename components
         "740 ms Comp" = "RC3",
         "124 ms Comp" = "RC5",
         "162 ms Comp" = "RC11",
         "Neg 259 ms Comp" = "RC12",
         "Pos 259 ms Comp" = "pos_RC12")

glimpse(cor_dat)

# create correlation matrix
cor_list_fun <- function(block_type) {
  tab <- cor_dat %>%
    filter(block == block_type) %>%
    select(`381 ms Comp`:valence,
           anticipating:reminiscing,
           pos_affectivity,
           neg_affectivity,
           erq_reappraisal,
           erq_suppression) %>%
    relocate(`124 ms Comp`,
             `162 ms Comp`,
             `Neg 259 ms Comp`,
             `Pos 259 ms Comp`,
             `381 ms Comp`,
             `740 ms Comp`,
             everything()) %>%
    rename("1. 124 ms Comp" = "124 ms Comp",
           "2. 162 ms Comp" = "162 ms Comp",
           "3. Neg 259 ms Comp" =  "Neg 259 ms Comp",
           "4. Pos 259 ms Comp" = "Pos 259 ms Comp",
           "5. 381 ms Comp" = "381 ms Comp",
           "6. 740 ms Comp" = "740 ms Comp",
           "7. Arousal Ratings" = arousal,
           "8. Valence Ratings" = valence,
           "9. Difficulty Ratings" = difficulty,
           "10. PA" = pos_affectivity,
           "11. NA" = neg_affectivity,
           "12. Ant" = anticipating,
           "13. StM" = savoring_moment,
           "14. Rem" = reminiscing,
           "15. Reapp" = erq_reappraisal,
           "16. Supp" = erq_suppression)

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
    mutate(" " = c("1. 124 ms Comp",
                   "2. 162 ms Comp",
                   "3. Neg 259 ms Comp",
                   "4. Pos 259 ms Comp",
                   "5. 381 ms Comp",
                   "6. 740 ms Comp",
                   "7. Arousal Ratings",
                   "8. Valence Ratings",
                   "9. Difficulty Ratings",
                   "10. PA",
                   "11. NA",
                   "12. Ant",
                   "13. StM",
                   "14. Rem",
                   "15. Reapp",
                   "16. Supp")) %>%
    select(" ", everything())
}

corr_list <- map(unique(cor_dat$block), ~ cor_list_fun(.x))

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_correlation_table-data", ".RData"))

## create tables
blocks <- c("Negative Increase", "Negative Decrease", "Negative Watch",
            "Neutral Watch",
            "Positive Watch", "Positive Decrease", "Positive Increase")

map2(corr_list, blocks, ~ {
  .x %>%
    kable(linesep = "", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 16)), caption = .y) %>%
    kable_styling(latex_options = "scale_down", font_size = 12) %>%
    column_spec(1:16, width = "3em") %>%
    footnote(general_title = "Note.",
             escape = FALSE,
             general = "*$p < .05$. PA = positive affectivity, NA = negative affectivity,
           Ant = anticipating, StM = savoring the moment, Rem = reminiscing,
           Reapp = reappraisal, Supp = suppression.",
             threeparttable = TRUE,
             footnote_as_chunk = TRUE)
})
