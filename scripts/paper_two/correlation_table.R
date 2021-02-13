# Correlation table

# load packages
library(tidyverse)
library(here)
library(magick)
library(kableExtra)

# load data
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv")

# rename components
cor_dat <- dat %>%
  rename("Early LPP" = "RC2",
         "Late LPP" = "RC3",
         "P125" = "RC5",
         "N170" = "RC11",
         "EPN" = "RC12",
         "EPP" = "pos_RC12")

# create correlation matrix
tab <- cor_dat %>%
  select(`Early LPP`:valence,
         pos_affectivity,
         savoring_moment,
         depression) %>%
  relocate(P125, N170, EPN, EPP, `Early LPP`, `Late LPP`, everything()) %>%
  rename("1. P125" = P125,
         "2. N170" = N170,
         "3. EPN" = EPN,
         "4. EPP" = EPP,
         "5. Early LPP" = `Early LPP`,
         "6. Late LPP" = `Late LPP`,
         "7. Arousal Ratings" = arousal,
         "8. Valence Ratings" = valence,
         "9. Difficulty Ratings" = difficulty,
         "10. PA" = pos_affectivity,
         "11. StM" = savoring_moment,
         "12. Depression" = depression)

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
colnames(tab) <- paste0(c(1:11), ".")
save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_correlation_table-data", ".RData"))
tab %>%
  mutate(" " = c("1. P125",
                 "2. N170",
                 "3. EPN",
                 "4. EPP",
                 "5. Early LPP",
                 "6. Late LPP",
                 "7. Arousal Ratings",
                 "8. Valence Ratings",
                 "9. Difficulty Ratings",
                 "10. PA",
                 "11. StM",
                 "12. Depression")) %>%
  select(" ", everything()) %>%
  ## create table
  kable(., linesep = "", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 13)), caption = "(ref:corr-table)") %>%
  kable_styling(latex_options = "scale_down", font_size = 12) %>%
  column_spec(1:12, width = "3em") %>%
  footnote(general_title = "Note.",
           escape = FALSE,
           general = "*$p < .05$. PA = positive affectivity, StM = savoring the moment. EEG components and behavioral ratings (i.e., arousal, valence, and difficulty) averaged across all regulation blocks.",
           threeparttable = TRUE,
           footnote_as_chunk = TRUE)

# table with EEG components and behavioral ratings conditioned on regulation block
cor_block_tab <- cor_dat %>%
  select(block,
         `Early LPP`:valence,
         pos_affectivity,
         savoring_moment,
         depression) %>%
  relocate(P125, N170, EPN, EPP, `Early LPP`, `Late LPP`, everything()) %>%
  pivot_wider(names_from = block,
              values_from = c(P125:`Late LPP`, arousal:valence),
              names_sep = "_") %>%
  relocate(P125_Neg_Inc:valence_Pos_Inc, savoring_moment, pos_affectivity, depression)

colnames(cor_block_tab) <- tools::toTitleCase(str_replace_all(colnames(cor_block_tab), "_", " "))

colnames(cor_block_tab) <- paste0(c(1:length(colnames(cor_block_tab))), ". ", colnames(cor_block_tab))

colnames(cor_block_tab) <- str_replace(colnames(cor_block_tab), "64. Savoring Moment", "64. StM")

colnames(cor_block_tab) <- str_replace(colnames(cor_block_tab), "65. Pos Affectivity", "65. PA")

col_names_for_table <- colnames(cor_block_tab)

#Compute correlation matrix
x <- as.matrix(cor_block_tab)
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
cor_block_tab <- cbind(Rnew[1:length(Rnew)-1])

## clean it up and prepare for latex rendering
cor_block_tab <- cor_block_tab %>%
  as.data.frame() %>%
  map_df(., ~ {
    as.character(.x) %>%
      str_trim() %>%
      str_remove("^0+")
  }) %>%
  map_df(., ~ {
    if_else(str_detect(.x, "-"), paste0("-", str_remove(.x, "-0")), .x)
  })

colnames(cor_block_tab) <- paste0(c(1:65), ".")

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "block_correlation_table-data", ".RData"))

cor_block_tab %>%
  mutate(" " = col_names_for_table) %>%
  select(" ", everything()) %>%
  ## create table
  kable(., format = "html",  linesep = "", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 13)), caption = "(ref:corr-block-table)") %>%
  kable_styling(latex_options = "scale_down", font_size = 12) %>%
  column_spec(2:66, width = "3em") %>%
  footnote(general_title = "Note.",
           escape = FALSE,
           general = "*$p < .05$. PA = positive affectivity, StM = savoring the moment",
           threeparttable = TRUE,
           footnote_as_chunk = TRUE)

# create correlation plot just highlighting the strongest relations among variables
