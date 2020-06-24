# Correlation table

# load packages
library(tidyverse)
library(here)
library(magick)

# load data
cor_dat <- read_csv("~/tmp/Kahrilas_Dissertation/data/paper_two/created_data/per_data_analyses_2020_5_19.csv")

# create correlation matrix
tab <- cor_dat %>%
  select(N170,
         EPN,
         LPP,
         LPP_front,
         arousal,
         valence,
         difficulty,
         pos_affectivity,
         anticipating,
         savoring_moment,
         reminiscing,
         depression,
         sex) %>%
  rename("1. N170" = N170,
         "2. EPN" = EPN,
         "3. LPP" = LPP,
         "4. Frontal LPP" = LPP_front,
         "5. Arousal Ratings" = arousal,
         "6. Valence Ratings" = valence,
         "7. Difficulty Ratings" = difficulty,
         "8. PA" = pos_affectivity,
         "9. Anticipating" = anticipating,
         "10. StM" = savoring_moment,
         "11. Reminiscing" = reminiscing,
         "12. Depression" = depression,
         "13. Sex" = sex)

#Compute correlation matrix
x <- as.matrix(tab)
correlation_matrix<-Hmisc::rcorr(x)
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
colnames(tab) <- paste0(c(1:12), ".")
save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_correlation_table-data", ".RData"))
# tab %>%
#   mutate(" " = c("1. N170",
#                  "2. EPN",
#                  "3. LPP",
#                  "4. Frontal LPP",
#                  "5. Arousal Ratings",
#                  "6. Valence Ratings",
#                  "7. Difficulty Ratings",
#                  "8. PA",
#                  "9. Anticipating",
#                  "10. StM",
#                  "11. Reminiscing",
#                  "12. Depression",
#                  "13. Sex")) %>%
#   select(" ", everything()) %>%
#   ## create table
#   kable(., linesep = "", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 13)), caption = "(ref:corr-table)") %>%
#   kable_styling(latex_options = "scale_down", font_size = 12) %>%
#   column_spec(2:13, width = "3em") %>%
#   footnote(general_title = "Note.",
#            escape = FALSE,
#            general = "*$p < .05$. PA = positive affectivity, StM = savoring the moment. EEG components and behavioral ratings (i.e., arousal, valence, and difficulty) averaged across all regulation blocks.",
#            threeparttable = TRUE,
#            footnote_as_chunk = TRUE)
