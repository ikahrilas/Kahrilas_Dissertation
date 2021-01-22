# 2-step PCA
## this script performs a temporospatial (i.e., two step) PCA on average-corrected ERP data.
## the end of the script provides a dataset containing factor scores that can be submitted to analyses.

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

# define modified principal function from psych package that can perform an oblique infomax rotation for spatial PCA
principal_info <- function(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA, covar=FALSE,scores=TRUE,missing=FALSE,impute="median",oblique.scores=TRUE,method="regression",use="pairwise",cor="cor",correct=.5,weight=NULL,...) {
  cl <- match.call()
  n <- dim(r)[2]

  if (!isCorrelation(r)  && (!isCovariance(r))){ #added  (isCovariance)  April 9, 2019 to handle the case of a covariance matrix
    raw <- TRUE
    n.obs <- dim(r)[1]
    if(scores) {x.matrix <- as.matrix(r)  #matrices are required for the substitution to work

    if(missing) {
      miss <- which(is.na(x.matrix),arr.ind=TRUE)
      if(impute=="mean") {
        item.means <- colMeans(x.matrix,na.rm=TRUE)   #replace missing values with means
        x.matrix[miss]<- item.means[miss[,2]]} else {
          item.med   <- apply(x.matrix,2,median,na.rm=TRUE) #replace missing with medians
          x.matrix[miss]<- item.med[miss[,2]]}
    }}
    # 2011.12.21  added the covar option

    switch(cor,
           cor = {if(!is.null(weight))  {r <- cor.wt(r,w=weight)$r} else  {
             r <- cor(r,use=use)}
           },
           cov = {r <- cov(r,use=use)
           covar <- TRUE},
           wtd = { r <- cor.wt(r,w=weight)$r},
           spearman = {r <- cor(r,use=use,method="spearman")},
           kendall = {r <- cor(r,use=use,method="kendall")},
           tet = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
           poly = {r <- polychoric(r,correct=correct,weight=weight)$rho},
           tetrachoric = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
           polychoric = {r <- polychoric(r,correct=correct,weight=weight)$rho},
           mixed = {r <- mixedCor(r,use=use,correct=correct)$rho},
           Yuleb = {r <- YuleCor(r,bonett=TRUE)$rho},
           YuleQ = {r <- YuleCor(r,1)$rho},
           YuleY = {r <- YuleCor(r,.5)$rho }
    )
    #   if(!covar) {r <- cor(r,use="pairwise")} else r <- cov(r,use="pairwise")  # if given a rectangular matrix, then find the correlations or covariances first
  } else {
    raw <- FALSE
    if(!is.matrix(r)) {  r <- as.matrix(r)}
    sds <- sqrt(diag(r))    #convert covariance matrices to correlation matrices
    if(!covar) r <- r/(sds %o% sds)  }  #added June 9, 2008
  if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)}
  #added 24/4/15 to stop with bad data and give more meaningful help
  if(any(is.na(r))) {
    bad <- TRUE
    tempr <-r
    wcl <-NULL
    while(bad) {
      wc <- table(which(is.na(tempr), arr.ind=TRUE))  #find the correlations that are NA
      wcl <- c(wcl,as.numeric(names(which(wc==max(wc)))))
      tempr <- r[-wcl,-wcl]
      if(any(is.na(tempr))) {bad <- TRUE} else {bad <- FALSE}
    }

    cat('\nLikely variables with missing values are ',colnames(r)[wcl],' \n')
    stop("I am sorry: missing values (NAs) in the correlation matrix do not allow me to continue.\nPlease drop those variables and try again." )
  }

  eigens <- eigen(r)    #call the eigen value decomposition routine
  result$values <- eigens$values
  eigens$values[ eigens$values < .Machine$double.eps] <-  .Machine$double.eps  #added May 14, 2009 to fix case of singular matrices
  loadings <- eigens$vectors %*% sqrt(diag(eigens$values,nrow=length(eigens$values))) #added May 2, 2016 for the weird case of a single variable with covariance > 1

  if(nfactors > 0) {loadings <- loadings[,1:nfactors]} else {nfactors <- n}
  if (nfactors > 1) {communalities <- rowSums(loadings^2)} else {communalities <- loadings^2 }
  uniquenesses <- diag(r) - communalities # 2011.12.21 uniqueness is now found if covar is true
  names(communalities) <- colnames(r)    # 2009.02.10 Make sure this is a named vector -- correction by Gumundur Arnkelsson



  #added January 19, 2009 to flip based upon colSums of loadings
  if (nfactors > 1) {sign.tot <- vector(mode="numeric",length=nfactors)
  sign.tot <- sign(colSums(loadings))
  sign.tot[sign.tot==0] <- 1
  loadings <- loadings %*% diag(sign.tot)
  } else { if (sum(loadings) < 0) {loadings <- -as.matrix(loadings)} else {loadings <- as.matrix(loadings)}
    colnames(loadings) <- "PC1" }


  colnames(loadings) <- paste("PC",1:nfactors,sep='')
  rownames(loadings) <- rownames(r)
  Phi <- NULL

  rot.mat <- NULL
  if(rotate != "none") {if (nfactors > 1) {
    if (rotate=="varimax" |rotate=="Varimax" | rotate=="quartimax" | rotate =="bentlerT" | rotate =="geominT" | rotate =="targetT" | rotate =="bifactor"   | rotate =="TargetT"|
        rotate =="equamax"| rotate =="varimin"|rotate =="specialT" | rotate =="Promax"  | rotate =="promax"| rotate =="cluster" |rotate == "biquartimin"  |rotate =="specialQ" ) {
      Phi <- NULL
      colnames(loadings) <- paste("RC",1:nfactors,sep='')   #for rotated component
      switch(rotate,  #The orthogonal cases  for GPArotation + ones developed for psych
             varimax = {rotated <- stats::varimax(loadings,...)  #varimax is from stats, the others are from GPArotation
             loadings <- rotated$loadings
             rot.mat <- rotated$rotmat},
             Varimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
               #varimax is from the stats package, Varimax is from GPArotations
               #rotated <- do.call(rotate,list(loadings,...))
               #rotated <- do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...))
               rotated <- GPArotation::Varimax(loadings,...)
               loadings <- rotated$loadings
               rot.mat <- t(solve(rotated$Th))} ,
             quartimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}

               #rotated <- do.call(rotate,list(loadings))
               rotated <- GPArotation::quartimax(loadings,...)
               loadings <- rotated$loadings
               rot.mat <- t(solve(rotated$Th))} ,
             bentlerT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}

               #rotated <- do.call(rotate,list(loadings,...))
               rotated <- GPArotation::bentlerT(loadings,...)
               loadings <- rotated$loadings
               rot.mat <- t(solve(rotated$Th))} ,
             geominT	= {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}

               #rotated <- do.call(rotate,list(loadings,...))
               rotated <- GPArotation::geominT(loadings,...)
               loadings <- rotated$loadings
               rot.mat <- t(solve(rotated$Th))} ,
             targetT = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
               rotated <- GPArotation::targetT(loadings,Tmat=diag(ncol(loadings)),...)
               loadings <- rotated$loadings
               rot.mat <- t(solve(rotated$Th))} ,

             bifactor = {rot <- bifactor(loadings,...)
             loadings <- rot$loadings
             rot.mat <- t(solve(rot$Th))},
             TargetT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
               rot <- GPArotation::targetT(loadings,Tmat=diag(ncol(loadings)),...)
               loadings <- rot$loadings
               rot.mat <- t(solve(rot$Th))},
             equamax =  {rot <- equamax(loadings,...)
             loadings <- rot$loadings
             rot.mat <- t(solve(rot$Th))},
             varimin = {rot <- varimin(loadings,...)
             loadings <- rot$loadings
             rot.mat <- t(solve(rot$Th))},
             specialT =  {rot <- specialT(loadings,...)
             loadings <- rot$loadings
             rot.mat <- t(solve(rot$Th))},
             Promax =   {pro <- Promax(loadings,...)
             loadings <- pro$loadings
             Phi <- pro$Phi
             rot.mat <- pro$rotmat},
             promax =   {pro <- stats::promax(loadings,...)   #from stats
             loadings <- pro$loadings
             rot.mat <- pro$rotmat
             ui <- solve(rot.mat)
             Phi <-  cov2cor(ui %*% t(ui))},
             cluster = 	 {loadings <- varimax(loadings,...)$loadings
             pro <- target.rot(loadings)
             loadings <- pro$loadings
             Phi <- pro$Phi
             rot.mat <- pro$rotmat},
             biquartimin =    {ob <- biquartimin(loadings,...)
             loadings <- ob$loadings
             Phi <- ob$Phi
             rot.mat <- t(solve(ob$Th))},
             #  TargetQ  =  {ob <- TargetQ(loadings,...)
             #                 loadings <- ob$loadings
             # 				 Phi <- ob$Phi
             # 				 rot.mat <- t(solve(ob$Th))},
             specialQ = {ob <- specialQ(loadings,...)
             loadings <- ob$loadings
             Phi <- ob$Phi
             rot.mat <- t(solve(pro$Th))})
    } else {
      colnames(loadings) <- paste("TC",1:nfactors,sep='') #for transformed components
      #The following oblique cases all use GPArotation
      if (rotate == "infomaxQ"|rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ"  |rotate == "targetQ"  ) {
        if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
          Phi <- NULL} else {

            ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...)))
            if(inherits(ob, as.character("try-error")))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
              ob <- Promax(loadings)}

            loadings <- ob$loadings
            Phi <- ob$Phi
            rot.mat <- t(solve(ob$Th))}
      } else {message("Specified rotation not found, rotate='none' used")
        colnames(loadings) <- paste("PC",1:nfactors,sep='')  }  #not rotated
    }
  }
  }

  #just in case the rotation changes the order of the components, sort them by size of eigen value
  if(nfactors >1) {
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated,decreasing=TRUE)
    loadings <- loadings[,ev.order]}
  if(!is.null(Phi)) {Phi <- Phi[ev.order,ev.order] } #January 20, 2009 but, then, we also need to change the order of the rotation matrix!
  signed <- sign(colSums(loadings))
  c.names <- colnames(loadings)
  signed[signed==0] <- 1
  loadings <- loadings %*% diag(signed)  #flips factors to be in positive direction but loses the colnames
  colnames(loadings) <- c.names
  if(!is.null(Phi)) {Phi <- diag(signed) %*% Phi %*% diag(signed)
  colnames(Phi) <- rownames(Phi) <- c.names}



  class(loadings) <- "loadings"
  #Find the summary statistics of Variance accounted for
  #normally just found in the print function  (added 4/22/17)
  #from  the print function
  if(is.null(Phi)) {if(nfactors > 1)  {vx <- colSums(loadings^2) } else {vx <- sum(loadings^2)
  }} else {vx <- diag(Phi %*% t(loadings) %*% loadings)
  }

  vtotal <- sum(communalities + uniquenesses)
  names(vx) <- colnames(loadings)
  varex <- rbind("SS loadings" =   vx)
  varex <- rbind(varex, "Proportion Var" =  vx/vtotal)
  if (nfactors > 1) {
    varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/vtotal))
    varex <- rbind(varex, "Proportion Explained"=  vx/sum(vx))
    varex <- rbind(varex, "Cumulative Proportion"=  cumsum(vx/sum(vx)))
  }

  result$n.obs <- n.obs
  stats <- factor.stats(r,loadings,Phi,n.obs,fm="pc")
  class(result) <- c("psych", "principal")
  result$fn <- "principal"
  result$loadings <- loadings
  result$Phi <- Phi
  result$Call <- cl
  result$communality <- communalities
  result$uniquenesses <- uniquenesses
  result$complexity <- stats$complexity
  #result$stats <- stats
  result$chi <- stats$chi
  result$EPVAL <- stats$EPVAL
  result$R2 <- stats$R2
  result$objective <- stats$objective
  result$residual <- stats$residual
  result$rms <- stats$rms
  result$fit <-  stats$fit
  result$fit.off <-  stats$fit.off
  result$factors <-  stats$factors
  result$dof <-  stats$dof
  result$null.dof <- stats$null.dof
  result$null.model <- stats$null.model
  result$criteria <-  stats$criteria
  result$STATISTIC <-  stats$STATISTIC
  result$PVAL <-  stats$PVAL
  result$weights <-  stats$weights
  result$r.scores <-  stats$r.scores
  result$rot.mat <- rot.mat
  result$Vaccounted <-varex
  if(!is.null(Phi) && oblique.scores) {
    result$Structure <- loadings %*% Phi} else {result$Structure <- loadings
    }

  if(scores && raw) {


    result$weights <- try(solve(r,result$Structure),silent=TRUE)
    if(inherits(result$weights, "try-error"))  {warning("The matrix is not positive semi-definite, scores found from Structure loadings")
      result$weights <- result$Structure} # else {
    result$scores <- scale(x.matrix,scale=!covar) %*% result$weights #}
  }

  #   result$scores <- factor.scores(scale(x.matrix,scale=!covar),result$Structure,Phi=Phi,method=method,rho=r)  # method = method added Nov 20, 2011
  #  result$weights<- result$scores$weights
  #   result$scores <- result$scores$scores}
  return(result)
}


# read in average referenced data set
dat_2000 <- read_csv(here("data", "paper_two", "pre_pca_dat.csv"))

# this is the code used to impute missing values for raw data set. it is commented out to save computing time
# # make dataset with timepoints as variables and 2000 ms cutoff
# ## beyond 2000 ms is not of interest and was also not subjected to artifact correction in the preprocessing stage
# dat_2000_missing <- dat %>%
#   filter(ms < 2000) %>%
#   pivot_longer(cols = c(A1:EXG2),
#                names_to = "elec",
#                values_to = "mv")
#   pivot_wider(names_from = ms,       # create variables that represent a combination of each electrode
#               values_from = mv) # at each time point
#
# dat_2000 <- map_df(unique(dat_2000$block), ~ {# mean impute each missing value by each block, such that averages for
#     tmp <- filter(dat_2000, block == .x)      ## each electrode/time variable are computed using only values from the
#     map_df(tmp, ~ {                           ## same block condition
#         impute(.x, fun = mean)
#       })
#   }) %>%
#   pivot_longer(-c(pid:prop_trials), names_to = "elec_ms", values_to = "mv") %>%
#   separate(col = "elec_ms", into = c("elec", "ms"), sep = "_") %>%
#   pivot_wider(names_from = ms, values_from = mv)

write_csv(dat_2000, file = here("data", "paper_two", "pre_pca_dat.csv"))

# The parallel code below is commented out since it is computationally intensive. Uncomment and run
# to see results of parallel analysis, though results are below in the the comments.
# conduct parallel analyses on temporal data to determine number of components to retain for temporal PCA
dat_pa <- dat_2000 %>% select(-c(pid, block, elec, n_trials, prop_trials)) # filter out variables
## run the parallel analysis - this takes a long time.
# paran(dat_pa,
#       centile = 95,
#       iterations = 100, # the default is 30p (p = # of columns) which takes really long
#       status = TRUE,    # i've seen little variation in using 1 - 100, so settled on 100
#       graph = TRUE,
#       cfa = FALSE
#       )
## results suggest that 22 components or 43 factors should be retained

# perform temporal PCA with covariance matrix and promax rotation
# with 22 components, which is informed by the parallel analysis
## promax rotation with kappa = 3, tends to give best results for ERPs and is the default for SAS
## covariance matrix (mean corrected)
## derive factor scores using "Harman" method, which finds weights based upon so-called "idealized" variables
dat_pca_promax <- principal(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
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
  relocate(ms, paste0("RC", 1:22))

# long form for plotting
cov_loadings_long <- pivot_longer(cov_loadings_df, cols = RC1:RC22, names_to = "component", values_to = "mv")
cov_loadings_long$component <- factor(cov_loadings_long$component, levels = paste0("RC", 1:22))

# plot it to analyze time course
temp_loadings_plot <- ggplot(cov_loadings_long, aes(ms, mv)) +
  geom_line() +
  facet_wrap(~ component, nrow = 5)
temp_loadings_plot

# components to retain based on time course: 2-5, 7, 9 ,11, & 12
comp_to_retain <- paste0("RC", c(2:5, 7, 9, 11, 12))

cov_loadings_df <- cov_loadings_df %>% select(ms, all_of(comp_to_retain))

# The topography of each factor is encoded by the mean amplitude of its factor scores at each site.
# One can use this information to reproduce the portion of an observation's waveform represented by
# a given factor by multiplying the time point factor loadings by the observation's factor score and
# then multiplying each time point by its standard deviation (Dien, 1998a).

# extract factor scores from PCA
factor_scores_df <- data.frame(dat_pca_promax$scores)
names(factor_scores_df) <- comp_vector

# merge with original data that has block type and electrode variables
dat_2000_fac_scores <- bind_cols(dat_2000 %>%
                                   select(pid:elec),
                                 factor_scores_df)

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
topo_dat <- dat_2000_fac_scores %>%
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
  drop_na() %>%
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
# save the plot
ggsave(here("images", "paper_2", "component_topos", paste0(component, ".png")),
       plot = p,
       device = "png",
       width = 10,
       dpi = "retina")
}
# iterate the function over each component
map(comp_to_retain, ~ topo_facet(.x))

# plot ERPs for components of interest
## iterate over each component and multiply covariance loadings by factor scores for each
## observation.
## create "comp" variable as well.
temp_raw_df <- map_df(1:length(comp_vector), ~ {
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
  mutate(comp = comp_vector[.x]) %>%
  bind_cols(dat_2000 %>%
              select(pid:elec)) %>%
  relocate(pid:elec, comp)
})

# derive ERP plots in mV units for each temporal component
## define list of electrodes of interest for each component based on visual
## inspection of topo plots

elec_selections <- list(c("A29", "B26"),
                        c("A29", "B26"),
                        "B28",
                        c("A29", "B26"),
                        c("A29", "A31", "B30", "B26"),
                        c("A29", "B26"),
                        "B28", # site of negative activity for RC12
                        c("A29", "B26"))

# list of component sites for map function
component_list <- list("RC2",
                       "RC3",
                       "RC4",
                       "RC5",
                       "RC7",
                       "RC11",
                       "RC12",
                       "RC9")

# iterate over each component and electrode selection to create ERP plots
map2(component_list, elec_selections, ~ {
temp_raw_df %>%
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
  ggsave(here("images", "paper_2", "temporal_component_ERPs", paste0(.x, ".png")),
         device = "png",
         width = 12,
         dpi = "retina")
})

# prepare factor score data from temporal PCA for 2nd step spatial PCA with electrodes as columns
pre_spatial_pca_dat <- bind_cols(select(dat_2000, pid:elec),
                                 select(factor_scores_df, all_of(comp_to_retain))) %>%
  pivot_longer(cols = RC2:RC12,
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
  select(contains(.x))

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
                        c("B29", "B27", "B23", "B19", "B22", "B28"),
                        "B26",
                        "B32",
                        c("B30", "B31", "B32"),
                        c("A27, A29"),
                        "A31",
                        "B32",
                        "B20",
                        "A30")

# list of component sites for map function
component_list_temp_spat <- list("TC2-SC1",
                       "TC2-SC2",
                       "TC2-SC6",
                       "TC3-SC1",
                       "TC5-SC1",
                       "TC5-SC4",
                       "TC5-SC6",
                       "TC9-SC1",
                       "TC9-SC2",
                       "TC11-SC2")

# iterate over each component and electrode selection to create ERP plots
map2(component_list_temp_spat, elec_selections_temp_spat, ~ {
  temp_spat_erp_df %>%
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
  ggsave(here("images", "paper_2", "temporospatial_components_ERPs", paste0(.x, ".png")),
         device = "png",
         width = 12,
         dpi = "retina")
})








plot_butterfly(rename(long_loading_loc, "time" = "ms", "electrode" = "elec", "amplitude" = "mv"))

long_loading_loc %>%
  rename("time" = "ms", "electrode" = "elec", "amplitude" = "mv") %>%
  filter(block %in% c("Pos_Watch")) %>%
  erp_scalp(.,
            montage = "biosemi64alpha")
# dat_pca_promax_cor <- principal(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
#                                 nfactors = 43,
#                                 rotate = "Promax",
#                                 m = 3,
#                                 normalize = TRUE,
#                                 cor = "cor",
#                                 missing = TRUE)
# # perform promax rotation on kaiser-normalized factor loadings
# dat_pca_promax <- Promax(dat_pca_nr,
#                          normalize = TRUE,
#                          m = 3)

# code below is for spatial PCA, which might be helpful for P3-like and positive slow wave components.
test <- principal(weighted_dfs_lst[[2]] %>% select(-c(pid, block, ms, n_trials, prop_trials, component)),
                  nfactors = 10,
                  rotate = "infomaxQ",
                  cor = "cov",
                  missing = TRUE)

GPFoblq(test$loadings, method = "infomax", normalize = TRUE, maxit = 100000)



save.image(file = paste0("scripts/paper_two/", Sys.Date(), "_two_step_pca", ".RData"))

# old unused code:
# # this code extracts the covariance loadings for each rotated factor along with its corresponding timepoint
# # in a data frame, such that it can be merged with the raw data to derive factor scores.
# ## name of electrodes to be used in tidyeval
# elec_vec <- c(paste0("A", 1:32), paste0("B", 1:32), "EXG1", "EXG2")
# ## covariance loading by component extraction
# loadings_lst <- map(paste0("RC", 1:43), ~ {
#   select(cov_loadings_df, all_of(.x), ms) %>%
#     rename("component" = all_of(.x))
# })
# ## 43 data frames with weighted ERPs for each component
# weighted_dfs_lst <- map(loadings_lst, ~ {
#   full_join(filter(dat, ms < 2000), all_of(.x), by = "ms") %>%
#     mutate(across(.cols = all_of(elec_vec), .fns = ~ all_of(.x) * component))
# })


# ## unused alternative parallel analysis code
# # conduct parallel analysis
# parallel <- fa.parallel(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
#                         fm = "ml",
#                         fa = "both",
#                         n.iter = 50,
#                         quant = .95,
#                         SMC = TRUE)
# ## results of parallel analysis suggest 43 factors and 22 components
# # cleaner ggplot style scree plot
# ## create dataframe
# obs <- data.frame(parallel$fa.values)
# obs$type <- c('Observed Data')
# obs$num <- c(row.names(obs))
# obs$num <- as.numeric(obs$num)
# colnames(obs) <- c('eigenvalue', 'type', 'num')
#
# ## create quantiles for simulated eigenvalues
# percentile <- apply(parallel$values,2,function(x) quantile(x,.95))
# min <- as.numeric(nrow(obs))
# min <- (4*min) - (min-1)
# max <- as.numeric(nrow(obs))
# max <- 4*max
# percentile1 <- percentile[min:max]
#
# ## create data frame for simulated data
# sim <-  data.frame(percentile)
# sim$type <-  c('Simulated Data (95th %ile)')
# sim$num <-  c(row.names(obs))
# sim$num <-  as.numeric(sim$num)
# colnames(sim) <-  c('eigenvalue', 'type', 'num')
#
# ## merge the dataframes
# eigendat <- rbind(obs, sim)
#
# ## apa theme for scree plot
# apatheme <- theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_blank(),
#         text=element_text(family='Arial'),
#         legend.title=element_blank(),
#         legend.position=c(.7,.8),
#         axis.line.x = element_line(color='black'),
#         axis.line.y = element_line(color='black'))
#
# ## now plot it
# p <-  ggplot(filter(eigendat, num < 43), aes(x=num, y=eigenvalue, shape=type)) +
#   #Add lines connecting data points
#   geom_line()+
#   #Add the data points.
#   geom_point(size=4)+
#   #Label the y-axis 'Eigenvalue'
#   scale_y_continuous(name='Eigenvalue')+
#   #Label the x-axis 'Factor Number', and ensure that it ranges from 1-max # of factors, increasing by one with each 'tick' mark.
#   scale_x_continuous(name='Factor Number', breaks=min(eigendat$num):max(eigendat$num))+
#   #Manually specify the different shapes to use for actual and simulated data, in this case, white and black circles.
#   scale_shape_manual(values=c(16,1)) +
#   #Add vertical line indicating parallel analysis suggested max # of factors to retain
#   geom_vline(xintercept = parallel$nfact, linetype = 'dashed')+
#   #Apply our apa-formatting theme
#   apatheme
# p
#
#
#
#
#
#
#
#
#
#
