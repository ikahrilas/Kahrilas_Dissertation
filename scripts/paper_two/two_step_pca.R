# 2-step PCA

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
## average referenced data
dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
## mastoid referenced data
dat_mast <- read_csv(here("data", "paper_two", "created_data", "erp_mast.csv"))

# make dataset with timepoints as variables and 2000 ms cutoff for early and later components
dat_2000 <- dat %>%
  filter(ms < 2000) %>%
  pivot_longer(c(A1:EXG2), names_to = "elec") %>%
  pivot_wider(names_from = ms, values_from = value)

dat_pa <- dat_2000 %>% select(-c(pid, block, elec, n_trials, prop_trials))
dat_pa <- dat_pa[complete.cases(dat_pa), ]

paran(dat_pa,
      iterations = 100,
      centile = 95)

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

# perform temporal PCA with covariance matrix and promax rotation
# with 43 factors, which is informed by the parallel analysis
## promax rotation with kappa = 3, tends to give best results for ERPs and is the default for SAS
## covariance matrix (mean corrected)
## derive factor scores using "Harman" method, which finds weights based upon so-called "idealized" variables
dat_pca_promax <- principal(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
                            nfactors = 43,
                            rotate = "promax", # SPSS seems to do a Kaiser normalization before doing
                            m = 3,             ## Promax, this is done here by the call to "promax"
                            cor = "cov",       ## which does the normalization before calling Promax in GPArotation.
                            missing = TRUE,
                            scores = TRUE,
                            method = "Harman")

# dat_pca_nr <- principal(select(dat_2000, -c(pid, block, elec, n_trials, prop_trials)),
#                         nfactors = 43,
#                         rotate = "none",
#                         cor = "cov",
#                         missing = TRUE,
#                         scores = TRUE,
#                         method = "Harman")
#
# dat_pca_promax <- kaiser(dat_pca_nr, rotate = "Promax")
# extract factor scores from PCA
factor_scores_df <- data.frame(dat_pca_promax$scores)

# merge with raw data that has timepoints as variables - this can be used for topo plots
dat_2000_fac_scores <- bind_cols(dat_2000, factor_scores_df)

# turn data into long form for ERP raw waveform plotting
dat_2000_fac_scores_long <- dat_2000_fac_scores %>%
  pivot_longer(cols = c("-200":"1999.14395738572"),
               names_to = "ms",
               values_to = "mv") %>%
  mutate(ms = as.numeric(ms))

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

# The topography of each factor is encoded by the mean amplitude of its factor scores at each site.
# One can use this information to reproduce the portion of an observation's waveform represented by
# a given factor by multiplying the time point factor loadings by the observation's factor score and
# then multiplying each time point by its standard deviation (Dien, 1998a).

# merge covariance loading and factor score data

comp_raw_dat <- full_join(cov_loadings_df %>%
                           rename_with(.cols = contains("RC"), .fn = ~ paste0(.x, "_cov_loading")),
                         dat_2000_fac_scores_long %>%
                           rename_with(.cols = contains("RC"), .fn = ~ paste0(.x, "_fac_score")),
                         by = "ms") %>%
                    # this code is so horrendous, forgive me god, but trying to put this in a map function
                    # keeps throwing an error due to running out of memory. could rewrite this code in the
                    # future by multiplying matrices of covariance loadings and factor scores together rather
                    # than employing algebra notation here.
                           mutate(RC1_raw = RC1_cov_loading * RC1_fac_score,
                                  RC2_raw = RC2_cov_loading * RC2_fac_score,
                                  RC3_raw = RC3_cov_loading * RC3_fac_score,
                                  RC4_raw = RC4_cov_loading * RC4_fac_score,
                                  RC5_raw = RC5_cov_loading * RC5_fac_score,
                                  RC6_raw = RC6_cov_loading * RC6_fac_score,
                                  RC7_raw = RC7_cov_loading * RC7_fac_score,
                                  RC8_raw = RC8_cov_loading * RC8_fac_score,
                                  RC9_raw = RC9_cov_loading * RC9_fac_score,
                                  RC10_raw = RC10_cov_loading * RC10_fac_score,
                                  RC11_raw = RC11_cov_loading * RC11_fac_score,
                                  RC12_raw = RC12_cov_loading * RC12_fac_score,
                                  RC13_raw = RC13_cov_loading * RC13_fac_score,
                                  RC14_raw = RC14_cov_loading * RC14_fac_score,
                                  RC15_raw = RC15_cov_loading * RC15_fac_score,
                                  RC16_raw = RC16_cov_loading * RC16_fac_score,
                                  RC17_raw = RC17_cov_loading * RC17_fac_score,
                                  RC18_raw = RC18_cov_loading * RC18_fac_score,
                                  RC19_raw = RC19_cov_loading * RC19_fac_score,
                                  RC20_raw = RC20_cov_loading * RC20_fac_score,
                                  RC21_raw = RC21_cov_loading * RC21_fac_score,
                                  RC22_raw = RC22_cov_loading * RC22_fac_score,
                                  RC23_raw = RC23_cov_loading * RC23_fac_score,
                                  RC24_raw = RC24_cov_loading * RC24_fac_score,
                                  RC25_raw = RC25_cov_loading * RC25_fac_score,
                                  RC26_raw = RC26_cov_loading * RC26_fac_score,
                                  RC27_raw = RC27_cov_loading * RC27_fac_score,
                                  RC28_raw = RC28_cov_loading * RC28_fac_score,
                                  RC29_raw = RC29_cov_loading * RC29_fac_score,
                                  RC30_raw = RC30_cov_loading * RC30_fac_score,
                                  RC31_raw = RC31_cov_loading * RC31_fac_score,
                                  RC32_raw = RC32_cov_loading * RC32_fac_score,
                                  RC33_raw = RC33_cov_loading * RC33_fac_score,
                                  RC34_raw = RC34_cov_loading * RC34_fac_score,
                                  RC35_raw = RC35_cov_loading * RC35_fac_score,
                                  RC36_raw = RC36_cov_loading * RC36_fac_score,
                                  RC37_raw = RC37_cov_loading * RC37_fac_score,
                                  RC38_raw = RC38_cov_loading * RC38_fac_score,
                                  RC39_raw = RC39_cov_loading * RC39_fac_score,
                                  RC40_raw = RC40_cov_loading * RC40_fac_score,
                                  RC41_raw = RC41_cov_loading * RC41_fac_score,
                                  RC42_raw = RC42_cov_loading * RC42_fac_score,
                                  RC43_raw = RC43_cov_loading * RC43_fac_score)
# clean up the dataset a bit and only retain raw component variables plus other vars
comp_raw_dat <- comp_raw_dat %>%
                  select(paste0("RC", 1:43, "_raw"),
                         mv,
                         ms,
                         pid,
                         block,
                         n_trials,
                         prop_trials,
                         elec,
                         )

# create data frame with valence and regulation variables and merge with electrode
# coordinate data and merge with covariance loadings

topo_dat <- dat_2000_fac_scores %>%
  group_by(block, elec) %>%
  summarize(across(paste0("RC", 1:43), ~ mean(.x, na.rm = TRUE))) %>%
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
p <- ggplot(topo_dat, aes(x = x, y = y, fill = topo_dat[[component]])) +
  stat_scalpmap() +
  geom_mask(scale_fac = 1.7) +
  geom_head() +
  geom_channels(size = 0.125) +
  scale_fill_viridis_c(limits = c(min(topo_dat[[component]]), max(topo_dat[[component]])), oob = scales::squish) +
  scale_color_manual(breaks = c("black", "white"),
                     values = c("black", "white"),
                     guide = FALSE) +
  labs(fill = "Average Mv") +
  coord_equal() +
  theme_void() +
  facet_grid(regulation ~ valence, switch = "both") +
  theme(strip.text.x = element_text(size = 12, vjust = 1),
        strip.text.y = element_text(size = 12, angle = 90)) +
  ggtitle(component) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 16))
p
# save the plot
ggsave(here("images", "paper_2", "component_topos", paste0(component, ".png")),
       plot = p,
       device = "png",
       width = 14)
}
# iterate the function over each component
map(paste0("RC", 1:43), ~ topo_facet(.x))




plot_butterfly(rename(long_loading_loc, "time" = "ms", "electrode" = "elec", "amplitude" = "mv"))

long_loading_loc %>%
  rename("time" = "ms", "electrode" = "elec", "amplitude" = "mv") %>%
  filter(block %in% c("Pos_Watch")) %>%
erp_scalp(.,
          montage = "biosemi64alpha")



# prepare factor score data from temporal PCA for 2nd step spatial PCA
pre_spatial_pca_dat <- bind_cols(select(dat_2000, pid:elec), factor_scores_df) %>%
  relocate(pid:elec, paste0("RC", 1:43)) %>%
  pivot_longer(cols = c(RC1:RC43),
               names_to = "comp",
               values_to = "mv") %>%
  pivot_wider(names_from = elec,
              values_from = mv) %>%
  drop_na()

# parallel analysis

parallel_fun <- function(component) {
test <- paran(pre_spatial_pca_dat %>%
              filter(comp == component) %>%
              select(-c(pid:comp)),
              centile = 95,
              status = FALSE,
              graph = FALSE)
return(test$Retained)
}

component_vector <- map_dbl(c(paste0("RC", 1:43)), ~ parallel_fun(.x))
names(component_vector) <- paste0("RC", 1:43)

peepee <- principal_info(pre_spatial_pca_dat %>%
            filter(comp == names(component_vector)[1]) %>%
            select(-c(pid:comp)),
          nfactors = component_vector[1],
          rotate = "none",
          cor = "cov",
          method = "Harman")

peepee

test <- GPFoblq(peepee$loadings, method = "infomax")
test$Gq
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
                  rotate = "none",
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












