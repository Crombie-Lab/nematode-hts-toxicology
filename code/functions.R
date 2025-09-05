#!/usr/bin/env Rscript
library(tidyverse)

today <- format(Sys.Date(), "%Y%m%d")
#==================================================#
# Prompt Toxcast configuration
#==================================================#
tcplConfPrompt <- function(db   = "invitrodb_v2", user = "root", host = "localhost", drvr = "MySQL"){
  # prompt for password
  pw <- readline(prompt="what's the magic word?: ")
  # message
  message(glue::glue("great, configuring tcpl:\ntcpl::tcplConf(db = {db},\nuser = {user},\nhost = {host},\ndrvr = {drvr},\npass = {pw})"))
  # pass to config function
  tcpl::tcplConf(db = db, user = user, host = host, drvr = drvr, pass = pw)
}

#==================================================#
# Geometric mean function from P. McMurdie 
#==================================================#
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#==================================================#
# 3D distance 
#==================================================#
# write a #D distance function
dist.fun <- function(x1, y1, z1, x2, y2, z2) {
  dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2) # slope, intercept, rsq
  return(dist)
}

# Wrap the distance function in with data scaling and ploting
distFromIdeal <- function(data, plot = T) {
  # add in ideal data and scale all
  dat <- data %>%
    tibble::add_row(x = "ideal", y = "ideal", 
                    orth.reg.n.observations = min(.$orth.reg.n.observations),
                    orth.reg.slope = 1,
                    orth.reg.intercept = 0,
                    orth.reg.ssq = 0,
                    orth.reg.mse = 0,
                    orth.reg.r.squared = 1) %>%
    dplyr::mutate(scaled.slope = as.double(scale(orth.reg.slope)),
                  scaled.int = as.double(scale(orth.reg.intercept)),
                  scaled.rsq = as.double(scale(orth.reg.r.squared)))
  
  # get the scaled ideal values
  ideal <- dat %>%
    filter(x == "ideal" & y == "ideal")
  
  # now calculate distance from scaled ideal values
  dat.dist.from.ideal <- dat %>%
    dplyr::mutate(type = ifelse((x == "ideal" & y == "ideal"), "ideal", "data")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dist = dist.fun(x1 = scaled.slope,
                                  y1 = scaled.int,
                                  z1 = scaled.rsq,
                                  x2 = ideal$scaled.slope,
                                  y2 = ideal$scaled.int,
                                  z2 = ideal$scaled.rsq)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dist) %>%
    dplyr::mutate(dist_rank = 1:n())
  
  if(plot == T){
    
    p <- plotly::plot_ly(data = dat.dist.from.ideal, x = ~scaled.slope, y = ~scaled.int, z = ~scaled.rsq,
                            type = "scatter3d",
                            mode = "markers",
                            color = ~dist,
                            colors = viridis::magma(5, alpha = 1, begin = 0, end = 1, direction = -1),
                            text = ~glue::glue("dist={dist}\nn={orth.reg.n.observations}\n{x}\n{y}"),
                            showlegend = F)
    
    # return
    return(list(data = dat.dist.from.ideal, plot = p))
  }
  else {
  #return
  return(dat.dist.from.ideal)
  }
}

#======================================================================#
# Pull out mean squared error and r-squared from orthReg model
#======================================================================#
# r squared - https://stackoverflow.com/questions/75067630/orthogonal-linear-regression-total-least-squares-fit-get-rmse-and-r-squared-i
r_squared_odreg <- function(object, y) {
  denom <- sum((y - mean(y))^2)
  1 - object$ssq/denom
}
# meas squeared error - https://stackoverflow.com/questions/75067630/orthogonal-linear-regression-total-least-squares-fit-get-rmse-and-r-squared-i
mse_odreg <- function(object){
  mean(object$resid^2)
}

#======================================================================#
# Run Pearson's product-moment correlation
#======================================================================#
# Testing - run zebrafish.R etox data
#data <- or4$plots[[1]]$data 

normReg <- function(data, x, y, min.cases = 3, QC = "ignore") {
  # make a list to hold it all
  out <- NULL
  
  # QC filter?
  if(QC == "filter") {
    data <- data %>%
      dplyr::filter(QC == "PASS" | is.na(QC))
  }
  if(QC == "ignore") {
    data <- data
  }
  
  # filter data to focal taxa and and geometric mean for group
  all_focal_dat <- data %>%
    dplyr::mutate(duration_d = ifelse(is.na(duration_d), "NA", duration_d)) %>% # THIS STEP HANDELS NAs in duration data
    dplyr::mutate(endpoint = ifelse(is.na(endpoint), "NA", endpoint)) %>% # THIS STEP HANDELS NAs in endpoint data???????
    dplyr::mutate(pair = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ x[1],
                                          (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ y[1],
                                          TRUE ~ NA_character_),
                  pair_gen = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ "x",
                                              (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ "y",
                                              TRUE ~ NA_character_)) %>% # label pairs, should handle giving a group or a latin_name since they are unique
    dplyr::filter(!is.na(pair_gen)) %>% # filter to pairs with labels NEW pair_gen OLD pair
    dplyr::group_by(cas, pair_gen) %>% # NEW pair_gen OLD pair
    dplyr::mutate(gm_mean = gm_mean(effect_value), # take the mean TAC 20250303
                  min = min(effect_value),
                  max = max(effect_value)) %>% # get geometric mean for chemical and pair
    dplyr::ungroup()
  
  
  
  # reshape for plotting
  plot_dat <- all_focal_dat %>%
    dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
    tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
    dplyr::filter(complete.cases(.)) # keep only complete cases
  
  # handle case number filter
  if((nrow(plot_dat) < min.cases)) {
    # build output
    regdf <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                            y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                            reg.n.observations = nrow(plot_dat),
                            reg.slope = NA_real_,
                            reg.intercept = NA_real_,
                            reg.r.squared = NA_real_,
                            reg.r.pearson = NA_real_,
                            reg.r.pvalue = NA_real_,
                            reg.df = NA_real_)
    # return Reg df
    return(reg_df)
  } 
  else {
    # perform regression
    lm_fit <- lm(log10(plot_dat$gm_mean_y) ~ log10(plot_dat$gm_mean_x))
    lm_s <- summary(lm_fit)
    # get pearson
    p_cor <- cor.test(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
    
    # build output
    reg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                            y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                            reg.n.observations = nrow(plot_dat),
                            reg.slope = round(lm_s$coefficients[2], digits = 3),
                            reg.intercept = round(lm_s$coefficients[1], digits = 3),
                            reg.r.squared = round(lm_s$r.squared, digits = 3),
                            reg.r.pearson = round(unname(p_cor$estimate), digits = 3),
                            reg.r.pvalue = round(unname(p_cor$p.value), digits = 3),
                            reg.df = round(unname(p_cor$parameter), digits = 3))
    # return Reg df
    return(reg_df)
  }
}
#======================================================================#
# function for data selection, shaping, orthogonal regression, plotting 
#======================================================================#
# x is a vector of specific latin_name/group, test_statistic, duration_d, endpoint
# y is as x, code NAs as "NA"
# min.cases is the minimum number of observations required to perform the regression?
# plot is true to output a plot. If false, the default, no plot is made.
# QC is either "filter" or NULL. The default = NULL, no action.
# "filter" will remove the failed QC from the dataframe. 
orthReg <- function(data, x, y, min.cases = 3, QC = "ignore", plot = F){

  # make a list to hold it all
  out <- NULL
  
  # QC filter?
  if(QC == "filter") {
    data <- data %>%
      dplyr::filter(QC == "PASS" | is.na(QC))
  }
  if(QC == "ignore") {
    data <- data
  }
  
  # filter data to focal taxa and and geometric mean for group
  all_focal_dat <- data %>%
    dplyr::mutate(duration_d = ifelse(is.na(duration_d), "NA", duration_d)) %>% # THIS STEP HANDELS NAs in duration data
    dplyr::mutate(endpoint = ifelse(is.na(endpoint), "NA", endpoint)) %>% # THIS STEP HANDELS NAs in endpoint data???????
    dplyr::mutate(pair = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ x[1],
                                          (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ y[1],
                                          TRUE ~ NA_character_),
                  pair_gen = dplyr::case_when((latin_name == x[1] | group == x[1]) & test_statistic == x[2] & duration_d == x[3] & endpoint == x[4] ~ "x",
                                              (latin_name == y[1] | group == y[1]) & test_statistic == y[2] & duration_d == y[3] & endpoint == y[4] ~ "y",
                                              TRUE ~ NA_character_)) %>% # label pairs, should handle giving a group or a latin_name since they are unique
    dplyr::filter(!is.na(pair_gen)) %>% # filter to pairs with labels NEW pair_gen OLD pair
    dplyr::group_by(cas, pair_gen) %>% # NEW pair_gen OLD pair
    dplyr::mutate(gm_mean = gm_mean(effect_value),
                  min = min(effect_value),
                  max = max(effect_value)) %>% # get geometric mean for chemical and pair
    dplyr::ungroup()
  
  
  
  # reshape for plotting
  plot_dat <- all_focal_dat %>%
    dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
    tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
    dplyr::filter(complete.cases(.)) # keep only complete cases
  
  # handle case number filter
  if((nrow(plot_dat) < min.cases)) {
    # build output
    orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                 y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                 orth.reg.n.observations = nrow(plot_dat),
                                 orth.reg.slope = NA_real_,
                                 orth.reg.intercept = NA_real_,
                                 orth.reg.ssq = NA_real_,
                                 orth.reg.mse = NA_real_,
                                 orth.reg.r.squared = NA_real_)
    if(plot == T){
      # make dummy plot
      plot <- ggplot2::ggplot(plot_dat) +
        #ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
        #ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
        #ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
        #ggplot2::geom_point(shape = 21, fill = "red", color = "black") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(μg/L)"),
                      y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(μg/L)"),
                      subtitle = glue::glue("LESS THAN {min.cases} CASES\nx={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
        ggplot2::scale_x_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::annotation_logticks()
      
      # setup output list
      out$orthreg_df <- orthreg_df
      out$plot <- plot
      
      # return list with data and plot
      return(out)
    }
    else {
      # return orthReg df
      return(orthreg_df)
    }
  }
  else {
  # ORIGINAL REGRESSION BASED ON gm_mean
  orthog_reg_model_log10 <- pracma::odregress(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
  
  # # old orthogonal regression
  oldOR_df <- tibble::tibble(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
  pcObject <- princomp(oldOR_df)
  myLoadings <- unclass(loadings(pcObject))[,1]
  OR.slope <- myLoadings["y"]/myLoadings["x"]
  OR.int <- mean(log10(plot_dat$gm_mean_y))-OR.slope*mean(log10(plot_dat$gm_mean_x))
  ###Rsq is empirical in this setting.  The first principal component will always explain at least half
  ###of the total variation, since the first two components will explain 100% of it in this simple setting,
  ###and the first must explain at least as much as the second.
  orthog_reg_model_r_squared <- as.double(((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5)
  # corrCoef <- cor(log10(plot_dat$gm_mean_x),log10(plot_dat$gm_mean_y))
  # corr.rsq <- corrCoef^2
  
  # run the extraction functions
  orthog_reg_model_mse <- mse_odreg(orthog_reg_model_log10)
  # NEW Rsq!!!!! DIFFERNET THAN OLD - using old for now
  #orthog_reg_model_r_squared <- r_squared_odreg(orthog_reg_model_log10, log10(plot_dat$gm_mean_y))
  
  #----------------------------------------------------------------------------#
  # Add deming function
  #----------------------------------------------------------------------------#
  library(deming)
  deming_model_log10 <- deming::deming(log10(plot_dat$gm_mean_y) ~ log10(plot_dat$gm_mean_x))
  
  #----------------------------------------------------------------------------#
  # End deming function test
  #----------------------------------------------------------------------------#
  
  # build output
  orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                   y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                   orth.reg.n.observations = nrow(plot_dat),
                                   orth.reg.slope = orthog_reg_model_log10$coeff[1],
                                   orth.reg.intercept = orthog_reg_model_log10$coeff[2],
                                   orth.reg.ssq = orthog_reg_model_log10$ssq[1],
                                   orth.reg.mse = orthog_reg_model_mse,
                                   orth.reg.r.squared = orthog_reg_model_r_squared,
                               deming.reg.slope = deming_model_log10$coefficients[[2]],
                               deming.reg.slope.lower.ci = deming_model_log10$ci[[2]],
                               deming.reg.slope.upper.ci = deming_model_log10$ci[[2,2]],
                               deming.reg.intercept = deming_model_log10$coefficients[[1]],
                               deming.reg.intercept.lower.ci = deming_model_log10$ci[[1]],
                               deming.reg.intercept.upper.ci = deming_model_log10$ci[[1,2]],)
  
  if(plot == T){
  # plot it with log ticks 1group, 2test_stat, 3durationd, 4endpoint
  plot <- ggplot2::ggplot(plot_dat) +
    ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
    ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
    ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
    ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
    ggplot2::geom_point(shape = 21, fill = "red", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(mg/L)"),
                  y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(mg/L)"),
                  subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::annotation_logticks()
  
  # setup output list
  out$orthreg_df <- orthreg_df
  out$plot <- plot
  
  # return list with data and plot
  return(out)
  }
  else{
    # return orthReg df
    return(orthreg_df)
  }
  }
} 

#============================================================#
# Pairwise orthogonal Regression function
#============================================================#
pwOrthReg <- function(data, group, limit.comp = NULL, min.n = 5, message = F, QC = "ignore", plot = F){
  # fix error if group is not "group"
  if(group != "group"){
    data = data %>%
      dplyr::select(-group)
  }
  # lets get an id for the pairs and filter if observations are less than min.n
  dat.id <- data %>%
    dplyr::rename_with(., ~ str_replace(.x, pattern = group, replacement = "group"), matches(group)) %>%
    #dplyr::rename_at(vars(matches(group)), ~ "group") %>%
    dplyr::mutate(id = paste(group, test_statistic, duration_d, endpoint, sep = ":")) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(n.test = length(unique(cas))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n.test >= min.n) %>%
    dplyr::select(-n.test)
  
  # setup unique pairs of group, test_statistic, duration, endpoint
  unique.pairs <- dat.id %>%
    dplyr::distinct(id)
  
  pairs.list <- combn(x = unique.pairs$id, m = 2, simplify = F)
  
  if(!is.null(limit.comp)){
    limits <- NULL
    # find unwanted comps
    for(j in 1:length(pairs.list)){
      if((grepl(pairs.list[[j]][1], pattern = limit.comp) == T & grepl(pairs.list[[j]][2], pattern = limit.comp) == T) |
         (grepl(pairs.list[[j]][1], pattern = limit.comp) == F & grepl(pairs.list[[j]][2], pattern = limit.comp) == F) == T){
        limits <- rlist::list.append(limits, j)
      }
      if((grepl(pairs.list[[j]][1], pattern = limit.comp) == F & grepl(pairs.list[[j]][2], pattern = limit.comp) == T)){
        pairs.list[[j]] <- rev(pairs.list[[j]])
      }
    }
    # dump unwanted comps
    if(!is.null(limits)) {
      pairs.list <- pairs.list[-limits]
    }
  }
  
  # STOP function if filtered pair list is length 0
  if(length(pairs.list) == 0){
    stop("all pairs filtered, check input data")
  }
  
  # further prune the pair list to make sure similar drugs tested in both groups
  limits2 <- NULL
  for(j in 1:length(pairs.list)){
    dx <- dat.id %>%
      dplyr::filter(id %in% pairs.list[[j]][1]) %>%
      dplyr::distinct(cas)
    dy <- dat.id %>%
      dplyr::filter(id %in% pairs.list[[j]][2]) %>%
      dplyr::distinct(cas)
    d.shared <- sum(dx$cas %in% dy$cas)
    # collect unwanted comps
    if(d.shared < min.n) {
      limits2 <- rlist::list.append(limits2, j)
    }
  }
  # dump unwanted comps b/c of too few common drugs
  if(!is.null(limits2)) {
    pairs.list <- pairs.list[-limits2]
  }
  
  # Create the progress bar
  pb <- progress::progress_bar$new(total = length(pairs.list))
  
  # make a safe funciton for looping
  orthReg.safe <- purrr::safely(orthReg)
  normReg.safe <- purrr::safely(normReg)
  
  # setup a list to hold outputs
  orthReg.list <- NULL
  normReg.list <- NULL
  
  # add plot list if needed
  if(plot == TRUE) {
  plot.list <- NULL
  }
  
  # fix error if group is not based on group - ugly fix
  if(group != "group") {
    # add fix column
    dat.id <- dat.id %>%
      dplyr::mutate(fix = NA_character_)
    # rename it to group by position
    rn.num <- ncol(dat.id)
    names(dat.id)[13] <- group
  }
    
  # loop over pairs
  for(i in 1:length(pairs.list)){
    # get the pair data we care about
    pair.dat <- dat.id %>%
      dplyr::filter(id %in% pairs.list[[i]])
    
    # setup x and y, set NA's to 0 - NEED TO BE CERTAIN THIS IS KOSHER 
    x <- stringr::str_split_fixed(pairs.list[[i]][1], n = 4, pattern = ":")
    y <- stringr::str_split_fixed(pairs.list[[i]][2], n = 4, pattern = ":")
    
    # make a message for running
    if(message == T){
      message(glue::glue("running orthReg and normReg for {i} of {length(pairs.list)} pairs - {paste0(pairs.list[[i]], collapse =' ')}"))
    }
    if(message == F){
      # Update the progress bar
      pb$tick()
    }
    
    if(plot == F) {
    # perform orthogonal regression without plotting
    orthReg.out <- orthReg.safe(data = pair.dat, x = x, y = y, QC = QC, plot = F)
    normReg.out <- normReg.safe(data = pair.dat, x = x, y = y, QC = QC)
    
    # handle orthReg.safe output with errors
    if(is.null(orthReg.out$result)){
      orthReg.out$result <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                           y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                           orth.reg.n.observations = NA_real_,
                                           orth.reg.slope = NA_real_,
                                           orth.reg.intercept = NA_real_,
                                           orth.reg.ssq = NA_real_,
                                           orth.reg.mse = NA_real_,
                                           orth.reg.r.squared = NA_real_)
    }
    # handle normReg.safe output with errors
    if(is.null(normReg.out$result)){
      normReg.out$result <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                           y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                           reg.n.observations = NA_real_,
                                           reg.slope = NA_real_,
                                           reg.intercept = NA_real_,
                                           reg.r.squared = NA_real_,
                                           reg.r.pearson = NA_real_,
                                           reg.r.pvalue = NA_real_,
                                           reg.df = NA_real_)
    }
    # add to the output list
    orthReg.list[[i]] <- orthReg.out$result
    normReg.list[[i]] <- normReg.out$result
    }
    if(plot == T) {
      # perform orthogonal regression with plotting
      orthReg.out <- orthReg.safe(data = pair.dat, x = x, y = y, QC = QC, plot = T)
      normReg.out <- normReg.safe(data = pair.dat, x = x, y = y, QC = QC)
      
      # handle orthReg.safe output with errors
      if(is.null(orthReg.out$result)){
        orthReg.out$result$orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                             y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                             orth.reg.n.observations = NA_real_,
                                             orth.reg.slope = NA_real_,
                                             orth.reg.intercept = NA_real_,
                                             orth.reg.ssq = NA_real_,
                                             orth.reg.mse = NA_real_,
                                             orth.reg.r.squared = NA_real_)
        orthReg.out$result$plot <- ggplot2::ggplot()
      }
      # handle normReg.safe output with errors
      if(is.null(normReg.out$result)){
        normReg.out$result <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                             y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                             reg.n.observations = NA_real_,
                                             reg.slope = NA_real_,
                                             reg.intercept = NA_real_,
                                             reg.r.squared = NA_real_,
                                             reg.r.pearson = NA_real_,
                                             reg.r.pvalue = NA_real_,
                                             reg.df = NA_real_)
      }
      # add to the output list
      orthReg.list[[i]] <- orthReg.out$result$orthreg_df
      normReg.list[[i]] <- normReg.out$result
      plot.list[[i]] <- orthReg.out$result$plot
    }
  }
  if(plot == FALSE){
  # return
    out <- list(orthregs = orthReg.list, normRegs = normReg.list)
  return(out)
  }
  if(plot == TRUE){
    # return
    out <- list(orthregs = orthReg.list, normRegs = normReg.list, plots = plot.list)
    return(out)
  }
}
