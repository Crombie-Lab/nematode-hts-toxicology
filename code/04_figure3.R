#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")

#==============================================================================#
# step 1: pull Imager vs boyd EC10 vs AC50 w/ Carbryl and shape for plotting
#==============================================================================#
# Imager vs boyd EC10 vs AC50, leave in Carbaryl for plotting
df2 <- dat %>%
  dplyr::filter(latin_name == "Caenorhabditis elegans" &
                  test_statistic %in% c("EC10", "AC50"))

# set args
# x is a vector of specific latin_name/group, test_statistic, duration_d, endpoint, so is y
x = c("Caenorhabditis elegans", "EC10", "2", "Growth")
y = c("Caenorhabditis elegans", "AC50", "NA", "Growth")

or2_proc <- df2 %>%
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

# reshap for plotting
plot_dat <- or2_proc %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) %>% # keep only complete cases
  dplyr::mutate(fill = ifelse(chem_name == "Carbaryl", "grey", "red"))
# 
# # ORIGINAL REGRESSION BASED ON gm_mean
# orthog_reg_model_log10 <- pracma::odregress(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
# 
# # old orthogonal regression
# wc_OR_df <- tibble::tibble(x = log10(plot_dat$gm_mean_x), y = log10(plot_dat$gm_mean_y))
# pcObject <- princomp(wc_OR_df)
# myLoadings <- unclass(loadings(pcObject))[,1]
# OR.slope <- myLoadings["y"]/myLoadings["x"]
# OR.int <- mean(log10(plot_dat$gm_mean_y))-OR.slope*mean(log10(plot_dat$gm_mean_x))
# orthog_reg_model_r_squared <- as.double(((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5)
# 
# # run the extraction functions
# orthog_reg_model_mse <- mse_odreg(orthog_reg_model_log10)
# 
# # build output for orthogonal regression
# wc_orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
#                                     y = paste(y[1], y[2], y[3], y[4], sep = ":"),
#                                     orth.reg.n.observations = nrow(plot_dat),
#                                     orth.reg.slope = orthog_reg_model_log10$coeff[1],
#                                     orth.reg.intercept = orthog_reg_model_log10$coeff[2],
#                                     orth.reg.ssq = orthog_reg_model_log10$ssq[1],
#                                     orth.reg.mse = orthog_reg_model_mse,
#                                     orth.reg.r.squared = orthog_reg_model_r_squared)
# 
# # plot all data with full orthogonal regression - leave in Carbaryl
# pcomp <- ggplot2::ggplot(plot_dat) +
#   ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
#   ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
#   ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
#   ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
#   ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
#   ggplot2::geom_point(aes(fill = fill), show.legend = F, shape = 21, color = "black") +
#   ggplot2::scale_fill_manual(values = c("grey" = "red", "red" = "red")) +
#   ggplot2::theme_bw() +
#   ggplot2::labs(x = bquote(~"Nematode Imager"~.(x[2])~"(mg/L)"),
#                 y = bquote(~"Nematode COPAS"~.(y[2])~"(mg/L)")) +
#   #ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(mg/L)"),
#   #              y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(mg/L)")) +
#   #subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
#   ggplot2::scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x)),
#     limits = c(0.00075, 500)
#   ) +
#   ggplot2::scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x)),
#     limits = c(0.00075, 500)
#   ) +
#   #annotate(geom = "text", x = 30, y = 0.02, label = glue::glue("carbaryl"), size = 4) +
#   #annotate(geom = "text", x = 3, y = 0.02, label = glue::glue("r = {corr.all}"), size = 3) +
#   ggplot2::annotation_logticks()
# pcomp
#==============================================================================#
# Step 2: take out Carbaryl to run the analysis
#==============================================================================#
df2_noCarb <- dat %>%
  dplyr::filter(latin_name == "Caenorhabditis elegans" &
                  test_statistic %in% c("EC10", "AC50") &
                  chem_name != "Carbaryl")

# # set args
# x is a vector of specific latin_name/group, test_statistic, duration_d, endpoint, so is y
x = c("Caenorhabditis elegans", "EC10", "2", "Growth")
y = c("Caenorhabditis elegans", "AC50", "NA", "Growth")

# shape data and calc geometric mean
or2_noCarb_proc <- df2_noCarb %>%
  dplyr::mutate(duration_d = ifelse(is.na(duration_d), "NA", duration_d)) %>% # THIS STEP HANDELS NAs in duration data
  dplyr::mutate(endpoint = ifelse(is.na(endpoint), "NA", endpoint)) %>% # THIS STEP HANDELS NAs in endpoint data
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

# reshape for plotting and orthogonal regression
noCarb_plot_dat <- or2_noCarb_proc %>%
  dplyr::distinct(chem_name, pair_gen, gm_mean, min, max) %>% # just get geom mean and ranges
  tidyr::pivot_wider(names_from = pair_gen, values_from = c(gm_mean, min, max)) %>% # give us x and a y vars
  dplyr::filter(complete.cases(.)) # keep only complete cases

# ORIGINAL REGRESSION BASED ON gm_mean
orthog_reg_model_log10 <- pracma::odregress(x = log10(noCarb_plot_dat$gm_mean_x), y = log10(noCarb_plot_dat$gm_mean_y))

# # old orthogonal regression
oldOR_df <- tibble::tibble(x = log10(noCarb_plot_dat$gm_mean_x), y = log10(noCarb_plot_dat$gm_mean_y))
pcObject <- princomp(oldOR_df)
myLoadings <- unclass(loadings(pcObject))[,1]
OR.slope <- myLoadings["y"]/myLoadings["x"]
OR.int <- mean(log10(noCarb_plot_dat$gm_mean_y))-OR.slope*mean(log10(noCarb_plot_dat$gm_mean_x))
orthog_reg_model_r_squared <- as.double(((summary(pcObject)$sdev[1]^2)/sum(summary(pcObject)$sdev^2)-.5)/.5)

# run the extraction functions
orthog_reg_model_mse <- mse_odreg(orthog_reg_model_log10)

# build output for orthogonal regression
noCarb_orthreg_df <- tibble::tibble(x = paste(x[1], x[2], x[3], x[4], sep = ":"),
                                    y = paste(y[1], y[2], y[3], y[4], sep = ":"),
                                    orth.reg.n.observations = nrow(noCarb_plot_dat),
                                    orth.reg.slope = orthog_reg_model_log10$coeff[1],
                                    orth.reg.intercept = orthog_reg_model_log10$coeff[2],
                                    orth.reg.ssq = orthog_reg_model_log10$ssq[1],
                                    orth.reg.mse = orthog_reg_model_mse,
                                    orth.reg.r.squared = orthog_reg_model_r_squared)

# plot all data with orthogonal regression from Carbaryl removed
pcomp <- ggplot2::ggplot(plot_dat) +
  ggplot2::aes(x = gm_mean_x, y = gm_mean_y) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.5) +
  ggplot2::geom_abline(slope = orthog_reg_model_log10$coeff[1], intercept = orthog_reg_model_log10$coeff[2], size = 0.5, color = "red") +
  ggplot2::geom_errorbar(aes(ymin = min_y, ymax = max_y), width = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_errorbarh(aes(xmin = min_x, xmax = max_x), height = 0, size = 0.25, color = "grey70") +
  ggplot2::geom_point(aes(fill = fill), show.legend = F, shape = 21, color = "black") +
  ggplot2::scale_fill_manual(values = c("grey" = "grey", "red" = "red")) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = bquote(~"Nematode Imager"~.(x[2])~"(mg/L)"),
                y = bquote(~"Nematode COPAS"~.(y[2])~"(mg/L)")) +
  #ggplot2::labs(x = bquote(~italic(.(x[1]))~"toxicity"~.(x[2])~"(mg/L)"),
  #              y = bquote(~italic(.(y[1]))~"toxicity"~.(y[2])~"(mg/L)")) +
  #subtitle = glue::glue("x={x[1]}_{x[2]}_{x[3]}d_{x[4]}\ny={y[1]}_{y[2]}_{y[3]}d_{y[4]}")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  ggplot2::scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(0.00075, 500)
  ) +
  annotate(geom = "text", x = 30, y = 0.02, label = glue::glue("carbaryl"), size = 4) +
  #annotate(geom = "text", x = 3, y = 0.02, label = glue::glue("r = {corr.all}"), size = 3) +
  ggplot2::annotation_logticks()
pcomp

#===========================================================#
# Figure 3 - B. elegans toxicity data compared with data from
#  the "EnviroTox DB": fish 4-day mortality
#===========================================================#
# Filter to Boyd, Widmayer, and day 4 fish mortality from the EnviroTox DB
# Merge is identical to Table 2 for LC50
df4a <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "FISH" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & duration_d == 4 & test_statistic %in% c("LC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
                          TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                             "Mortality, Mortality",
                                                             "Mortality, Survival") ~ "Mortality",
                                             endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                             endpoint == "Intoxication, Immobile" ~ "Immobility",
                                             TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2)

# # run the orth regressions across all pairs to NEMATODE with QC filter
# or4a_QC <- pwOrthReg(data = df4a, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "filter", plot = T)
# or4adf_QC <- data.table::rbindlist(or4a_QC$orthregs)
# print(glue::glue("slope = {round(or4adf_QC[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf_QC[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf_QC[1]$orth.reg.r.squared, digits = 2)}"))
# or4ap_QC <- or4a_QC$plots[[1]]

# run the orth regressions across all pairs to NEMATODE without QC filter
or4a <- pwOrthReg(data = df4a, group = "group",  limit.comp = "NEMATODE", min.n = 5, QC = "ignore", message = T, plot = T)
or4adf <- data.table::rbindlist(or4a$orthregs)
print(glue::glue("slope = {round(or4adf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf[1]$orth.reg.r.squared, digits = 2)}"))
or4ap <- or4a$plots[[1]]
or4a$plots[[1]]

# # Look at Boyd chems shared with Widmayer and compare to FISH
# wb_overlap_fish <- df4a %>%
#   dplyr::filter(chem_name %in% or4a$plots[[1]]$data$chem_name)
# or4a2 <- pwOrthReg(data = wb_overlap_fish, group = "group",  limit.comp = "NEMATODE", min.n = 5, QC = "ignore", message = T, plot = T)
# or4adf2 <- data.table::rbindlist(or4a2$orthregs)
# or4a2$plots[[3]]
# print(glue::glue("slope = {round(or4adf2[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4adf2[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4adf2[3]$orth.reg.r.squared, digits = 2)}"))
# # need to add option to color QC fail data or filter QC fail data to orthReg function. - DONE, QC filter makes fit worse. Went with best fit and least data filtering
#==============================================================================#
# Figure 3 - C Imager vs INVERTEBRATES 
# The Invert mortality data can sometimes be labelled as immobility or
# immobilization or survival select LC50 and EC50 could both be used b/c of
# immobilization endppoint which is essentially mortality
#==============================================================================#
# Merge is identical to Table 2 for LC50
df4b <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "INVERT" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & test_statistic %in% c("LC50", "EC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
                          TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                             "Mortality, Mortality",
                                                             "Mortality, Survival") ~ "Mortality",
                                             endpoint == "Immobilization: Change in the failure to respond or lack of movement after mechanical stimulation." ~ "Immobility",
                                             endpoint == "Intoxication, Immobile" ~ "Immobility",
                                             TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2) %>%
  dplyr::filter(case_when(source == "EnviroTox DB" & !(endpoint %in% c("Mortality", "Immobility"))  ~ F, # take either Mort or immobilization
                          TRUE ~ T)) %>%
  dplyr::mutate(endpoint = case_when(source == "EnviroTox DB" & endpoint %in% c("Mortality", "Immobility")  ~ "Mortality",
                                     TRUE ~ endpoint), # set both to mort
                test_statistic = case_when(source == "EnviroTox DB"  ~ "LC50",
                                           TRUE ~ test_statistic)) # set test_stat to LC50 for both
# Some improvement to merge EC50 Immobilization to LC50 Mort - 13 to 17 compounds max at LC50 Mort 2day, with improved R^2

# run the orth regressions across all pairs without QC filter
or4b <- pwOrthReg(data = df4b, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, QC = "ignore", plot = T)
or4bdf <- data.table::rbindlist(or4b$orthregs)
or4bp <- or4b$plots[[1]]
# Widmayer vs. LC50 2d INVERT
or4b$plots[[1]]
print(glue::glue("slope = {round(or4bdf[1]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[1]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[1]$orth.reg.r.squared, digits = 2)}"))
# # Widmayer vs. LC50 1d INVERT 
# or4b$plots[[3]]
# print(glue::glue("slope = {round(or4bdf[3]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[3]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[3]$orth.reg.r.squared, digits = 2)}"))
# # Boyd vs. LC50 2d INVERT
# or4b$plots[[9]]
# print(glue::glue("slope = {round(or4bdf[9]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[9]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[9]$orth.reg.r.squared, digits = 2)}"))
# # Boyd vs. LC50 1d INVERT 
# or4b$plots[[11]]
# print(glue::glue("slope = {round(or4bdf[11]$orth.reg.slope, digits = 2)}, y-intercept = {round(or4bdf[11]$orth.reg.intercept, digits = 2)}, r^2 = {round(or4bdf[11]$orth.reg.r.squared, digits = 2)}"))
# # run the orth regressions across all pairs with QC filter, makes it worse for Boyd again, set to ignore.
#==============================================================================#
# Figure 3 - D Widmayer vs Daphnia magna:LC50:2:Mortality 
# The Daphnia magna data are compelling 
#==============================================================================#
# New for this figure
daph <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016") | latin_name == "Daphnia magna") %>%
  dplyr::filter(group == "NEMATODE" | latin_name == "Daphnia magna" | group == "NEMATODE_COPAS1")

# run the orth regressions across all pairs without QC filter
ordaph <- pwOrthReg(data = daph, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 5, message = F, QC = "ignore", plot = T)
ordaph.df <- data.table::rbindlist(ordaph$orthregs)
ordaphp <- ordaph$plots[[5]]
print(glue::glue("slope = {round(ordaph.df[5]$orth.reg.slope, digits = 2)}, y-intercept = {round(ordaph.df[5]$orth.reg.intercept, digits = 2)}, r^2 = {round(ordaph.df[5]$orth.reg.r.squared, digits = 2)}"))
#==============================================================================#
# Put them all together
#==============================================================================#
# put them together
fig3 <- cowplot::plot_grid(pcomp + labs(subtitle = "", x = bquote("Nematode Imager" ~ EC[10] ~ "(mg/L)"), y = bquote("Nematode COPAS" ~ AC[50] ~ "(mg/L)")),
                           or4ap + labs(subtitle = "", x = bquote("Nematode Imager" ~ EC[10] ~ "(mg/L)"), y = bquote("Fish" ~ LC[50] ~ "(mg/L)")), 
                           or4bp + labs(subtitle = "", x = bquote("Nematode Imager" ~ EC[10] ~ "(mg/L)"), y = bquote("Invertebrate" ~ LC[50] ~ "(mg/L)")),
                           ordaphp + labs(subtitle = "", x = bquote("Nematode Imager" ~ EC[10] ~ "(mg/L)"), y = bquote(~italic("Daphnia magna") ~ LC[50] ~ "(mg/L)")),
                           labels = c("A", "B", "C", "D"), align = "vh", ncol = 2)
cowplot::ggsave2(fig3, filename = glue::glue("figures/figure3.png"), width = 7.5, height = 7.5)

#==============================================================================#
# Export data for figure
#==============================================================================#
# make it
fig3_df <- noCarb_orthreg_df %>%
  dplyr::select(x:orth.reg.r.squared) %>%
  dplyr::bind_rows(or4a$orthregs[[1]] %>% dplyr::select(x:orth.reg.r.squared),
                   or4b$orthregs[[1]] %>% dplyr::select(x:orth.reg.r.squared),
                   ordaph$orthregs[[5]] %>% dplyr::select(x:orth.reg.r.squared)) %>%
  dplyr::mutate(figure = c("fig3a", "fig3b", "fig3c", "fig3d"), .before =x)
# export it                   
rio::export(fig3_df, file = "data/processed/fig3_data.csv")                   
