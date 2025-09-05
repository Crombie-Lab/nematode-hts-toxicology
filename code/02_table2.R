#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")

# make the comparisons in table 2
#===========================================================#
# Get FISH data, prioritize Fish 96 hr LC50
#===========================================================#
FISH <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "FISH" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & duration_d == 4 & test_statistic %in% c("LC50") ~ T, # removing EC50s has a big filtering effect, even though some are coded as mortality
                          TRUE ~ F)) %>%
  dplyr::mutate(endpoint2 = dplyr::case_when(endpoint %in% c("Mortality/Growth",
                                                             "Mortality, Mortality",
                                                             "Mortality, Survival") ~ "Mortality",
                                             TRUE ~ endpoint)) %>%
  dplyr::mutate(endpoint = endpoint2) %>%
  dplyr::select(-endpoint2)

# run the orth regressions across all pairs to NEMATODE with QC filter. QC filter makes orth reg worse. # QC makes it worse for reg regression too SO DROP QC filter, set to ignore
FISHreg <- pwOrthReg(data = FISH, group = "group",  limit.comp = "NEMATODE", min.n = 5, message = T, QC = "ignore", plot = T)
FISHreg_df <- data.table::rbindlist(FISHreg$orthregs) %>%
  cbind(data.table::rbindlist(FISHreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 3),
                or.int = round(orth.reg.intercept, digits = 3)) %>%
  dplyr::select(x, y,
                or.slope,
                or.int,
                or.r2 = orth.reg.r.squared,
                lr.slope = reg.slope,
                lr.int = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df)
#===========================================================#
# Get INVERTs, prioritize 48 hr LC50
#===========================================================#
INVERT <- dat %>%
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
INVERTreg <- pwOrthReg(data = INVERT, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, QC = "ignore", plot = T)
INVERTreg_df <- data.table::rbindlist(INVERTreg$orthregs) %>%
  cbind(data.table::rbindlist(INVERTreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 3),
                or.int = round(orth.reg.intercept, digits = 3)) %>%
  dplyr::select(x, y,
                or.slope,
                or.int,
                or.r2 = orth.reg.r.squared,
                lr.slope = reg.slope,
                lr.int = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df)

#===========================================================#
# Get ALGAE
#===========================================================#
# OECD test 201 is 3day biomass https://www.oecd.org/en/publications/test-no-201-alga-growth-inhibition-test_9789264069923-en.html
ALGAE <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "EnviroTox DB")) %>%
  dplyr::filter(group == "NEMATODE" | group == "ALGAE" | group == "NEMATODE_COPAS1") %>%
  dplyr::filter(case_when(group == "NEMATODE" ~ T,
                          group == "NEMATODE_COPAS1" ~ T,
                          source == "EnviroTox DB" & test_statistic == "EC50" & endpoint == "Population"  ~ T, 
                          TRUE ~ F))

# run the orth regressions across all pairs without QC filter
ALGAEreg <- pwOrthReg(data = ALGAE, group = "group", limit.comp = "NEMATODE", min.n = 5, message = F, QC = "ignore", plot = T)
ALGAEreg_df <- data.table::rbindlist(ALGAEreg$orthregs) %>%
  cbind(data.table::rbindlist(ALGAEreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 3),
                or.int = round(orth.reg.intercept, digits = 3)) %>%
  dplyr::select(x, y,
                or.slope,
                or.int,
                or.r2 = orth.reg.r.squared,
                lr.slope = reg.slope,
                lr.int = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df)

#------------------------------------------------------------------------------#
# ICE - Mammal data
#------------------------------------------------------------------------------#
mammal <- dat %>%
  dplyr::filter(source %in% c("Widmayer et al. 2022", "Boyd et al. 2016", "NIEHS_ICE"))

mammalreg <- pwOrthReg(data = mammal, group = "group",  limit.comp = "NEMATODE", min.n = 3, message = F, QC = "ignore", plot = T)
mammalreg_df <- data.table::rbindlist(mammalreg$orthregs) %>%
  cbind(data.table::rbindlist(mammalreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 3),
                or.int = round(orth.reg.intercept, digits = 3)) %>%
  dplyr::select(x, y,
                or.slope,
                or.int,
                or.r2 = orth.reg.r.squared,
                lr.slope = reg.slope,
                lr.int = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df)

#------------------------------------------------------------------------------#
# Add Padilla Zebrafish embryo data
#------------------------------------------------------------------------------#
zf <- dat %>%
  dplyr::filter(group %in% c("NEMATODE", "TC_ZF_Padilla", "NEMATODE_COPAS1")) %>%
  dplyr::filter(case_when(group == "TC_ZF_Padilla" & QC == "PASS" ~ T,
                          group == "NEMATODE" | group == "NEMATODE_COPAS1"~ T,
                          TRUE ~ F))

zfreg <- pwOrthReg(data = zf, group = "group",  limit.comp = "NEMATODE", min.n = 3, message = F, QC = "ignore", plot = T)
zfreg_df <- data.table::rbindlist(zfreg$orthregs) %>%
  cbind(data.table::rbindlist(zfreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 3),
                or.int = round(orth.reg.intercept, digits = 3)) %>%
  dplyr::select(x, y,
                or.slope,
                or.int,
                or.r2 = orth.reg.r.squared,
                lr.slope = reg.slope,
                lr.int = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df)

#===========================================================#
# Join groups and do full pairwise
#===========================================================#
# setup all
ALL <- dplyr::bind_rows(FISH, INVERT, ALGAE, mammal, zf) %>%
  dplyr::distinct(.keep_all = T) # get distinct rows

# run all regressions
ALLreg <- pwOrthReg(data = ALL, group = "group", limit.comp = NULL, min.n = 5, message = F, QC = "ignore", plot = T)

# make table 2 data
ALLreg_df <- data.table::rbindlist(ALLreg$orthregs) %>%
  cbind(data.table::rbindlist(ALLreg$normRegs)[,3:9]) %>%
  dplyr::arrange(desc(x), desc(y), desc(orth.reg.n.observations)) %>%
  dplyr::mutate(or.slope = round(orth.reg.slope, digits = 2),
                or.slope.95ci = glue::glue("({round(deming.reg.slope.lower.ci, digits = 2)} to {round(deming.reg.slope.upper.ci, digits = 2)})"),
                or.int = round(orth.reg.intercept, digits = 2),
                or.int.95ci = glue::glue("({round(deming.reg.intercept.lower.ci, digits = 2)} to {round(deming.reg.intercept.upper.ci, digits = 2)})"),
                or.r2 = round(orth.reg.r.squared,digits = 2),
                reg.r.pearson = round(reg.r.pearson, digits = 2),
                reg.r.pvalue = round(reg.r.pvalue, digits = 2)) %>%
  dplyr::select(x, y,
                `OR Slope` = or.slope,
                or.slope.95ci,
                `OR Intercept` = or.int,
                or.int.95ci,
                `OR R2` = or.r2,
                Slope = reg.slope,
                Intercept = reg.intercept,
                `Pearson's r` = reg.r.pearson,
                `p value` = reg.r.pvalue,
                `Deg. freedom` = reg.df) %>%
  dplyr::mutate(N = `Deg. freedom` + 2) # calculating slope and intercept so add 2 to get N

#==================================================#
# Summarize range of groups
#==================================================#
sumy <- ALL %>%
  dplyr::group_by(group, test_statistic, effect_unit) %>%
  dplyr::summarise(min_value = min(effect_value),
                   max_value = max(effect_value))
#===========================================================#
# Create Table 2
#===========================================================#
# Make table 2
tab2_df <- ALLreg_df %>%
  dplyr::mutate(keep = dplyr::case_when(x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "NEMATODE:EC10:2:Growth" ~ T,
                                        x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "FISH:LC50:4:Mortality" ~ T,
                                        x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "INVERT:LC50:2:Mortality" ~ T,
                                        x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "ALGAE:EC50:4:Population" ~ T,
                                        x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        x == "NEMATODE_COPAS1:AC50:NA:Growth" & y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "NEMATODE_COPAS1:AC50:NA:Growth" ~ T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "FISH:LC50:4:Mortality" ~ T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "INVERT:LC50:2:Mortality" ~ T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "ALGAE:EC50:4:Population" ~ T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        x == "NEMATODE:EC10:2:Growth" & y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~T,
                                        x == "FISH:LC50:4:Mortality" & y == "ALGAE:EC50:4:Population" ~ T,
                                        x == "FISH:LC50:4:Mortality" & y == "INVERT:LC50:2:Mortality" ~ T,
                                        x == "FISH:LC50:4:Mortality" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        x == "FISH:LC50:4:Mortality" & y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~ T,
                                        x == "INVERT:LC50:2:Mortality" & y == "ALGAE:EC50:4:Population" ~ T,
                                        x == "INVERT:LC50:2:Mortality" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        x == "INVERT:LC50:2:Mortality" & y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~ T,
                                        x == "ALGAE:EC50:4:Population" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        x == "ALGAE:EC50:4:Population" & y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~ T,
                                        x == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" & y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ T,
                                        TRUE ~ F)) %>%
  dplyr::mutate(`Species group 1` = case_when(x == "NEMATODE_COPAS1:AC50:NA:Growth" ~ "Nematode COPAS",
                                              x == "NEMATODE:EC10:2:Growth" ~ "Nematode Imager",
                                              x == "FISH:LC50:4:Mortality" ~ "Fish",
                                              x == "INVERT:LC50:2:Mortality" ~ "Invertebrate",
                                              x == "ALGAE:EC50:4:Population" ~ "Algae",
                                              x == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~ "Rat",
                                              x == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ "ZF embryo"),
                `Species group 2` = case_when(y == "NEMATODE_COPAS1:AC50:NA:Growth" ~ "Nematode COPAS",
                                              y == "NEMATODE:EC10:2:Growth" ~ "Nematode Imager",
                                              y == "FISH:LC50:4:Mortality" ~ "Fish",
                                              y == "INVERT:LC50:2:Mortality" ~ "Invertebrate",
                                              y == "ALGAE:EC50:4:Population" ~ "Algae",
                                              y == "NIEHS_RAT:LD50:NA:Acute Oral Toxicity" ~ "Rat",
                                              y == "TC_ZF_Padilla:AC50:6:TERATOSCORE" ~ "ZF embryo")) %>%
  dplyr::filter(keep == T) %>%
  dplyr::mutate(`Species group 1` = factor(`Species group 1`, levels = c("Nematode Imager",
                                                                         "Nematode COPAS",
                                                                         "Fish",
                                                                         "Invertebrate",
                                                                         "Algae",
                                                                         "Rat")),
                `Species group 2` = factor(`Species group 2`, levels = c("Nematode COPAS",
                                                                         "Algae",
                                                                         "Fish",
                                                                         "Invertebrate",
                                                                         "Rat",
                                                                         "ZF embryo"))) %>%
  dplyr::arrange(`Species group 1`, `Species group 2`) %>%
  dplyr::select(`Species group 1`, `Species group 2`, `OR Slope`:N, -`Deg. freedom`)

# make option with 95CI
option <- tab2_df %>%
  dplyr::mutate(`OR Slope 95% CI` = paste(`OR Slope`, or.slope.95ci, sep = " "),
                `OR Intercept 95% CI` = paste(`OR Intercept`, or.int.95ci, sep = " ")) %>%
  dplyr::select(`Species group 1`, `Species group 2`,
                `OR Slope 95% CI`, `OR Intercept 95% CI`,
                `OR R2`:N)

# output table 2
rio::export(option, file = glue::glue("data/processed/table2.csv"))

