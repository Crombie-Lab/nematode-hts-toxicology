#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data and drop group for now - causes eror with pwOrthReg function when assigning latin_name to group
dat <- data.table::fread("data/processed/00_data.csv")
#==================================================#
# Part 1: Summarize range of groups
#==================================================#
sumy <- dat %>%
  dplyr::filter(group == "NEMATODE") %>%
  dplyr::group_by(test_statistic, effect_unit) %>%
  dplyr::summarise(min_value = min(effect_value),
                   max_value = max(effect_value))
#==================================================#
# Part 2: Setup the priority Orth Regs
#==================================================#
# Widmayer vs Oncorhynchus mykiss 96 hr LC50
d1 <- dat %>%
  dplyr::filter((latin_name == "Oncorhynchus mykiss" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs - NEED TO DEBUG
# to debug - debug(pwOrthReg) debug(orthReg)
d1l <- pwOrthReg(data = d1, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis elegans", message = T, plot = T)
d1df <- data.table::rbindlist(d1l$orthregs)
d1plots <- cowplot::plot_grid(plotlist = d1l$plots)

# Widmayer vs Pimephales promelas 96 hr LC50
d2 <- dat %>%
  dplyr::filter((latin_name == "Pimephales promelas" &
                   duration_d == 4 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d2l <- pwOrthReg(data = d2, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d2df <- data.table::rbindlist(d2l$orthregs)
d2plots <- cowplot::plot_grid(plotlist = d2l$plots)


# Widmayer vs RAT NIEHS_RAT:LD50:NA:Acute Oral Toxicity
d3 <- dat %>%
  dplyr::filter((group == "NIEHS_RAT" &
                   test_statistic == "LD50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2))

# run the orth regressions across all pairs
d3l <- pwOrthReg(data = d3, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d3df <- data.table::rbindlist(d3l$orthregs)
d3plots <- cowplot::plot_grid(plotlist = d3l$plots)

# Widmayer vs Zebrafish TC_ZF_Padilla:AC50:6:TERATOSCORE
d4 <- dat %>%
  dplyr::filter((group == "TC_ZF_Padilla" &
                   duration_d == 6 &
                   test_statistic == "AC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) %>%
  dplyr::filter(case_when(group == "TC_ZF_Padilla" & QC == "PASS" ~ T,
                          group == "NEMATODE" | group == "NEMATODE_COPAS1"~ T,
                          TRUE ~ F))

# run the orth regressions across all pairs
d4l <- pwOrthReg(data = d4, group = "latin_name", min.n = 5, limit.comp = "Caenorhabditis", message = T, plot = T)
d4df <- data.table::rbindlist(d4l$orthregs)
d4plots <- cowplot::plot_grid(plotlist = d4l$plots)

# Widmayer vs Daphnia manga Daphnia magna:LC50:2:Mortality
d5 <- dat %>%
  dplyr::filter((latin_name == "Daphnia magna" &
                   duration_d == 2 &
                   test_statistic == "LC50") |
                  (latin_name == "Caenorhabditis elegans" &
                     duration_d == 2)) 

#==================================================#
# Part 3: Package the priority regressions
#==================================================#
proc_df <- dplyr::bind_rows(d1df, d2df, d3df, d4df) %>%
  dplyr::filter(x == "Caenorhabditis elegans:EC10:2:Growth")
