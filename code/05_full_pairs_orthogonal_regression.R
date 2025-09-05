#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

# load cleaned data
dat <- data.table::fread("data/processed/00_data.csv")
#==============================================================================#
# Part 1: loop over group pairs with the pairwise orthogonal regression
# function. This function will loop over all pairwise comparisons found in the
# data. There are over 2 million unique pairs if we consider species.
# If we just consider groups there are ~160,000 pairs. 
#==============================================================================#
pw.orth.reg.list <- pwOrthReg(data = dat, group = "group", min.n = 5, message = F, QC = "ignore")
pw.orth.reg.df <- data.table::rbindlist(pw.orth.reg.list$orthregs) %>%
  cbind(data.table::rbindlist(pw.orth.reg.list$normRegs)[,3:9])
# Save orth reg data b/c it takes a long time to generate - 1min per 1K pairs
save(pw.orth.reg.df, file = glue::glue("data/processed/{today}.pw.orth.reg.df.rda"))

#load orth reg data (UPDATE TO MOST RECENT IF NEEDED)
#load("data/processed/20250315.pw.orth.reg.df.rda") # initial data
#load("data/processed/20250519.pw.orth.reg.df.rda") # second iteration
load("data/processed/20250723.pw.orth.reg.df.rda") # final iteration)
#==============================================================================#
# Part 2: Summarize the results from pwOrthReg with over 10 observations?
#==============================================================================#
# filter the pairwise regressions to those with >= 5 observations
pw.orth.reg.df.proc <- pw.orth.reg.df %>%
  dplyr::filter(orth.reg.n.observations >= 10) 

# reshape for plotting orth reg fit parameter distributions
pw.orth.reg.df.proc.long <- pw.orth.reg.df.proc %>%
  dplyr::select(x:orth.reg.r.squared) %>%
  tidyr::pivot_longer(cols = c(orth.reg.n.observations:orth.reg.r.squared))

# now plot it
p1 <- ggplot(pw.orth.reg.df.proc.long) +
  aes(x = value) +
  geom_histogram(bins = 50) +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  labs(x = "", subtitle = glue::glue("{nrow(pw.orth.reg.df.proc)} of {nrow(pw.orth.reg.df)} possible orthogonal regressions\nwith >= 10 observations"))
p1  

# save the plot
cowplot::ggsave2(p1, filename = glue::glue("plots/{today}_orth.reg.model.stat.dist.png"), width = 7.5, height = 7.5)

#==============================================================================#
# Part 3: Prioritize elegans comps
#==============================================================================#
n_orthreg <- pw.orth.reg.df %>%
  dplyr::filter(grepl(x, pattern = "NEMATODE") | grepl(y, pattern = "NEMATODE")) %>%
  dplyr::arrange(orth.reg.slope, orth.reg.n.observations)

#==============================================================================#
# Part 4: Prioritize elegans comps with species
#==============================================================================#
pwol <- pwOrthReg(data = dat, group = "latin_name", limit.comp = "Caenorhabditis elegans", min.n = 5, message = F, QC = "ignore")
pwol.df <- data.table::rbindlist(pwol$orthregs) %>%
  cbind(data.table::rbindlist(pwol$normRegs)[,3:9])
# Save orth reg data b/c it takes a long time to generate - 1min per 1K pairs
save(pwol.df, file = glue::glue("data/processed/{today}.pw.orth.reg.species.df.rda"))

#load orth reg data (UPDATE TO MOST RECENT IF NEEDED)
#load("data/processed/20250315.pw.orth.reg.species.df.rda") # initial data
#load("data/processed/20250519.pw.orth.reg.species.df.rda") # second iteration
load("data/processed/20250723.pw.orth.reg.species.df.rda") # final iteration)

#==============================================================================#
# Part 5: Export all comparisons
#==============================================================================#
rio::export(pw.orth.reg.df, file = "data/processed/pw_orth_reg_df_by_group.csv")
rio::export(pwol.df, file = "data/processed/pw_orth_reg_df_by_species.csv")
