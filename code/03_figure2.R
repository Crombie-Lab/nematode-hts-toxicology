#!/usr/bin/env Rscript
library(tidyverse)

# set working dir to base directory of repository
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load functions
source("code/functions.R")

#==============================================================================#
# Step 1: Source data from Widmayer 2022, fig 3 and shape
#==============================================================================#
# get the strain colors from widmayer 2022
strain_colors       <- c("blue",   "orange","#5A0C13","#C51B29",  "#a37000","#627264","#67697C","purple")
names(strain_colors) <- c("CB4856","N2", "ECA36",  "ECA396"  ,"CB4855", "RC301",   "MY16", "XZ1516")

# source the widmayer fig3 data
raw <- readr::read_csv(url("https://raw.githubusercontent.com/AndersenLab/toxin_dose_responses/master/manuscript_tables/supp.table.3.csv"))

# shape it
shaped <- raw %>%
  tidyr::pivot_longer(cols = -Toxicant, values_to = "EC10", names_to = "strain") %>%
  tidyr::separate(EC10, into = c("EC10", "se"), sep = " ± ") %>%
  dplyr::mutate(EC10 = as.numeric(EC10),
                se = as.numeric(se),
                upper = EC10 + se,
                lower = EC10 - se,
                Toxicant = ifelse(Toxicant == "2_4-D", "2,4-D", Toxicant),
                class = dplyr::case_when(Toxicant == "Triphenyl phosphate" ~ "Flame Retardant",
                                         Toxicant %in% c("Pyraclostrobin",
                                                         "Mancozeb",
                                                         "Chlorothalonil",
                                                         "Carboxin") ~ "Fungicides",
                                         Toxicant %in% c("Paraquat",
                                                         "Atrazine",
                                                         "2,4-D") ~ "Herbicides",
                                         Toxicant %in% c("Propoxur",
                                                         "Methomyl",
                                                         "Chlorpyrifos",
                                                         "Carbaryl",
                                                         "Aldicarb") ~ "Insecticides",
                                         Toxicant %in% c("Zinc chloride",
                                                         "Silver nitrate",
                                                         "Nickel chloride",
                                                         "Methylmercury chloride",
                                                         "Lead(II) nitrate",
                                                         "Copper(II) chloride",
                                                         "Cadmium chloride",
                                                         "Arsenic trioxide") ~ "Metals"))
# make relative to N2
rel <- shaped %>%
  dplyr::group_by(Toxicant) %>%
  dplyr::mutate(N2 = ifelse(strain == "N2", EC10, NA_real_)) %>%
  tidyr::fill(N2, .direction = "updown") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rel.ec10 = (EC10 / N2) - 1,
                rel.upper = (upper / N2) - 1,
                rel.lower = (lower / N2) - 1)

# plot relative
rel.plot <- rel %>%
  ggplot(., mapping = aes(y = Toxicant, x = rel.ec10, 
                          xmin = rel.lower, 
                          xmax = rel.upper,
                          color = strain)) + 
  theme_bw(base_size = 11) +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.5) + 
  #scale_color_manual(values = strain_colors[c(1,3:8)], name = "Strain") +
  #scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "right") + 
  labs(x = "Relative EC10 (uM)")
rel.plot

# plot actual
abs.plot <- rel %>%
  ggplot(., mapping = aes(y = Toxicant, x = EC10, 
                          xmin = lower, 
                          xmax = upper,
                          color = factor(strain, levels = c("N2", "CB4855", "CB4856", "ECA36", "ECA396", "MY16", "RC301", "XZ1516")))) + 
  theme_bw(base_size = 11) +
  #geom_vline(xintercept = 0, linetype = 3, colour = "orange") + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  size = 0.15) + 
  scale_color_manual(values = strain_colors, name = "Strain") +
  #scale_alpha_manual(values = c(0.15,1), guide = "none") +
  facet_grid(class~., scales = "free", space = "free") + 
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "right") + 
  labs(x = bquote(~ EC[10] ~ "(μM)")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
abs.plot

# save them
#ggsave(rel.plot, filename = "plots/fig3_naturalVaritaionRel.png", width = 7, height = 7)
#ggsave(abs.plot, filename = "plots/fig3_naturalVaritaionAbs.png", width = 7, height = 7)
ggsave(abs.plot, filename = "figures/figure2.png", width = 6.5, height = 6.5)


