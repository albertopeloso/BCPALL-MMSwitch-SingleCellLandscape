
setwd('/media/data/lab/bcpall_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)
library(ggplot2)
library(patchwork)

set.seed(1998)

# ---- load filtered and merged object ----
dx_blasts <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_blasts.rds")
dx_blasts <- JoinLayers(dx_blasts)

sw_objects <- SplitObject(dx_blasts, split.by = "Switch")

sw_pos <- sw_objects$`sw+`
sw_neg <- sw_objects$`sw-`

# --- subset 1000 cells to apply signature ----
sw_pos_subset1000 <- sw_pos[, sample(colnames(sw_pos), size = 1000, replace = F)]
sw_neg_subset1000 <- sw_neg[, sample(colnames(sw_neg), size = 1000, replace = F)]

# ---- set a nice theme for plotting ----
nice_theme3 <- theme(plot.title=element_text(face = "bold",
                                             family = "Karla",
                                             size=18,
                                             hjust=0.5),
                     axis.text.y = element_text(family = "Karla", size = 12),
                     axis.text.x = element_text(family = "Karla", size = 12),
                     axis.line.x = element_line(color = "black"),
                     axis.line.y = element_line(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     axis.title.x = element_text(face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     axis.title.y = element_text(angle=90, face = "bold",
                                                 family = "Karla", size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "right",
                     # legend.justification.right = "bottom",
                     legend.title = element_blank())

################################ SIGNATURES ####################################

# --- hsc signature from Azimuth ----
azimuth_hsc <- c("CRHBP","AVP","MYCT1","BEX1",
                 "NPR3","CRYGD","MSRB3","CD34",
                 "NPDC1","MLLT3")

azimuth_hsc_markers <- list(azimuth_hsc)
sw_pos_subset1000 <- AddModuleScore(sw_pos_subset1000,
                                    features = azimuth_hsc_markers,
                                    name = "hsc_priming")

FeaturePlot(sw_pos_subset1000, features = "hsc_priming1")

sw_neg_subset1000 <- AddModuleScore(sw_neg_subset1000,
                                    features = azimuth_hsc_markers,
                                    name = "azimuth_hsc_priming")

FeaturePlot(sw_neg_subset1000, features = "azimuth_hsc_priming1")

# ---- myeloid lineage signature from ----
# Pellin D, Loperfido M, Baricordi C, et al. A comprehensive single cell 
# transcriptional landscape of human hematopoietic progenitors. 
# Nat Commun. 2019;10(1):2395. Published 2019 Jun 3. 
# doi:10.1038/s41467-019-10291-0
myelo_pellin <- c("LYZ","MS4A6A","ANXA2")
myelo_pellin_markers <- list(myelo_pellin)
sw_pos_subset1000 <- AddModuleScore(sw_pos_subset1000,
                                    features = myelo_pellin_markers,
                                    name = "myelo_priming")
sw_neg_subset1000 <- AddModuleScore(sw_neg_subset1000,
                                    features = myelo_pellin_markers,
                                    name = "myelo_priming")

################################### PRIMING ####################################

# ---- priming calculation on swpos ----
df_hsc_myelo_swpos <- data.frame(hsc = sw_pos_subset1000@meta.data$hsc_priming1,
                                myelo = sw_pos_subset1000@meta.data$myelo_priming1)

# plotting
p1=ggplot(df_hsc_myelo_swpos, aes(x = hsc, y = myelo, fill = ..level..)) +
  stat_density2d(geom = "polygon") +
  scale_fill_viridis_c(option = "mako", begin = 0.2, limits = c(1,20))+
  ggtitle("mmSWpos") +
  theme_classic() +
  ylim(-0.5,0.6) +
  xlim(-0.2,0.3) +
  labs(x = "HSC", y = "Myeloid") +
  nice_theme3 +
  theme(legend.position = "none")

# ---- priming calculation on swneg ----
df_hsc_myelo_swneg <- data.frame(hsc = sw_neg_subset1000@meta.data$hsc_priming1,
                                 myelo = sw_neg_subset1000@meta.data$myelo_priming1)

# plotting
p2=ggplot(df_hsc_myelo_swneg, aes(x = hsc, y = myelo, fill = ..level..)) +
  stat_density2d(geom = "polygon") +
  scale_fill_viridis_c(option = "mako", begin = 0.2, limits = c(1,20),
                       breaks = c(1,5,10,15,20),
                       labels = c("Low", "","","","High")) + 
  ggtitle("mmSWneg") +
  theme_classic() +
  labs(x = "HSC", y = "Myeloid") +
  nice_theme3 +
  ylim(-0.5,0.6) +
  xlim(-0.2,0.3)

png("./priming_on_blasts/priming_azimuth_hsc_pellin_myelo.png", width = 3400, height = 1400, res = 300, units = "px")
cowplot::plot_grid(p1, NULL, p2, ncol = 3, rel_widths = c(1, 0.5, 1))
dev.off()


