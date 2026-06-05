

setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(Nebulosa)
library(data.table)
library(findPC)
library(harmony)
library(ggtext)

set.seed(1998)

# ---- load filtered object ----
dx_blasts <- readRDS(file = "./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_blasts.rds")

DefaultAssay(dx_blasts) <- "AB"
dx_blasts <- NormalizeData(dx_blasts, normalization.method = "CLR")

############################# UMAP w/ density ##################################

# ---- set a nice theme for plots ----
nice_theme3 <- theme(plot.title=element_markdown(family = "Karla",
                                                 size = 30, hjust = 0.5),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     plot.title.position = "plot",
                     axis.title.x = element_text(hjust = 0.0475, vjust = 4.5,
                                                 face = "bold", family = "Karla",
                                                 size = 16, color = "black"),
                     axis.title.y = element_text(hjust = 0.05, vjust = -3,
                                                 angle=90, face = "bold",
                                                 family = "Karla", size = 16,
                                                 color = "black"),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "right",
                     # legend.justification.right = "bottom",
                     legend.title =  element_text(family = "Karla", size = 22,
                                                  lineheight = 0.2),
                     legend.title.position = "top")

# ---- tiff format ----

# CD10
tiff(file="abseq_plots_blasts_dx_density_tiff/CD10-MME-AHS0051-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD10** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD10-MME-AHS0051-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.0001,0.0015-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD11b
tiff(file="abseq_plots_blasts_dx_density_tiff/CD11b:ICRF44-ITGAM-AHS0184-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD11b** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD11b:ICRF44-ITGAM-AHS0184-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.00015,0.0015-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD14
tiff(file="abseq_plots_blasts_dx_density_tiff/CD14:MPHIP9-CD14-AHS0037-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD14** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD14:MPHIP9-CD14-AHS0037-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.00015,0.0015-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD15
tiff(file="abseq_plots_blasts_dx_density_tiff/CD15-FUT4-AHS0196-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD15** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD15-FUT4-AHS0196-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.0001,0.0015-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD19
tiff(file="abseq_plots_blasts_dx_density_tiff/CD19:SJ25C1-CD19-AHS0030-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD19** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD19:SJ25C1-CD19-AHS0030-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.0001,0.0015-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD20
tiff(file="abseq_plots_blasts_dx_density_tiff/CD20:L27-MS4A1-AHS0143-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD20** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD20:L27-MS4A1-AHS0143-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.00015,0.0015-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD2
tiff(file="abseq_plots_blasts_dx_density_tiff/CD2-CD2-AHS0029-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD2** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD2-CD2-AHS0029-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0025), breaks=c(0+0.00015,0.0025-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33
tiff(file="abseq_plots_blasts_dx_density_tiff/CD33:P67.6-CD33-AHS0171-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD33** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD33:P67.6-CD33-AHS0171-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.002), breaks=c(0+0.00015,0.002-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33 + CD45
tiff(file="abseq_plots_blasts_dx_density_tiff/CD33_CD45_joint_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD33<sup>+</sup> CD45<sup>+</sup>** (pAbO)'
p=Nebulosa::plot_density(dx_blasts, features = c("CD33:P67.6-CD33-AHS0171-pAbO", 
                                                 "CD45-PTPRC-AHS0040-pAbO"),
                         joint = T, method = "wkde",
                         shape=20, size=0.2, combine = F) 
p[[3]]+ 
  scale_color_viridis_c("Joint density\n",option = "mako", direction = 1,
                        begin = 0.2, limits=c(0,0.0000011), 
                        breaks=c(0+0.0000001,0.0000011-0.0000001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD34
tiff(file="abseq_plots_blasts_dx_density_tiff/CD34:8G12-CD34-AHS0182-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD34** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD34:8G12-CD34-AHS0182-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.000075,0.0015-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD371
tiff(file="abseq_plots_blasts_dx_density_tiff/CD371-CLEC12A-AHS0097-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD371** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD371-CLEC12A-AHS0097-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.000075,0.0015-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD38
tiff(file="abseq_plots_blasts_dx_density_tiff/CD38:HB7-CD38-AHS0189-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD38** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD38:HB7-CD38-AHS0189-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.000075,0.0015-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD45
tiff(file="abseq_plots_blasts_dx_density_tiff/CD45-PTPRC-AHS0040-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD45** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD45-PTPRC-AHS0040-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0015), breaks=c(0+0.000075,0.0015-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD58
tiff(file="abseq_plots_blasts_dx_density_tiff_tiff/CD58:1C3-CD58-AHS0237-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD58** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD58:1C3-CD58-AHS0237-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0010), breaks=c(0+0.000075,0.0010-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD3
tiff(file="abseq_plots_blasts_dx_density_tiff/CD3:SK7-CD3E-AHS0033-pAbO_density_wkde.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD3** (pAbO)'
Nebulosa::plot_density(dx_blasts, features = c("CD3:SK7-CD3E-AHS0033-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.008), breaks=c(0+0.0015,0.008-0.0015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# ---- abseq featureplots ----

# CD10
tiff(file="abseq_plots_blasts_dx_tiff/CD10-MME-AHS0051-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD10** (pAbO)'
FeaturePlot(dx_blasts, features = "CD10-MME-AHS0051-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD11b
tiff(file="abseq_plots_blasts_dx_tiff/CD11b:ICRF44-ITGAM-AHS0184-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD11b** (pAbO)'
FeaturePlot(dx_blasts, features = "CD11b:ICRF44-ITGAM-AHS0184-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()


# CD14
tiff(file="abseq_plots_blasts_dx_tiff/CD14:MPHIP9-CD14-AHS0037-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD14** (pAbO)'  
FeaturePlot(dx_blasts, features = "CD14:MPHIP9-CD14-AHS0037-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD15
tiff(file="abseq_plots_blasts_dx_tiff/CD15-FUT4-AHS0196-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD15** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD15-FUT4-AHS0196-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD19
tiff(file="abseq_plots_blasts_dx_tiff/CD19:SJ25C1-CD19-AHS0030-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD19** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD19:SJ25C1-CD19-AHS0030-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD20
tiff(file="abseq_plots_blasts_dx_tiff/CD20:L27-MS4A1-AHS0143-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD20** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD20:L27-MS4A1-AHS0143-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD2
tiff(file="abseq_plots_blasts_dx_tiff/CD2-CD2-AHS0029-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD2** (pAbO)'
FeaturePlot(dx_blasts, features = "CD2-CD2-AHS0029-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33
tiff(file="abseq_plots_blasts_dx_tiff/CD33:P67.6-CD33-AHS0171-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD33** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD33:P67.6-CD33-AHS0171-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD34
tiff(file="abseq_plots_blasts_dx_tiff/CD34:8G12-CD34-AHS0182-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD34** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD34:8G12-CD34-AHS0182-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD371
tiff(file="abseq_plots_blasts_dx_tiff/CD371-CLEC12A-AHS0097-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD371** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD371-CLEC12A-AHS0097-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD38
tiff(file="abseq_plots_blasts_dx_tiff/CD38:HB7-CD38-AHS0189.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD38** (pAbO)'
FeaturePlot(dx_blasts, features = "CD38:HB7-CD38-AHS0189-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD45
tiff(file="abseq_plots_blasts_dx_tiff/CD45-PTPRC-AHS0040-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD45** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD45-PTPRC-AHS0040-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD58
tiff(file="abseq_plots_blasts_dx_tiff/CD58:1C3-CD58-AHS0237-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD58** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD58:1C3-CD58-AHS0237-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD3
tiff(file="abseq_plots_blasts_dx_tiff/CD3:SK7-CD3E-AHS0033-pAbO.tiff",
     height = 2000, width = 2000, units = "px", res=300)
title <- '**CD3** (pAbO)' 
FeaturePlot(dx_blasts, features = "CD3:SK7-CD3E-AHS0033-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-9, y=-12, xend=-6, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-9, y=-12, xend=-9, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()
