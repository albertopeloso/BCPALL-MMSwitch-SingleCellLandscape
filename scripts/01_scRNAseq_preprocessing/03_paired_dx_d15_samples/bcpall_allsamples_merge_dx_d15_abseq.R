

setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)
library(Nebulosa)
library(ggtext)

set.seed(1998)

dx_d15_combined <- readRDS(file = "./seurat_objects_merge_dx_d15/bcpall_merge_dx_d15_seurat_filtered_with_ab_step1.rds")

DefaultAssay(dx_d15_combined) <- "RNA"

dx_d15_combined <- NormalizeData(dx_d15_combined,
                                 normalization.method = "LogNormalize")

dx_d15_combined <- FindVariableFeatures(dx_d15_combined,
                                    selection.method = "vst",
                                    nfeatures = 2000)
dx_d15_combined <- ScaleData(dx_d15_combined,
                             features = rownames(dx_d15_combined))

dx_d15_combined <- RunPCA(dx_d15_combined,npcs = 60)

DefaultAssay(dx_d15_combined) = "RNA"
ElbowPlot(dx_d15_combined, ndims = 50)
DimHeatmap(dx_d15_combined, dims = 35:60, cells = 500, balanced = TRUE)
dx_d15_combined <- RunUMAP(dx_d15_combined, dims = 1:46)
dx_d15_combined <- FindNeighbors(dx_d15_combined, dims = 1:46)
dx_d15_combined <- FindClusters(dx_d15_combined, resolution = 0.5)
DimPlot(dx_d15_combined, group.by = "Sample_Name")
DimPlot(dx_d15_combined, group.by = "orig.ident")
DimPlot(dx_d15_combined, label = T)

################################################################################

# plots for abseq expression in projection umap
output_dir <- "/media/data/lab/CLL1_sc_final_analyses_2024/abseq_plots_merge_dx_d15/"
abseq_names <- as.vector(dx_d15_combined@assays$AB@counts@Dimnames[[1]])
for(ab in abseq_names) {
  print(ab)
  file_path <- file.path(output_dir, paste0("umap_", ab, ".png"))
  png(file=file_path, height = 2000, width = 2000, units = "px", res=300)
  print(FeaturePlot(dx_d15_combined, features = ab) + 
          ggtitle(ab) +
          scale_color_gradientn(colours = c("grey99","#7ed7afff","#3aaeadff","#38639D")))
  dev.off()
}


########################### UMAP w/ density ###########################

nice_theme3 <- theme(plot.title=element_markdown(family = "Karla",size = 30, hjust = 0.5),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     plot.title.position = "plot",
                     axis.title.x = element_text(hjust = 0.0475, vjust = 4.5,
                                                 face = "bold", family = "Karla",
                                                 size = 16, color = "black"),
                     axis.title.y = element_text(hjust = 0.05, vjust = -3, angle=90,
                                                 face = "bold", family = "Karla",
                                                 size = 16, color = "black"),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "right",
                     # legend.justification.right = "bottom",
                     legend.title =  element_text(family = "Karla",size = 22, lineheight = 0.2),
                     legend.title.position = "top")



DefaultAssay(dx_d15_combined) <- "AB"
Nebulosa::plot_density(dx_d15_combined, features = c("CD14:MPHIP9-CD14-AHS0037-pAbO"), joint = F) + 
  # scale_color_gradientn(colours = rev(c("#c7e8e7","#75C6C5","#3aaeadff","#213B5E")))
  # scale_color_viridis_c(option = "mako", direction = 1, begin = 0.2, end = 0.9)
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.1, limits=c(0,0.06))
# + scale_color_gradientn(colours = c("grey80","#7ed7afff","#3aaeadff","#38639D","#213B5E"))
# + scale_color_viridis_c(option = "mako", direction = -1)


# CD10
png(file="abseq_plots_merge_dx_d15_density/CD10-MME-AHS0051-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD10** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD10-MME-AHS0051-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0020), breaks=c(0+0.0001,0.0020-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD11b
png(file="abseq_plots_merge_dx_d15_density/CD11b:ICRF44-ITGAM-AHS0184-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD11b** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD11b:ICRF44-ITGAM-AHS0184-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.003), breaks=c(0+0.00015,0.003-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()


# CD14
png(file="abseq_plots_merge_dx_d15_density/CD14:MPHIP9-CD14-AHS0037-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD14** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD14:MPHIP9-CD14-AHS0037-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0025), breaks=c(0+0.00015,0.0025-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD15
png(file="abseq_plots_merge_dx_d15_density/CD15-FUT4-AHS0196-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD15** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD15-FUT4-AHS0196-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0020), breaks=c(0+0.0001,0.0020-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD19
png(file="abseq_plots_merge_dx_d15_density/CD19:SJ25C1-CD19-AHS0030-pAbO_density_wkde.png", 
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD19** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD19:SJ25C1-CD19-AHS0030-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2, 
                        limits=c(0,0.0015), breaks=c(0+0.0001,0.0015-0.0001),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD20
png(file="abseq_plots_merge_dx_d15_density/CD20:L27-MS4A1-AHS0143-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD20** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD20:L27-MS4A1-AHS0143-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0025), breaks=c(0+0.00015,0.0025-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD20 + CD45
png(file="abseq_plots_merge_dx_d15_density/CD20_CD45_joint_density_wkde.png",
    height = 2000, width = 2400, 
    units = "px", res=300)
title <- '**CD20<sup>+</sup> CD45<sup>+</sup>** (pAbO)'
p=Nebulosa::plot_density(dx_d15_combined,
                         features = c("CD20:L27-MS4A1-AHS0143-pAbO","CD45-PTPRC-AHS0040-pAbO"),
                         joint = T, method = "wkde",
                         shape=20, size=0.2, combine = F) 

p[[3]]+ 
  scale_color_viridis_c("Joint density\n",option = "mako", direction = 1, 
                        begin = 0.2, limits=c(0,0.0000021), 
                        breaks=c(0+0.00000015,0.0000021-0.00000015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD2
png(file="abseq_plots_merge_dx_d15_density/CD2-CD2-AHS0029-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD2** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD2-CD2-AHS0029-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0025), breaks=c(0+0.00015,0.0025-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD2 + CD45
png(file="abseq_plots_merge_dx_d15_density/CD2_CD45_joint_density_wkde.png",
    height = 2000, width = 2400, 
    units = "px", res=300)
title <- '**CD2<sup>+</sup> CD45<sup>+</sup>** (pAbO)'
p=Nebulosa::plot_density(dx_d15_combined, features = c("CD2-CD2-AHS0029-pAbO","CD45-PTPRC-AHS0040-pAbO"),
                         joint = T, method = "wkde",
                         shape=20, size=0.2, combine = F) 

p[[3]]+ 
  scale_color_viridis_c("Joint density\n",option = "mako", direction = 1,
                        begin = 0.2, limits=c(0,0.000002), 
                        breaks=c(0+0.0000001,0.000002-0.0000001), labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33 + CD45
png(file="./CD33_CD45_joint_density_wkde.png", height = 2000, width = 2400, 
    units = "px", res=300)
title <- '**CD33<sup>+</sup> CD45<sup>+</sup>** (pAbO)'
p=Nebulosa::plot_density(dx_d15_combined, 
                         features = c("CD33:P67.6-CD33-AHS0171-pAbO","CD45-PTPRC-AHS0040-pAbO"),
                         joint = T, method = "wkde",
                         shape=20, size=0.2, combine = F) 

p[[3]]+ 
  scale_color_viridis_c("Joint density\n",option = "mako", direction = 1, 
                        begin = 0.2, limits=c(0,0.000002), 
                        breaks=c(0+0.0000001,0.000002-0.0000001), 
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()


# CD33 + CD45
tiff(file="./CD33_CD45_joint_density_wkde.tiff", height = 2000, width = 2400, 
    units = "px", res=300)
title <- '**CD33<sup>+</sup> CD45<sup>+</sup>** (pAbO)'
p=Nebulosa::plot_density(dx_d15_combined,
                         features = c("CD33:P67.6-CD33-AHS0171-pAbO","CD45-PTPRC-AHS0040-pAbO"),
                         joint = T, method = "wkde",
                         shape=20, size=0.2, combine = F) 

p[[3]]+ 
  scale_color_viridis_c("Joint density\n",option = "mako", direction = 1, 
                        begin = 0.2, limits=c(0,0.000002), 
                        breaks=c(0+0.0000001,0.000002-0.0000001), labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33
png(file="abseq_plots_merge_dx_d15_density/CD33:P67.6-CD33-AHS0171-pAbO_density_wkde.png", 
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD33** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD33:P67.6-CD33-AHS0171-pAbO"),
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.002), breaks=c(0+0.00015,0.002-0.00015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD34
png(file="abseq_plots_merge_dx_d15_density/CD34:8G12-CD34-AHS0182-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD34** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD34:8G12-CD34-AHS0182-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.001), breaks=c(0+0.000075,0.001-0.000075),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD371
png(file="abseq_plots_merge_dx_d15_density/CD371-CLEC12A-AHS0097-pAbO_density_wkde.png",
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD371** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD371-CLEC12A-AHS0097-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2, 
                        limits=c(0,0.0012), breaks=c(0+0.00009,0.0012-0.00009),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD38
png(file="abseq_plots_merge_dx_d15_density/CD38:HB7-CD38-AHS0189-pAbO_density_wkde.png", 
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD38** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD38:HB7-CD38-AHS0189-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2,
                        limits=c(0,0.0011), breaks=c(0+0.00009,0.0011-0.00009),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD45
png(file="abseq_plots_merge_dx_d15_density/CD45-PTPRC-AHS0040-pAbO_density_wkde.png", 
    height = 2000, width = 2200, 
    units = "px", res=300)
title <- '**CD45** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD45-PTPRC-AHS0040-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2, 
                        limits=c(0,0.001), breaks=c(0+0.00009,0.001-0.00009),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-9.7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD58
png(file="abseq_plots_merge_dx_d15_density/CD58:1C3-CD58-AHS0237-pAbO_density_wkde.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD58** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD58:1C3-CD58-AHS0237-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2, 
                        limits=c(0,0.0007), breaks=c(0+0.00005,0.0007-0.00005),
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
png(file="abseq_plots_merge_dx_d15_density/CD3:SK7-CD3E-AHS0033-pAbO_density_wkde.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD3** (pAbO)'
Nebulosa::plot_density(dx_d15_combined, features = c("CD3:SK7-CD3E-AHS0033-pAbO"), 
                       joint = F, method = "wkde",shape=20, size=0.2) + 
  scale_color_viridis_c("Density\n",option = "mako", direction = 1, begin = 0.2, 
                        limits=c(0,0.015), breaks=c(0+0.0015,0.015-0.0015),
                        labels=c("Min", "Max")) +
  nice_theme3 +
  geom_segment(aes(x=-13.5, y=-13.5, xend=-10, yend=-13.5), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-13.5, xend=-13.5, yend=-10), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()


################################################################################

# abseq plots expression values

# CD10
png(file="abseq_plots_merge_dx_d15/CD10-MME-AHS0051-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD10** (pAbO)'
FeaturePlot(dx_d15_combined, features = "CD10-MME-AHS0051-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD11b
png(file="abseq_plots_merge_dx_d15/CD11b:ICRF44-ITGAM-AHS0184-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD11b** (pAbO)'
FeaturePlot(dx_d15_combined, features = "CD11b:ICRF44-ITGAM-AHS0184-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()


# CD14
png(file="abseq_plots_merge_dx_d15/CD14:MPHIP9-CD14-AHS0037-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD14** (pAbO)'  
FeaturePlot(dx_d15_combined, features = "CD14:MPHIP9-CD14-AHS0037-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD15
png(file="abseq_plots_merge_dx_d15/CD15-FUT4-AHS0196-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD15** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD15-FUT4-AHS0196-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD19
png(file="abseq_plots_merge_dx_d15/CD19:SJ25C1-CD19-AHS0030-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD19** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD19:SJ25C1-CD19-AHS0030-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD20
png(file="abseq_plots_merge_dx_d15/CD20:L27-MS4A1-AHS0143-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD20** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD20:L27-MS4A1-AHS0143-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), a
               rrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD2
png(file="abseq_plots_merge_dx_d15/CD2-CD2-AHS0029-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD2** (pAbO)'
FeaturePlot(dx_d15_combined, features = "CD2-CD2-AHS0029-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD33
png(file="abseq_plots_merge_dx_d15/CD33:P67.6-CD33-AHS0171-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD33** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD33:P67.6-CD33-AHS0171-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD34
png(file="abseq_plots_merge_dx_d15/CD34:8G12-CD34-AHS0182-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD34** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD34:8G12-CD34-AHS0182-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD371
png(file="abseq_plots_merge_dx_d15/CD371-CLEC12A-AHS0097-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD371** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD371-CLEC12A-AHS0097-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD38
png(file="abseq_plots_merge_dx_d15/CD38:HB7-CD38-AHS0189.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD38** (pAbO)'
FeaturePlot(dx_d15_combined, features = "CD38:HB7-CD38-AHS0189-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD45
png(file="abseq_plots_merge_dx_d15/CD45-PTPRC-AHS0040-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD45** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD45-PTPRC-AHS0040-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD58
png(file="abseq_plots_merge_dx_d15/CD58:1C3-CD58-AHS0237-pAbO.png",
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD58** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD58:1C3-CD58-AHS0237-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()

# CD3
png(file="abseq_plots_merge_dx_d15/CD3:SK7-CD3E-AHS0033-pAbO.png", 
    height = 2000, width = 2000, 
    units = "px", res=300)
title <- '**CD3** (pAbO)' 
FeaturePlot(dx_d15_combined, features = "CD3:SK7-CD3E-AHS0033-pAbO") + 
  scale_color_viridis_c(option = "mako", direction = 1, begin = 0.25) +
  nice_theme3 +
  geom_segment(aes(x=-12, y=-13, xend=-8.5, yend=-13),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-9),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle(title)
dev.off()
