
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)
library(monocle3)

set.seed(1998)

dx_d15_blasts <- readRDS("./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_seurat_filtered_with_ab_blasts_202507.rds")

##################### Split couples and re-compute DimRed ######################

# PT06
dx_d15_blasts_pt06 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT06_DX","PT06_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt06$Sample_Name)),
             c("PT06_DX","PT06_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt06) <- "RNA"

dx_d15_blasts_pt06 <- FindVariableFeatures(dx_d15_blasts_pt06,
                                      selection.method = "vst",
                                      nfeatures = 2000)

dx_d15_blasts_pt06 <- ScaleData(dx_d15_blasts_pt06, 
                                features = rownames(dx_d15_blasts_pt06))

dx_d15_blasts_pt06 <- RunPCA(dx_d15_blasts_pt06,npcs = 60)

ElbowPlot(dx_d15_blasts_pt06, ndims = 50)
DimHeatmap(dx_d15_blasts_pt06, dims = 1:20, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt06) = "RNA"
dx_d15_blasts_pt06 <- RunUMAP(dx_d15_blasts_pt06, dims = 1:10, n.neighbors = 20)
dx_d15_blasts_pt06 <- FindNeighbors(dx_d15_blasts_pt06, dims = 1:10)
dx_d15_blasts_pt06 <- FindClusters(dx_d15_blasts_pt06, resolution = 1)

DimPlot(dx_d15_blasts_pt06, group.by = "Sample_Name")
DimPlot(dx_d15_blasts_pt06)
FeaturePlot(dx_d15_blasts_pt06, features = "SRGN")

# final colors
dx_d15_blasts_pt06_joined <- JoinLayers(dx_d15_blasts_pt06)
dx_d15_blasts_pt06.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt06_joined)
nice_theme2 <- theme(plot.title=element_text(face = "bold", 
                                             family = "Helvetica Light", 
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, vjust = 5.5, 
                                                 face = "bold", family = "Karla", 
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5, angle=90, 
                                                 face = "bold", family = "Karla", 
                                                 size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "bottom",
                     # legend.justification.right = "bottom",
                     legend.title = element_blank())

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt06_202507.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt06.cds, color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 1.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT06_DX" = "#d992eb","PT06_D15" = "#AC92EB"),
                     labels = c("PT06_DX" = "BCPALL#06_DX","PT06_D15" = "BCPALL#06_D15")) +
  geom_segment(aes(x=-6, y=-11, xend=-4.1, yend=-11), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-6, y=-11, xend=-6, yend=-9),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

# PT07
dx_d15_blasts_pt07 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT07_DX","PT07_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt07$Sample_Name)), 
             c("PT07_DX","PT07_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt07) <- "RNA"

dx_d15_blasts_pt07 <- FindVariableFeatures(dx_d15_blasts_pt07,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt07 <- ScaleData(dx_d15_blasts_pt07, 
                                features = rownames(dx_d15_blasts_pt07))

dx_d15_blasts_pt07 <- RunPCA(dx_d15_blasts_pt07,npcs = 60)

ElbowPlot(dx_d15_blasts_pt07, ndims = 50)
DimHeatmap(dx_d15_blasts_pt07, dims = 1:10, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt07) = "RNA"
dx_d15_blasts_pt07 <- RunUMAP(dx_d15_blasts_pt07, dims = 1:5, n.neighbors = 25)
dx_d15_blasts_pt07 <- FindNeighbors(dx_d15_blasts_pt07, dims = 1:5)
dx_d15_blasts_pt07 <- FindClusters(dx_d15_blasts_pt07, resolution = 0.5)

DimPlot(dx_d15_blasts_pt07, group.by = "Sample_Name")

dx_d15_blasts_pt07_joined <- JoinLayers(dx_d15_blasts_pt07)
dx_d15_blasts_pt07.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt07_joined)

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt07_202507.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt07.cds, color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 1.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT07_DX" = "#4FC1E8","PT07_D15" = "#4f75e8"),
                     labels = c("PT07_DX" = "BCPALL#07_DX","PT07_D15" = "BCPALL#07_D15")) +
  geom_segment(aes(x=-8.5, y=-11, xend=-6.5, yend=-11), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-8.5, y=-11, xend=-8.5, yend=-9.1), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

# PT08
dx_d15_blasts_pt08 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT08_DX","PT08_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt08$Sample_Name)), 
             c("PT08_DX","PT08_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt08) <- "RNA"

dx_d15_blasts_pt08 <- FindVariableFeatures(dx_d15_blasts_pt08,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt08 <- ScaleData(dx_d15_blasts_pt08, 
                                features = rownames(dx_d15_blasts_pt08))

dx_d15_blasts_pt08 <- RunPCA(dx_d15_blasts_pt08,npcs = 60)

ElbowPlot(dx_d15_blasts_pt08, ndims = 50)
DimHeatmap(dx_d15_blasts_pt08, dims = 1:35, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt08) = "RNA"
dx_d15_blasts_pt08 <- RunUMAP(dx_d15_blasts_pt08, dims = 1:10, n.neighbors = 20)
dx_d15_blasts_pt08 <- FindNeighbors(dx_d15_blasts_pt08, dims = 1:10)
dx_d15_blasts_pt08 <- FindClusters(dx_d15_blasts_pt08, resolution = 0.5)

DimPlot(dx_d15_blasts_pt08, group.by = "Sample_Name")

dx_d15_blasts_pt08_joined <- JoinLayers(dx_d15_blasts_pt08)
dx_d15_blasts_pt08.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt08_joined)

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt08_2020507.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt08.cds, color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 2) + 
  nice_theme2 +
  scale_color_manual(values = c("PT08_DX" = "#94b76c","PT08_D15" = "#607f3e"),
                     labels = c("PT08_DX" = "BCPALL#08_DX","PT08_D15" = "BCPALL#08_D15")) +
  geom_segment(aes(x=-6.5, y=-7, xend=-5.25, yend=-7), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-6.5, y=-7, xend=-6.5, yend=-5.3), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

# PT13
dx_d15_blasts_pt13 <- subset(dx_d15_blasts
                             , subset = Sample_Name == c("PT13_DX","PT13_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt13$Sample_Name)), 
             c("PT13_DX","PT13_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt13) <- "RNA"

dx_d15_blasts_pt13 <- FindVariableFeatures(dx_d15_blasts_pt13,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt13 <- ScaleData(dx_d15_blasts_pt13, 
                                features = rownames(dx_d15_blasts_pt13))

dx_d15_blasts_pt13 <- RunPCA(dx_d15_blasts_pt13,npcs = 60)

ElbowPlot(dx_d15_blasts_pt13, ndims = 50)
DimHeatmap(dx_d15_blasts_pt13, dims = 1:35, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt13) = "RNA"
dx_d15_blasts_pt13 <- RunUMAP(dx_d15_blasts_pt13, dims = 1:20, n.neighbors = 10)
dx_d15_blasts_pt13 <- FindNeighbors(dx_d15_blasts_pt13, dims = 1:20)
dx_d15_blasts_pt13 <- FindClusters(dx_d15_blasts_pt13, resolution = 0.5)

DimPlot(dx_d15_blasts_pt13, group.by = "Sample_Name")

dx_d15_blasts_pt13_joined <- JoinLayers(dx_d15_blasts_pt13)
dx_d15_blasts_pt13.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt13_joined)

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt13_202507.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt13.cds, color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 1.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT13_DX" = "#FFCE54","PT13_D15" = "#d39700"),
                     labels = c("PT13_DX" = "BCPALL#09_DX","PT13_D15" = "BCPALL#09_D15")) +
  geom_segment(aes(x=-11, y=-6, xend=-9, yend=-6), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-11, y=-6, xend=-11, yend=-4.65), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

# PT14
dx_d15_blasts_pt14 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT14_DX","PT14_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt14$Sample_Name)),
             c("PT14_DX","PT14_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt14) <- "RNA"

dx_d15_blasts_pt14 <- FindVariableFeatures(dx_d15_blasts_pt14,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt14 <- ScaleData(dx_d15_blasts_pt14, 
                                features = rownames(dx_d15_blasts_pt14))

dx_d15_blasts_pt14 <- RunPCA(dx_d15_blasts_pt14,npcs = 60)

ElbowPlot(dx_d15_blasts_pt14, ndims = 50)
DimHeatmap(dx_d15_blasts_pt14, dims = 1:35, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt14) = "RNA"
dx_d15_blasts_pt14 <- RunUMAP(dx_d15_blasts_pt14, dims = 1:10, n.neighbors = 15)
dx_d15_blasts_pt14 <- FindNeighbors(dx_d15_blasts_pt14, dims = 1:10)
dx_d15_blasts_pt14 <- FindClusters(dx_d15_blasts_pt14, resolution = 0.5)

DimPlot(dx_d15_blasts_pt14, group.by = "Sample_Name")

dx_d15_blasts_pt14_joined <- JoinLayers(dx_d15_blasts_pt14)
dx_d15_blasts_pt14.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt14_joined)

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt14_202507.png",
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt14.cds, color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 1.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT14_DX" = "#ED5564","PT14_D15" = "#c71627"),
                     labels = c("PT14_DX" = "BCPALL#10_DX","PT14_D15" = "BCPALL#10_D15")) +
  geom_segment(aes(x=-7, y=-8, xend=-4.85, yend=-8), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-7, y=-8, xend=-7, yend=-6.1), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

# PT15
dx_d15_blasts_pt15 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT15_DX","PT15_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt15$Sample_Name)), 
             c("PT15_DX","PT15_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt15) <- "RNA"

dx_d15_blasts_pt15 <- FindVariableFeatures(dx_d15_blasts_pt15,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt15 <- ScaleData(dx_d15_blasts_pt15, 
                                features = rownames(dx_d15_blasts_pt15))

dx_d15_blasts_pt15 <- RunPCA(dx_d15_blasts_pt15,npcs = 60)

ElbowPlot(dx_d15_blasts_pt15, ndims = 50)
DimHeatmap(dx_d15_blasts_pt15, dims = 1:35, cells = 500, balanced = TRUE)
DefaultAssay(dx_d15_blasts_pt15) = "RNA"
dx_d15_blasts_pt15 <- RunUMAP(dx_d15_blasts_pt15, dims = 1:10, n.neighbors = 15)
dx_d15_blasts_pt15 <- FindNeighbors(dx_d15_blasts_pt15, dims = 1:10)
dx_d15_blasts_pt15 <- FindClusters(dx_d15_blasts_pt15, resolution = 0.5)

DimPlot(dx_d15_blasts_pt15, group.by = "Sample_Name")

dx_d15_blasts_pt15_joined <- JoinLayers(dx_d15_blasts_pt15)
dx_d15_blasts_pt15.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt15_joined)

png(file="./merge_plots_dx_d15_blasts_couples/umap_pt15_202507.png",
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt15.cds, color_cells_by = "Sample_Name",
                     label_cell_groups = F, cell_size = 1.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT15_DX" = "#ED9255","PT15_D15" = "#985726"),
                     labels = c("PT15_DX" = "BCPALL#11_DX","PT15_D15" = "BCPALL#11_D15")) +
  geom_segment(aes(x=-14.5, y=-8.5, xend=-12, yend=-8.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-14.5, y=-8.5, xend=-14.5, yend=-6.5),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

################################################################################

nice_theme3 <- theme(plot.title=element_text(face = "bold", 
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, vjust = 5.5,
                                                 face = "bold", family = "Karla", 
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5, 
                                                 angle=90, face = "bold", 
                                                 family = "Karla", size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "bottom",
                     legend.title = element_blank())

dx_d15_blasts_joined <- JoinLayers(dx_d15_blasts)
dx_d15_blasts.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_joined)

# color switch
png(file="./merge_plots_dx_d15_blasts/umap_split_switch_monocle_blasts.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts.cds, color_cells_by = "Switch",
                     label_cell_groups = F, cell_size = 0.5) + 
  nice_theme3 +
  scale_color_manual(values =  c("sw-"="#46acc8", "sw+"="#fa9c32"),
                     labels = c("mmSWneg", "mmSWpos")) + 
  geom_segment(aes(x=-13.5, y=-12, xend=-11.2, yend=-12),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + 
  geom_segment(aes(x=-13.5, y=-12, xend=-13.5, yend=-9), 
               arrow = arrow(length=unit(.4, 'cm')), color = "black") +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()
