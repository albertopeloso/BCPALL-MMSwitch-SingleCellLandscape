
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ---- load filtered object ----
dx_d15_combined <- readRDS("./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_seurat_filtered_with_ab_dimred.rds")

# ---- add metadata ----
md <- dx_d15_combined@meta.data %>% as.data.frame
md <- md %>%
  mutate(Switch = case_when(
    Sample_Name == "PT06_DX" ~ "sw+",
    Sample_Name == "PT07_DX" ~ "sw+",
    Sample_Name == "PT08_DX" ~ "sw-",
    Sample_Name == "PT13_DX" ~ "sw+",
    Sample_Name == "PT14_DX" ~ "sw+",
    Sample_Name == "PT15_DX" ~ "sw+",
    Sample_Name == "PT06_D15" ~ "sw+",
    Sample_Name == "PT07_D15" ~ "sw+",
    Sample_Name == "PT08_D15" ~ "sw-",
    Sample_Name == "PT13_D15" ~ "sw+",
    Sample_Name == "PT14_D15" ~ "sw+",
    Sample_Name == "PT15_D15" ~ "sw+"
  )) %>%
  mutate(Timepoint = case_when(
    grepl("DX", Sample_Name) ~ "DX",
    grepl("D15", Sample_Name) ~ "D15",
  ))

dx_d15_combined <- AddMetaData(dx_d15_combined, metadata = md)

# ---- load markers ----

cluster3_T.markers <- read.table(file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_T_markers.txt", 
                                 sep = "\t",
                                 header = TRUE)

cluster14_B.markers <- read.table(file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_B_markers.txt",
                                  sep = "\t",
                                  header = TRUE)

cluster16_myelo_ery.markers <- read.table(file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_myelo_ery_markers.txt",
                                          sep = "\t",
                                          header = TRUE)

# ---- plotting ----

# check precomputed umap
dx_d15_combined <- JoinLayers(dx_d15_combined)
DimPlot(dx_d15_combined, group.by = "Sample_Name")

# set nice themes for plotting
nice_theme2 <- theme(plot.title=element_text(face = "bold",
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, vjust = 6.5,
                                                 face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5,
                                                 angle=90, face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     legend.text = element_text(family = "Karla",size = 14),
                     legend.position = "right",
                     legend.title = element_text(face = "bold", 
                                                 family = "Karla",
                                                 size = 16))

# B-cell score
B_cell <- cluster14_B.markers$Gene.Name[1:10]
B_cell_markers <- list(B_cell)
dx_d15_combined <- AddModuleScore(object = dx_d15_combined, 
                                  features = B_cell_markers,
                                  name = "B_cell_score")

# myelo/ery score
myelo_ery <- cluster16_myelo_ery.markers$Gene.Name[1:10]
myelo_ery_markers <- list(myelo_ery)
dx_d15_combined <- AddModuleScore(object = dx_d15_combined, 
                                  features = myelo_ery_markers, 
                                  name = "myelo_ery_score")

# T-cell score
T_cell <- cluster3_T.markers$Gene.Name[1:10]
T_cell_markers <- list(T_cell)
dx_d15_combined <- AddModuleScore(object = dx_d15_combined, 
                                  features = T_cell_markers, 
                                  name = "T_cell_score")

dx_d15_combined.cds <- SeuratWrappers::as.cell_data_set(dx_d15_combined)

png(file="./merge_plots_dx_d15/umap_b_cell_score_monocle.png", 
    height = 2500, width = 2700, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined.cds, color_cells_by = "B_cell_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "B-cell \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5,
                                    face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5, 
                                    angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-13, xend=-9.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-10.25), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

png(file="./merge_plots_dx_d15/umap_myelo_ery_score_monocle.png", 
    height = 2500, width = 2800, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined.cds, color_cells_by = "myelo_ery_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "Myelo/ery \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, 
                                    face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5, 
                                    angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-13, xend=-9.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-10.25), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

png(file="./merge_plots_dx_d15/umap_t_cell_score_monocle.png", 
    height = 2500, width = 2700, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined.cds, color_cells_by = "T_cell_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "T-cell \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, 
                                    face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5,
                                    angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-13, xend=-9.5, yend=-13), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-13, xend=-12, yend=-10.25), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()
