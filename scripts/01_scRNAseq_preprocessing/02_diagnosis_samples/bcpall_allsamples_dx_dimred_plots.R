
setwd('/media/data/lab/bcpall_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ---- read in filtered and merged seurat object ----
dx_combined <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_dimred.rds")

# ---- check precomputed umap
DimPlot(dx_combined, group.by = "Sample_Name")
DimPlot(dx_combined, group.by = "orig.ident")
DimPlot(dx_combined, label = T)

# add metadata
md <- dx_combined@meta.data %>% as.data.frame
md <- md %>%
  mutate(Switch = case_when(
    Sample_Name == "PT01" ~ "sw+",
    Sample_Name == "PT02" ~ "sw+",
    Sample_Name == "PT03" ~ "sw-",
    Sample_Name == "PT04" ~ "sw-",
    Sample_Name == "PT05" ~ "sw+",
    Sample_Name == "PT06_DX" ~ "sw+",
    Sample_Name == "PT07_DX" ~ "sw+",
    Sample_Name == "PT08_DX" ~ "sw-",
    Sample_Name == "PT13_DX" ~ "sw+",
    Sample_Name == "PT14_DX" ~ "sw+",
    Sample_Name == "PT15_DX" ~ "sw+",
  ))

dx_combined <- AddMetaData(dx_combined, metadata = md)

# set nice themes for plotting
nice_theme1 <- theme(plot.title=element_text(face = "bold",
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_text(family = "Helvetica Light"),
                     axis.text.x=element_text(family = "Helvetica Light"),
                     axis.title.x = element_text(family = "Helvetica Light"),
                     axis.title.y = element_text(family = "Helvetica Light"),
                     legend.text = element_text(family = "Helvetica Light"),
                     legend.title = element_text(family = "Helvetica Light",
                                                 face = "bold"))

nice_theme2 <- theme(plot.title=element_text(face = "bold",
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, vjust = 5.5,
                                                 face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5,
                                                 angle=90, face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "bottom",
                     # legend.justification.right = "bottom",
                     legend.title = element_blank())

# ---- UMAP split switch ----
png(file="./merge_plots_dx/umap_split_switch.png",
    height = 2500, width = 2500, units = "px", res = 300)
DimPlot(dx_combined, group.by = "Switch") + 
  nice_theme1 +
  scale_color_manual(values = c("sw-"="#46acc8", "sw+"="#fa9c32"))
dev.off()

# monocle plots
dx_combined_joined <- JoinLayers(dx_combined)
dx_combined.cds <- SeuratWrappers::as.cell_data_set(dx_combined_joined)

png(file="./merge_plots_dx/umap_split_switch_monocle.png", height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_combined.cds, color_cells_by = "Switch", label_cell_groups = F, cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values =  c("sw-"="#46acc8", "sw+"="#fa9c32"), labels = c("mmSWneg", "mmSWpos")) + 
  geom_segment(aes(x=-12, y=-12, xend=-9.7, yend=-12), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.1), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

# UMAP pts
png(file="./merge_plots_dx/umap_color_pts_monocle.png", height = 2800, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_combined.cds, color_cells_by = "Sample_Name", label_cell_groups = F, cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT01" = "#b1e222",
                                "PT02" = "#3066be",
                                "PT03" = "#119da4",
                                "PT04" = "#921961",
                                "PT05" = "#efc8d5",
                                "PT06_DX" = "#d992eb",
                                "PT07_DX" = "#4FC1E8",
                                "PT08_DX" = "#94b76c",
                                "PT13_DX" = "#FFCE54",
                                "PT14_DX" = "#ED5564",
                                "PT15_DX" = "#ED9255"),
                     labels = c("PT01" = "BCPALL#01",
                                "PT02" = "BCPALL#02",
                                "PT03" = "BCPALL#03",
                                "PT04" = "BCPALL#04",
                                "PT05" = "BCPALL#05",
                                "PT06_DX" = "BCPALL#06",
                                "PT07_DX" = "BCPALL#07",
                                "PT08_DX" = "BCPALL#08",
                                "PT13_DX" = "BCPALL#09",
                                "PT14_DX" = "BCPALL#10",
                                "PT15_DX" = "BCPALL#11")) +
  
  geom_segment(aes(x=-12, y=-12, xend=-9.7, yend=-12), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.1), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()
