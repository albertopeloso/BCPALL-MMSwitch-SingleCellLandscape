
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

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
                                                 face = "bold", family = "Karla", 
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5,
                                                 angle=90, face = "bold", 
                                                 family = "Karla", size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "bottom",
                     # legend.justification.right = "bottom",
                     legend.title = element_blank())

dx_d15_combined <- readRDS(file = "./seurat_objects_merge_dx_d15/bcpall_merge_dx_d15_seurat_filtered_with_ab_step1.rds")

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
    Sample_Name == "PT15_D15" ~ "sw+",
  )) %>%
  mutate(Timepoint = case_when(
    grepl("DX", Sample_Name) ~ "DX",
    grepl("D15", Sample_Name) ~ "D15",
  ))

DefaultAssay(dx_d15_combined) <- "RNA"
dx_d15_combined <- AddMetaData(dx_d15_combined, metadata = md)

# log normalizaion of UMI counts
dx_d15_combined <- NormalizeData(dx_d15_combined, normalization.method = "LogNormalize")

# PCA
# identify the most variable genes
dx_d15_combined <- FindVariableFeatures(dx_d15_combined,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# scale the data
dx_d15_combined <- ScaleData(dx_d15_combined, features = rownames(dx_d15_combined))

# perform PCA
dx_d15_combined <- RunPCA(dx_d15_combined,npcs = 60)

# evaluating the effects of mitochondrial percentage
# check quartile values
summary(dx_d15_combined$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0478  7.9830 11.3402 13.1906 17.0645 37.9893 
# create a table based on cutoffs:
#  - Low --> below 1st quartile (7.9830)
#  - Medium --> between 1st and median (11.3402)
#  - Medium High --> between median and 3rd quantile (13.1906)
#  - High --> above 3rd quartile (17.0645)

dx_d15_combined$quartile.mt <- 
  cut(dx_d15_combined$percent.mt,
      breaks = c(-Inf, 7.9830, 11.3402, 13.1906, Inf),
      labels = c("Low", "Medium", "Medium High", "High"))
# check if it has been created
head(dx_d15_combined@meta.data)

# check whether the mitochondrial percentage is a source of variation
DimPlot(dx_d15_combined,
        reduction = "pca",
        group.by = "quartile.mt",
        split.by = "quartile.mt")


# first UMAP, with random dims
DefaultAssay(dx_d15_combined) = "RNA"
ElbowPlot(dx_d15_combined, ndims = 50)
DimHeatmap(dx_d15_combined, dims = 35:60, cells = 500, balanced = TRUE)
dx_d15_combined <- RunUMAP(dx_d15_combined, dims = 1:46)
dx_d15_combined <- FindNeighbors(dx_d15_combined, dims = 1:46)
dx_d15_combined <- FindClusters(dx_d15_combined, resolution = 2)
DimPlot(dx_d15_combined, group.by = "Sample_Name") +
  scale_color_manual(values = c("PT06_DX" = "#AC92EB","PT06_D15" = "#5b27d6",
                                "PT07_DX" = "#4FC1E8","PT07_D15" = "#126c8b",
                                "PT08_DX" = "#607f3e","PT08_D15" = "#1f2a14",
                                "PT13_DX" = "#FFCE54","PT13_D15" = "#d39700",
                                "PT14_DX" = "#ED5564","PT14_D15" = "#c71627",
                                "PT15_DX" = "#ED9255","PT15_D15" = "#b05213"),
                     labels = c("PT06_DX" = "BCPALL#06DX","PT06_D15" = "BCPALL#06D15",
                                "PT07_DX" = "BCPALL#07DX","PT07_D15" = "BCPALL#07D15",
                                "PT08_DX" = "BCPALL#08DX","PT08_D15" = "BCPALL#08D15",
                                "PT13_DX" = "BCPALL#09DX","PT13_D15" = "BCPALL#09D15",
                                "PT14_DX" = "BCPALL#10DX","PT14_D15" = "BCPALL#10D15",
                                "PT15_DX" = "BCPALL#11DX","PT15_D15" = "BCPALL#11D15")) +
  nice_theme1

# monocle
dx_d15_combined_joined <- JoinLayers(dx_d15_combined)
dx_d15_combined_joined.cds <- SeuratWrappers::as.cell_data_set(dx_d15_combined_joined)
#
png(file="./merge_plots_dx_d15/umap_color_pts_dx_d15_monocle.png", 
    height = 2700, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined_joined.cds, 
                     color_cells_by = "Sample_Name", 
                     label_cell_groups = F, cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values = c("PT06_DX" = "#d992eb","PT06_D15" = "#AC92EB",
                                "PT07_DX" = "#4FC1E8","PT07_D15" = "#4f75e8",
                                "PT08_DX" = "#94b76c","PT08_D15" = "#607f3e",
                                "PT13_DX" = "#FFCE54","PT13_D15" = "#d39700",
                                "PT14_DX" = "#ED5564","PT14_D15" = "#c71627",
                                "PT15_DX" = "#ED9255","PT15_D15" = "#985726"),
                     labels = c("PT06_DX" = "BCPALL#06_DX","PT06_D15" = "BCPALL#06_D15",
                                "PT07_DX" = "BCPALL#07_DX","PT07_D15" = "BCPALL#07_D15",
                                "PT08_DX" = "BCPALL#08_DX","PT08_D15" = "BCPALL#08_D15",
                                "PT13_DX" = "BCPALL#09_DX","PT13_D15" = "BCPALL#09_D15",
                                "PT14_DX" = "BCPALL#10_DX","PT14_D15" = "BCPALL#10_D15",
                                "PT15_DX" = "BCPALL#11_DX","PT15_D15" = "BCPALL#11_D15")) +
  geom_segment(aes(x=-12, y=-14, xend=-9.7, yend=-14), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-14, xend=-12, yend=-10.9), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

####################### UMAP with color on DX/D+15 #############################

png(file="./merge_plots_dx_d15/umap_color_pts_dx_d15_color_timepoint_monocle.png",
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined_joined.cds, color_cells_by = "Timepoint",
                     label_cell_groups = F, cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values = c("DX" = "#9999ff", "D15" = "#0000e5"),
                     labels = c("DX" = "Diagnosis", "D15" = "Day+15"),
                     breaks = c("DX", "D15")) +
  geom_segment(aes(x=-12, y=-14, xend=-9.7, yend=-14), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-14, xend=-12, yend=-10.9), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

####################### UMAP with color on sw+/sw- #############################

png(file="./merge_plots_dx_d15/umap_color_pts_dx_d15_color_switching_monocle.png", 
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_combined_joined.cds, color_cells_by = "Switch", 
                     label_cell_groups = F, cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values =  c("sw+"="#fa9c32", "sw-"="#46acc8"), 
                     labels = c("mmSWpos", "mmSWneg")) + 
  geom_segment(aes(x=-12, y=-14, xend=-9.7, yend=-14), 
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-14, xend=-12, yend=-10.9), 
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()
