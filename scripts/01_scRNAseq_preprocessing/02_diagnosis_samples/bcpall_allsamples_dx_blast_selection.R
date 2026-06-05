
setwd('/media/data/lab/bcpall_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ---- read in filtered and merged seurat object ----
dx_combined <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_dimred.rds")

# ---- assign metadata for group ----
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
  )) %>%
  mutate(Timepoint = case_when(
    grepl("DX", Sample_Name) ~ "DX",
    grepl("D15", Sample_Name) ~ "D15",
  ))


dx_combined <- AddMetaData(dx_combined, metadata = md)

# ---- check precomputed umap ----
DimPlot(dx_combined, label = T)

# --- remove cluster 3, 14, 16 - healthy BM cells ----
dx_blasts <- subset(dx_combined, idents = c(3, 14, 16), invert = T)

DefaultAssay(dx_blasts) <- "RNA"

# ---- log normalizaion of UMI counts ----
dx_blasts <- NormalizeData(dx_blasts, normalization.method = "LogNormalize")

# ---- PCA ----
# identify the most variable genes
dx_blasts <- FindVariableFeatures(dx_blasts,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# scale the data
dx_blasts <- ScaleData(dx_blasts, features = rownames(dx_blasts))

# perform PCA
dx_blasts <- RunPCA(dx_blasts,npcs = 60)

# ---- perform umap ----
DefaultAssay(dx_blasts) = "RNA"
dx_blasts <- RunUMAP(dx_blasts, dims = 1:44)
dx_blasts <- FindNeighbors(dx_blasts, dims = 1:44)
dx_blasts <- FindClusters(dx_blasts, resolution = 0.5)

saveRDS(dx_blasts, file = "./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_blasts.rds")

################################## plotting ####################################

# ---- check blasts selection with dimreds ----
DimPlot(dx_blasts, group.by = "Sample_Name")
DimPlot(dx_blasts, group.by = "orig.ident")
DimPlot(dx_blasts, label = T)

# ---- set a nice theme for UMAPs ----
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

# ---- plotting using monocle3 ----
dx_blasts_joined <- JoinLayers(dx_blasts)
dx_blasts.cds <- SeuratWrappers::as.cell_data_set(dx_blasts_joined)

# ---- color switch ----
png(file="./merge_plots_dx_blasts/umap_split_switch_monocle_blasts.png",
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_blasts.cds, color_cells_by = "Switch", label_cell_groups = F,
                     cell_size = 0.5) + 
  nice_theme2 +
  scale_color_manual(values =  c("sw-"="#46acc8", "sw+"="#fa9c32"),
                     labels = c("mmSWneg", "mmSWpos")) + 
  geom_segment(aes(x=-12, y=-12, xend=-9.7, yend=-12),
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.1),
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

# ---- color sample ----- 
png(file="./merge_plots_dx_blasts/umap_color_patients_monocle_blasts.png",
    height = 2500, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_blasts.cds, color_cells_by = "Sample_Name",
                     label_cell_groups = F, cell_size = 0.5) + 
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
  
  geom_segment(aes(x=-12, y=-12, xend=-9.7, yend=-12),
               arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.1),
               arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()
