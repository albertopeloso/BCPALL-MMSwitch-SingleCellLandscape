
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')
library(data.table)
library(Seurat)
library(Matrix)
library(harmony)
library(dplyr)
library(ggplot2)

dx_d15_blasts_pos <- readRDS('./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_swpos_seurat_filtered_with_ab_blasts_harmony_projection_202509.rds')

# import from scenic analysis
auc_mtx <- as.matrix(fread("./scenic_v2/merge_dx_d15_blasts_sw_pos_harmony_proj_auc_mtx.csv"))
rownames(auc_mtx) <- auc_mtx[,1]
auc_mtx <- auc_mtx[,-1]
dimnames <- dimnames(auc_mtx)
auc_mtx <- matrix(as.numeric(auc_mtx), dimnames = dimnames, ncol = ncol(auc_mtx))
auc_mtx <- as(auc_mtx, "sparseMatrix")
dx_d15_blasts_pos[["scenic"]] <- CreateAssayObject(data =  auc_mtx)
DefaultAssay(dx_d15_blasts_pos) <- "scenic"

dx_d15_blasts_pos <- FindVariableFeatures(dx_d15_blasts_pos,
                                  selection.method = "vst",
                                  nfeatures = 500)
dx_d15_blasts_pos <- ScaleData(dx_d15_blasts_pos,
                               features = rownames(dx_d15_blasts_pos))
dx_d15_blasts_pos <- RunPCA(dx_d15_blasts_pos, reduction.name = "pca.scenic")
dx_d15_blasts_pos <- FindNeighbors(dx_d15_blasts_pos, reduction = "pca.scenic",
                                   dims = 1:20,
                                   graph.name = c("scenic_nn", "scenic_snn"))
dx_d15_blasts_pos <- FindClusters(dx_d15_blasts_pos, resolution = 0.5, 
                                  graph.name = "scenic_snn")
dx_d15_blasts_pos <- RunUMAP(dx_d15_blasts_pos, reduction = "pca.scenic",
                             dims = 1:20, reduction.name = "umap.scenic")

#----------- top regulons by cluster ---------------#

# Get top regulons
DefaultAssay(dx_d15_blasts_pos) <- "scenic"
Idents(dx_d15_blasts_pos) <- "scenic_snn_res.0.5"

# Find all markers
all_markers <- FindAllMarkers(
  dx_d15_blasts_pos,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Get top 10 regulons per cluster
top_regulons <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

patient_colors_subset <- c(
  "BCPALL#06_DX" = "#d992eb",
  "BCPALL#06_D15" = "#AC92EB",
  "BCPALL#07_DX" = "#4FC1E8",
  "BCPALL#07_D15" = "#4f75e8",
  "BCPALL#09_DX" = "#FFCE54",
  "BCPALL#09_D15" = "#d39700",
  "BCPALL#10_DX" = "#ED5564",
  "BCPALL#10_D15" = "#c71627",
  "BCPALL#11_DX" = "#ED9255",
  "BCPALL#11_D15" = "#985726"
)

dx_d15_blasts_pos_subset <- subset(
  dx_d15_blasts_pos,
  subset = Sample_Name %in% c(
    "PT06_DX","PT06_D15",
    "PT07_DX","PT07_D15",
    "PT13_DX","PT13_D15",
    "PT14_DX","PT14_D15",
    "PT15_DX","PT15_D15"
  )
)

labels_subset <- c(
  "PT06_DX" = "BCPALL#06_DX",
  "PT06_D15" = "BCPALL#06_D15",
  "PT07_DX" = "BCPALL#07_DX",
  "PT07_D15" = "BCPALL#07_D15",
  "PT13_DX" = "BCPALL#09_DX",
  "PT13_D15" = "BCPALL#09_D15",
  "PT14_DX" = "BCPALL#10_DX",
  "PT14_D15" = "BCPALL#10_D15",
  "PT15_DX" = "BCPALL#11_DX",
  "PT15_D15" = "BCPALL#11_D15"
)

dx_d15_blasts_pos_subset$Sample_Name <- factor(
  dx_d15_blasts_pos_subset$Sample_Name,
  levels = names(labels_subset),
  labels = labels_subset
)

png("./scenic_v2/heatmap_all_samples.png",
    width = 3600, height = 3600, res = 300, units = "px")
DoHeatmap(
  dx_d15_blasts_pos_subset, 
  features = unique(top_regulons$gene),
  group.by = "Sample_Name",
  group.colors = patient_colors_subset,
  size = 3,
  angle = 45
) + NoLegend() +
theme(plot.margin = unit(c(1,3,1,1), "cm"))
dev.off()
