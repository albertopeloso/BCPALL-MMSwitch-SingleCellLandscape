
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ------------------------------ load files ---------------------------------- #
dx_combined <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_dimred.rds")

# ---- check precomputed umap ----
# ---- healthy bm cells cluster are 3, 14, 16 ----
DimPlot(dx_combined)

# ---- prepare files for infercnv ----
sample_names <- unique(dx_combined@meta.data$Sample_Name)

for (sample in sample_names) {
  print(paste("Processing sample",sample))
  dx_subset <- subset(dx_combined, Sample_Name == sample)
  # prepare files for infercnv 
  dx_subset <- JoinLayers(dx_subset)
  # Extract cell barcodes and annotations
  annotations_subset <- data.frame(
    Cell = colnames(dx_subset),
    Cluster = case_when(
      dx_subset@meta.data$seurat_clusters == 3 ~ "T cells",
      dx_subset@meta.data$seurat_clusters == 14 ~ "B cells",
      dx_subset@meta.data$seurat_clusters == 16 ~ "Myelo/ery",
      dx_subset@meta.data$seurat_clusters %in% c(0:2, 4:13, 15, 17) ~ 
        paste0("Blasts seu_c",dx_subset@meta.data$seurat_clusters)
      )
    )
  print(paste("Saving",sample))
  annotations_filename <- paste0("dx_", tolower(sample),
                                 "_seurat_cluster_annotation_infercnv.txt")
  write.table(annotations_subset,
              file = paste0("./infercnv/input_single_patients_v2/",
                            annotations_filename),
              sep = "\t", row.names = F, col.names = F)
  print(paste("Processed",sample))
}

# ---- run infercnv on hpc ----

