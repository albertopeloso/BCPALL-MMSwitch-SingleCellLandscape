
setwd('/media/data/lab/bcpall_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ---- read in filtered and merged seurat object ----
dx_combined <- readRDS(file = "./seurat_objects_merge_dx/bcpall_merge_dx_seurat_filtered_with_ab_step1.rds")

DefaultAssay(dx_combined) <- "RNA"

# log normalizaion of UMI counts
dx_combined <- NormalizeData(dx_combined, normalization.method = "LogNormalize")

# ---- PCA ----
dx_combined <- FindVariableFeatures(dx_combined,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# scale the data
dx_combined <- ScaleData(dx_combined, features = rownames(dx_combined))

# perform PCA
dx_combined <- RunPCA(dx_combined,npcs = 60)

# evaluating the effects of mitochondrial percentage
# check quartile values
summary(dx_combined$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2801 15.8304 19.7312 20.5184 24.9597 34.9854
# create a table based on cutoffs:
#  - Low --> below 1st quartile (9.5940)
#  - Medium --> between 1st and median (14.0696)
#  - Medium High --> between median and 3rd quantile (15.4431)
#  - High --> above 3rd quartile (19.9639)

dx_combined$quartile.mt <- 
  cut(dx_combined$percent.mt,
      breaks = c(-Inf, 9.5940, 15.4431, 19.9639, Inf),
      labels = c("Low", "Medium", "Medium High", "High"))

# check if it has been created
head(dx_combined@meta.data)

# check whether the mitochondrial percentage is a source of variation
DimPlot(dx_combined,
        reduction = "pca",
        group.by = "quartile.mt",
        split.by = "quartile.mt")

# ---- perform umap ----
DefaultAssay(dx_combined) = "RNA"
ElbowPlot(dx_combined, ndims = 50)
DimHeatmap(dx_combined, dims = 35:60, cells = 500, balanced = TRUE)
dx_combined <- RunUMAP(dx_combined, dims = 1:44)
dx_combined <- FindNeighbors(dx_combined, dims = 1:44)
dx_combined <- FindClusters(dx_combined, resolution = 2, )
DimPlot(dx_combined, group.by = "Sample_Name")
DimPlot(dx_combined, group.by = "orig.ident")
DimPlot(dx_combined, label = T)

################################## plotting ####################################

# ---- test umaps with different dimensions and resolutions ----

# grouped by cluster
dimensions <- c(25:50)
output_dir <- "/media/data/lab/bcpall_sc_final_analyses_2024/test_umap_dims/"
for(dimension in dimensions) {
  print(dimension)
  dx_combined <- RunUMAP(dx_combined, dims = 1:dimension)
  dx_combined <- FindNeighbors(dx_combined, dims = 1:dimension)
  dx_combined <- FindClusters(dx_combined, resolution = 0.5)
  file_path <- file.path(output_dir, paste0("umap_", dimension, "_dims", ".png"))
  png(file=file_path, height = 800, width = 800, units = "px")
  print(DimPlot(dx_combined)) 
  dev.off()
}

# grouped by sample name
dimensions <- c(25:50)
output_dir <- "/media/data/lab/bcpall_sc_final_analyses_2024/test_umap_dims_sn/"
for(dimension in dimensions) {
  print(dimension)
  dx_combined <- RunUMAP(dx_combined, dims = 1:dimension)
  dx_combined <- FindNeighbors(dx_combined, dims = 1:dimension)
  dx_combined <- FindClusters(dx_combined, resolution = 0.5)
  file_path <- file.path(output_dir, paste0("umap_", dimension, "_dims", ".png"))
  png(file=file_path, height = 800, width = 800, units = "px")
  print(DimPlot(dx_combined, group.by = "Sample_Name")) 
  dev.off()
}

# different res
resolutions = seq(0.5,2, by = 0.25)
output_dir <- "/media/data/lab/bcpall_sc_final_analyses_2024/test_umap_44_res/"
for(resolution in resolutions) {
  dx_combined <- RunUMAP(dx_combined, dims = 1:44)
  dx_combined <- FindNeighbors(dx_combined, dims = 1:44)
  dx_combined <- FindClusters(dx_combined, resolution = resolution)
  file_path <- file.path(output_dir, paste0("umap_44_dims_", "_res",
                                            resolution, ".png"))
  png(file=file_path, height = 800, width = 800, units = "px")
  print(DimPlot(dx_combined, label = T)) 
  dev.off()
}
