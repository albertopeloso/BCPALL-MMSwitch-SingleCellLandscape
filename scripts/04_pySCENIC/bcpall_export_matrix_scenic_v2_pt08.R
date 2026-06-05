
setwd("/media/data/lab/CLL1_sc_final_analyses_2024/")

library(data.table)
library(dplyr)
library(Seurat)
library(convert2anndata)
library(anndata)
library(S4Vectors)

sobj_pt08 <- readRDS('./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_seurat_filtered_with_ab_blasts_proj_ann_pt08.rds')
sobj_pt08

sobj_pt08@meta.data$cells = rownames(sobj_pt08@meta.data)
sobj_pt08@commands <- list()
sce <- convert_seurat_to_sce(sobj_pt08)

# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "RNA")

# Save the AnnData object
write_h5ad(ad, "./palantir_on_couples_202501/h5ad_objects/pt08.h5ad")
