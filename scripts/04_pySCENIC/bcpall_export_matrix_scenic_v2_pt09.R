
setwd("/media/data/lab/CLL1_sc_final_analyses_2024/")

library(data.table)
library(dplyr)
library(Seurat)
library(convert2anndata)
library(anndata)
library(S4Vectors)

sobj_pt13 <- readRDS('./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_seurat_filtered_with_ab_blasts_proj_ann_pt13.rds')
sobj_pt13

sobj_pt13@meta.data$cells = rownames(sobj_pt13@meta.data)
sobj_pt13@commands <- list()
sce <- convert_seurat_to_sce(sobj_pt13)

# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "RNA")

# Save the AnnData object
write_h5ad(ad, "./palantir_on_couples_202501/h5ad_objects/pt09.h5ad")
