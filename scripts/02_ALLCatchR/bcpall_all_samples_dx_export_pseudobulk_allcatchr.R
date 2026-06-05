
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)

set.seed(1998)

# ------------------------------ load files ---------------------------------- #
dx_combined <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_dimred.rds")

# ---- add information on patient and group ----
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

# --------------------- export pseudobulk expression ------------------------- #

sample_names <- as.vector(unique(dx_combined$Sample_Name))

pseudobulk_list <- list()

for (sample in sample_names) {
  print(sample)
  dx_subset <- subset(dx_combined, Sample_Name == sample)
  
  if (identical(as.character(unique(dx_subset$Sample_Name)), sample)) {
    print("OK")
  }
  
  pb_subset_dx <- AggregateExpression(dx_subset, assay = "RNA",
                                      group.by = "Sample_Name")
  pb_subset_dx_df <- as.data.frame(pb_subset_dx$RNA)
  
  # Rename the column to the sample name
  colnames(pb_subset_dx_df) <- sample
  pseudobulk_list[[sample]] <- pb_subset_dx_df
}

pseudobulk_combined <- Reduce(function(x, y) {
  merged <- merge(x, y, by = "row.names", all = TRUE)
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  return(merged)
}, pseudobulk_list)

write.table(pseudobulk_combined,
            file = "./pseudobulk_dx/pb_counts_allcells_combined.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)
