
setwd('/media/data/lab/CLL1_sc_final_analyses_2024/')

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)
library(monocle3)
library(SingleR)
library(celldex)
library(symphony)
library(ggpubr)
library(patchwork)
library(BoneMarrowMap)
library(jcolors)

set.seed(1998)

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

dx_d15_blasts <- readRDS("./seurat_objects_merge_dx_d15/bcp_all_merge_dx_d15_seurat_filtered_with_ab_blasts.rds")

##################### Split couples and re-compute DimRed ######################

# PT08
dx_d15_blasts_pt08 <- subset(dx_d15_blasts, 
                             subset = Sample_Name == c("PT08_DX","PT08_D15"))
if(identical(as.character(unique(dx_d15_blasts_pt08$Sample_Name)),
             c("PT08_DX","PT08_D15"))) {
  print("OK")
}
DefaultAssay(dx_d15_blasts_pt08) <- "RNA"

dx_d15_blasts_pt08 <- FindVariableFeatures(dx_d15_blasts_pt08,
                                           selection.method = "vst",
                                           nfeatures = 2000)

dx_d15_blasts_pt08 <- ScaleData(dx_d15_blasts_pt08, 
                                features = rownames(dx_d15_blasts_pt08))

dx_d15_blasts_pt08 <- RunPCA(dx_d15_blasts_pt08,npcs = 60)

ElbowPlot(dx_d15_blasts_pt08, ndims = 50)
DimHeatmap(dx_d15_blasts_pt08, dims = 1:20, cells = 500, balanced = TRUE)

DefaultAssay(dx_d15_blasts_pt08) = "RNA"
dx_d15_blasts_pt08 <- RunUMAP(dx_d15_blasts_pt08, dims = 1:10, n.neighbors = 20)
dx_d15_blasts_pt08 <- FindNeighbors(dx_d15_blasts_pt08, dims = 1:10)
dx_d15_blasts_pt08 <- FindClusters(dx_d15_blasts_pt08, resolution = 0.5)

DimPlot(dx_d15_blasts_pt08, group.by = "Sample_Name")
DimPlot(dx_d15_blasts_pt08)

############################ SingleR annotation ################################

# Loading reference data with Ensembl annotations
ref.data <- HumanPrimaryCellAtlasData(ensembl=F)

# Performing predictions
predictions_PT08 <- SingleR(test=dx_d15_blasts_pt08@assays$RNA$counts, assay.type.test=1, 
                            ref=ref.data, labels=ref.data$label.main, de.method = 'classic')

table(predictions_PT08$labels)

dx_d15_blasts_pt08$ann_celltype <- predictions_PT08$labels

dx_d15_blasts_pt08_ann_joined <- JoinLayers(dx_d15_blasts_pt08)
dx_d15_blasts_pt08_ann.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt08_ann_joined)

# Define the updated cell types
cell_types <- c(
  "DC", "Epithelial_cells", "B_cell", "Neutrophils", 
  "Monocyte", "Erythroblast", "BM & Prog.", "Endothelial_cells", 
  "Gametocytes", "HSC_-G-CSF", "Macrophage", "NK_cell", 
  "Embryonic_stem_cells", "Tissue_stem_cells", "T_cells", "Osteoblasts", 
  "BM", "iPS_cells", "MSC", "HSC_CD34+", "CMP", "GMP", 
  "MEP", "Myelocyte", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", 
  "Pro-Myelocyte"
)

# Generate Monocle3-style palette for the updated list
n_colors <- length(cell_types)
monocle3_colors <- scales::hue_pal()(n_colors)

# Create a named vector for use in ggplot2
cell_colors <- setNames(monocle3_colors, cell_types)

# Print the resulting color map
print(cell_colors)

DimPlot(dx_d15_blasts_pt08, reduction = "umap", group.by = "ann_celltype", label = F, pt.size = 1.5, label.size = 5) +
  theme(plot.title=element_text(face = "bold", family = "Karla", size=20),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -4.5, angle=90, face = "bold", family = "Karla"))+ 
  ggtitle('Cell type annotation - immune cells') +
  geom_segment(aes(x=-11, y=-9.5, xend=-8.5, yend=-9.5), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-11, y=-9.5, xend=-11, yend=-7), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")

labels <- c("DC" = "DC", 
            "Epithelial_cells" = "Epithelial cells", 
            "B_cell"= "B cells", 
            "Neutrophils" = "Neutrophils", 
            "T_cells" = "T cells", 
            "Monocyte" = "Monocytes", 
            "Erythroblast" = "Erythroblasts", 
            "BM & Prog." = "BM & Prog.", 
            "Endothelial_cells" = "Endothelial cells", 
            "Gametocytes" = "Gametocytes", 
            "HSC_-G-CSF" = "HSC G-CSF", 
            "Macrophage" = "Macrophages", 
            "NK_cell" = "NK cells", 
            "Embryonic_stem_cells" = "Embryonic stem cells", 
            "Tissue_stem_cells" = "Tissue stem cells", 
            "Osteoblasts" = "Osteoblasts", 
            "BM" = "BM", 
            "iPS_cells" = "iPS cells", 
            "MSC" = "MSCs", 
            "HSC_CD34+" = bquote('HSCs CD'*34^'+'), 
            "CMP" = "CMPs", 
            "GMP" = "GMPs", 
            "MEP" = "MEPs", 
            "Myelocyte" = "Myelocytes", 
            "Pre-B_cell_CD34-" = bquote('Pre-B cells CD'*34^'-'), 
            "Pro-B_cell_CD34+" = bquote('Pro-B cells CD'*34^'+'), 
            "Pro-Myelocyte" = "Pro-Myelocytes")

png(file="./merge_plots_dx_d15_blasts_couples_annotation/umap_pt08_annotation.png", height = 2600, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt08_ann.cds, color_cells_by = "ann_celltype",
                     label_cell_groups = F,  cell_size = 2) + 
  scale_color_manual(values = cell_colors, labels = labels) +
  nice_theme2 +
  geom_segment(aes(x=-6.5, y=-7, xend=-5.25, yend=-7), arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-6.5, y=-7, xend=-6.5, yend=-5.3), arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2") +
  guides(color = guide_legend(ncol=4, override.aes = list(size = 5)))
dev.off()

################################## Projection ##################################

# Set directory to store projection reference files
projection_path = './zeng_reference/'

# Load Symphony reference
BM_ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Set uwot path for UMAP projection
BM_ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

BM_ref_obj <- create_ReferenceObject(BM_ref)

DimPlot(BM_ref_obj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', raster=FALSE, label=TRUE, label.size = 4) 

DimPlot(BM_ref_obj, reduction = 'umap', group.by = 'CellType_Broad', raster=FALSE, label=TRUE, label.size = 4) +
  scale_color_manual(values = proj_colors)

# our data
dx_d15_blasts_pt08 <- JoinLayers(dx_d15_blasts_pt08)
abseq <- GetAssay(dx_d15_blasts_pt08, "AB")

# MAPPING #

batchvar <- c('Sample_Name')

# Map query dataset using Symphony (Kang et al 2021)
dx_d15_blasts_pt08 <- map_Query(
  exp_query = dx_d15_blasts_pt08[["RNA"]]$counts, 
  metadata_query = dx_d15_blasts_pt08@meta.data,
  ref_obj = BM_ref,
  vars = batchvar
)

donor_key <- "Sample_Name" 
dx_d15_blasts_pt08 <- dx_d15_blasts_pt08 %>% calculate_MappingError(., reference = BM_ref, MAD_threshold = 3, 
                                                                    threshold_by_donor = TRUE, donor_key = donor_key)


# Get QC Plots
QC_plots <- plot_MappingErrorQC(dx_d15_blasts_pt08)

# CELL TYPE ASSIGNMENT #

# Predict Hematopoietic Cell Types by KNN classification 
dx_d15_blasts_pt08 <- predict_CellTypes(
  query_obj = dx_d15_blasts_pt08, 
  ref_obj = BM_ref, 
  initial_label = 'initial_CellType_BoneMarrowMap',
  final_label = 'predicted_CellType_BoneMarrowMap' 
) 

DimPlot(subset(dx_d15_blasts_pt08, mapping_error_QC == 'Pass'), 
        reduction = 'umap_projected',
        group.by = c('predicted_CellType_BoneMarrowMap'), 
        raster=FALSE, label=T, label.size = 4, split.by = "Sample_Name",
        ncol = 4) + 
  theme(legend.position = "none") +
  xlim(-15,10) + ylim(-15,15)

DimPlot(subset(dx_d15_blasts_pt08, mapping_error_QC == 'Pass'),
        reduction = 'umap_projected',
        group.by = c('predicted_CellType_BoneMarrowMap_Broad'), 
        raster=FALSE, label=F, label.size = 4, split.by = "Sample_Name", 
        ncol = 4) + 
  theme(legend.position = "none") +
  xlim(-15,10) + ylim(-15,15)

################################# PLOTS ########################################

#define colors
celltypes <- unique(BM_ref_obj@meta.data$CellType_Broad)
colors <- rev(scales::hue_pal(h = c(20,330))(24))
proj_colors <- setNames(colors, celltypes)

nice_theme2 <- theme(plot.title=element_text(face = "bold", 
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, 
                                                 vjust = 5.5, face = "bold", 
                                                 family = "Karla", size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5, 
                                                 angle=90, face = "bold", 
                                                 family = "Karla", size = 16),
                     legend.text = element_text(family = "Karla",size = 16),
                     legend.position = "bottom",
                     # legend.justification.right = "bottom",
                     legend.title = element_blank())

dx_d15_blasts_pt08_joined <- JoinLayers(subset(dx_d15_blasts_pt08, 
                                               mapping_error_QC == 'Pass'))
dx_d15_blasts_pt08.cds <- SeuratWrappers::as.cell_data_set(dx_d15_blasts_pt08_joined)

png(file="./merge_plots_dx_d15_blasts_couples_projection_labels/umap_pt08_projection_labels.png", 
    height = 2600, width = 2500, units = "px", res = 300)
monocle3::plot_cells(dx_d15_blasts_pt08.cds, color_cells_by = "predicted_CellType_BoneMarrowMap_Broad",
                     label_cell_groups = F,  cell_size = 2) + 
  scale_color_manual(values = proj_colors) +
  nice_theme2 +
  geom_segment(aes(x=-6.5, y=-7, xend=-5.25, yend=-7),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #x
  geom_segment(aes(x=-6.5, y=-7, xend=-6.5, yend=-5.3),
               arrow = arrow(length=unit(.4, 'cm')), color = "black") + #y
  labs(x = "UMAP1", y = "UMAP2") +
  guides(color = guide_legend(ncol=4, override.aes = list(size = 5)))
dev.off()
