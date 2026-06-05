
setwd("/media/data/lab/bcpall_sc_final_analyses_2024/blasts_projection_zeng/")

library(Seurat)
library(tidyverse)
library(symphony)
library(ggpubr)
library(patchwork)
library(BoneMarrowMap)

### The script for the projection is  adapted from
### Zeng AGX, Iacobucci I, Shah S, et al. Single-cell Transcriptional Atlas of 
### Human Hematopoiesis Reveals Genetic and Hierarchy-Based Determinants of 
### Aberrant AML Differentiation. 
### Blood Cancer Discov. 2025;6(4):307-324. doi:10.1158/2643-3230.BCD-24-0342
### in the BoneMarrowMap R package 

#################### preprocess and setting references #########################

# Set directory to store projection reference files
projection_path = '/media/data/lab/bcpall_sc_final_analyses_2024/zeng_reference/'

# Download Bone Marrow Reference 
# curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_SymphonyRef.rds', 
# destfile = paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# Download uwot model file
# curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
# destfile = paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot'))

# Load Symphony reference
BM_ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))

# Set uwot path for UMAP projection
BM_ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

BM_ref_obj <- create_ReferenceObject(BM_ref)
# check  reference UMAP
DimPlot(BM_ref_obj, reduction = 'umap', group.by = 'CellType_Annotation_formatted',
        raster=FALSE, label=TRUE, label.size = 4)

# reference UMAP with cell cycle phase and lineage pseudotime estimates
p1 <- DimPlot(BM_ref_obj, reduction = 'umap',
              group.by = 'CyclePhase', raster=FALSE)
p2 <- FeaturePlot(BM_ref_obj, reduction = 'umap', 
                  features = 'Pseudotime', raster=FALSE) 

# load our scRNAseq dataset
merge_dx <- readRDS('../seurat_objects_merge_dx/bcpall_merge_dx_seurat_filtered_with_ab_blasts.rds')
merge_dx
merge_dx <- JoinLayers(merge_dx)
abseq <- GetAssay(merge_dx, "AB")

################################# mapping ######################################

# save original and integrated embeddings
original_reductions <- list(
  pca = merge_dx[["pca"]],
  umap = merge_dx[["umap"]]
)

# batch variable to correct in the merge_dx data
batchvar <- c('Sample_Name', 'orig.ident')

# Map merge_dx dataset using Symphony (Kang et al 2021)
merge_dx <- map_Query(
  exp_query = merge_dx[["RNA"]]$counts, 
  metadata_query = merge_dx@meta.data,
  ref_obj = BM_ref,
  vars = batchvar
)

# Run QC based on mapping error score
donor_key <- "Sample_Name"
merge_dx <- merge_dx %>% calculate_MappingError(., reference = BM_ref,
                                                MAD_threshold = 3, 
                                                threshold_by_donor = TRUE,
                                                donor_key = donor_key) 

# Plot distribution by patient to ensure you are catching the tail
png(file="../projection_plots_zeng_ref_dx/reference_complete_mapping_error.png",
    height = 800, width = 1600, units = "px")
merge_dx@meta.data %>% 
  ggplot(aes(x = mapping_error_score, fill = mapping_error_QC)) + 
  geom_histogram(bins = 200) + facet_wrap(.~get(donor_key))
dev.off()

# Get QC Plots
QC_plots <- plot_MappingErrorQC(merge_dx)

# Plot together
png(file="../projection_plots_zeng_ref_dx/reference_complete_qc_plots.png",
    height = 800, width = 1600, units = "px")
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

############################ cell type assignment ##############################

# Predict Hematopoietic Cell Types by KNN classification 
merge_dx <- predict_CellTypes(
  query_obj = merge_dx, 
  ref_obj = BM_ref, 
  initial_label = 'initial_CellType_BoneMarrowMap',
  final_label = 'predicted_CellType_BoneMarrowMap'  #
) 

DimPlot(subset(merge_dx, mapping_error_QC == 'Pass'), reduction = 'umap_projected',
        group.by = c('predicted_CellType_BoneMarrowMap'), 
        raster=FALSE, label=T,
        label.size = 4, split.by = "Sample_Name", ncol = 4) + 
  theme(legend.position = "none") +
  xlim(-15,10) + ylim(-15,15)

DimPlot(subset(merge_dx, mapping_error_QC == 'Pass'), reduction = 'umap_projected',
        group.by = c('predicted_CellType_BoneMarrowMap_Broad'),
        raster=FALSE, label=F,
        label.size = 4, split.by = "Sample_Name", ncol = 4) + 
  theme(legend.position = "none") +
  xlim(-15,10) + ylim(-15,15)

################################# plots ########################################

# add back original embeddings
merge_dx[["pca"]] <- original_reductions$pca
Key(object = merge_dx[["umap_projected"]]) <- "umapprojected_"
merge_dx[["umap"]] <- original_reductions$umap

sample_names <- as.vector(unique(merge_dx@meta.data$Sample_Name))
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected"
for(sample_name in sample_names) {
  print(sample_name)
  file_path <- file.path(output_dir, paste0("projection_", sample_name, ".png"))
  png(file=file_path, height = 800, width = 800, units = "px")
  print(DimPlot(subset(merge_dx, mapping_error_QC == 'Pass' & 
                         Sample_Name == sample_name), reduction = 'umap_projected',
                group.by = c('predicted_CellType_BoneMarrowMap'), 
                raster=FALSE, label=T, label.size = 4, split.by = "Sample_Name",
                ncol = 4) + theme(legend.position = "none") +
          xlim(-15,10) + ylim(-15,15)) 
  dev.off()
}

############################ projection density ################################

sample_names <- as.vector(unique(merge_dx@meta.data$Sample_Name))
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected_density"
batch_key <- "Sample_Name"
for(sample_name in sample_names) {
  print(sample_name)
  file_path <- file.path(output_dir, paste0("projection_",
                                            sample_name,
                                            "_complete_density.png"))
  png(file=file_path, height = 600, width = 600, units = "px")
  print(plot_Projection_BM(
    query_obj = merge_dx, 
    sample_name = sample_name, 
    batch_key = batch_key, 
    ref_obj = BM_ref, 
    Hierarchy_only = FALSE,
    downsample_reference = TRUE, 
    downsample_frac = 0.25,  
    query_point_size = 0.2,
  ))
  dev.off()
}

############################ merge_dx composition ############################## 

merge_dx_composition <- get_Composition(
  merge_dx_obj = merge_dx, 
  donor_key = 'Sample_Name', 
  celltype_label = 'predicted_CellType_BoneMarrowMap', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = 0.5, 
  return_type = 'long')

merge_dx_composition 

merge_dx_composition <- get_Composition(
  merge_dx_obj = merge_dx, 
  donor_key = 'Sample_Name', 
  celltype_label = 'predicted_CellType_BoneMarrowMap', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = 0.5, 
  return_type = 'proportion')

merge_dx_composition 

# Simple heatmap to visualize composition of projected samples
p <- merge_dx_composition %>% 
  # show celltypes present in >1% of total cells
  select(Sample_Name, colnames(merge_dx_composition)[-1][
    colSums(merge_dx_composition[-1]) > 0.01]) %>% 
  # convert to matrix and display heatmap
  column_to_rownames('Sample_Name') %>% data.matrix() %>%
  ComplexHeatmap::Heatmap()
p


############################# B-cell development ###############################

# Set directory to store projection reference files
projection_path = '/media/data/lab/bcpall_sc_final_analyses_2024/zeng_reference/'

# Download Bone Marrow Reference - 187 Mb
# curl::curl_download('https://bdevelopmentmap.s3.us-east-2.amazonaws.com/BDevelopment_RefMap_SymphonyRef.rds',
# destfile = paste0(projection_path, 'BDevelopment_RefMap_SymphonyRef.rds'))
# Download uwot model file - 99 Mb
# curl::curl_download('https://bdevelopmentmap.s3.us-east-2.amazonaws.com/BDevelopment_RefMap_uwot_model.uwot',
# destfile = paste0(projection_path, 'BDevelopment_RefMap_uwot_model.uwot'))

# Load Symphony reference
bdev_ref <- readRDS(paste0(projection_path,
                           'BDevelopment_RefMap_SymphonyRef.rds'))
bdev_ref$umap$embedding[,1] = -bdev_ref$umap$embedding[,1] 

# Set uwot path for UMAP projection
bdev_ref$save_uwot_path <- paste0(projection_path,
                                  'BDevelopment_RefMap_uwot_model.uwot')

bdev_ref_obj <- create_ReferenceObject(bdev_ref)

# reference UMAP
DimPlot(bdev_ref_obj, reduction = 'umap', group.by = 'BDevelopment_CellType',
        raster=FALSE, label=TRUE, label.size = 4)

# filter for celltypes along B cell development
B_development_celltypes <- c('HSC', 'HSC/MPP', 'MPP-MyLy', 'MPP-LMPP', 'LMPP',
                             'Early GMP', 'MLP', 'MLP-II', 'Pre-pDC',
                             'Pre-pDC Cycling', 'pDC','CLP', 'EarlyProB',
                             'Pre-ProB', 'Pro-B VDJ', 'Pro-B Cycling',
                             'Large Pre-B', 'Small Pre-B',
                             'Immature B', 'Mature B')

# subset only B cell development cells
merge_dx <- subset(merge_dx, predicted_CellType_BoneMarrowMap %in%
                     B_development_celltypes)
merge_dx

# regenerate a new abseq assay, with subset cells
subset_abseq <- merge_dx[["AB"]] 

################################ mapping #######################################

# batch variable to correct in the merge_dx data
batchvar <- c('Sample_Name', 'orig.ident')

# Map merge_dx dataset using Symphony (Kang et al 2021)
merge_dx <- map_Query(
  exp_query = merge_dx[["RNA"]]$counts, 
  metadata_query = merge_dx@meta.data,
  ref_obj = bdev_ref,
  vars = batchvar
)

# Flip UMAP1 to go from left to right
merge_dx[['umap_projected']]@cell.embeddings[,1] = -merge_dx[['umap_projected']]@cell.embeddings[,1]
merge_dx[['umap_projected']]
bdev_ref$meta_data

# Predict Hematopoietic Cell Types by KNN classification
merge_dx <- predict_CellTypes(
  query_obj = merge_dx, 
  ref_obj = bdev_ref, 
  ref_label = 'BDevelopment_CellType',   
  initial_label = 'initial_CellType_BDevelopment',
  final_label = 'predicted_CellType_BDevelopment'
) 

# Set batch/condition to be visualized individually
batch_key <- 'Sample_Name'
# batch_key <- 'orig.ident'

# returns a list of plots for each donor from a pre-specified batch variable
projection_plots_bdev <- plot_Projection_byDonor(
  query_obj = merge_dx, 
  batch_key = batch_key, 
  ref_obj = bdev_ref, 
  downsample_reference = TRUE, 
  downsample_frac = 0.25,  
  query_point_size = 0.2,
  saveplot = F, 
  save_folder = 'projectionFigures/'
)

# show plots together with patchwork
patchwork::wrap_plots(projection_plots_bdev, ncol = 3)

projection_plots_bdev <- projection_plots_bdev[c("PT01", "PT02", "PT03", "PT04",
                                                 "PT05", "PT06_DX", "PT07_DX",
                                                 "PT08_DX", "PT13_DX", "PT14_DX",
                                                 "PT15_DX")]

# show plots together with patchwork
png(file="../projection_plots_zeng_ref_dx/umap_projected_density/projection_bdev_density.png",
    height = 1200, width = 1400, units = "px")
patchwork::wrap_plots(projection_plots_bdev, ncol = 3)
dev.off()

# plot density
sample_names <- as.vector(unique(merge_dx@meta.data$Sample_Name))
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected_density_bdev"
for(sample_name in sample_names) {
  print(sample_name)
  file_path <- file.path(output_dir, paste0("projection_", sample_name,
                                            "_bdev_density.png"))
  png(file=file_path, height = 500, width = 500, units = "px")
  print(plot_Projection_byDonor(
    query_obj = subset(merge_dx, Sample_Name == sample_name), 
    batch_key = batch_key, 
    ref_obj = bdev_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25,
    query_point_size = 0.2,
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ))
  dev.off()
  
}

# plots umap projection b development
sample_names <- as.vector(unique(merge_dx@meta.data$Sample_Name))
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected_density_bdev"
for(sample_name in sample_names) {
  print(sample_name)
  file_path <- file.path(output_dir, paste0("projection_", sample_name,
                                            "_bdev.png"))
  png(file=file_path, height = 800, width = 800, units = "px")
  print(DimPlot(subset(merge_dx, mapping_error_QC == 'Pass' & 
                         Sample_Name == sample_name),
                reduction = 'umap_projected',
                group.by = c('predicted_CellType_BDevelopment'), 
                raster=FALSE, label=TRUE, label.size = 4,
                split.by = 'Sample_Name') +
          xlim(-10,12) + ylim(-10,8))
  dev.off()
}

# re-add Abseq
merge_dx[["AB"]] <- subset_abseq
DefaultAssay(merge_dx) <- "AB"
merge_dx <- NormalizeData(merge_dx, normalization.method = "CLR", margin = 2)
Idents(merge_dx) <- "predicted_CellType_BDevelopment"


# plots for abseq expression in projection umap
abseq_names <- as.vector(merge_dx@assays$AB@counts@Dimnames[[1]])
for(ab in abseq_names) {
  print(ab)
  file_path <- file.path(output_dir, paste0("projection_bdev", ab, ".png"))
  png(file=file_path, height = 600, width = 800, units = "px")
  print(FeaturePlot(merge_dx, features = ab) + 
          ggtitle(ab) +
          scale_colour_gradientn(colours = brewer.pal(n = 9, name = "OrRd")))
  dev.off()
}

############# splitted by switching status and change color palette ############

# modified function from BoneMarrowMap
plot_Projection_BM_mod <- function (query_obj, batch_key, sample_name, ref_obj,
                                    Hierarchy_only = FALSE, 
                                    downsample_reference = TRUE,
                                    downsample_frac = 0.25,
                                    query_point_size = 0.2, 
                                    query_contour_size = 0.3,
                                    saveplot = TRUE,
                                    device = "pdf", 
                                    save_folder = "projectionFigures/") 
{
  if (!batch_key %in% colnames(query_obj@meta.data)) {
    stop("Label \"{batch_key}\" is not available in the query metadata.")
  }
  if (!sample_name %in% unique(query_obj@meta.data[[batch_key]])) {
    stop("Label \"{sample_name}\" not found in the \"{batch_key}\" variable within the query metadata")
  }
  dat <- query_obj@meta.data %>% tibble::rownames_to_column("Cell") %>% 
    dplyr::filter(.data[[batch_key]] == sample_name) %>% 
    dplyr::left_join(query_obj@reductions$umap@cell.embeddings %>% 
                       data.frame() %>% tibble::rownames_to_column("Cell"), 
                     by = "Cell") %>% dplyr::filter(mapping_error_QC == "Pass") %>% 
    dplyr::mutate(ref_query = "query")
  if (Hierarchy_only) {
    dat <- dat %>% dplyr::filter(umap_1 > -5)
  }
  background <- data.frame(ref_obj$umap$embedding) %>% dplyr::rename(umap_1 = X1, 
                                                                     umap_2 = X2)
  if (Hierarchy_only) {
    background <- background %>% dplyr::filter(umap_1 > -5)
  }
  if (downsample_reference) {
    set.seed(123)
    background <- background %>% dplyr::sample_frac(downsample_frac)
  }
  heatpalette <- grDevices::heat.colors(12)
  p <- dat %>% ggplot2::ggplot(aes(x = umap_1, y = umap_2)) + 
    ggplot2::geom_point(data = background, color = "grey78", 
                        size = 0.05, alpha = 0.5) + 
    ggpointdensity::geom_pointdensity(size = query_point_size) + 
    jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, 
                                        bias = 1.75) + 
    ggplot2::geom_density_2d(alpha = 0.4, color = "black", h = 1.5,
                             size = query_contour_size) + 
    ggplot2::theme_void() + ggplot2::ggtitle(sample_name) + 
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                   legend.position = "none")
  if (saveplot) {
    if (!file.exists(save_folder)) {
      dir.create(file.path(paste0("./", save_folder)))
    }
    if (Hierarchy_only) {
      if (device %in% c("pdf", "PDF", "Pdf")) {
        ggsave(paste0(save_folder, "density_", sample_name, 
                      "_projectedUMAP_HierarchyOnly.pdf"), height = 4, 
               width = 4.5, device = "pdf")
      }
      else if (device %in% c("png", "PNG", "Png")) {
        ggsave(paste0(save_folder, "density_", sample_name, 
                      "_projectedUMAP_HierarchyOnly.png"), height = 4, 
               width = 4.5, device = "png", dpi = 240)
      }
      else {
        stop("Please specify device as either a pdf or png.")
      }
    }
    else {
      if (device %in% c("pdf", "PDF", "Pdf")) {
        ggsave(paste0(save_folder, "density_", sample_name, 
                      "_projectedUMAP.pdf"), height = 4, width = 6, 
               device = "pdf")
      }
      else if (device %in% c("png", "PNG", "Png")) {
        ggsave(paste0(save_folder, "density_", sample_name, 
                      "_projectedUMAP.pdf"), height = 4, width = 6, 
               device = "png", dpi = 240)
      }
      else {
        stop("Please specify device as either a pdf or png.")
      }
    }
  }
  return(p)
}

###################### all hematopoietic subpopulations ########################

library(paletteer)

switching <- as.vector(unique(merge_dx@meta.data$Switch))
switch_key <- "Switch"
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected_density_allemato"
for(switching_state in switching) {
  print(switching_state)
  file_path <- file.path(output_dir, paste0("projection_", switching_state,
                                            "_allemato_density.png"))
  png(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = switching_state,
    batch_key = switch_key, 
    ref_obj = BM_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25, 
    query_point_size = 0.2,
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

for(switching_state in switching) {
  print(switching_state)
  file_path <- file.path(output_dir, paste0("projection_", switching_state,
                                            "_allemato_density.tiff"))
  tiff(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = switching_state,
    batch_key = switch_key, 
    ref_obj = BM_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25
    query_point_size = 0.2,
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

############################## split PATIENT ###################################

for(sample in sample_names) {
  print(sample)
  file_path <- file.path(output_dir, paste0("projection_split_", sample,
                                            "_allemato_density_newcol.png"))
  png(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = sample,
    batch_key = batch_key, 
    ref_obj = BM_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25, 
    query_point_size = 0.2,  
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

for(sample in sample_names) {
  print(sample)
  file_path <- file.path(output_dir, paste0("projection_split_", sample,
                                            "_allemato_density_newcol.tiff"))
  tiff(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = sample,
    batch_key = batch_key, 
    ref_obj = BM_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25,   
    query_point_size = 0.2,   
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

############################ bdev subpopulations ###############################

library(paletteer)

switching <- as.vector(unique(merge_dx@meta.data$Switch))
switch_key <- "Switch"
output_dir <- "../projection_plots_zeng_ref_dx/umap_projected_density_bdev"
for(switching_state in switching) {
  print(switching_state)
  file_path <- file.path(output_dir, paste0("projection_", switching_state,
                                            "_bdev_density.png"))
  png(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = switching_state,
    batch_key = switch_key, 
    ref_obj = bdev_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25, 
    query_point_size = 0.2,  
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

for(switching_state in switching) {
  print(switching_state)
  file_path <- file.path(output_dir, paste0("projection_", switching_state,
                                            "_bdev_density.tiff"))
  tiff(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = switching_state,
    batch_key = switch_key, 
    ref_obj = bdev_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25,   
    query_point_size = 0.2,   
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

############################## split PATIENT ###################################

for(sample in sample_names) {
  print(sample)
  file_path <- file.path(output_dir, paste0("projection_split_", sample,
                                            "_bdev_density_newcol.png"))
  png(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = sample,
    batch_key = batch_key, 
    ref_obj = bdev_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25,  
    query_point_size = 0.2, 
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}

for(sample in sample_names) {
  print(sample)
  file_path <- file.path(output_dir, paste0("projection_split_", sample,
                                            "_bdev_density_newcol.tiff"))
  tiff(file=file_path, width = 2600, height = 2000, units = "px", res = 300)
  p <- plot_Projection_BM_mod(
    query_obj = merge_dx, 
    sample_name = sample,
    batch_key = batch_key, 
    ref_obj = bdev_ref, 
    downsample_reference = TRUE, 
    downsample_frac = 0.25,  
    query_point_size = 0.2, 
    saveplot = F, 
    save_folder = 'projectionFigures/'
  ) +
    scale_color_gradientn(colours = paletteer_c("grDevices::BluYl", 30)) +
    theme_void() +
    theme(legend.position = "none") 
  print(p)
  dev.off()
}
