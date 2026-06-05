
setwd("/media/data/lab/bcpall_sc_final_analyses_2024/")

library(tidyverse)
library(Seurat)
library(data.table)
library(findPC)
library(harmony)
library(tibble)
library(viridis) 

set.seed(1998)

# ---- load object and metadata ----

dx_combined <- readRDS("./seurat_objects_merge_dx/bcp_all_merge_dx_seurat_filtered_with_ab_dimred.rds")

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
dx_combined <- JoinLayers(dx_combined)

# ---- marker extraction from non-malignant clusters ----

# find all markers of cluster 3 (normal T cells)
cluster3_T.markers <- FindMarkers(dx_combined, 
                                  ident.1 = 3,
                                  min.pct = 0.5)
cluster3_T.markers <- cluster3_T.markers %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column(var = "Gene Name") %>%
  mutate("Non-malignant population" = "T-cells")
write.table(cluster3_T.markers, 
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_T_markers.txt", 
            sep = "\t",
            row.names = FALSE)

# find all markers of cluster 14 (normal B cells)
cluster14_B.markers <- FindMarkers(dx_combined, 
                                   ident.1 = 14,
                                   min.pct = 0.5)
cluster14_B.markers <- cluster14_B.markers %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column(var = "Gene Name") %>%
  mutate("Non-malignant population" = "B-cells")
write.table(cluster14_B.markers, 
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_B_markers.txt", 
            sep = "\t",
            row.names = FALSE)

# find all markers of cluster 16 (normal myelo/ery)
cluster16_myelo_ery.markers <- FindMarkers(dx_combined, 
                                           ident.1 = 16,
                                           min.pct = 0.5)
cluster16_myelo_ery.markers <- cluster16_myelo_ery.markers %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column(var = "Gene Name") %>%
  mutate("Non-malignant population" = "myelo-ery")
write.table(cluster16_myelo_ery.markers,
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/single_cell/normal_myelo_ery_markers.txt",
            sep = "\t",
            row.names = FALSE)

# ---- plotting ----

# dotplot for marker genes
cluster3_allmarkers <- FindMarkers(dx_combined, 
                                  ident.1 = 3,
                                  logfc.threshold = 0.01)
cluster14_allmarkers <- FindMarkers(dx_combined, 
                                   ident.1 = 14,
                                   logfc.threshold = 0.01)
cluster16_allmarkers <- FindMarkers(dx_combined, 
                                    ident.1 = 16,
                                    logfc.threshold = 0.01)

non_malignant_genes <- c(cluster3_T.markers$`Gene Name`[1:10],
                         cluster16_myelo_ery.markers$`Gene Name`[1:10],
                         cluster14_B.markers$`Gene Name`[1:10])

# non_malignant <- rbind(
#   cbind(cluster3_allmarkers[non_malignant_genes, ],
#         data.frame(celltype = "T-cells")),
#   cbind(cluster16_allmarkers[non_malignant_genes, ],
#         data.frame(celltype = "myelo-ery")),
#   cbind(cluster14_allmarkers[non_malignant_genes, ],
#         data.frame(celltype = "B-cells"))
# ) %>% 
#   rownames_to_column(var = "gene") %>%
#   mutate(gene = gsub("[0-9]+$", "", gene))
# 

non_malignant <- bind_rows(
  cluster3_allmarkers %>% rownames_to_column("gene") %>% mutate(celltype = "T-cells"),
  cluster16_allmarkers %>% rownames_to_column("gene") %>% mutate(celltype = "myelo-ery"),
  cluster14_allmarkers %>% rownames_to_column("gene") %>% mutate(celltype = "B-cells")
) %>%
  filter(gene %in% non_malignant_genes)

non_malignant$celltype <- factor(non_malignant$celltype,
                                 levels = c("B-cells","myelo-ery", "T-cells"))
non_malignant$gene <- factor(non_malignant$gene,
                             levels = non_malignant_genes)
non_malignant$pct.1_100 <- non_malignant$pct.1 * 100

png(file="./merge_plots_dx/dotplot_non_malignant_markers.png", height = 4000, width = 4000, units = "px", res = 300)
ggplot(non_malignant, aes(x = celltype,
                          y = gene)) +
  geom_point(aes(size = pct.1_100, color = avg_log2FC)) +
  scale_size(range = c(1, 12), name = "Percentage \nexpressing") +
  theme_minimal(base_family = "karla") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Log2(FC)",
                     breaks = c(-4, -2, 0, 2, 4, 6))+
  scale_x_discrete(position = "top") +
  theme(text = element_text(size=40),
        axis.text.y = element_text(face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.key.size = unit(1.2, "cm")) +
  labs(
    x = "Cell type",
    y = "")
dev.off()

# check precomputed umap
DimPlot(dx_combined, group.by = "Sample_Name")

# set nice themes for plotting
nice_theme2 <- theme(plot.title=element_text(face = "bold",
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_text(hjust = 0.0475, vjust = 6.5,
                                                 face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     axis.title.y = element_text(hjust = 0.05, vjust = -5,
                                                 angle=90, face = "bold",
                                                 family = "Karla",
                                                 size = 16),
                     legend.text = element_text(family = "Karla",size = 14),
                     legend.position = "right",
                     legend.title = element_text(face = "bold", 
                                                 family = "Karla",
                                                 size = 16))

# B-cell score
B_cell <- cluster14_B.markers$`Gene Name`[1:10]
B_cell_markers <- list(B_cell)
dx_combined <- AddModuleScore(object = dx_combined, 
                              features = B_cell_markers,
                              name = "B_cell_score")

# myelo/ery score
myelo_ery <- cluster16_myelo_ery.markers$`Gene Name`[1:10]
myelo_ery_markers <- list(myelo_ery)
dx_combined <- AddModuleScore(object = dx_combined, 
                              features = myelo_ery_markers, 
                              name = "myelo_ery_score")

# T-cell score
T_cell <- cluster3_T.markers$`Gene Name`[1:10]
T_cell_markers <- list(T_cell)
dx_combined <- AddModuleScore(object = dx_combined, 
                              features = T_cell_markers, 
                              name = "T_cell_score")

dx_combined.cds <- SeuratWrappers::as.cell_data_set(dx_combined)

png(file="./merge_plots_dx/umap_b_cell_score_monocle.png", height = 2500, width = 2700, units = "px", res = 300)
monocle3::plot_cells(dx_combined.cds, color_cells_by = "B_cell_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "B-cell \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5, angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-12, xend=-9.5, yend=-12), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.25), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

png(file="./merge_plots_dx/umap_myelo_ery_score_monocle.png", height = 2500, width = 2800, units = "px", res = 300)
monocle3::plot_cells(dx_combined.cds, color_cells_by = "myelo_ery_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "Myelo/ery \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5, angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-12, xend=-9.5, yend=-12), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.25), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()

png(file="./merge_plots_dx/umap_t_cell_score_monocle.png", height = 2500, width = 2700, units = "px", res = 300)
monocle3::plot_cells(dx_combined.cds, color_cells_by = "T_cell_score1",
                     label_cell_groups = F, cell_size = 0.5) + 
  scale_color_gradientn(colours = viridis::magma(100), name = "T-cell \nscore",
                        breaks = c(0, 1, 2)) +
  theme(plot.title=element_text(face = "bold.italic", family = "Karla", size=20)) +
  theme(plot.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0.05, vjust = 4.5, face = "bold", family = "Karla"),
        axis.title.y = element_text(hjust = 0.05, vjust = -3.5, angle=90, face = "bold", family = "Karla"),
        legend.title = element_text(face = "bold", family = "Karla"),
        legend.text = element_text(face = "bold", family = "Karla")) + 
  nice_theme2 +
  geom_segment(aes(x=-12, y=-12, xend=-9.5, yend=-12), arrow = arrow(length=unit(.4, 'cm'))) + 
  geom_segment(aes(x=-12, y=-12, xend=-12, yend=-9.25), arrow = arrow(length=unit(.4, 'cm'))) +
  labs(x = "UMAP1", y = "UMAP2")
dev.off()



