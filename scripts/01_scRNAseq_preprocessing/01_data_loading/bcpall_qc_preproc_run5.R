
setwd("/media/data/lab/CLL1_sc_final_analyses_2024/")

######## setup the seurat object using UMI-DBEC normalized from 7bridges #######

library(tidyverse)
library(Seurat)
library(data.table)
library(ggplot2)
library(ggpubr)

set.seed(1998)

# ---- read in the  matrix ----
expMat <- ReadMtx(mtx="./matrixes_run5/BCP-ALL-POOL1-POOL2-MAY24_DBEC_MolsPerCell_MEX/matrix.mtx",
                  cells="./matrixes_run5/BCP-ALL-POOL1-POOL2-MAY24_DBEC_MolsPerCell_MEX/barcodes.tsv",
                  features="./matrixes_run5/BCP-ALL-POOL1-POOL2-MAY24_DBEC_MolsPerCell_MEX/features.tsv") 
expMat <- as.data.frame(t(as.matrix(expMat)))

# ---- separate abseq counts from gene counts, using ending with "pAbO" ----
ab <- expMat %>% dplyr::select(ends_with("pAbO"))

# ---- create seurat object ----
ab <- as.sparse(t(ab))
rna <- as.sparse(t(rna))
seuratObj <- CreateSeuratObject(counts = rna, min.cells = 1,
                                min.features = 100, assay = "RNA")

# ---- add cell assignment ----
cell_assign <- fread(file = "./matrixes_run5/BCP-ALL-POOL1-POOL2-MAY24_Sample_Tag_Calls.csv", sep = ",", header = T) %>%
  data.frame(row.names = 1)

# ---- add to the metadata ----
seuratObj <- AddMetaData(object = seuratObj, metadata = cell_assign)

#################################### QC ########################################

# ---- calculate mito percentage ----
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-") 

# ---- set theme and colors for plotting ----
nice_theme1 <- theme(plot.title=element_text(face = "bold",
                                             family = "Helvetica Light",
                                             size=14),
                     axis.text.y=element_text(family = "Helvetica Light"),
                     axis.text.x=element_text(family = "Helvetica Light"),
                     axis.title.x = element_text(family = "Helvetica Light"),
                     axis.title.y = element_text(family = "Helvetica Light"),
                     legend.text = element_text(family = "Helvetica Light"),
                     legend.title = element_text(family = "Helvetica Light",
                                                 face = "bold"))

sample_colors = c("PT15_DX" = "#0ccbee",
                  "PT15_D15" = "#07788d",
                  "Multiplet" = "grey",
                  "Undetermined" ="grey")

################################## plotting ####################################

# ---- number of cells for each sample ----
png(file="./qc_plots_run5/n_cells_x_sample.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(x = Sample_Name, fill = Sample_Name)) +
  ylab("Cells count") +
  xlab("Sample Name") +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of cells per sample") +  
  labs(color = "Sample Name", fill="Sample Name") +
  scale_fill_manual(values = sample_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  nice_theme1
dev.off()

# ---- remove multiplets and undetermined ----
`%notin%` <- Negate(`%in%`) # generate a function to exclude rows  
seuratObj <- subset(seuratObj , subset = Sample_Name %notin% c("Multiplet", "Undetermined"))

# ---- verify number of cells for each sample ----
png(file="./qc_plots_run5/n_cells_x_sample_no_multiplets_undet.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(x = Sample_Name, fill = Sample_Name)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of cells per sample") +
  ylab("Cells count") +
  xlab("Sample Name") +
  labs(color = "Sample Name", fill="Sample Name") +
  scale_fill_manual(values = sample_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  nice_theme1
dev.off()

# ---- visualize the number of UMI counts per cell ----
# (valid cell at least 500 UMI --> 800 UMI?)
png(file="./qc_plots_run5/umi_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(color = Sample_Name, x = nCount_RNA, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800) +
  geom_text(aes(x=800, label="\n Minimum (x=800)", y=1.5), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  ggtitle("Number of UMIs per cell") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- visualize the number of genes per cell ----
# (valid cell at least 300 genes)
png(file="./qc_plots_run5/genes_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(color = Sample_Name, x = nFeature_RNA, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 5000) +
  ggtitle("Number of genes per cell") +
  geom_text(aes(x=500, label="\n Minimum (x=500)", y=1.12), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  geom_text(aes(x=5000, label="\n Maximum (x=5000)", y=1.12), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- visualize the mitochondrial percentage per cell ----
png(file="./qc_plots_run5/mito_perc_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(color = Sample_Name, x = percent.mt, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab("Mito reads percentage") +
  geom_vline(xintercept = 25, linetype = "dashed",colour = "grey60") +
  geom_text(aes(x=25, label="\n25% (suggested)", y=1.5), colour="grey60",
            angle=90, size = 3.5, family= "Helvetica Light") +
  geom_vline(xintercept = 38, colour = "grey60") +
  geom_text(aes(x=38, label="\n38%", y=1.5), colour="grey60",
            angle=90, size = 3.5, family= "Helvetica Light") +
  ggtitle("Mitochondrial percentage") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- vlnplots for recap ----
png(file="./qc_plots_run5/vlnplot_recap.png",
    height = 800, width = 1400, units = "px", res = 150)
v1=VlnPlot(seuratObj, features = c("nFeature_RNA")) + nice_theme1 + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45)) +
  scale_x_discrete(labels=c("SeuratProject"="all cells")) +
  scale_y_continuous(breaks = c(0,500,1000,2000,4000,6000,8000)) +
  geom_hline(yintercept = 500, color="grey60", linetype = "dashed") +
  geom_hline(yintercept = 5000, color="grey60", linetype = "dashed") +
  ggtitle("Number of genes per cell")
v2=VlnPlot(seuratObj, features = c("nCount_RNA")) + nice_theme1 + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45)) +
  scale_x_discrete(labels=c("SeuratProject"="all cells")) +
  scale_y_continuous(breaks = c(0,2500,5000,10000,20000,30000), limits = c(0,30000)) +
  geom_hline(yintercept = 800, color="grey60", linetype = "dashed") +
  ggtitle("Number of UMIs per cell") 
v3=VlnPlot(seuratObj, features = c("percent.mt")) + nice_theme1 + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45)) +
  scale_x_discrete(labels=c("SeuratProject"="all cells")) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  ylim(0,100) +
  geom_hline(yintercept = 38, color="grey60", linetype = "dashed") +
  ggtitle("Mitochondrial percentage") 
ggarrange(v1,v2,v3, ncol = 3, legend = "none")  
dev.off()

# ---- Add number of genes per UMI for each cell to metadata ----
seuratObj$log10GenesPerUMI <- log10(seuratObj$nFeature_RNA) / log10(seuratObj$nCount_RNA)

# ---- visualize overall complexity ----
png(file="./qc_plots_run5/complexity.png", height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(x = log10GenesPerUMI, color = Sample_Name, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("Overall complexity") +
  ylab("Cell density") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- combine quality control metrics ----
png(file="./qc_plots_run5/qc_metrics_combo.png",
    height = 800, width = 1000, units = "px", res = 150)
seuratObj@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point() +
  scale_color_gradient(low = "gray90", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500, colour = "red") +
  geom_hline(yintercept = 300, colour = "red") +
  facet_wrap(~Sample_Name) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(color = "Mito reads \npercentage") +
  nice_theme1 
dev.off()

saveRDS(seuratObj,
        file = "./seurat_objects_run5/bcpall_run5_seurat_unfiltered_step1.rds")

################################# FILTERING ####################################

# ---- filter on cells ----
filtered.seuratObj <- subset(x = seuratObj,
                                 subset = (nCount_RNA >= 800) &
                                          (nFeature_RNA >= 500) &
                                          (nFeature_RNA <= 5000) &
                                          (percent.mt < 38))

# ---- filter on genes ----
counts <- GetAssayData(object = filtered.seuratObj, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 5
filtered_counts <- counts[keep_genes, ] 
final_filtered.seuratObj <- CreateSeuratObject(counts = filtered_counts,
                                               meta.data = filtered.seuratObj@meta.data)

############################# re-assess QC metrics #############################

# ---- verify number of cells for each sample ----
png(file="./qc_plots_run5/n_cells_x_sample_filtered.png",
    height = 800, width = 1000, units = "px", res = 150)
final_filtered.seuratObj@meta.data %>%
  ggplot(aes(x = Sample_Name, fill = Sample_Name)) +
  ylab("Cells count") +
  xlab("Sample Name") +
  geom_bar() +
  theme_classic() +
  ggtitle("Number of cells per sample") +  
  labs(color = "Sample Name", fill="Sample Name") +
  scale_fill_manual(values = sample_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  nice_theme1
dev.off()

# ---- visualize UMI counts per cell ----
png(file="./qc_plots_run5/umi_per_cell_filtered.png",
    height = 800, width = 1000, units = "px", res = 150)
final_filtered.seuratObj@meta.data %>%
  ggplot(aes(color = Sample_Name, x = nCount_RNA, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  scale_x_log10(limits=c(100,1e5)) +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 800) +
  geom_text(aes(x=800, label="\n Minimum (x=800)", y=1.5), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  ggtitle("Number of UMIs per cell") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- visualize the number of genes per cell ----
# (valid cell at least 300 genes)
png(file="./qc_plots_run5/genes_per_cell_filtered.png",
    height = 800, width = 1000, units = "px", res = 150)
final_filtered.seuratObj@meta.data %>%
  ggplot(aes(color = Sample_Name, x = nFeature_RNA, fill = Sample_Name)) +
  geom_density(alpha = 0.2) +
  scale_x_log10(limits=c(500,6000)) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 5000) 
  ggtitle("Number of genes per cell") +
  geom_text(aes(x=500, label="\n Minimum (x=500)", y=1.12), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  geom_text(aes(x=5000, label="\n Maximum (x=5000)", y=1.12), colour="black",
            angle=90, size = 4.5, family= "Helvetica Light") +
  labs(color = "Sample Name", fill="Sample Name") +
  nice_theme1
dev.off()

# ---- add ADT ----
cell_index <- colnames(final_filtered.seuratObj)
ab <- as.data.frame(ab)
ab <- ab %>%
  tibble::rownames_to_column('protein') %>%
  dplyr::select(protein, cell_index) %>%
  column_to_rownames("protein")
head(ab)[c(1,2), c(1:5)] 

# ---- create a new assay for ADT ----
ab_assay <- CreateAssayObject(counts = ab)
final_filtered.seuratObj[["AB"]] <- ab_assay
final_filtered.seuratObj

saveRDS(final_filtered.seuratObj,
        file = "./seurat_objects_run5/bcpall_run5_seurat_filtered_with_ab_step1.rds")

saveRDS(ab, file = "./seurat_objects_run5/bcpall_run5_ab_step_1.rds")


