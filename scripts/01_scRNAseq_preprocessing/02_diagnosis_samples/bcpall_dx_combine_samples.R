
setwd("/media/data/lab/bcpall_sc_final_analyses_2024/")

library(tidyverse)
library(Seurat)
library(data.table)
library(ggplot2)
library(ggpubr)

# ---- read in filtered seurat objects from different runs ----
seuratObj_run1 <- readRDS("./seurat_objects_run1/bcpall_run1_seurat_filtered_with_ab_step1.rds") 
seuratObj_run2 <- readRDS("./seurat_objects_run2/bcpall_run2_seurat_filtered_with_ab_step1.rds") 
seuratObj_run4 <- readRDS("./seurat_objects_run4/bcpall_run4_seurat_filtered_with_ab_step1.rds") 
seuratObj_run5 <- readRDS("./seurat_objects_run5/bcpall_run5_seurat_filtered_with_ab_step1.rds") 

# re-name run1
seuratObj_run1_dx <- seuratObj_run1 

# ---- filter dx on runs ----
`%notin%` <- Negate(`%in%`) # generate a function to exclude rows  
# filter dx on run 2
seuratObj_run2_dx <- subset(seuratObj_run2, subset = Sample_Name %notin% 
                              c("PT06_D15", "PT07_D15", "PT08_D15"))
# filter dx on run 4
seuratObj_run4_dx <- subset(seuratObj_run4, subset = Sample_Name %notin%
                              c("PT13_D15", "PT14_D15"))
# filter dx on run 5
seuratObj_run5_dx <- subset(seuratObj_run5, subset = Sample_Name %notin%
                              c("PT15_D15"))

# add ids for run1
seuratObj_run1_dx@meta.data$orig.ident <- rep("run1",
                                              length(seuratObj_run1_dx@meta.data$Sample_Name))
# add ids for run2
seuratObj_run2_dx@meta.data$orig.ident <- rep("run2",
                                              length(seuratObj_run2_dx@meta.data$Sample_Name))
# add ids for run4
seuratObj_run4_dx@meta.data$orig.ident <- rep("run4",
                                              length(seuratObj_run4_dx@meta.data$Sample_Name))
# add ids for run5
seuratObj_run5_dx@meta.data$orig.ident <- rep("run5",
                                              length(seuratObj_run5_dx@meta.data$Sample_Name))

# ---- generate new combined seurat object for dx ----
dx_combined <- merge(seuratObj_run1_dx, y = c(seuratObj_run2_dx, seuratObj_run4_dx, seuratObj_run5_dx),
                     add.cell.ids = c("r1", "r2", "r4", "r5"), project = "dx_combined")


########################## re-assess QC metrics ################################


# ---- set theme and colors for plotting ----

nice_theme1 <- theme(plot.title=element_text(face = "bold", family = "Helvetica Light", size=14),
                     axis.text.y=element_text(family = "Helvetica Light"),
                     axis.text.x=element_text(family = "Helvetica Light"),
                     axis.title.x = element_text(family = "Helvetica Light"),
                     axis.title.y = element_text(family = "Helvetica Light"),
                     legend.text = element_text(family = "Helvetica Light"),
                     legend.title = element_text(family = "Helvetica Light", face = "bold"))

sample_colors <- c("PT01" = "#AC92EB",
                   "PT02" = "#4FC1E8",
                   "PT03" = "#607f3e",
                   "PT04" = "#FFCE54",
                   "PT05" = "#ED5564",
                   "PT06_DX" = "#92a5eb",
                   "PT07_DX" = "#A0D568",
                   "PT08_DX" = "#ED9255",
                   "PT13_DX" = "#769de4",
                   "PT14_DX" = "#e5c07c", 
                   "PT15_DX" = "#0ccbee")

# ---- verify number of cells for each sample ----
png(file="./qc_plots_merge_dx/n_cells_x_sample.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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
# (valid cell at least 800 UMI)
png(file="./qc_plots_merge_dx/umi_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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
png(file="./qc_plots_merge_dx/genes_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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
png(file="./qc_plots_merge_dx/mito_perc_per_cell.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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
png(file="./qc_plots_merge_dx/vlnplot_recap.png",
    height = 800, width = 1400, units = "px", res = 150)
v1=VlnPlot(dx_combined, features = c("nFeature_RNA")) + nice_theme1 + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45)) +
  scale_x_discrete(labels=c("SeuratProject"="all cells")) +
  scale_y_continuous(breaks = c(0,500,1000,2000,4000,6000,8000)) +
  geom_hline(yintercept = 500, color="grey60", linetype = "dashed") +
  geom_hline(yintercept = 5000, color="grey60", linetype = "dashed") +
  ggtitle("Number of genes per cell")
v2=VlnPlot(dx_combined, features = c("nCount_RNA")) + nice_theme1 + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 45)) +
  scale_x_discrete(labels=c("SeuratProject"="all cells")) +
  scale_y_continuous(breaks = c(0,2500,5000,10000,20000,30000),
                     limits = c(0,30000)) +
  geom_hline(yintercept = 800, color="grey60", linetype = "dashed") +
  ggtitle("Number of UMIs per cell") 
v3=VlnPlot(dx_combined, features = c("percent.mt")) + nice_theme1 + 
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
dx_combined$log10GenesPerUMI <- log10(dx_combined$nFeature_RNA) / log10(dx_combined$nCount_RNA)

# ---- visualize overall complexity ----
png(file="./qc_plots_merge_dx/complexity.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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
png(file="./qc_plots_merge_dx/qc_metrics_combo.png",
    height = 800, width = 1000, units = "px", res = 150)
dx_combined@meta.data %>%
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

saveRDS(dx_combined,
        file = "./seurat_objects_merge_dx/bcpall_merge_dx_seurat_filtered_with_ab_step1.rds")




