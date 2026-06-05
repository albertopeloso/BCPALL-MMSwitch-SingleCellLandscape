
setwd("/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/methylation/")

library(ChAMP)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)

# ---- load IDAT files ----
myLoad <- champ.load(directory = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/idat_files/", 
                     arraytype="EPIC")
myLoad$pd

champ.QC(resultsDir = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/CHAMP_QCimages")
saveRDS(myLoad, file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/myLoad_20260410.rds")

# ---- perform normalization ----
myNorm <- champ.norm(beta = myLoad$beta, 
                     arraytype = "EPIC", 
                     cores = 24,
                     resultsDir = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/CHAMP_Normalization")
saveRDS(myNorm, file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/myNorm_20260410.rds")

# ---- check batch effect and perform correction ----
myLoad$pd$Slide <- as.factor(myLoad$pd$Slide) #use it as factor to avoid bugs
champ.SVD(beta = myNorm %>% as.data.frame(), 
          pd = myLoad$pd, 
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/CHAMP_SVDimages/")

pca <- prcomp(t(myNorm), scale. = TRUE)
pca_df <- data.frame(pca$x, myLoad$pd)
ggplot(pca_df, aes(x = PC1, y = PC2, 
                   color = Sample_Group, 
                   shape = Slide)) +
  geom_point(size = 3) +
  theme_minimal()

myCombat <- champ.runCombat(beta = myNorm,
                            pd = myLoad$pd,
                            batchname = c("Slide"))

# ---- extraction of top1000 variable probes ----
top1000_probes <- names(sort(apply(myCombat, 1, sd),
                                      decreasing = TRUE))[1:1000]
top1000_probes_table <- myCombat[top1000_probes, ]
write.table(top1000_probes_table,
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/top1000_probes.txt",
            sep = "\t")

# ---- mds calculation and plotting ----
dist_matrix <- dist(t(top1000_probes_table))
mds_result <- cmdscale(dist_matrix, k = 2)
mds_df <- as.data.frame(mds_result)

group <- c(rep("mmSWneg", 6), rep("mmSWpos", 13))

samples <- c("BCPALL#15", "BCPALL#20", "BCPALL#03", "BCPALL#04",
             "BCPALL#21", "BCPALL#08", "BCPALL#09", "BCPALL#12",
             "BCPALL#01", "BCPALL#13", "BCPALL#14", "BCPALL#02",
             "BCPALL#16", "BCPALL#07", "BCPALL#17", "BCPALL#18",
             "BCPALL#06", "BCPALL#05", "BCPALL#19")

subtype <- c("DUX4r", "KMT2Ar", "DUX4r", "DUX4r", "DUX4r", "ETV6::RUNX1",
             "DUX4r", "DUX4r", "DUX4r", "DUX4r", "ZNF384r", "DUX4r", "ZNF384r",
             "DUX4r", "DUX4r", "DUX4r", "DUX4r", "ZNF384r", "Hyperdiploid")

mds_result_bind <- cbind(mds_result, group, samples, subtype)
mds_df_bind <- as.data.frame(mds_result_bind)
mds_df_bind$V1 <- as.numeric(mds_df_bind$V1)
mds_df_bind$V2 <- as.numeric(mds_df_bind$V2)

###############################################

mdsplot <- ggplot() +
  
###-------------------------------
# Group ellipses (no legend)
###-------------------------------

stat_ellipse(
  data = mds_df_bind,
  aes(V1, V2, fill = group),
  geom = "polygon",
  type = "norm",
  level = 0.95,
  alpha = 0.15,
  colour = NA,
  show.legend = FALSE
) +
  stat_ellipse(
    data = mds_df_bind,
    aes(V1, V2, colour = group),
    geom = "path",
    type = "norm",
    level = 0.95,
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c("mmSWpos" = "#FA9C32", "mmSWneg" = "#46ACC8"),
    breaks = NULL
  ) +
  
  # -------------------------------
# FAKE GROUP LEGEND (dummy geom)
# -------------------------------
geom_point(
  data = dplyr::distinct(mds_df_bind, group),
  aes(x = -Inf, y = -Inf, fill = group),
  shape = 22,
  size = 4,
  colour = "black",
  alpha = 0,
  stroke = 0.3,
  inherit.aes = FALSE,
  show.legend = TRUE
) +
  
  scale_fill_manual(
    name = "Group",
    values = c("mmSWpos" = "#FA9C32", "mmSWneg" = "#46ACC8"),
    guide = guide_legend(
      order = 1,
      override.aes = list(
        alpha  = 1,
        shape  = 22,
        colour = "black",
        stroke = 0.3,
        size   = 4
      )
    )
  ) +
  new_scale_fill() +
  
  # -------------------------------
# Subtype points + legend
# -------------------------------
geom_point(
  data = mds_df_bind,
  aes(V1, V2, fill = subtype),
  shape = 21,
  size = 3,
  stroke = 0.3,
  color = "black",
  alpha = 0.9
) +
  scale_fill_manual(
    name = "Subtype",
    values = c(
      "DUX4r" = "#FF0000",
      "KMT2Ar" = "#00A08A",
      "ETV6::RUNX1" = "#F7CD5C",
      "ZNF384r" = "#9966CC",
      "Hyperdiploid" = "#8FC5D6"
    ),
    labels = c(
      expression(paste(italic("DUX4-"), "r")),
      expression(paste(italic("KMT2A-"), "r")),
      expression(italic("ETV6::RUNX1")),
      expression(paste(italic("ZNF384-"), "r")),
      "Hyperdiploid"
    ),
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size  = 4,
        stroke = 0.3,
        colour = "black"
      )
    )
  ) +
  
  coord_cartesian(expand = TRUE) +
  labs(
    title = "",
    x = "MDS dimension 1",
    y = "MDS dimension 2"
  ) +
  theme_bw(base_family = "Karla") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm")
  )

ggsave("./plots/mds_plot_20260412.png",
       plot = mdsplot, width = 7, height = 5, dpi = 600)

saveRDS(myCombat,
        file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/myCombat_20260410.rds")

write.table(as.data.frame(myCombat), 
            row.names = TRUE,
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/MatrixSignalGEO.txt",
            sep = "\t")

top1000_probes_samplenames <- top1000_probes_table
names(top1000_probes_samplenames) <- samples
write.table(top1000_probes_samplenames,
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/top1000_probes_samplenames.txt",
            sep = "\t")

