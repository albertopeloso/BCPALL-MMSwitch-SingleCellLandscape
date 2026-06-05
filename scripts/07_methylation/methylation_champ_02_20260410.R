
setwd("/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/")

library(ChAMP)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(vioplot)
library(ggdendro)
library(patchwork)

myCombat <- readRDS("/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/myCombat_20260410.rds")
myLoad <- readRDS(file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/myLoad_20260410.rds")
myDMP <- champ.DMP(beta = myCombat, 
                   pheno = myLoad$pd$Sample_Group, 
                   arraytype = "EPIC")

# ---- barplot with mean beta value ----

row_means <- rowMeans(myCombat)
col_means <- colMeans(myCombat)

subgroup <- c(rep("mmSWneg", 6), rep("mmSWpos", 13))
samples <- c("BCPALL#15", "BCPALL#20", "BCPALL#03", "BCPALL#04", "BCPALL#21",
             "BCPALL#08", "BCPALL#09", "BCPALL#12", "BCPALL#01", "BCPALL#13",
             "BCPALL#14", "BCPALL#02", "BCPALL#16", "BCPALL#07","BCPALL#17",
             "BCPALL#18", "BCPALL#06", "BCPALL#05", "BCPALL#19")

mean_beta <- col_means
means_df <- data.frame(samples, mean_beta, subgroup)
means_df$samples <- factor(means_df$samples, levels = samples)

means_df <- means_df %>%
  mutate(sample_num = as.numeric(str_extract(samples, "\\d+"))) %>%
  group_by(subgroup) %>%
  arrange(sample_num, .by_group = TRUE) %>%
  mutate(samples_ordered = factor(samples, levels = samples)) %>%
  ungroup()

barmeans <-  ggplot(means_df, aes(x = samples_ordered,
                                  y = mean_beta,
                                  fill = subgroup)) +
  geom_bar(stat = "identity") +
  ylab(expression(paste("Mean ", beta, "  -value"))) +
  scale_fill_manual(name = "Group",
                    values = c("mmSWneg" = "#57a0d3",
                               "mmSWpos" = "#f4874c"),
                    guide = guide_legend(
                      override.aes = list(
                        shape  = 22,
                        size   = 2)
                    )
  ) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y = element_text(size = 8, 
                                    family = "Karla", 
                                    face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 7,
                                   color = "black", 
                                   family = "Karla"),
        axis.text.y = element_text( size = 7,
                                    color = "black",
                                    family = "Karla"),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 8,
                                   family = "Karla"),
        legend.title = element_text(size = 10,
                                    family = "Karla",
                                    face = "bold"),
        panel.border = element_rect(colour = "black")
  )

barmeans

ggsave("./methylation/plots/meansplot.png",
plot = barmeans, dpi = 300, width = 5, height = 3)

# pvalue between groups
t_test_res <- t.test(mean_beta ~ subgroup, data = means_df)
t_test_res$p.value
# 
# Welch Two Sample t-test
# 
# data:  mean_beta by subgroup
# t = -0.88453, df = 7.9772, p-value = 0.4023
# alternative hypothesis: true difference in means between group mmSWneg and group mmSWpos is not equal to 0
# 95 percent confidence interval:
#   -0.02758816  0.01229680
# sample estimates:
#   mean in group mmSWneg mean in group mmSWpos 
# 0.6021057             0.6097514 


# ---- DMP calculation ----

DMPdf <- myDMP[[1]]
class(DMPdf)

DMPsubset <- subset(DMPdf, deltaBeta > 0.3 | deltaBeta < -0.3)

feature <- DMPsubset %>% 
  group_by(feature) %>% 
  summarise(count = n())
class(feature)
feature <- as.data.frame(feature)


DMPsubsetSWNEG <- subset(DMPsubset, deltaBeta<0)
featureSWNEG <- DMPsubsetSWNEG %>% count(feature)
featureSWNEG <- DMPsubsetSWNEG %>%
  unnest(feature) %>%
  group_by(feature) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

DMPsubsetSWPOS <- subset(DMPsubset, deltaBeta>0)
featureSWPOS <- DMPsubsetSWPOS %>% 
  unnest(feature) %>%
  group_by(feature) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))


featureTOT <- feature %>%
  left_join(featureSWNEG, by = "feature") %>%
  left_join(featureSWPOS, by = "feature")

colnames(featureTOT) <- c("Features", "All_samples", "mmSWneg", "mmSWpos")

featureTOT_long <- featureTOT %>%
  pivot_longer(
    cols = -Features, 
    names_to = "Group",
    values_to = "Value"
  )

# ---- donutplots for genomic regions ----

# ---- swneg ----

featureSWNEG$fraction <- featureSWNEG$count / sum (featureSWNEG$count)
featureSWNEG$ymax <- cumsum(featureSWNEG$fraction)
featureSWNEG$ymin <- c(0, head(featureSWNEG$ymax, n=-1))
featureSWNEG$labelPosition <- (featureSWNEG$ymax + featureSWNEG$ymin) / 2
featureSWNEG$label <- paste0(featureSWNEG$feature, "\n value: ", featureSWNEG$count)

neg_palette <- c(
  "#e8f7fe",
  "#beeaf8",
  "#94dcf5",
  "#6ec8ed",
  "#43b4e3",
  "#1a9bd3",
  "#007dc0",
  "#0065af"
)


donut_neg <- ggplot(featureSWNEG,
                    aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                        fill = feature)) +
  geom_rect() +
  geom_label_repel(
    aes(x = 4, y = labelPosition, label = feature),
    size = 4,
    nudge_x = 1.5,
    direction = "y",
    force = 0.5,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.color = "gray50",
    family = "Karla",
    fontface = "bold",
    box.padding = 0.5
  ) +
  scale_fill_manual(values = neg_palette) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) +
  annotate(
    "text",
    x = 1,
    y = 0.5,
    label = "mmSWneg",
    family = "Karla",
    fontface = "bold",
    size = 4,
    color = "black"
  ) +
  theme_void() +
  theme(legend.position = "none")

donut_neg

ggsave("./methylation/plots/donutneg.png",
       plot = donut_neg, width= 5, height= 4, dpi=600, bg = "white")


# ---- swpos ----

featureSWPOS_ordered <- featureSWPOS[match(featureSWPOS$feature,
                                           featureSWNEG$feature), ]
featureSWPOS_ordered$fraction <- featureSWPOS_ordered$count / sum (featureSWPOS_ordered$count)
featureSWPOS_ordered$ymax <- cumsum(featureSWPOS_ordered$fraction)
featureSWPOS_ordered$ymin <- c(0, head(featureSWPOS_ordered$ymax, n = -1))
featureSWPOS_ordered$labelPosition <- (featureSWPOS_ordered$ymax +
                                         featureSWPOS_ordered$ymin) / 2
featureSWPOS_ordered$label <- paste0(featureSWPOS_ordered$feature,
                                     "\n value: ",
                                     featureSWPOS_ordered$count)


pos_palette <- c(
  "#ffffd4",
  "#fee391",
  "#fec44f",
  "#ffac26",
  "#ff9a2f",
  "#ff7a36",
  "#ef6c00",
  "#e65100"
)


donut_pos <- ggplot(featureSWPOS_ordered, 
                    aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                        fill = feature)) +
  geom_rect() +
  geom_label_repel(
    aes(x = 4, y = labelPosition, label = feature),
    size = 4,
    nudge_x = 1.5,
    direction = "y",
    force = 0.5,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.color = "gray50",
    family = "Karla",
    fontface = "bold",
    box.padding = 0.5
  ) +
  scale_fill_manual(values = pos_palette) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) +
  annotate(
    "text",
    x = 1,
    y = 0.5,
    label = "mmSWpos",
    family = "Karla",
    fontface = "bold",
    size = 4,
    color = "black"
  ) +
  theme_void() +
  theme(legend.position = "none")

donut_pos

ggsave("./methylation/plots/donutpos.png",
       plot = donut_pos, width= 5, height= 4, dpi=600, bg = "white")


# ---- barplots for genomic regions swneg ----

features_layout <- list(
  c("TSS200", "TSS1500"),
  c("1stExon", "ExonBnd"),
  c("5'UTR", "3'UTR"),
  c("Body", "IGR")
)

cgi_levels <- c("island", "opensea", "shelf", "shore")

bar_color_neg <- "#57a0d3"

# ---- function for generating barplots ----

make_cgi_barplot <- function(data, feature_name, bar_color) {
  
  df_cgi <- data %>%
    filter(feature == feature_name) %>%
    count(cgi, name = "n_DMPs") %>%
    rename(cgi_position = cgi) %>%
    mutate(cgi_position = factor(cgi_position, levels = cgi_levels)) %>%
    complete(cgi_position, fill = list(n_DMPs = 0))
  
  ggplot(df_cgi, aes(x = cgi_position, y = n_DMPs)) +
    geom_bar(stat = "identity", fill = bar_color) +
    labs(title = feature_name) +
    theme_bw() +
    theme(
      text = element_text(family = "karla", face = "bold"),
      plot.title = element_text(hjust = 0.5, color = "black"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      legend.position = "none"
    )
}


plots_swneg <- lapply(
  unlist(features_layout),
  function(f) make_cgi_barplot(DMPsubsetSWNEG, f, bar_color_neg)
)

names(plots_swneg) <- unlist(features_layout)

bar_tot_swneg <-
  (plots_swneg$TSS200 | plots_swneg$TSS1500) /
  (plots_swneg$`1stExon` | plots_swneg$ExonBnd) /
  (plots_swneg$`5'UTR` | plots_swneg$`3'UTR`) /
  (plots_swneg$Body | plots_swneg$IGR)

bar_tot_swneg

ggsave("./methylation/plots/bar_tot_swneg.png",
       plot = bar_tot_swneg, width = 5, height = 7, dpi = 600)


# ---- barplots for genomic regions swpos ----

bar_color_pos <- "#f4874c"

plots_swpos <- lapply(
  unlist(features_layout),
  function(f) make_cgi_barplot(DMPsubsetSWPOS, f, bar_color_pos)
)

names(plots_swpos) <- unlist(features_layout)

bar_tot_swpos <-
  (plots_swpos$TSS200 | plots_swpos$TSS1500) /
  (plots_swpos$`1stExon` | plots_swpos$ExonBnd) /
  (plots_swpos$`5'UTR` | plots_swpos$`3'UTR`) /
  (plots_swpos$Body | plots_swpos$IGR)

bar_tot_swpos

ggsave("./methylation/plots/bar_tot_swpos.png",
       plot = bar_tot_swpos, width = 5, height = 7, dpi = 600)


# ---- GSEA ----

# ---- TSS200 CGIs swpos ----

set.seed(1998)
myDMR <- champ.DMR(beta = myCombat,
                   pheno = myLoad$pd$Sample_Group,
                   cores = 24)
DMP_SWPOS_TSS200 <- DMPsubsetSWPOS %>%
  filter(feature == "TSS200")

DMP_SWPOS_TSS200_CGI <- DMPsubsetSWPOS %>%
  filter(feature == "TSS200", cgi == "island")

myGSEA_swPOS_TSS200_CGI <- champ.GSEA(beta = myCombat, 
                                      DMP = DMP_SWPOS_TSS200_CGI,
                                      DMR = myDMR,
                                      arraytype = "EPIC",
                                      adjPval = 0.05,
                                      method = "fisher",
                                      cores = 24)

myGSEA_swPOS_TSS200_CGI_DMP <- myGSEA_swPOS_TSS200_CGI$DMP
myGSEA_swPOS_TSS200_CGI_DMP$log10adjpval <- -log10(myGSEA_swPOS_TSS200_CGI_DMP$adjPval)

myGSEA_swPOS_TSS200_CGI_DMP$Gene_List <-
  gsub("_", " ", myGSEA_swPOS_TSS200_CGI_DMP$Gene_List)

gsea_pos_tss200_cgi <- ggplot(myGSEA_swPOS_TSS200_CGI_DMP[1:10,], aes(
  x = reorder(Gene_List, log10adjpval), 
  y = log10adjpval)) +
  geom_bar(stat = "identity", fill = "#f4874c") +
  coord_flip() +
  labs(title = "",
       x = "Pathway",
       y = "-log10(adjusted p-value)") +
  theme_void() +
  theme(text = element_text(family = "Karla"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8, hjust = 1),
        axis.title.x = element_text(size = 8),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
  )

gsea_pos_tss200_cgi

ggsave("./methylation/plots/gsea_tss200_cgi.png",
       plot = gsea_pos_tss200_cgi, width = 6, height = 3,
       dpi = 600, bg = "white")

write.table(DMP_SWPOS_TSS200_CGI, 
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/DMP_TSS200_CGI.txt",
            sep = "\t")

write.table(myGSEA_swPOS_TSS200_CGI_DMP,
            file = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape_data/methylation/GSEA_on_DMP_TSS200_CGI.txt",
            sep = "\t")

