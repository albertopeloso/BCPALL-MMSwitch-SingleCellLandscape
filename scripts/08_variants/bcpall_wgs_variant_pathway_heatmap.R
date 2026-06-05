################################################################################
# HEATMAP
################################################################################


library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(showtext)
library(grid)
library(gtable)

showtext_auto()

################################################################################
# Data Loading And Preparation
################################################################################

data <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "WGS")
gene_list <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "GENE_Pathway")

# ---- Patient Counts ----
num_patients <- length(unique(data$ID))
patients_pos <- data %>% filter(SWITCH == "yes") %>% distinct(ID) %>% nrow()
patients_neg <- data %>% filter(SWITCH == "no") %>% distinct(ID) %>% nrow()

# ---- Calculate Percentages Per Gene ----
gene_subgroup <- data %>%
  group_by(SYMBOL, SWITCH) %>%
  summarise(unique_patients = n_distinct(ID), .groups = "drop") %>% 
  pivot_wider(names_from = SWITCH, values_from = unique_patients, values_fill = 0) %>%
  mutate(
    perc_pos_subgroup = round((yes / patients_pos) * 100, 1),
    perc_neg_subgroup = round((no / patients_neg) * 100, 1)
  )

# ---- Add Pathway Information ----
gene_subgroup_path <- gene_subgroup %>%
  mutate(SYMBOL = toupper(trimws(SYMBOL))) %>%
  left_join(
    gene_list %>% mutate(SYMBOL = toupper(trimws(SYMBOL))),
    by = "SYMBOL"
  ) %>%
  mutate(Pathway = ifelse(is.na(Pathway), "Unknown", Pathway))

################################################################################
# Gene Ordering
################################################################################

desired_order <- c(
  "CYTOSKELETON",
  "DNA DAMAGE AND CELL CYCLE",
  "EPIGENETIC AND CHROMATIN REMODELLING",
  "RNA MACHINERY",
  "RAS PATHWAY",
  "TRANSCRIPTION",
  "OTHER PATHWAYS"
)

# ---- Pathway Colors ----
palette <- c(
  "#61B470", "#E4DD8B", "#D47E57",
  "#EB3636", "#F48FB1", "#CE93D9", "#BBAAA5"
)
names(palette) <- desired_order

# ---- Calculate Maximum Percentage For Sorting ----
gene_subgroup_path$max_perc <- pmax(
  gene_subgroup_path$perc_pos_subgroup,
  gene_subgroup_path$perc_neg_subgroup
)

gene_subgroup_path$Pathway <- factor(gene_subgroup_path$Pathway, levels = desired_order)

# ---- Sort By Pathway And Percentage ----
gene_subgroup_path_sorted <- gene_subgroup_path %>%
  arrange(Pathway, desc(max_perc))

################################################################################
# Create Heatmap 
################################################################################

# ---- Create Matrix Of Values ----
heatmap_data <- t(gene_subgroup_path_sorted[, c("perc_pos_subgroup", "perc_neg_subgroup")])
colnames(heatmap_data) <- gene_subgroup_path_sorted$SYMBOL
rownames(heatmap_data) <- c("mmSWpos", "mmSWneg")

# ---- Column Annotations ----
annotation_col <- data.frame(Pathway = gene_subgroup_path_sorted$Pathway)
rownames(annotation_col) <- gene_subgroup_path_sorted$SYMBOL

annotation_colors <- list(Pathway = palette)

# ---- Calculate Gaps Between Pathways ----
gaps_col <- which(diff(as.numeric(gene_subgroup_path_sorted$Pathway)) != 0)

# ---- Color Scale ----
my_colors <- colorRampPalette(c("white", "#6baed6", "#08306b"))(50)

# ---- Create Heatmap ----

p <- pheatmap(
  heatmap_data,
  color = my_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = gaps_col,
  fontsize_col = 19,
  fontsize_row = 17,
  fontsize = 17,
  angle_col = 90,
  show_colnames = TRUE,
  cellheight = 25,
  cellwidth = 15,
  border_color = NA,
  main = "",
  fontfamily = "Karla",
  fontfamily_row = "Karla",
  fontfamily_col = "Karla"
)

################################################################################
# Modify Gtable
################################################################################

g <- p$gtable

# Apply italic formatting to gene names
for (i in 1:length(g$grobs)) {
  if (inherits(g$grobs[[i]], "text") && !is.null(g$grobs[[i]]$label)) {
    if (any(g$grobs[[i]]$label %in% colnames(heatmap_data))) {
      g$grobs[[i]]$gp$fontface <- "italic"
    }
  }
}

# ---- Add Title To Legend ----
legend_index <- which(g$layout$name == "legend")

if(length(legend_index) > 0) {
  legend_top <- g$layout[legend_index, "t"]
  g <- gtable_add_rows(g, unit(1.5, "cm"), pos = legend_top - 1)
  
  title_grob <- textGrob(
    "Percentage (%)", 
    gp = gpar(fontsize = 17, fontface = "bold", fontfamily = "Karla"),
    hjust = 0,
    x = unit(0.10, "npc")
  )
  
  left_col <- max(1, g$layout[legend_index, "l"] - 1)
  right_col <- min(ncol(g), g$layout[legend_index, "r"] + 1)
  
  g <- gtable_add_grob(
    g, title_grob, 
    t = legend_top, 
    l = left_col,
    r = right_col
  )
}

################################################################################
# Display And Save
################################################################################

# ---- Display ----
grid.newpage()
grid.draw(g)

# ---- Save With Optimized Dimensions ----
png(
  "heatmap_final.png", 
  width = 5000,
  height = 1000,
  res = 300,
  type = "cairo"
)


showtext_begin()
grid.newpage()
grid.draw(g)
showtext_end()

dev.off()

cat("Heatmap saved as 'heatmap_final.png'\n")