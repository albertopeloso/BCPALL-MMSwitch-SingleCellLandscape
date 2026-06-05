################################################################################
# ONCOPLOT
################################################################################

# ---- Font Configuration ----

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(readxl)
library(showtext)

showtext_auto()

################################################################################
# Data Loading and Preparation
################################################################################

rna <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "RNA")
wgs <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "WGS")
info <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "Metadata")
gene_pathway <- read_excel("./all_wgs_variants_filtered_pass.xlsx", sheet = "GENE_Pathway")

# ---- mutations ----
mutations <- bind_rows(
  rna %>% select(ID, SYMBOL, Consequence) %>% mutate(Source="RNA"),
  wgs %>% select(ID, SYMBOL, Consequence) %>% mutate(Source="WGS"))
mutations <- mutations %>%
  mutate(Gene_ID = paste(SYMBOL, ID, sep="_"))

# ---- alterations ----
alterations <- mutations %>%
  group_by(ID, SYMBOL) %>%
  summarise(Alteration = paste(unique(Consequence), collapse = ";"), .groups = 'drop') %>%
  pivot_wider(names_from = ID, values_from = Alteration) %>%
  column_to_rownames("SYMBOL")

# ---- source matrix ----
source_data <- mutations %>%
  group_by(ID, SYMBOL) %>%
  summarise(Sources = paste(sort(unique(Source)), collapse = "+"), .groups = 'drop') %>%
  pivot_wider(names_from = ID, values_from = Sources) %>%
  column_to_rownames("SYMBOL")

mat <- as.matrix(alterations)
mat[is.na(mat)] <- ""
source_mat <- as.matrix(source_data)
source_mat[is.na(source_mat)] <- ""

################################################################################
# switch=yes [nras,ptpn11,kras] + other switch=yes + switch=no
################################################################################

# ---- NRAS ----
nras_rna_ids <- character(0)
if("NRAS" %in% rownames(mat)) {
  nras_rna_ids <- rna %>%
    filter(SYMBOL == "NRAS") %>%
    pull(ID) %>%
    unique()
}
nras_wgs_ids <- character(0)
if("NRAS" %in% rownames(mat)) {
  nras_wgs_ids <- wgs %>%
    filter(SYMBOL == "NRAS") %>%
    pull(ID) %>%
    unique()
}
# ---- PTPN11 ----
ptpn11_rna_ids <- character(0)
if ("PTPN11" %in% rownames(mat)) {
  ptpn11_rna_ids <- rna %>%
    filter(SYMBOL == "PTPN11") %>%
    pull(ID) %>%
    unique()
}
ptpn11_wgs_ids <- character(0)
if ("PTPN11" %in% rownames(mat)) {
  ptpn11_wgs_ids <- wgs %>%
    filter(SYMBOL == "PTPN11") %>%
    pull(ID) %>%
    unique()
}
# ---- KRAS ----
kras_rna_ids <- character(0)
if ("KRAS" %in% rownames(mat)) {
  kras_rna_ids <- rna %>%
    filter(SYMBOL == "KRAS") %>%
    pull(ID) %>%
    unique()
}
kras_wgs_ids <- character(0)
if ("KRAS" %in% rownames(mat)) {
  kras_wgs_ids <- wgs %>%
    filter(SYMBOL == "KRAS") %>%
    pull(ID) %>%
    unique()
}
# ---- RNA & WGS ----
nras_any_ids <- unique(c(nras_rna_ids, nras_wgs_ids))
ptpn11_any_ids <- unique(c(ptpn11_rna_ids, ptpn11_wgs_ids))
kras_any_ids <- unique(c(kras_rna_ids, kras_wgs_ids))

# ---- switch=yes ----
switch_yes_nras <- info %>%
  filter(SWITCH == "yes", ID %in% nras_any_ids) %>%
  pull(ID)
switch_yes_ptpn11 <- info %>%
  filter(SWITCH == "yes", ID %in% ptpn11_any_ids, !ID %in% nras_any_ids) %>%
  pull(ID)
switch_yes_kras <- info %>%
  filter(SWITCH == "yes", ID %in% kras_any_ids, !ID %in% nras_any_ids, !ID %in% ptpn11_any_ids) %>%
  pull(ID)
patients_with_mut <- colnames(mat)[colSums(mat != "") > 0]
switch_yes_other <- info %>%
  filter(SWITCH == "yes",
         ID %in% patients_with_mut,
         !ID %in% c(nras_any_ids, ptpn11_any_ids, kras_any_ids)) %>%
  pull(ID)
switch_yes_none <- info %>%
  filter(SWITCH == "yes",
         !ID %in% patients_with_mut) %>%
  pull(ID)

# ---- switch=no ----
switch_no <- info %>%
  filter(SWITCH == "no") %>%
  pull(ID)

# ---- final order ----
ordered_ids <- c(switch_yes_nras, switch_yes_ptpn11, switch_yes_kras, switch_yes_other, switch_yes_none, switch_no)

# ---- no mutations ----
patient_no_mut <- setdiff(ordered_ids, colnames(mat))
if(length(patient_no_mut) > 0) {
  empty_cols <- matrix("", nrow = nrow(mat), ncol = length(patient_no_mut))
  colnames(empty_cols) <- patient_no_mut
  mat <- cbind(mat, empty_cols)
}

# ---- ordering info and matrix ----
info <- info[match(ordered_ids, info$ID), ]
mat <- mat[, ordered_ids]


################################################################################
# Pathway
################################################################################

# ---- gene freq ----
gene_freq <- data.frame(
  SYMBOL = rownames(mat),
  freq = rowSums(mat != ""),
  stringsAsFactors = FALSE
)
# ---- pathway order ----
gene_order <- gene_freq %>%
  left_join(gene_pathway, by = "SYMBOL") %>%
  mutate(Pathway = ifelse(is.na(Pathway), "OTHER PATHWAYS", Pathway))

pathway_levels <- unique(gene_order$Pathway)
pathway_levels <- c(
  "Ras Pathway",
  setdiff(pathway_levels, c("RAS PATHWAY", "OTHER PATHWAYS")),
  "Other"
)
gene_order <- gene_order %>%
  mutate(Pathway = factor(Pathway, levels = pathway_levels)) %>%
  arrange(Pathway, desc(freq))

mat <- mat[gene_order$SYMBOL, ]

pathway_split <- gene_order$Pathway
names(pathway_split) <- gene_order$SYMBOL

# ---- RAS Pathway Order
gene_freq <- data.frame(
  SYMBOL = rownames(mat),
  freq = rowSums(mat != ""),
  stringsAsFactors = FALSE
)
gene_order <- gene_freq %>%
  left_join(gene_pathway, by = "SYMBOL") %>%
  mutate(Pathway = ifelse(is.na(Pathway), "Other", Pathway))

ras_order <- c("NRAS", "PTPN11", "KRAS", "RANBP2")

gene_order <- gene_order %>%
  mutate(
    Pathway = factor(Pathway, levels = c("RAS PATHWAY", setdiff(unique(Pathway), "RAS PATHWAY"))),
    SYMBOL = factor(SYMBOL, levels = c(
      ras_order,
      setdiff(SYMBOL[Pathway == "RAS PATHWAY"], ras_order),
      SYMBOL[Pathway != "RAS PATHWAY"]
    ))
  ) %>%
  arrange(Pathway, SYMBOL)

mat <- mat[as.character(gene_order$SYMBOL), ]

pathway_split <- gene_order$Pathway
names(pathway_split) <- gene_order$SYMBOL

################################################################################
# Define Color
################################################################################

# ---- Pathway Color ----
unique_pathways <- unique(pathway_split)
pathway_colors_base <- c(
  "#4D4D4D","#4D4D4D","#4D4D4D","#4D4D4D",
  "#4D4D4D","#4D4D4D","#4D4D4D")
if(length(unique_pathways) > length(pathway_colors_base)) {
  pathway_colors_base <- rep(pathway_colors_base, ceiling(length(unique_pathways)/length(pathway_colors_base)))
}
pathway_colors <- setNames(pathway_colors_base[1:length(unique_pathways)], unique_pathways)

# ---- Alteration Color ----
alterations_present <- unique(unlist(strsplit(mat[mat != ""], ";")))
alterations_present <- alterations_present[alterations_present != ""]
all_col <- c(
  "Missense Mutation" = "#648B2B",
  "Nonsense Mutation" = "#D2963A",
  "Frameshift Insertion" = "#A20033",
  "Frameshift Deletion" = "#92604E",
  "Inframe Insertion" = "#CB7292",
  "Inframe Deletion" = "#367E7F",
  "Splice Site" = "#D8B8EA"
)

col <- all_col[names(all_col) %in% alterations_present]

################################################################################
# Define Alter Fun
################################################################################

# ---- Data Preparation ----
source_mat_full <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
rownames(source_mat_full) <- rownames(mat)
colnames(source_mat_full) <- colnames(mat)


common_rows <- intersect(rownames(mat), rownames(source_mat))
common_cols <- intersect(colnames(mat), colnames(source_mat))
source_mat_full[common_rows, common_cols] <- source_mat[common_rows, common_cols]


source_mat <- source_mat_full

combined_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
rownames(combined_mat) <- rownames(mat)
colnames(combined_mat) <- colnames(mat)

for(i in 1:nrow(mat)) {
  for(j in 1:ncol(mat)) {
    if(mat[i,j] != "") {
      variant_type <- mat[i,j]
      source_type <- if(rownames(mat)[i] %in% rownames(source_mat) && 
                        colnames(mat)[j] %in% colnames(source_mat)) {
        source_mat[i,j]
      } else {
        ""
      }
      if(source_type != "") {
        combined_mat[i,j] <- paste(variant_type, source_type, sep=";")
      } else {
        combined_mat[i,j] <- variant_type
      }
    }
  }
}


combined_mat[] <- gsub("RNA\\+WGS", "RNA;WGS", combined_mat)
alter_fun_combined = list()

alter_fun_combined <- function(x, y, w, h, v) {
  
  h <- 0.9 * h
  alts <- names(which(v))
  n <- length(alts)
  
  # ---- Background ----
  if (n == 0) {
    grid.rect(x, y, w*0.9, h,
              gp = gpar(fill = "#CCCCCC", col = NA))
    return()
  }
  
  # ---- Source And Variant Types ----
  variants <- intersect(alts, names(col))
  sources  <- intersect(alts, c("RNA", "WGS"))
  
  # ---- Variant Types ----
  if (
    length(variants) == 2 &&
    all(c("Inframe Insertion", "Missense Mutation") %in% variants)
  ) {
    
    grid.polygon(
      unit.c(x-0.45*w, x-0.45*w, x+0.45*w),
      unit.c(y - h*0.45, y + h*0.45, y - h*0.45),
      gp = gpar(fill = col["Inframe Insertion"], col = NA)
    )
    grid.polygon(
      unit.c(x+0.46*w, x+0.46*w, x-0.46*w),
      unit.c(y + h*0.49, y - h*0.49, y + h*0.49),
      gp = gpar(fill = col["Missense Mutation"], col = "white")
    )
  } else if (length(variants) >= 1) {
    grid.rect(x, y, w*0.9, h,
              gp = gpar(fill = col[variants[1]], col = NA))
  }
  
  # ---- Source ----
  if ("RNA" %in% sources) {
    grid.rect(x, y, w*0.9, h,
              gp = gpar(fill = NA, col = "black", lwd = 1))
  }
  
  if ("WGS" %in% sources) {
    grid.points(x, y, pch = 20, size = unit(2, "mm"),
                gp = gpar(col = "black"))
  }
}

################################################################################
# Create Oncoplot
################################################################################

library(ComplexHeatmap)

ht_opt(
  legend_title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Karla"),
  legend_labels_gp = gpar(fontsize = 10, fontfamily = "Karla")
)

annotation_data <- info %>%
  select(SWITCH, WGS, RNA) %>%
  rename(
    Group = SWITCH,
    WGS = WGS,
    RNAseq = RNA
  ) %>%
  mutate(
    Group = ifelse(Group == "yes", "mmSWpos", "mmSWneg"),
  )

showtext_auto()

annotation_colors <- list(
  Group = c("mmSWpos" = "#FA9C32", "mmSWneg" = "#46ACC8"),
  WGS = c("yes" = "orchid4", "no" = "#CCCCCC"),
  RNAseq = c("yes" = "gold1", "no" = "#CCCCCC")
)

op <- oncoPrint(
  combined_mat,  
  alter_fun = alter_fun_combined,
  col = col,
  alter_fun_is_vectorized = FALSE,
  show_heatmap_legend = FALSE,  
  column_order = 1:ncol(mat),
  row_split = pathway_split,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Karla"),
  row_gap = unit(2, "mm"),
  heatmap_legend_param = list(
    title = "Variant Type"
  ),
  top_annotation = HeatmapAnnotation(
    column_barplot = anno_oncoprint_barplot(
      border = FALSE,
      height = unit(3, "cm"),
    ),
    Group = annotation_data$Group,
    WGS = annotation_data$WGS,
    RNAseq = annotation_data$RNAseq,
    col = annotation_colors,
    show_legend = c(Group = TRUE, WGS = FALSE, RNAseq = FALSE),
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Karla"),
    annotation_legend_param = list(
      Group = list(
        title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Karla"),
        labels_gp = gpar(fontsize = 10, fontfamily = "Karla")
      )
    )
  ),
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 10, fontface = "italic", fontfamily = "Karla"),
  column_names_gp = gpar(fontsize = 10, fontfamily = "Karla"),
  pct_gp = gpar(fontsize = 8, fontfamily = "Karla"),
  left_annotation = rowAnnotation(
    Pathway = pathway_split,
    col = list(Pathway = pathway_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(0.5, "mm"),
    simple_anno_size = unit(0.5, "mm")
  ),
  right_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot(
      width = unit(4, "cm")
    )
  )
)

################################################################################
# Legend
################################################################################

data_avail_legend <- Legend(
  labels = c("WGS", "RNAseq", "N/A"),
  legend_gp = gpar(fill = c("orchid4", "gold1", "#CCCCCC")),
  title = "Data Availability",
  title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Karla"),
  labels_gp = gpar(fontsize = 10, fontfamily = "Karla")
)
variant_legend <- Legend(
  labels = names(all_col),
  legend_gp = gpar(fill = unname(all_col)),
  title = "Variant Type",
  title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Karla"),
  labels_gp = gpar(fontsize = 10, fontfamily = "Karla")
)

variant_source_legend <- Legend(
  labels = c("RNAseq", "WGS"),
  graphics = list(
    function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, 
                                   gp = gpar(fill = NA, col = "black", lwd = 1.5)),
    function(x, y, w, h) grid.points(x, y, pch = 20, 
                                     size = unit(2, "mm"), 
                                     gp = gpar(col = "black"))
  ),
  title = "Variant Call Source",
  title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Karla"),
  labels_gp = gpar(fontsize = 10, fontfamily = "Karla")
)
################################################################################
# Save
################################################################################

pdf("oncoprint.pdf", width = 14, height = 10, useDingbats = FALSE)

draw(op,
     annotation_legend_side = "right",
     annotation_legend_list = list(data_avail_legend,variant_legend,variant_source_legend),
     merge_legend = TRUE)
dev.off()

ht_opt(RESET = TRUE)

library(magick)
img <- image_read_pdf("oncoprint.pdf", density = 600)
image_write(img, "oncoprint.png")


