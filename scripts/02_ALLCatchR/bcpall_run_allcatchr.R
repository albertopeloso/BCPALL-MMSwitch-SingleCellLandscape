################################################################################
# ALLCatchR
################################################################################
library(ALLCatchR)

# ---- Read the data ----
gene_data <- read.table("./bcpall_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_data <- na.omit(gene_data)

# ---- Get gene IDs from rownames ----
gene_ids <- rownames(gene_data)

# ---- Filter: keep only rows that do NOT start with "ENS" ----
non_ensembl_rows <- !grepl("^ENS", gene_ids)
gene_data <- gene_data[non_ensembl_rows, ]

# ---- Remove any duplicates in gene names ----
gene_data <- gene_data[!duplicated(rownames(gene_data)), ]

# ---- Prepare for ALLCatchR: move rownames to GeneID column ----
gene_data_out <- gene_data
gene_data_out$GeneID <- rownames(gene_data_out)
gene_data_out <- gene_data_out[, c("GeneID", setdiff(colnames(gene_data_out), "GeneID"))]

# ---- Save the file ----
write.table(gene_data_out, 
            file = "./gene_bcpall_counts_input_allcatchr.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# ---- Run ALLCatchR ----
allcatch("./gene_bcpall_counts_input_allcatchr.txt", ID_class = "symbol", sep = "\t")
