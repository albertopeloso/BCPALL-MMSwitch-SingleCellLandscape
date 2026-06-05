setwd("/home/apeloso/bcpall_infercnv_single_patient_v2")

library(infercnv)

sample_names = c("PT05", "PT02", "PT04", "PT03", "PT01", "PT08_DX",
                 "PT07_DX", "PT06_DX", "PT14_DX", "PT13_DX", "PT15_DX")
ref_groups <- c("T cells","B cells","Myelo/ery")

for (sample in sample_names) {

	print(paste("Processing sample", sample))

   	counts_file <- paste0("./counts/dx_", tolower(sample),
   	                      "_counts_infercnv.txt")
    	ann_file    <- paste0("./annotations/dx_", tolower(sample),
    	                      "_seurat_cluster_annotation_infercnv.txt")
    	out_dir     <- paste0(sample, "_output_dir_v2")
    
    	ann <- read.table(ann_file, header = FALSE, stringsAsFactors = FALSE)
    	group_counts <- table(ann$V2)
    	valid_groups <- names(group_counts[group_counts >= 10])
    
    	if (length(valid_groups) == 0) {
        	warning("All groups have fewer than 10 cells for sample ",
        	        sample, ". Skipping sample.")
        	next
    	}
    
    	ann <- ann[ann$V2 %in% valid_groups, ]
    	ann_clean_file <- paste0("./annotations/",
    	                         sample, "_filtered_annotation.txt")
	write.table(ann, ann_clean_file, sep = "\t", quote = FALSE,
	            row.names = FALSE, col.names = FALSE)

	final_ref_groups <- intersect(ref_groups, valid_groups)
	if (length(final_ref_groups) == 0) {
    	warning("No valid reference groups in sample ",
    	        sample, ". Skipping sample.")
    	next
	}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_file,
                                    annotations_file=ann_clean_file,
                                    delim="\t",
                                    gene_order_file="./hg38_gencode_v27.txt",
                                    ref_group_names=final_ref_groups)


infercnv_obj = infercnv::run(infercnv_obj,
			    # analysis_mode='subclusters',
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=out_dir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
			     num_threads=44)

print(paste("Processed sample", sample))

}

