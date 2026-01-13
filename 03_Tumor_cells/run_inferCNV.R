#!/usr/bin/env Rscript

# run inside singularity image
library(infercnv)

projdir <- '/data/projects/p682_Single_Cell_RNA_Sequencing_of_MPM/'

eefc <- readRDS(paste0(projdir, "infercnv/eefc.rds"))

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=eefc,
                                               annotations_file=paste0(projdir, "infercnv/cell_annotation_infercnv.txt"),
                                               delim="\t",
                                               gene_order_file=paste0(projdir, "infercnv/gencode_v21_gen_pos.complete.symbols.txt"),
                                               ref_group_names=c("SC18_Mesothelial cells", 
                                                                 "SC20_Mesothelial cells", 
                                                                 "SC21_Mesothelial cells",
                                                                 "SC22_Mesothelial cells"))

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = paste0(projdir, "infercnv/output_infercnv_pdf"),  # dir is auto-created for storing outputs
                              # cluster_by_groups=T,   # cluster
                              denoise=T,
                              HMM = T,
                              # analysis_mode="subclusters",
                              num_threads=24, 
                              output_format = "pdf"
)


