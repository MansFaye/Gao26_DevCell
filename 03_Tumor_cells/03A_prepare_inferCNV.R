## This script runs locally but writes to a server
## ../../remote is a directory mounted from the server

sc <- readRDS("analysis/sc_integrated_annotated.rds")

eef <- subset(sc, subset = celltype_main %in% c(
  "Mesothelial cells",
  "Endothelial cells",
  "Fibroblasts",
  "Cycling mesothelial cells"
))
table(eef$orig.ident, eef$celltype_main)

eef <- subset(eef, subset = orig.ident != "SC10")
eef <- subset(eef, subset = celltype_main != "Fibroblasts")

table(eef$orig.ident, eef$celltype_main)

eefc <- eef@assays$RNA@counts

cells <- colnames(eefc)
annot <- paste(eef$orig.ident, eef$celltype_main, sep = "_")

annot[endsWith(annot, "_Endothelial cells")] <- "Endothelial cells"
annot <- gsub("Cycling mesothelial cells", "Mesothelial cells", annot)

# remove groups with n=1
t <- table(annot)
nis1 <- names(t)[t == 1]

annodf <- data.frame(cell = cells, annotation = annot)
if (length(nis1) > 0) {
  annodf <- annodf[-which(annodf$annot %in% nis1), ]
}


eefc <- eefc[, annodf$cell]
saveRDS(eefc, "../../remote/infercnv/eefc.rds")

write.table(annodf,
  "../../remote/infercnv/cell_annotation_infercnv.txt",
  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

gene_order <- read.delim("gencode_v21_gen_pos.complete.txt", header = FALSE)
gene_order$V1 <- do.call(rbind, strsplit(gene_order$V1, "|", fixed = TRUE))[, 1]
gene_order <- gene_order[!duplicated(gene_order$V1), ]
write.table(gene_order, "../../remote/infercnv/gencode_v21_gen_pos.complete.symbols.txt",
  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)
