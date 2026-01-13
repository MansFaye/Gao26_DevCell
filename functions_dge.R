library(scuttle)
library(Seurat)
library(scran)
library(edgeR)
library(clusterProfiler)
library(xlsx)
library(ggplot2)

perform_dge_pseudodbulk <- function(x, ident, contrast, gsea_gene_sets, outdir = "."){
  
  dir.create(paste0(outdir, "/analysis"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/plots"), recursive = TRUE, showWarnings = FALSE)

  x <- as.SingleCellExperiment(x)
  
  colData(x) <- colData(x)[,c("orig.ident", ident)]
  
  summed <- aggregateAcrossCells(x, 
                                 use.assay.type = 1,
                                 ids = colData(x)[,"orig.ident", ident])

  y <- DGEList(counts(summed), samples=colData(summed)$orig.ident)
  
  ##filter lowly expressed (recommanded for limma)
  keep <- filterByExpr(y, group=summed[[ident]])
  y <- y[keep,]
  
  ##see how many genes were kept 
  summary(keep)
  
  ## Create the design matrix and include the technology as a covariate:
  design <- model.matrix(~0 + summed[[ident]])
  
  # change column/rownames names to more simple group names: 
  colnames(design) <- gsub("summed\\[\\[ident\\]\\]", "", colnames(design))
  rownames(design)<-colData(summed)$orig.ident
  
  colnames(design)[colnames(design) == contrast[1]] <- "contrast1"
  colnames(design)[colnames(design) == contrast[2]] <- "contrast2"
  
  # Create contrasts, i.e. specify which groups we want to compare, here we want
  # to find genes differentially expressed between cluster 1 and cluster 2.
  contrast.mat <- limma::makeContrasts(contrast1 - contrast2,
                                       levels = design)
  
  
  #FROM HERE IT IS THE SAME AS BEFORE
  dge <- edgeR::calcNormFactors(y)  
  
  #Do limma
  vm <- limma::voom(dge, design = design, plot = TRUE)
  fit <- limma::lmFit(vm, design = design)
  fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
  fit.contrasts <- limma::eBayes(fit.contrasts)
  
  # Show the top differentially expressed genes:
  limma::topTable(fit.contrasts, number = 10, sort.by = "P")
  deg <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
  length(which(deg$adj.P.Val<0.05))
  
  saveRDS(deg, paste0(outdir, "/analysis/deg_", contrast[1], "_", contrast[2], ".rds"))
  write.xlsx(deg, paste0(outdir, "/analysis/deg_", contrast[1], "_", contrast[2], ".xlsx"))
  
  gene_list <- deg$t
  names(gene_list) <- rownames(deg)
  gene_list <- gene_list[order(gene_list, decreasing = TRUE)]
  
  msigdb <- GSEA(geneList = gene_list, 
                 TERM2GENE = gsea_gene_sets[, c("gs_name", "gene_symbol")],
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
  
  if(nrow(msigdb) > 0) {
    p <- dotplot(msigdb, 
                 # showCategory=ids, 
                 title="Enriched MSigDB hallmark gene sets",
                 split=".sign") + 
      facet_grid(.~.sign)
    
    write(paste0("\nWriting to: ", outdir, "/plots/dotplot_gsea_", contrast[1], "_",contrast[2],".pdf"), stdout())
    
    pdf(paste0(outdir, "/plots/dotplot_gsea_", contrast[1], "_",contrast[2],".pdf"), width = 10, height = 10)
    print(p)
    dev.off()
    
    res <- msigdb@result
    write.xlsx(res, paste0(outdir, "/analysis/gsea_", contrast[1],"_",contrast[2],".xlsx"))
  }
  

  
  entrez_ids <- AnnotationDbi::select(hs, 
                                      keys = names(gene_list),
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "SYMBOL")
  
  names(gene_list) <- entrez_ids$ENTREZID[match(names(gene_list), entrez_ids$SYMBOL)]
  
  kegg <- gseKEGG(geneList     = gene_list,
                  organism     = 'hsa',
                  keyType       = "ncbi-geneid",
                  pvalueCutoff = 0.1,
                  pAdjustMethod = "BH")
  
  if(nrow(kegg) > 0) {
    p <- dotplot(kegg, 
                 # showCategory=ids, 
                 title="Enriched MSigDB hallmark gene sets",
                 split=".sign") + 
      facet_grid(.~.sign)
    
    pdf(paste0(outdir, "/plots/dotplot_kegg_", contrast[1],"_",contrast[2],".pdf"), 
        width = 10, height = 10)
    print(p)
    dev.off()
    
    res <- kegg@result
    write.xlsx(res, paste0(outdir, "/analysis/kegg_", contrast[1],"_",contrast[2],".xlsx"))
  }

  
  diffgenes <- list()
  diffgenes[["down"]] <- rownames(deg)[deg$adj.P.Val < 0.05 & deg$logFC < 0]
  diffgenes[["up"]] <- rownames(deg)[deg$adj.P.Val < 0.05 & deg$logFC > 0]
  
  
  for(dir in names(diffgenes)){
    d2gs <- diffgenes[[dir]]
    d2GO <- clusterProfiler::enrichGO(d2gs, 
                                      "org.Hs.eg.db", 
                                      keyType = "SYMBOL", 
                                      universe = rownames(deg),
                                      ont = "BP", qvalueCutoff = 0.05)
    
    if(!is.null(d2GO)){
      d2GOsim <- simplify(d2GO)
      if(nrow(d2GOsim@result) > 1){
        write.xlsx(d2GOsim@result, 
                   paste0(outdir, "/analysis/GO_", contrast[1],"_",contrast[2], "_", dir, ".xlsx"), 
                   row.names = FALSE)
        pdf(paste0(outdir, "/plots/", contrast[1],"_",contrast[2], "_", dir, "_emapplotGO.pdf"))
        p <- enrichplot::emapplot(enrichplot::pairwise_termsim(d2GOsim),
                                  showCategory = 30, cex_label_category = 0.5)
        print(p)
        dev.off()
      }
    }
  }
  
}

differential_abundance <- function(cts, coldata, ident, contrast) {

  y.ab <- DGEList(counts = cts, group = coldata[,ident])
  
  keep <- filterByExpr(y.ab, group = y.ab$samples$group)
  
  y.ab <- y.ab[keep, ncol(y.ab) |> seq()]

  y.ab <- estimateDisp(y.ab, trend="none")
  
  res <- exactTest(y.ab, pair = contrast)
  
  tt <- topTags(res) |> as.data.frame()
  
  return(list(
    counts = cts,
    da_results = tt,
    coldata = coldata
  ))
}
