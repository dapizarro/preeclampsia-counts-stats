source("R/00_setup.R")
run_import_qc <- function(){
  cfg <- load_config(); make_dirs(); counts <- read_counts(cfg$paths$gene_counts); meta <- read_metadata(cfg$paths$metadata); rownames(meta) <- meta[[cfg$columns$sample_id]]
  common <- intersect(colnames(counts), rownames(meta)); if(length(common)<4) stop("Too few matching samples between counts and metadata")
  counts <- filter_low_counts(counts[,common,drop=FALSE], cfg$thresholds$min_count, cfg$thresholds$min_samples); meta <- meta[common,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(counts, meta, design=~1); vsd <- vst(dds, blind=TRUE)
  pca <- plotPCA(vsd, intgroup=c(cfg$columns$group,cfg$columns$placenta), returnData=TRUE); pct <- round(100*attr(pca,"percentVar"))
  p <- ggplot(pca, aes(PC1, PC2, color=.data[[cfg$columns$group]], shape=.data[[cfg$columns$placenta]])) + geom_point(size=3) + theme_bw() + labs(x=paste0("PC1: ",pct[1],"%"), y=paste0("PC2: ",pct[2],"%"))
  ggsave("figures/qc/Figure_1_PCA_group_placenta.pdf", p, width=7, height=5)
  sample_dists <- dist(t(assay(vsd))); pdf("figures/qc/Figure_2_sample_distance_heatmap.pdf", width=8, height=7); pheatmap(as.matrix(sample_dists)); dev.off()
  write.csv(meta,"results/tables/metadata_matched.csv",row.names=FALSE); saveRDS(list(counts=counts, meta=meta, vsd=vsd), "results/tables/imported_qc_objects.rds")
}
