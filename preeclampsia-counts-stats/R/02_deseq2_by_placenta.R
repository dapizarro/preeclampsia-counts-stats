source("R/00_setup.R")
run_deseq2_by_placenta <- function(){
  cfg <- load_config(); obj <- readRDS("results/tables/imported_qc_objects.rds"); counts <- obj$counts; meta <- obj$meta
  contrasts <- yaml::read_yaml("config/contrasts.yaml")$by_placenta; pc <- cfg$columns$placenta; gc <- cfg$columns$group; sc <- cfg$columns$sex; dc <- cfg$columns$delivery
  for(pl in unique(meta[[pc]])){
    samples <- rownames(meta)[meta[[pc]]==pl]; m <- droplevels(meta[samples,,drop=FALSE]); c <- counts[,samples,drop=FALSE]
    dds <- DESeqDataSetFromMatrix(c,m,as.formula(paste("~",sc,"+",dc,"+",gc))); dds[[gc]] <- relevel(factor(dds[[gc]]), ref="Control"); dds <- DESeq(dds); saveRDS(dds,paste0("results/deseq2/dds_genes_",pl,".rds"))
    for(ct in contrasts){
      res <- results(dds, contrast=c(gc,ct$numerator,ct$denominator)); tab <- as.data.frame(res) |> rownames_to_column("gene_id") |> arrange(padj)
      out <- paste0("results/deseq2/DEGs_",pl,"_",ct$name,".csv"); write.csv(tab,out,row.names=FALSE)
      sig <- tab |> filter(!is.na(padj), padj <= cfg$thresholds$fdr, abs(log2FoldChange) >= cfg$thresholds$abs_log2fc); write.csv(sig, sub(".csv","_FDR005_absLFC1.csv",out), row.names=FALSE)
      p <- ggplot(tab, aes(log2FoldChange, -log10(padj))) + geom_point(alpha=.35) + theme_bw() + geom_vline(xintercept=c(-1,1), linetype=2) + geom_hline(yintercept=-log10(.05), linetype=2) + labs(title=paste(pl,ct$name))
      ggsave(paste0("figures/volcano/Volcano_",pl,"_",ct$name,".pdf"), p, width=6, height=5)
    }
  }
}
