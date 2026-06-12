source("R/00_setup.R")
run_maternal_vs_fetal <- function(){
  cfg <- load_config(); obj <- readRDS("results/tables/imported_qc_objects.rds"); counts <- obj$counts; meta <- obj$meta; tests <- yaml::read_yaml("config/contrasts.yaml")$maternal_fetal
  gc <- cfg$columns$group; pc <- cfg$columns$placenta; sc <- cfg$columns$sex; dc <- cfg$columns$delivery
  for(ct in tests){ samples <- rownames(meta)[meta[[gc]]==ct$group]; m <- droplevels(meta[samples,,drop=FALSE]); c <- counts[,samples,drop=FALSE]; dds <- DESeqDataSetFromMatrix(c,m,as.formula(paste("~",sc,"+",dc,"+",pc))); dds[[pc]] <- relevel(factor(dds[[pc]]), ref="Fetal"); dds <- DESeq(dds); tab <- as.data.frame(results(dds, contrast=c(pc,"Materna","Fetal"))) |> rownames_to_column("gene_id") |> arrange(padj); write.csv(tab,paste0("results/deseq2/DEGs_",ct$name,".csv"),row.names=FALSE) }
}
