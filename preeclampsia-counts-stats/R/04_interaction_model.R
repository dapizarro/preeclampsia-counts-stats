source("R/00_setup.R")
run_interaction_model <- function(){
 cfg <- load_config(); obj <- readRDS("results/tables/imported_qc_objects.rds"); counts <- obj$counts; meta <- obj$meta
 full <- as.formula(paste("~",cfg$columns$sex,"+",cfg$columns$delivery,"+",cfg$columns$placenta,"+",cfg$columns$group,"+",paste0(cfg$columns$placenta,":",cfg$columns$group)))
 red <- as.formula(paste("~",cfg$columns$sex,"+",cfg$columns$delivery,"+",cfg$columns$placenta,"+",cfg$columns$group))
 dds <- DESeqDataSetFromMatrix(counts, meta, full); dds <- DESeq(dds, test="LRT", reduced=red); tab <- as.data.frame(results(dds)) |> rownames_to_column("gene_id") |> arrange(padj); write.csv(tab,"results/deseq2/DEGs_interaction_placenta_group_LRT.csv",row.names=FALSE)
}
