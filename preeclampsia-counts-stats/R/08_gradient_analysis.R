source("R/00_setup.R")
run_gradient_analysis <- function(){
 cfg <- load_config(); obj <- readRDS("results/tables/imported_qc_objects.rds"); counts <- obj$counts; meta <- obj$meta; meta$group_ordinal <- as.numeric(factor(meta[[cfg$columns$group]], levels=c("Control","Temprana","Tardia")))
 dds <- DESeqDataSetFromMatrix(counts, meta, as.formula(paste("~",cfg$columns$sex,"+",cfg$columns$delivery,"+",cfg$columns$placenta,"+ group_ordinal"))); dds <- DESeq(dds); tab <- as.data.frame(results(dds, name="group_ordinal")) |> rownames_to_column("gene_id") |> arrange(padj); write.csv(tab,"results/deseq2/DEGs_gradient_Control_Temprana_Tardia.csv",row.names=FALSE)
}
