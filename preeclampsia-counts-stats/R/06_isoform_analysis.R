source("R/00_setup.R")
run_isoform_analysis <- function(){
 cfg <- load_config(); if(!file.exists(cfg$paths$transcript_counts)){ message("No transcript_counts.csv found; skipping isoform module."); return(invisible(NULL)) }
 counts <- read_counts(cfg$paths$transcript_counts); meta <- read_metadata(cfg$paths$metadata); rownames(meta) <- meta[[cfg$columns$sample_id]]; common <- intersect(colnames(counts), rownames(meta)); counts <- filter_low_counts(counts[,common,drop=FALSE], cfg$thresholds$min_count, cfg$thresholds$min_samples); meta <- meta[common,,drop=FALSE]
 dds <- DESeqDataSetFromMatrix(counts, meta, as.formula(paste("~",cfg$columns$sex,"+",cfg$columns$delivery,"+",cfg$columns$placenta,"+",cfg$columns$group))); dds <- DESeq(dds); saveRDS(dds,"results/deseq2/dds_transcripts_global.rds")
}
